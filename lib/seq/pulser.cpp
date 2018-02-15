/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *                                                                       *
 *   This library is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   Lesser General Public License for more details.                     *
 *                                                                       *
 *   You should have received a copy of the GNU Lesser General Public    *
 *   License along with this library. If not,                            *
 *   see <http://www.gnu.org/licenses/>.                                 *
 *************************************************************************/

#include "pulser.h"

#include <nacs-utils/utils.h>
#include <nacs-utils/base64.h>

#include <type_traits>

namespace NaCs {
namespace Seq {

NACS_EXPORT() Sequence
Sequence::fromBinary(const uint32_t *bin, size_t len)
{
    auto exectx = IR::ExeContext::get();
    // [TTL default: 4B]
    // [n_non_ttl: 4B]
    // [[[chn_type: 4B][chn_id: 4B][defaults: 8B]] x n_non_ttl]
    // [n_pulses: 4B]
    // [[[chn_type: 4B][chn_id: 4B][t_start: 8B][t_len: 8B]
    //  [[0: 4B][val: 8B] / [code_len: 4B][code: code_len x 4B]]] x n_pulses]
    // Optional:
    // [[n_clocks: 4B][[[t_start_ns: 8B][t_len_ns: 8B][clock_div: 4B]] x n_clocks]]

    std::map<Channel,Val> defaults;

    uint32_t ttl_default = bin[0];
    defaults[{Channel::TTL, 0}] = Val::get<uint32_t>(ttl_default);
    uint32_t n_non_ttl = bin[1];
    size_t cursor = 2;
    for (uint32_t i = 0;i < n_non_ttl;i++) {
        auto chn_type = Channel::Type(bin[cursor]);
        int chn_id = int(bin[cursor + 1]);
        double val;
        memcpy(&val, &bin[cursor + 2], 8);
        defaults[{chn_type, chn_id}] = Val::get<double>(val);
        cursor += 4;
    }

    uint32_t n_pulses = bin[cursor];
    cursor++;
    std::vector<Pulse> seq(n_pulses);
    for (uint32_t i = 0;i < n_pulses;i++) {
        auto chn_type = Channel::Type(bin[cursor]);
        int chn_id = int(bin[cursor + 1]);
        Channel chn{chn_type, chn_id};
        double t_startf;
        double t_lenf;
        memcpy(&t_startf, &bin[cursor + 2], 8);
        memcpy(&t_lenf, &bin[cursor + 4], 8);
        uint64_t t_start = uint64_t(t_startf / 10e-9);
        uint64_t t_len = uint64_t(t_lenf / 10e-9);
        uint32_t code_len = bin[cursor + 6];
        if (code_len == 0) {
            double val;
            memcpy(&val, &bin[cursor + 7], 8);
            seq[i] = Pulse{t_start, t_len, chn,
                           PulseData(Val::get<double>(val))};
            cursor += 9;
            continue;
        }
        cursor += 7;
        auto func = exectx->getFunc<double(double, double)>(IR::Function(&bin[cursor], code_len));
        seq[i] = Pulse{t_start, t_len, chn, PulseData(std::move(func))};
        cursor += code_len;
    }
    if (cursor >= len)
        return Sequence(std::move(seq), std::move(defaults), {}, std::move(exectx));
    uint32_t n_clocks = bin[cursor];
    cursor++;
    std::vector<Clock> clocks(n_clocks);
    for (uint32_t i = 0;i < n_clocks;i++) {
        auto &clock = clocks[i];
        uint64_t t_start_ns;
        uint64_t t_len_ns;
        memcpy(&t_start_ns, &bin[cursor], 8);
        memcpy(&t_len_ns, &bin[cursor + 2], 8);
        uint32_t clock_div = bin[cursor + 4];
        cursor += 5;
        clock.t = t_start_ns / 10;
        clock.len = t_len_ns / 10;
        clock.div = clock_div;
    }
    return Sequence(std::move(seq), std::move(defaults), std::move(clocks), std::move(exectx));
}

NACS_EXPORT() Sequence
Sequence::fromBase64(const uint8_t *data, size_t len)
{
    size_t bin_len = Base64::decode_len(data, len);
    assert(bin_len % 4 == 0);
    std::vector<uint32_t> bin(bin_len / 4);
    Base64::decode((uint8_t*)bin.data(), data, len);
    return fromBinary(&bin[0], bin.size());
}

NACS_EXPORT() void
PulsesBuilder::schedule(Sequence &sequence, seq_cb_t seq_cb,
                        uint32_t *ttl_mask_out, Time::Constraints t_cons)
{
    auto &seq = sequence.pulses;
    auto &defaults = sequence.defaults;
    sort(seq);

    // Merge TTL pulses that happen at the same time into a single TTL pulse
    uint32_t ttl_val = 0;
    auto ttl_it = defaults.find({Channel::TTL, 0});
    if (ttl_it != defaults.end())
        ttl_val = ttl_it->second.val.i32;
    uint32_t used_ttl_mask = 0;
    uint64_t prev_ttl_t = 0;
    ssize_t prev_ttl_idx = -1;
    size_t to = 0, from = 0;
    for (;from < seq.size();from++, to++) {
        auto &pulse = seq[from];
        if (pulse.chn.typ != Channel::TTL) {
            if (from != to)
                seq[to] = std::move(pulse);
            continue;
        }
        assert(pulse.len == 0);
        assert(pulse.chn.id < 32);
        uint32_t mask = uint32_t(1) << pulse.chn.id;
        used_ttl_mask = used_ttl_mask | mask;
        bool val = pulse.cb(pulse.t, Val::get<double>((ttl_val & mask) != 0)).val.f64 != 0;
        uint32_t new_ttl_val;
        if (val) {
            new_ttl_val = ttl_val | mask;
        } else {
            new_ttl_val = ttl_val & ~mask;
        }
        if (new_ttl_val == ttl_val && prev_ttl_idx != -1) {
            to--;
            continue;
        }
        ttl_val = new_ttl_val;
        if (prev_ttl_idx != -1 && prev_ttl_t == pulse.t) {
            seq[prev_ttl_idx].cb = PulseData(Val::get<uint32_t>(new_ttl_val));
            to--;
        } else {
            prev_ttl_idx = to;
            prev_ttl_t = pulse.t;
            pulse.chn = {Channel::Type::TTL, 0};
            pulse.cb = PulseData(Val::get<uint32_t>(new_ttl_val));
            if (from != to) {
                seq[to] = std::move(pulse);
            }
        }
    }
    if (ttl_mask_out)
        *ttl_mask_out = used_ttl_mask;
    seq.resize(to);
    Channel clock_chn{Channel::Type::CLOCK, 0};
    Seq::schedule(*this, seq, t_cons, defaults, sequence.clocks, Filter{}, std::move(seq_cb),
                  [] (auto clock_div) { return Val::get<uint32_t>(clock_div); }, clock_chn);
}

} // Seq
} // NaCs
