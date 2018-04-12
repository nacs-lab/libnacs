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

#include "seq.h"

#include <nacs-utils/utils.h>

#include <type_traits>
#include <iomanip>

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
            seq[i] = Pulse(t_start, t_len, chn, Val::get<double>(val));
            cursor += 9;
            continue;
        }
        cursor += 7;
        auto func = exectx->getFunc<double(double, double)>(IR::Function(&bin[cursor], code_len));
        seq[i] = Pulse(t_start, t_len, chn, std::move(func));
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

static void print_channel(std::ostream &stm, Channel::Type chn_type, int chn_id,
                          bool allow_ttl=true, bool allow_clock=false)
{
    auto invalid_type = [&] { stm << "<invalid channel type" << (int)chn_type << ">"; };
    switch (chn_type) {
    case Channel::DDS_FREQ:
        stm << "dds_freq";
        break;
    case Channel::DDS_AMP:
        stm << "dds_amp";
        break;
    case Channel::DAC:
        stm << "dac";
        break;
    case Channel::TTL:
        if (!allow_ttl) {
            invalid_type();
        }
        else {
            stm << "ttl";
        }
        break;
    case Channel::CLOCK:
        if (!allow_clock) {
            invalid_type();
        }
        else {
            stm << "clock";
        }
        break;
    default:
        invalid_type();
        break;
    }
    stm << "(" << chn_id << ")";
}

NACS_EXPORT() void Sequence::dumpBinary(std::ostream &stm, const uint32_t *bin, size_t len)
{
    stm << "TTL default: 0x" << std::hex << std::setfill('0')
        << std::setw(8) << bin[0] << std::dec << std::endl;
    uint32_t n_non_ttl = bin[1];
    stm << "Channels(non TTL): " << n_non_ttl << std::endl;
    size_t cursor = 2;

    for (uint32_t i = 0;i < n_non_ttl;i++) {
        auto chn_type = Channel::Type(bin[cursor]);
        int chn_id = int(bin[cursor + 1]);
        double val;
        memcpy(&val, &bin[cursor + 2], 8);
        stm << "    " << i + 1 << ". ";
        print_channel(stm, chn_type, chn_id, false);
        stm << std::endl
            << "      default: " << val << std::endl;
        cursor += 4;
    }

    uint32_t n_pulses = bin[cursor];
    cursor++;
    stm << "Pulses: " << n_pulses << std::endl;
    for (uint32_t i = 0;i < n_pulses;i++) {
        auto chn_type = Channel::Type(bin[cursor]);
        int chn_id = int(bin[cursor + 1]);
        double t_startf;
        double t_lenf;
        memcpy(&t_startf, &bin[cursor + 2], 8);
        memcpy(&t_lenf, &bin[cursor + 4], 8);
        stm << "    " << i + 1 << ". ";
        print_channel(stm, chn_type, chn_id, true);
        stm << std::endl
            << "      time: " << t_startf << " + " << t_lenf << std::endl;
        uint32_t code_len = bin[cursor + 6];
        if (code_len == 0) {
            double val;
            memcpy(&val, &bin[cursor + 7], 8);
            stm << "      val: " << val << std::endl;
            cursor += 9;
            continue;
        }
        cursor += 7;
        stm << IR::Function(&bin[cursor], code_len) << std::endl;
        cursor += code_len;
    }
    if (cursor >= len)
        return;
    uint32_t n_clocks = bin[cursor];
    cursor++;
    stm << "Clocks: " << n_clocks << std::endl;
    for (uint32_t i = 0;i < n_clocks;i++) {
        uint64_t t_start_ns;
        uint64_t t_len_ns;
        memcpy(&t_start_ns, &bin[cursor], 8);
        memcpy(&t_len_ns, &bin[cursor + 2], 8);
        uint32_t clock_div = bin[cursor + 4];
        cursor += 5;
        stm << "    " << i + 1 << ". start (ns): " << t_start_ns << std::endl
            << "      len (ns): " << t_len_ns << std::endl
            << "      clock_div: " << clock_div << std::endl;
    }
}

} // Seq
} // NaCs
