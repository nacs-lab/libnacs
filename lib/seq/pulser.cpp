/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

namespace NaCs {
namespace Seq {

NACS_EXPORT void
PulsesBuilder::schedule(std::vector<Pulse> &seq,
                        const std::map<Channel,Val> &defaults,
                        seq_cb_t seq_cb, Time::Constraints t_cons)
{
    Filter filter;
    sort(seq);
    Seq::schedule(*this, seq, t_cons, defaults, filter, std::move(seq_cb));
}

struct IRPulse {
    IRPulse(IR::Function &&func)
        : m_func(std::move(func)),
          m_ctx(m_func)
    {}
    IRPulse(const IR::Function &func)
        : m_func(func),
          m_ctx(m_func)
    {}
    IRPulse(IRPulse &&other)
        : IRPulse(std::move(other.m_func))
    {}
    IRPulse(const IRPulse &other)
        : IRPulse(other.m_func)
    {}
    Val operator()(uint64_t t, Val start, uint64_t)
    {
        double t_start = double(t) * 10e-9;
        m_ctx.reset(0, IR::TagVal(t_start).val);
        m_ctx.reset(1, start.val);
        return m_ctx.eval().val;
    }
private:
    IR::Function m_func;
    IR::EvalContext m_ctx;
};

NACS_EXPORT std::pair<std::vector<Pulse>,std::map<Channel,Val>>
PulsesBuilder::fromBase64(const uint8_t *data, size_t len)
{
    size_t bin_len = Base64::decode_len(data, len);
    assert(bin_len % 4 == 0);
    std::vector<uint32_t> bin(bin_len / 4);
    Base64::decode((uint8_t*)bin.data(), data, len);

    // [TTL default: 4B]
    // [n_non_ttl: 4B]
    // [[[chn_type: 4B][chn_id: 4B][defaults: 8B]] x n_non_ttl]
    // [n_pulses: 4B]
    // [[[chn_type: 4B][chn_id: 4B][t_start: 8B][t_len: 8B]
    //  [[0: 4B][val: 8B] / [code_len: 4B][code: code_len x 4B]]] x n_pulses]

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
        assert(chn_type != Channel::Type::TTL);
        IR::Function func(&bin[cursor], code_len);
        seq[i] = Pulse{t_start, t_len, chn,
                       PulseData(IRPulse(std::move(func)))};
        cursor += code_len;
    }

    return std::make_pair(std::move(seq), std::move(defaults));
}

}
}
