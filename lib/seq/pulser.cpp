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

NACS_EXPORT() void
PulsesBuilder::schedule(Sequence &sequence, seq_cb_t seq_cb, Time::Constraints t_cons)
{
    auto &seq = sequence.pulses;
    auto &defaults = sequence.defaults;
    Filter filter;
    sort(seq);
    uint32_t ttl_val = 0;
    auto ttl_it = defaults.find({Channel::TTL, 0});
    if (ttl_it != defaults.end())
        ttl_val = ttl_it->second.val.i32;
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
        bool val = pulse.cb(pulse.t, Val::get<double>((ttl_val & mask) != 0),
                            pulse.len).val.f64 != 0;
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
    seq.resize(to);
    Channel clock_chn{Channel::Type::CLOCK, 0};
    Seq::schedule(*this, seq, t_cons, defaults, sequence.clocks, filter, std::move(seq_cb),
                  [] (auto clock_div) { return Val::get<uint32_t>(clock_div); }, clock_chn);
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

NACS_EXPORT() Sequence
PulsesBuilder::fromBinary(const uint32_t *bin, size_t len)
{
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
        IR::Function func(&bin[cursor], code_len);
        seq[i] = Pulse{t_start, t_len, chn,
                       PulseData(IRPulse(std::move(func)))};
        cursor += code_len;
    }
    if (cursor >= len)
        return Sequence(std::move(seq), std::move(defaults));
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
    return Sequence(std::move(seq), std::move(defaults), std::move(clocks));
}

NACS_EXPORT() Sequence
PulsesBuilder::fromBase64(const uint8_t *data, size_t len)
{
    size_t bin_len = Base64::decode_len(data, len);
    assert(bin_len % 4 == 0);
    std::vector<uint32_t> bin(bin_len / 4);
    Base64::decode((uint8_t*)bin.data(), data, len);
    return fromBinary(&bin[0], bin_len);
}

namespace {

struct DDSState {
    uint32_t freq;
    uint16_t amp;
};

struct PulserState {
    uint64_t cur_t = 0;
    uint32_t cur_ttl = 0;
    DDSState dds[22] = {};
    uint16_t dac[4] = {};
    size_t last_timed_inst = 0;
    uint8_t max_time_left = 0;
    std::vector<uint8_t> code{};

    static constexpr int start_ttl = 0;
    static constexpr int start_ttl_mask = (1 << start_ttl);
    static constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;

    template<typename Inst>
    size_t addInst(Inst inst)
    {
        auto len = code.size();
        code.resize(len + sizeof(inst));
        memcpy(&code[len], &inst, sizeof(inst));
        max_time_left = 0;
        return len;
    }

    void incLastTime(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t b = code[last_timed_inst];
        uint8_t op = b & 0xf;
        uint8_t tmask = 0xf0;
        switch (op) {
        case 0: // TTLAll
            break;
        case 1: // TTL2
            tmask = 0x30;
            break;
        case 3: // TTL5
            tmask = 0x70;
            break;
        default:
            abort();
        }
        t = uint8_t(t + ((b & tmask) >> 4));
        code[last_timed_inst] = uint8_t((b & ~tmask) | ((t << 4) & tmask));
    }

    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        cur_t += dt;
        if (dt <= max_time_left) {
            incLastTime(uint8_t(dt));
            return;
        }
        uint16_t times[16];
        memset(times, 0, sizeof(times));
        for (int i = 15; i >= 0; i--) {
            // The first number not representable by the next one is
            // 0x10000 * 2^(3 * i - 3)
            // or
            // 0x2000 << (3 * i)
            uint64_t tnext_max = uint64_t(0x2000) << (3 * i);
            if (dt < tnext_max)
                continue;
            times[i] = uint16_t(dt >> (3 * i));
            if (dt == 0) {
                break;
            }
            else if (dt <= max_time_left) {
                incLastTime(uint8_t(dt));
                break;
            }
        }
        for (int i = 15; i >= 0; i--) {
            if (!times[i])
                continue;
            addInst(ByteInst::Wait{4, uint8_t(i & 0xf), times[i]});
        }
    }

    uint64_t addTTLReal(uint64_t t, uint32_t ttl)
    {
        if (ttl == cur_ttl)
            return 1;
        addWait(t - cur_t);
        cur_t += 3;
        auto changes = ttl ^ cur_ttl;
        cur_ttl = ttl;
        auto nchgs = __builtin_popcountll(changes);
        assert(nchgs > 0);
        if (nchgs == 1) {
            auto bit = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(ByteInst::TTL2{1, 0, uint8_t(bit & 0x1f),
                        uint8_t(bit & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 2) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(ByteInst::TTL2{1, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 3) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            addInst(ByteInst::TTL4{2, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit3 & 0x1f)});
        }
        else if (nchgs == 4) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit3);
            auto bit4 = __builtin_ffs(changes) - 1;
            addInst(ByteInst::TTL4{2, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit4 & 0x1f)});
        }
        else if (nchgs == 5) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit3);
            auto bit4 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit4);
            auto bit5 = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(ByteInst::TTL5{3, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit5 & 0x1f)});
            max_time_left = 7;
        }
        else {
            last_timed_inst = addInst(ByteInst::TTLAll{0, 0, ttl});
            max_time_left = 15;
        }
        return 3;
    }

    uint64_t addTTLSingle(uint64_t t, uint8_t chn, bool val)
    {
        if (val) {
            return addTTLReal(t, cur_ttl | (1 << chn));
        }
        else {
            return addTTLReal(t, cur_ttl & ~uint32_t(1 << chn));
        }
    }

    uint64_t addTTL(uint64_t t, uint32_t ttl)
    {
        return addTTLReal(t, ttl & ~start_ttl_mask);
    }

    uint64_t addClock(uint64_t t, uint8_t period)
    {
        addWait(t - cur_t);
        cur_t += 5;
        addInst(ByteInst::Clock{5, 0, period});
        return 5;
    }

    uint64_t addDDSFreq(uint64_t t, uint8_t chn, uint32_t freq)
    {
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        if (freq == dds[chn].freq)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint32_t dfreq = freq - dds[chn].freq;
        dds[chn].freq = freq;
        if (dfreq <= 0x3f || dfreq >= 0xffffffc0) {
            addInst(ByteInst::DDSDetFreq2{7, uint8_t(chn & 0x1f), uint8_t(dfreq & 0x7f)});
        }
        else if (dfreq <= 0x3fff || dfreq >= 0xffffc000) {
            addInst(ByteInst::DDSDetFreq3{8, uint8_t(chn & 0x1f), uint16_t(dfreq & 0x7fff)});
        }
        else if (dfreq <= 0x003fffff || dfreq >= 0xffc00000) {
            addInst(ByteInst::DDSDetFreq4{9, uint8_t(chn & 0x1f), uint32_t(dfreq & 0x7fffff)});
        }
        else {
            addInst(ByteInst::DDSFreq{6, uint8_t(chn & 0x1f), uint32_t(freq & 0x7fffffff)});
        }
        return 50;
    }

    uint64_t addDDSAmp(uint64_t t, uint8_t chn, uint16_t amp)
    {
        if (amp > 4095)
            amp = 4095;
        if (amp == dds[chn].amp)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint16_t damp = uint16_t(amp - dds[chn].amp);
        dds[chn].amp = amp;
        if (damp <= 0x3f || damp >= 0xffc0) {
            addInst(ByteInst::DDSDetAmp{11, uint8_t(chn & 0x1f), uint8_t(damp & 0x7f)});
        }
        else {
            addInst(ByteInst::DDSAmp{10, 0, uint8_t(chn & 0x1f), uint16_t(amp & 0xfff)});
        }
        return 50;
    }

    uint64_t addDAC(uint64_t t, uint8_t chn, uint16_t V)
    {
        if (V == dac[chn])
            return 1;
        addWait(t - cur_t);
        cur_t += 45;
        uint16_t dV = uint16_t(V - dac[chn]);
        dac[chn] = V;
        if (dV <= 0x1ff || dV >= 0xfe00) {
            addInst(ByteInst::DACDet{13, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            addInst(ByteInst::DAC{12, 0, uint8_t(chn & 0x3), V});
        }
        return 45;
    }

    uint64_t start(uint64_t cur_t)
    {
        // wait 100us
        cur_t += 10000;
        addTTLSingle(cur_t, start_ttl, true);
        // 1us
        cur_t += 100;
        addTTLSingle(cur_t, start_ttl, false);
        // 5us
        return cur_t + 500;
    }

    uint64_t end(uint64_t cur_t)
    {
        // This is a hack that is believed to make the NI card happy.
        // 1us
        cur_t += 100;
        addClock(cur_t, 59);
        // 30ms
        cur_t += 3000000;
        // Turn off the clock even when it is not used just as a
        // place holder for the end of the sequence.
        return cur_t + addClock(cur_t, 255);
    }

    static constexpr uint16_t getDACVoltData(double volt)
    {
        if (volt >= 10)
            return uint16_t(0);
        if (volt <= -10)
            return uint16_t(0xffff);
        // this is for the DAC8814 chip in SPI0
        double scale = 65535 / 20.0;
        double offset = 10.0;
        return uint16_t(((offset - volt) * scale) + 0.5);
    }
};

}

NACS_EXPORT() std::vector<uint8_t> PulsesBuilder::toByteCode(const Sequence &seq)
{
    PulserState state;

    Seq::PulsesBuilder seq_builder =
        [&] (Seq::Channel chn, Seq::Val val, uint64_t t, uint64_t tlim) -> uint64_t {
        uint64_t mint = 50;
        if (chn.typ == Seq::Channel::TTL) {
            mint = 3;
        }
        else if (chn.typ == Seq::Channel::CLOCK) {
            mint = 5;
        }
        else if (chn.typ == Seq::Channel::DAC) {
            mint = 45;
        }
        if (t + mint > tlim)
            return 0;
        switch (chn.typ) {
        case Seq::Channel::TTL:
            return state.addTTL(t, val.val.i32);
        case Seq::Channel::DDS_FREQ:
            return state.addDDSFreq(t, uint8_t(chn.id),
                                    uint32_t(0.5 + val.val.f64 * state.freq_factor));
        case Seq::Channel::DDS_AMP:
            return state.addDDSAmp(t, uint8_t(chn.id),
                                   uint16_t(val.val.f64 * 4095.0 + 0.5));
        case Seq::Channel::DAC:
            return state.addDAC(t, uint8_t(chn.id), state.getDACVoltData(val.val.f64));
        case Seq::Channel::CLOCK:
            return state.addClock(t, uint8_t(val.val.i32 - 1));
        default:
            throw std::runtime_error("Invalid Pulse.");
        }
        return mint;
    };
    auto seq_cb = [&] (auto &, uint64_t cur_t, Seq::Event evt) {
        if (evt == Seq::Event::start) {
            return state.start(cur_t);
        } else {
            return state.end(cur_t);
        }
    };
    seq_builder.schedule(const_cast<Seq::Sequence&>(seq), seq_cb);

    return state.code;
}

}
}
