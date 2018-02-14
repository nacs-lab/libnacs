/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "bytecode.h"

namespace NaCs {
namespace Seq {

namespace ByteCode {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len)
{
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b & 0xf;
        assert(op < 14);
        auto len = inst_size[op];
        i += len;
        count++;
    }
    return count;
}

namespace {

struct Printer {
    void ttl(uint32_t ttl, uint64_t t)
    {
        stm << "TTL: val=" << std::hex << ttl << std::dec << " t=" << t << std::endl;
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        stm << "DDS Freq: chn=" << int(chn)
            << " freq=" << std::hex << freq << std::dec << std::endl;
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        stm << "DDS Amp: chn=" << int(chn)
            << " amp=" << std::hex << amp << std::dec << std::endl;
    }
    void dac(uint8_t chn, uint16_t V)
    {
        stm << "DAC: chn=" << int(chn) << " V=" << std::hex << V << std::dec << std::endl;
    }
    void wait(uint64_t t)
    {
        stm << "Wait: t=" << t << std::endl;
    }
    void clock(uint8_t period)
    {
        stm << "Clock: period=" << int(period) << std::endl;
    }
    std::ostream &stm;
};

struct TimeKeeper {
    void ttl(uint32_t, uint64_t t)
    {
        total_t += t;
    }
    void dds_freq(uint8_t, uint32_t)
    {
        total_t += 50;
    }
    void dds_amp(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dac(uint8_t, uint16_t)
    {
        total_t += 45;
    }
    void wait(uint64_t t)
    {
        total_t += t;
    }
    void clock(uint8_t)
    {
        total_t += 5;
    }
    uint64_t total_t = 0;
};

}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len)
{
    Printer printer{stm};
    ExeState state;
    state.run(printer, code, code_len);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len)
{
    TimeKeeper keeper;
    ExeState state;
    state.run(keeper, code, code_len);
    return keeper.total_t;
}

namespace {

template<typename Vec>
struct ScheduleState {
    struct DDS {
        uint32_t freq = 0;
        uint16_t amp = 0;
        bool freq_set = false;
        bool amp_set = false;
    };
    struct DAC {
        uint16_t V = 0;
        bool set = false;
    };

    ScheduleState(uint32_t ttl)
        : cur_ttl(ttl)
    {
    }

    uint64_t cur_t = 0;
    uint32_t cur_ttl = 0;
    uint32_t _all_ttl_mask = 0;
    bool ttl_set = false;
    DDS dds[22] = {};
    DAC dac[4] = {};
    size_t last_timed_inst = 0;
    uint8_t max_time_left = 0;
    Vec code;

    static constexpr int start_ttl = 0;
    static constexpr int start_ttl_mask = (1 << start_ttl);
    static constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;

    uint32_t all_ttl_mask() const
    {
        return _all_ttl_mask | start_ttl_mask;
    }

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
        uint16_t times[17];
        memset(times, 0, sizeof(times));
        for (int i = 15; i >= 0; i--) {
            // The first number not representable by the next one is
            // 0x10000 * 2^(3 * i - 3)
            // or
            // 0x2000 << (3 * i)
            if (i > 0) {
                uint64_t tnext_max = uint64_t(0x2000) << (3 * i);
                if (dt < tnext_max) {
                    continue;
                }
            }
            else if (dt <= uint16_t(2047 + max_time_left)) {
                if (max_time_left) {
                    dt -= max_time_left;
                    incLastTime(max_time_left);
                }
                times[16] = uint16_t(dt);
                dt = 0;
                break;
            }
            times[i] = uint16_t(dt >> (3 * i));
            dt -= uint64_t(times[i]) << (3 * i);
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
            addInst(Inst::Wait{OpCode::Wait, uint8_t(i & 0xf), times[i]});
        }
        if (times[16]) {
            addInst(Inst::Wait2{OpCode::Wait2, 1, uint16_t(times[16] & 2047)});
        }
    }

    uint64_t addTTLReal(uint64_t t, uint32_t ttl)
    {
        ttl = ttl & all_ttl_mask();
        cur_ttl = cur_ttl & all_ttl_mask();
        if (ttl == cur_ttl && ttl_set)
            return 1;
        addWait(t - cur_t);
        cur_t += 3;
        auto changes = ttl ^ cur_ttl;
        if (!ttl_set) {
            ttl_set = true;
            changes = ttl;
        }
        cur_ttl = ttl;
        auto nchgs = __builtin_popcountll(changes);
        assert(nchgs > 0);
        if (nchgs == 1) {
            auto bit = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit & 0x1f),
                        uint8_t(bit & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 2) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 3) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            addInst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
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
            addInst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
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
            last_timed_inst = addInst(Inst::TTL5{OpCode::TTL5, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit5 & 0x1f)});
            max_time_left = 7;
        }
        else {
            last_timed_inst = addInst(Inst::TTLAll{OpCode::TTLAll, 0, ttl});
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
        addInst(Inst::Clock{OpCode::Clock, 0, period});
        return 5;
    }

    uint64_t addDDSFreq(uint64_t t, uint8_t chn, double freqf)
    {
        uint32_t freq = uint32_t(0.5 + freqf * freq_factor);
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        if (dds[chn].freq_set && freq == dds[chn].freq)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint32_t dfreq = freq - dds[chn].freq;
        dds[chn].freq = freq;
        dds[chn].freq_set = true;
        if (dfreq <= 0x3f || dfreq >= 0xffffffc0) {
            addInst(Inst::DDSDetFreq2{OpCode::DDSDetFreq2, uint8_t(chn & 0x1f),
                        uint8_t(dfreq & 0x7f)});
        }
        else if (dfreq <= 0x3fff || dfreq >= 0xffffc000) {
            addInst(Inst::DDSDetFreq3{OpCode::DDSDetFreq3, uint8_t(chn & 0x1f),
                        uint16_t(dfreq & 0x7fff)});
        }
        else if (dfreq <= 0x003fffff || dfreq >= 0xffc00000) {
            addInst(Inst::DDSDetFreq4{OpCode::DDSDetFreq4, uint8_t(chn & 0x1f),
                        uint32_t(dfreq & 0x7fffff)});
        }
        else {
            addInst(Inst::DDSFreq{OpCode::DDSFreq, uint8_t(chn & 0x1f),
                        uint32_t(freq & 0x7fffffff)});
        }
        return 50;
    }

    uint64_t addDDSAmp(uint64_t t, uint8_t chn, double ampf)
    {
        uint16_t amp = uint16_t(ampf * 4095.0 + 0.5);
        if (amp > 4095)
            amp = 4095;
        if (dds[chn].amp_set && amp == dds[chn].amp)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint16_t damp = uint16_t(amp - dds[chn].amp);
        dds[chn].amp = amp;
        dds[chn].amp_set = true;
        if (damp <= 0x3f || damp >= 0xffc0) {
            addInst(Inst::DDSDetAmp{OpCode::DDSDetAmp, uint8_t(chn & 0x1f),
                        uint8_t(damp & 0x7f)});
        }
        else {
            addInst(Inst::DDSAmp{OpCode::DDSAmp, 0, uint8_t(chn & 0x1f),
                        uint16_t(amp & 0xfff)});
        }
        return 50;
    }

    uint64_t addDAC(uint64_t t, uint8_t chn, double Vf)
    {
        uint16_t V;
        if (Vf >= 10) {
            V = 0;
        }
        else if (Vf <= -10) {
            V = 0xffff;
        }
        else {
            // this is for the DAC8814 chip in SPI0
            double scale = 65535 / 20.0;
            double offset = 10.0;
            V = uint16_t(((offset - Vf) * scale) + 0.5);
        }
        if (dac[chn].set && V == dac[chn].V)
            return 1;
        addWait(t - cur_t);
        cur_t += 45;
        uint16_t dV = uint16_t(V - dac[chn].V);
        dac[chn].V = V;
        dac[chn].set = true;
        if (dV <= 0x1ff || dV >= 0xfe00) {
            addInst(Inst::DACDet{OpCode::DACDet, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            addInst(Inst::DAC{OpCode::DAC, 0, uint8_t(chn & 0x3), V});
        }
        return 45;
    }

    uint64_t start(uint64_t t)
    {
        // wait 20us
        t += 2000;
        addTTLSingle(t, start_ttl, true);
        // 1us
        t += 100;
        addTTLSingle(t, start_ttl, false);
        // 5us
        return t + 500;
    }

    uint64_t end(uint64_t t)
    {
        // 1us
        t += 100;
        addWait(t - cur_t);
        return t;
    }
};

} // anonymous

template<typename Vec>
Vec SeqToByteCode(Sequence &seq, uint32_t *ttl_mask)
{
    // The bytecode is guaranteed to not enable any channel that is not present in the mask.
    uint32_t ttl_default = 0;
    auto ttl_it = seq.defaults.find({Channel::TTL, 0});
    if (ttl_it != seq.defaults.end())
        ttl_default = ttl_it->second.val.i32;

    ScheduleState<Vec> state(ttl_default);

    PulsesBuilder seq_builder =
        [&] (Channel chn, Val val, uint64_t t, uint64_t tlim) -> uint64_t {
        uint64_t mint = 50;
        if (chn.typ == Channel::TTL) {
            mint = 3;
        }
        else if (chn.typ == Channel::CLOCK) {
            mint = 5;
        }
        else if (chn.typ == Channel::DAC) {
            mint = 45;
        }
        if (t + mint > tlim)
            return 0;
        switch (chn.typ) {
        case Channel::TTL:
            return state.addTTL(t, val.val.i32);
        case Channel::DDS_FREQ:
            return state.addDDSFreq(t, uint8_t(chn.id), val.val.f64);
        case Channel::DDS_AMP:
            return state.addDDSAmp(t, uint8_t(chn.id), val.val.f64);
        case Channel::DAC:
            return state.addDAC(t, uint8_t(chn.id), val.val.f64);
        case Channel::CLOCK:
            return state.addClock(t, uint8_t(val.val.i32 - 1));
        default:
            throw std::runtime_error("Invalid Pulse.");
        }
        return mint;
    };
    auto seq_cb = [&] (auto &, uint64_t cur_t, Event evt) {
        if (evt == Event::start) {
            return state.start(cur_t);
        } else {
            return state.end(cur_t);
        }
    };
    seq_builder.schedule(seq, seq_cb, &state._all_ttl_mask, {50, 40, 4096});

    if (ttl_mask)
        *ttl_mask = state.all_ttl_mask();

    return state.code;
}

} // ByteCode

NACS_EXPORT() std::vector<uint8_t> Sequence::toByteCode(uint32_t *ttl_mask)
{
    return ByteCode::SeqToByteCode<std::vector<uint8_t>>(*this, ttl_mask);
}

namespace {

// Simple wrapper around malloc for C interop.
// Implement an API close enough to `std::vector` for our purpose.
// Note that the destructor of this type does **NOT** free the memory.
template<typename ET>
struct MallocVector {
    MallocVector()
        : m_data(nullptr),
          m_len(0)
    {}
    static_assert(std::is_trivial<ET>::value, "");
    ET &operator[](size_t i)
    {
        return m_data[i];
    }
    const ET &operator[](size_t i) const
    {
        return m_data[i];
    }
    size_t size() const
    {
        return m_len;
    }
    void resize(size_t sz)
    {
        m_data = (ET*)realloc(m_data, sizeof(ET) * sz);
        m_len = sz;
    }
private:
    ET *m_data;
    size_t m_len;
};

}

NACS_EXPORT() uint8_t *Sequence::toByteCode(size_t *sz, uint32_t *ttl_mask)
{
    auto vec = ByteCode::SeqToByteCode<MallocVector<uint8_t>>(*this, ttl_mask);
    *sz = vec.size();
    return &vec[0];
}

}
}

using namespace NaCs::Seq;

extern "C" NACS_EXPORT() uint8_t *nacs_seq_bin_to_bytecode(const uint32_t *data,
                                                           size_t data_len,
                                                           size_t *code_len,
                                                           uint32_t *ttl_mask)
{
    return Sequence::fromBinary(data, data_len).toByteCode(code_len, ttl_mask);
}

extern "C" NACS_EXPORT() uint64_t nacs_seq_bytecode_total_time(const uint8_t *code,
                                                               size_t code_len)
{
    return ByteCode::total_time(code, code_len);
}