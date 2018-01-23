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

#include <nacs-utils/ir.h>
#include <nacs-seq/seq.h>

#include <type_traits>
#include <functional>
#include <memory>

#include <math.h>

#ifndef __NACS_SEQ_PULSER_H__
#define __NACS_SEQ_PULSER_H__

namespace NaCs {
namespace Seq {

struct Channel {
    enum Type {
        TTL = 1,
        DDS_FREQ = 2,
        DDS_AMP = 3,
        DAC = 4,
        CLOCK = 5
    };
    Type typ;
    int id;
};

bool operator<(const Channel &id1, const Channel id2)
{
    if (id1.typ < id2.typ) {
        return true;
    } else if (id1.typ > id2.typ) {
        return false;
    }
    return id1.id < id2.id;
}

struct Val {
    IR::GenVal val;
    Val(IR::GenVal v=getGenVal<double>(0)) : val(v)
    {}
    template<typename T, typename T2> static inline Val get(T2 val=0)
    {
        Val v;
        v.val = getGenVal<T>(val);
        return v;
    }
    template<typename T> static inline Val get()
    {
        return get<T>(0);
    }
private:
    template<typename T, typename T2>
    static inline std::enable_if_t<std::is_same<T,double>::value, IR::GenVal>
    getGenVal(T2 val)
    {
        IR::GenVal gv;
        gv.f64 = val;
        return gv;
    }
    template<typename T, typename T2>
    static inline std::enable_if_t<std::is_same<T,uint32_t>::value, IR::GenVal>
    getGenVal(T2 val)
    {
        IR::GenVal gv;
        gv.i32 = val;
        return gv;
    }
};

struct PulseData {
    typedef std::function<Val(uint64_t,Val,uint64_t)> func_t;
    PulseData(PulseData&&) = default;
    PulseData &operator=(PulseData&&) = default;
    PulseData()
        : m_val(),
          m_cb()
    {}
    PulseData(const Val &val)
        : m_val(val),
          m_cb()
    {}
    PulseData(Val &&val)
        : m_val(val),
          m_cb()
    {}
    template<typename T>
    PulseData(T &&v)
        : m_val(),
          m_cb(new func_t(std::forward<T>(v)))
    {}
    Val operator()(uint64_t t, Val start, uint64_t len) const
    {
        if (m_cb)
            return (*m_cb)(t, start, len);
        return m_val;
    }
private:
    Val m_val;
    mutable std::unique_ptr<func_t> m_cb;
};

typedef BasePulse<Channel, PulseData> Pulse;

struct Filter {
    bool operator()(Channel)
    {
        return true;
    }
    bool operator()(Channel cid, Val val1, Val val2)
    {
        switch (cid.typ) {
        case Channel::TTL:
            return val1.val.i32 != val2.val.i32;
        case Channel::DDS_FREQ:
            return fabs(val1.val.f64 - val2.val.f64) > 0.2;
        case Channel::DDS_AMP:
            return fabs(val1.val.f64 - val2.val.f64) > 0.00005;
        case Channel::DAC:
            return fabs(val1.val.f64 - val2.val.f64) > 0.0002;
        default:
            return val1.val.f64 != val2.val.f64;
        }
    }
};

struct Sequence {
    std::vector<Pulse> pulses;
    std::map<Channel,Val> defaults;
    std::vector<Clock> clocks;
    Sequence(std::vector<Pulse> &&_pulses, std::map<Channel,Val> &&_defaults)
        : pulses(std::move(_pulses)),
          defaults(std::move(_defaults)),
          clocks{}
    {
    }
    Sequence(std::vector<Pulse> &&_pulses, std::map<Channel,Val> &&_defaults,
             std::vector<Clock> &&_clocks)
        : pulses(std::move(_pulses)),
          defaults(std::move(_defaults)),
          clocks(std::move(_clocks))
    {
    }
    static Sequence fromBase64(const uint8_t *data, size_t len);
    static Sequence fromBinary(const uint32_t *data, size_t len);
};

namespace ByteCode {

/**
 * Byte code format:
 * TTL all: [#0: 4][t: 4][val: 32] (5 bytes)
 * TTL flip2: [#1: 4][t: 2][val1: 5][val2: 5] (2 bytes, val1 == val2 means single flip)
 * TTL flip4: [#2: 4][val: 5][val: 5][val: 5][val: 5] (3 bytes, len=3)
 * TTL flip5: [#3: 4][t: 3][val: 5][val: 5][val: 5][val: 5][val: 5] (4 bytes)
 * Wait: [#4: 4][exp: 4][t: 16] (3 bytes, len=t * 2^(3 * exp))
 * Wait2: [#4: 4][#1: 1][t: 11] (2 bytes)
 * Clock: [#5: 4][#0: 4][period: 8] (2 bytes)
 * DDS Freq: [#6: 4][chn: 5][freq: 31] (5 bytes)
 * DDS det Freq: [#7: 4][chn: 5][det_freq: 7] (2 bytes)
 * DDS det Freq: [#8: 4][chn: 5][det_freq: 15] (3 bytes)
 * DDS det Freq: [#9: 4][chn: 5][det_freq: 23] (4 bytes)
 * DDS Amp: [#10: 4][#0: 3][chn: 5][amp: 12] (3 bytes)
 * DDS det Amp: [#11: 4][chn: 5][det_amp: 7] (2 bytes)
 * DAC: [#12: 4][#0: 2][chn: 2][val: 16] (3 bytes)
 * DAC det: [#13: 4][chn: 2][val: 10] (2 bytes)
 **/

size_t count(const uint8_t *code, size_t code_len);
static inline size_t count(const std::vector<uint8_t> &code)
{
    return count(&code[0], code.size());
}

void print(std::ostream &stm, const uint8_t *code, size_t code_len);
static inline void print(std::ostream &stm, const std::vector<uint8_t> &code)
{
    print(stm, &code[0], code.size());
}

}

namespace ByteInst {

struct __attribute__((__packed__)) TTLAll {
    uint8_t op: 4; // 0
    uint8_t t: 4;
    uint32_t val;
};
static_assert(sizeof(TTLAll) == 5);

struct __attribute__((__packed__)) TTL2 {
    uint8_t op: 4; // 1
    uint8_t t: 2;
    uint16_t val1: 5;
    uint16_t val2: 5;
};
static_assert(sizeof(TTL2) == 2);

struct __attribute__((__packed__)) TTL4 {
    uint8_t op: 4; // 2
    uint16_t val1: 5;
    uint16_t val2: 5;
    uint16_t val3: 5;
    uint16_t val4: 5;
};
static_assert(sizeof(TTL4) == 3);

struct __attribute__((__packed__)) TTL5 {
    uint8_t op: 4; // 3
    uint8_t t: 3;
    uint16_t val1: 5;
    uint16_t val2: 5;
    uint16_t val3: 5;
    uint16_t val4: 5;
    uint16_t val5: 5;
};
static_assert(sizeof(TTL5) == 4);

struct __attribute__((__packed__)) Wait {
    uint8_t op: 4; // 4
    uint8_t exp: 4;
    uint16_t t;
};
static_assert(sizeof(Wait) == 3);

struct __attribute__((__packed__)) Wait2 {
    uint8_t op: 4; // 5
    uint8_t _0: 1; // 1
    uint16_t t: 11;
};
static_assert(sizeof(Wait2) == 2);

struct __attribute__((__packed__)) Clock {
    uint8_t op: 4; // 5
    uint8_t _0: 4; // 0
    uint8_t period;
};
static_assert(sizeof(Clock) == 2);

struct __attribute__((__packed__)) DDSFreq {
    uint8_t op: 4; // 6
    uint16_t chn: 5;
    uint32_t freq: 31;
};
static_assert(sizeof(DDSFreq) == 5);

struct __attribute__((__packed__)) DDSDetFreq2 {
    uint8_t op: 4; // 7
    uint16_t chn: 5;
    uint32_t freq: 7;
};
static_assert(sizeof(DDSDetFreq2) == 2);

struct __attribute__((__packed__)) DDSDetFreq3 {
    uint8_t op: 4; // 8
    uint16_t chn: 5;
    uint32_t freq: 15;
};
static_assert(sizeof(DDSDetFreq3) == 3);

struct __attribute__((__packed__)) DDSDetFreq4 {
    uint8_t op: 4; // 9
    uint16_t chn: 5;
    uint32_t freq: 23;
};
static_assert(sizeof(DDSDetFreq4) == 4);

struct __attribute__((__packed__)) DDSAmp {
    uint8_t op: 4; // 10
    uint8_t _0: 3;
    uint16_t chn: 5;
    uint16_t amp: 12;
};
static_assert(sizeof(DDSAmp) == 3);

struct __attribute__((__packed__)) DDSDetAmp {
    uint8_t op: 4; // 11
    uint16_t chn: 5;
    uint16_t amp: 7;
};
static_assert(sizeof(DDSDetAmp) == 2);

struct __attribute__((__packed__)) DAC {
    uint8_t op: 4; // 12
    uint8_t _0: 2;
    uint16_t chn: 2;
    uint16_t amp: 16;
};
static_assert(sizeof(DAC) == 3);

struct __attribute__((__packed__)) DACDet {
    uint8_t op: 4; // 13
    uint16_t chn: 2;
    uint16_t amp: 10;
};
static_assert(sizeof(DACDet) == 2);

static constexpr uint8_t codelen[14] = {
    5, 2, 3, 4, // TTL
    3, 2, // wait, clock
    5, 2, 3, 4, // DDS Freq
    3, 2, // DDS Amp
    3, 2, // DAC
};

}

struct PulsesBuilder {
    typedef std::function<uint64_t(Channel,Val,uint64_t,uint64_t)> cb_t;
    typedef std::function<uint64_t(PulsesBuilder&,uint64_t,Event)> seq_cb_t;
    template<typename T>
    PulsesBuilder(T &&_cb)
        : cb(std::forward<T>(_cb))
    {}
    uint64_t operator()(Channel chn, Val val, uint64_t t, uint64_t tlim)
    {
        return cb(chn, val, t, tlim);
    }
    __attribute__((deprecated)) static Sequence fromBase64(const uint8_t *data, size_t len);
    void schedule(Sequence &, seq_cb_t seq_cb,
                  Time::Constraints t_cons={50, 40, 4096});
    void schedule(Sequence &&seq, seq_cb_t seq_cb,
                  Time::Constraints t_cons={50, 40, 4096})
    {
        schedule(seq, seq_cb, t_cons);
    }
    static std::vector<uint8_t> toByteCode(const Sequence &seq);
private:
    cb_t cb;
};

struct ByteCodeExeState {
    template<typename T>
    void runByteCode(T &&cb, const uint8_t *code, size_t len);

private:
    template<typename T>
    T loadInst(const uint8_t *code, size_t idx=0)
    {
        T v;
        memcpy(&v, &code[idx], sizeof(T));
        return v;
    }
    struct DDS {
        uint32_t freq = 0;
        uint16_t amp = 0;
    };
    uint32_t m_ttl = 0;
    DDS m_dds[22] = {};
    uint16_t m_dac[4] = {};
};

template<typename T>
void ByteCodeExeState::runByteCode(T &&cb, const uint8_t *code, size_t code_len)
{
    for (size_t i = 0; i < code_len;) {
        auto *p = &code[i];
        uint8_t b = *p;
        uint8_t op = b & 0xf;
        auto inst_len = ByteInst::codelen[op];
        i += inst_len;
        auto consumeAllWait = [&] {
            uint64_t t = 0;
            while (true) {
                auto *p2 = &code[i];
                uint8_t b2 = *p2;
                uint8_t op2 = b2 & 0xf;
                if (op2 == 5) {
                    if (!(b2 & 0x10))
                        break;
                    i += sizeof(ByteInst::Wait2);
                    auto inst = loadInst<ByteInst::Wait2>(p2);
                    t += uint64_t(inst.t);
                }
                else if (op2 != 4) {
                    break;
                }
                else {
                    i += sizeof(ByteInst::Wait);
                    auto inst = loadInst<ByteInst::Wait>(p2);
                    t += uint64_t(inst.t) << (inst.exp * 3);
                }
            }
            return t;
        };
        auto runTTL = [&] (uint32_t ttl, uint64_t t) {
            t += consumeAllWait();
            cb.ttl(ttl, t);
            m_ttl = ttl;
        };
        auto runDDSFreq = [&] (uint8_t chn, uint32_t freq) {
            cb.dds_freq(chn, freq);
            m_dds[chn].freq = freq;
        };
        auto runDDSAmp = [&] (uint8_t chn, uint16_t amp) {
            cb.dds_amp(chn, amp & 4095);
            m_dds[chn].amp = amp;
        };
        auto runDAC = [&] (uint8_t chn, uint16_t V) {
            cb.dac(chn, V);
            m_dac[chn] = V;
        };
        switch (op) {
        case 0: {
            auto inst = loadInst<ByteInst::TTLAll>(p);
            runTTL(inst.val, inst.t + 3);
            break;
        }
        case 1: {
            auto inst = loadInst<ByteInst::TTL2>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            if (inst.val2 != inst.val1)
                ttl = ttl ^ (1 << inst.val2);
            runTTL(ttl, inst.t + 3);
            break;
        }
        case 2: {
            auto inst = loadInst<ByteInst::TTL4>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            ttl = ttl ^ (1 << inst.val2);
            ttl = ttl ^ (1 << inst.val3);
            if (inst.val4 != inst.val3)
                ttl = ttl ^ (1 << inst.val4);
            runTTL(ttl, 3);
            break;
        }
        case 3: {
            auto inst = loadInst<ByteInst::TTL5>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            ttl = ttl ^ (1 << inst.val2);
            ttl = ttl ^ (1 << inst.val3);
            ttl = ttl ^ (1 << inst.val4);
            ttl = ttl ^ (1 << inst.val5);
            runTTL(ttl, inst.t + 3);
            break;
        }
        case 4: {
            auto inst = loadInst<ByteInst::Wait>(p);
            uint64_t t = uint64_t(inst.t) << (inst.exp * 3);
            cb.wait(t + consumeAllWait());
            break;
        }
        case 5: {
            if (b & 0x10) {
                auto inst = loadInst<ByteInst::Wait2>(p);
                cb.wait(uint64_t(inst.t) + consumeAllWait());
            }
            else {
                cb.clock(loadInst<ByteInst::Clock>(p).period);
            }
            break;
        }
        case 6: {
            auto inst = loadInst<ByteInst::DDSFreq>(p);
            runDDSFreq(inst.chn, inst.freq);
            break;
        }
        case 7: {
            auto inst = loadInst<ByteInst::DDSDetFreq2>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x40)
                freq = freq | 0xffffff80;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case 8: {
            auto inst = loadInst<ByteInst::DDSDetFreq3>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x4000)
                freq = freq | 0xffff8000;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case 9: {
            auto inst = loadInst<ByteInst::DDSDetFreq4>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x400000)
                freq = freq | 0xff800000;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case 10: {
            auto inst = loadInst<ByteInst::DDSAmp>(p);
            runDDSAmp(inst.chn, inst.amp);
            break;
        }
        case 11: {
            auto inst = loadInst<ByteInst::DDSDetAmp>(p);
            uint8_t chn = inst.chn;
            uint16_t amp = inst.amp;
            if (amp & 0x40)
                amp = amp | 0xff80;
            runDDSAmp(chn, uint16_t(amp + m_dds[chn].amp));
            break;
        }
        case 12: {
            auto inst = loadInst<ByteInst::DAC>(p);
            runDAC(inst.chn, inst.amp);
            break;
        }
        case 13: {
            auto inst = loadInst<ByteInst::DACDet>(p);
            uint8_t chn = inst.chn;
            uint16_t amp = inst.amp;
            if (amp & 0x200)
                amp = amp | 0xfc00;
            runDAC(chn, uint16_t(amp + m_dac[chn]));
            break;
        }
        default:
            abort();
        }
    }
}

}
}

#endif
