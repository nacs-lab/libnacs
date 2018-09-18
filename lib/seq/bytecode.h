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

#ifndef __NACS_SEQ_BYTECODE_H__
#define __NACS_SEQ_BYTECODE_H__

/**
 * This file defines the compat format used for execution on the FPGA board.
 * This is currently sequential and being executed on the CPU but a future version
 * could be executed on the FPGA direction and exploid more parallelism.
 *
 * The main way to generate bytecode is from a `Sequence`.
 * There is no way provided to manually construct a bytecode since that is not
 * the (main) usecase.
 */

#include "seq.h"

namespace NaCs {
namespace Seq {
namespace ByteCode {

// Each bytecode instruction has a 4 bit opcode followed by a instruction length that's
// a function of the opcode.
enum OpCode : uint8_t {
    TTLAll = 0,
    TTL2 = 1,
    TTL4 = 2,
    TTL5 = 3,
    Wait = 4,
    Wait2 = 5,
    Clock = 5,
    DDSFreq = 6,
    DDSDetFreq2 = 7,
    DDSDetFreq3 = 8,
    DDSDetFreq4 = 9,
    DDSAmp = 10,
    DDSDetAmp = 11,
    DAC = 12,
    DACDet = 13,
};

namespace Inst {

/**
 * Byte code format:
 * TTL all: [#0: 4][t: 4][val: 32] (5 bytes)
 * TTL flip2: [#1: 4][t: 2][val1: 5][val2: 5] (2 bytes, val1 == val2 means single flip)
 * TTL flip4: [#2: 4][val: 5][val: 5][val: 5][val: 5] (3 bytes, len=3)
 * TTL flip5: [#3: 4][t: 3][val: 5][val: 5][val: 5][val: 5][val: 5] (4 bytes)
 * Wait: [#4: 4][exp: 4][t: 16] (3 bytes, len=t * 2^(3 * exp))
 * Wait2: [#5: 4][#1: 1][t: 11] (2 bytes)
 * Clock: [#5: 4][#0: 4][period: 8] (2 bytes)
 * DDS Freq: [#6: 4][chn: 5][freq: 31] (5 bytes)
 * DDS det Freq: [#7: 4][chn: 5][det_freq: 7] (2 bytes)
 * DDS det Freq: [#8: 4][chn: 5][det_freq: 15] (3 bytes)
 * DDS det Freq: [#9: 4][chn: 5][det_freq: 23] (4 bytes)
 * DDS Amp: [#10: 4][#0: 3][chn: 5][amp: 12] (3 bytes)
 * DDS det Amp: [#11: 4][chn: 5][det_amp: 7] (2 bytes)
 * DAC: [#12: 4][#0: 2][chn: 2][val: 16] (3 bytes)
 * DAC det: [#13: 4][chn: 2][val: 10] (2 bytes)
 *
 * The TTL flip and Freq/Amp/DAC det instructions are for compression and
 * to save memory/network bandwidth. They have exactly the same semantics as the full version
 * and uses a user maintained state to generate the full instruction. See `ExeState`.
 * `Wait2` is also for saving space for short waits
 **/

struct NACS_PACKED TTLAll {
    uint8_t op: 4; // 0
    uint8_t t: 4;
    uint32_t val;
};
static_assert(sizeof(TTLAll) == 5, "");

struct NACS_PACKED TTL2 {
    uint8_t op: 4; // 1
    uint8_t t: 2;
    uint16_t val1: 5;
    uint16_t val2: 5;
};
static_assert(sizeof(TTL2) == 2, "");

struct NACS_PACKED TTL4 {
    uint8_t op: 4; // 2
    uint16_t val1: 5;
    uint16_t val2: 5;
    uint16_t val3: 5;
    uint16_t val4: 5;
};
static_assert(sizeof(TTL4) == 3, "");

struct NACS_PACKED TTL5 {
    uint8_t op: 4; // 3
    uint8_t t: 3;
    uint16_t val1: 5;
    uint16_t val2: 5;
    uint16_t val3: 5;
    uint16_t val4: 5;
    uint16_t val5: 5;
};
static_assert(sizeof(TTL5) == 4, "");

struct NACS_PACKED Wait {
    uint8_t op: 4; // 4
    uint8_t exp: 4;
    uint16_t t;
};
static_assert(sizeof(Wait) == 3, "");

struct NACS_PACKED Wait2 {
    uint8_t op: 4; // 5
    uint8_t _0: 1; // 1
    uint16_t t: 11;
};
static_assert(sizeof(Wait2) == 2, "");

struct NACS_PACKED Clock {
    uint8_t op: 4; // 5
    uint8_t _0: 4; // 0
    uint8_t period;
};
static_assert(sizeof(Clock) == 2, "");

struct NACS_PACKED DDSFreq {
    uint8_t op: 4; // 6
    uint16_t chn: 5;
    uint32_t freq: 31;
};
static_assert(sizeof(DDSFreq) == 5, "");

struct NACS_PACKED DDSDetFreq2 {
    uint8_t op: 4; // 7
    uint16_t chn: 5;
    uint32_t freq: 7;
};
static_assert(sizeof(DDSDetFreq2) == 2, "");

struct NACS_PACKED DDSDetFreq3 {
    uint8_t op: 4; // 8
    uint16_t chn: 5;
    uint32_t freq: 15;
};
static_assert(sizeof(DDSDetFreq3) == 3, "");

struct NACS_PACKED DDSDetFreq4 {
    uint8_t op: 4; // 9
    uint16_t chn: 5;
    uint32_t freq: 23;
};
static_assert(sizeof(DDSDetFreq4) == 4, "");

struct NACS_PACKED DDSAmp {
    uint8_t op: 4; // 10
    uint8_t _0: 3;
    uint16_t chn: 5;
    uint16_t amp: 12;
};
static_assert(sizeof(DDSAmp) == 3, "");

struct NACS_PACKED DDSDetAmp {
    uint8_t op: 4; // 11
    uint16_t chn: 5;
    uint16_t amp: 7;
};
static_assert(sizeof(DDSDetAmp) == 2, "");

struct NACS_PACKED DAC {
    uint8_t op: 4; // 12
    uint8_t _0: 2;
    uint16_t chn: 2;
    uint16_t amp: 16;
};
static_assert(sizeof(DAC) == 3, "");

struct NACS_PACKED DACDet {
    uint8_t op: 4; // 13
    uint16_t chn: 2;
    uint16_t amp: 10;
};
static_assert(sizeof(DACDet) == 2, "");

}

static constexpr uint8_t inst_size[14] = {
    5, 2, 3, 4, // TTL
    3, 2, // wait, clock
    5, 2, 3, 4, // DDS Freq
    3, 2, // DDS Amp
    3, 2, // DAC
};

/**
 * Returns the number of instructions.
 *
 * Note that this does not merge wait instructions.
 */
size_t count(const uint8_t *code, size_t code_len);
static inline size_t count(const std::vector<uint8_t> &code)
{
    return count(&code[0], code.size());
}

/**
 * Print all instructions in a human readable format.
 *
 * Similar to most other functions on the bytecode, this merge wait instructions together
 * and also merge wait into TTL instructions.
 */
void print(std::ostream &stm, const uint8_t *code, size_t code_len, uint32_t ttl_mask=0);
static inline void print(std::ostream &stm, const std::vector<uint8_t> &code,
                         uint32_t ttl_mask=0)
{
    print(stm, &code[0], code.size(), ttl_mask);
}

/**
 * Total time it takes to execute the bytecode.
 */
uint64_t total_time(const uint8_t *code, size_t code_len);
static inline uint64_t total_time(const std::vector<uint8_t> &code)
{
    return total_time(&code[0], code.size());
}

// Keep track of user state during execution of bytecode
struct ExeState {
    /**
     * Execute the bytecode and use the user provided `cb` to generate action.
     * Wait will be merged together and into TTL pulses when possible.
     *
     * The `cb` must have the following field/method.
     *
     * * `ttl(uint32_t ttl, uint64_t t)`:
     *
     *    Generate a TTL pulse with the value `ttl` and wait a total of time `t` cycles.
     *
     * * `dds_freq(uint8_t chn, uint32_t freq)`:
     *
     *    Generate a DDS frequency pulse. Should take `50` cycles.
     *
     * * `dds_amp(uint8_t chn, uint16_t amp)`:
     *
     *    Generate a DDS amplitude pulse. Should take `50` cycles.
     *
     * * `dac(uint8_t chn, uint16_t V)`:
     *
     *    Generate a DAC pulse. Should take `45` cycles.
     *
     * * `wait(uint64_t t)`:
     *
     *    Wait a total of time `t` cycles.
     *
     * * `clock(uint8_t period)`:
     *
     *    Generate a clock pulse (`255` is off). Should take 5 cycles.
     */
    template<typename T>
    void run(T &&cb, const uint8_t *code, size_t len);

private:
    // Unaligned and typed load.
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
void ExeState::run(T &&cb, const uint8_t *code, size_t code_len)
{
    for (size_t i = 0; i < code_len;) {
        auto *p = &code[i];
        uint8_t b = *p;
        uint8_t op = b & 0xf;
        auto inst_len = inst_size[op];
        i += inst_len;
        // Look forward until the end of the code or a non-wait pulse is find.
        // Return the total time consumed.
        auto consumeAllWait = [&] {
            uint64_t t = 0;
            while (i < code_len) {
                auto *p2 = &code[i];
                uint8_t b2 = *p2;
                uint8_t op2 = b2 & 0xf;
                if (op2 == OpCode::Wait2) {
                    if (!(b2 & 0x10))
                        break;
                    i += sizeof(Inst::Wait2);
                    auto inst = loadInst<Inst::Wait2>(p2);
                    t += uint64_t(inst.t);
                }
                else if (op2 != OpCode::Wait) {
                    break;
                }
                else {
                    i += sizeof(Inst::Wait);
                    auto inst = loadInst<Inst::Wait>(p2);
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
        case OpCode::TTLAll: {
            auto inst = loadInst<Inst::TTLAll>(p);
            runTTL(inst.val, inst.t + 3);
            break;
        }
        case OpCode::TTL2: {
            auto inst = loadInst<Inst::TTL2>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            if (inst.val2 != inst.val1)
                ttl = ttl ^ (1 << inst.val2);
            runTTL(ttl, inst.t + 3);
            break;
        }
        case OpCode::TTL4: {
            auto inst = loadInst<Inst::TTL4>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            ttl = ttl ^ (1 << inst.val2);
            ttl = ttl ^ (1 << inst.val3);
            if (inst.val4 != inst.val3)
                ttl = ttl ^ (1 << inst.val4);
            runTTL(ttl, 3);
            break;
        }
        case OpCode::TTL5: {
            auto inst = loadInst<Inst::TTL5>(p);
            uint32_t ttl = m_ttl;
            ttl = ttl ^ (1 << inst.val1);
            ttl = ttl ^ (1 << inst.val2);
            ttl = ttl ^ (1 << inst.val3);
            ttl = ttl ^ (1 << inst.val4);
            ttl = ttl ^ (1 << inst.val5);
            runTTL(ttl, inst.t + 3);
            break;
        }
        case OpCode::Wait: {
            auto inst = loadInst<Inst::Wait>(p);
            uint64_t t = uint64_t(inst.t) << (inst.exp * 3);
            cb.wait(t + consumeAllWait());
            break;
        }
        case OpCode::Clock: {
            if (b & 0x10) {
                auto inst = loadInst<Inst::Wait2>(p);
                cb.wait(uint64_t(inst.t) + consumeAllWait());
            }
            else {
                cb.clock(loadInst<Inst::Clock>(p).period);
            }
            break;
        }
        case OpCode::DDSFreq: {
            auto inst = loadInst<Inst::DDSFreq>(p);
            runDDSFreq(inst.chn, inst.freq);
            break;
        }
        case OpCode::DDSDetFreq2: {
            auto inst = loadInst<Inst::DDSDetFreq2>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x40)
                freq = freq | 0xffffff80;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case OpCode::DDSDetFreq3: {
            auto inst = loadInst<Inst::DDSDetFreq3>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x4000)
                freq = freq | 0xffff8000;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case OpCode::DDSDetFreq4: {
            auto inst = loadInst<Inst::DDSDetFreq4>(p);
            uint8_t chn = inst.chn;
            uint32_t freq = inst.freq;
            if (freq & 0x400000)
                freq = freq | 0xff800000;
            runDDSFreq(chn, freq + m_dds[chn].freq);
            break;
        }
        case OpCode::DDSAmp: {
            auto inst = loadInst<Inst::DDSAmp>(p);
            runDDSAmp(inst.chn, inst.amp);
            break;
        }
        case OpCode::DDSDetAmp: {
            auto inst = loadInst<Inst::DDSDetAmp>(p);
            uint8_t chn = inst.chn;
            uint16_t amp = inst.amp;
            if (amp & 0x40)
                amp = amp | 0xff80;
            runDDSAmp(chn, uint16_t(amp + m_dds[chn].amp));
            break;
        }
        case OpCode::DAC: {
            auto inst = loadInst<Inst::DAC>(p);
            runDAC(inst.chn, inst.amp);
            break;
        }
        case OpCode::DACDet: {
            auto inst = loadInst<Inst::DACDet>(p);
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
}

extern "C" {

uint8_t *nacs_seq_bin_to_bytecode(const uint32_t *data, size_t data_len,
                                  size_t *code_len, uint32_t *ttl_mask);
uint64_t nacs_seq_bytecode_total_time(const uint8_t *code, size_t code_len);

}

#endif
