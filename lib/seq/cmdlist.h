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

#ifndef __NACS_SEQ_CMDLIST_H__
#define __NACS_SEQ_CMDLIST_H__

/**
 * This file defines the format used for running simple sequences on the FPGA board.
 * This is not for experiment sequences, which uses the bytecode.
 *
 * Comparing to bytecode, this is not as compact and is not optimized for
 * small changes. It's simpler for decoding and encoding
 * (though the bytecode isn't much harder) and support all operations instead of
 * only the ones that can be used in an experimental sequence.
 */

#include "seq.h"

#include <istream>

namespace NaCs {

class buff_ostream;

namespace Seq {
namespace CmdList {

// Each cmdlist instruction has a 4 bit opcode followed by a instruction length that's
// a function of the opcode.
enum OpCode : uint8_t {
    TTLAll = 0,
    TTL1 = 1,
    Wait = 2,
    Clock = 3,
    DDSFreq = 4,
    DDSAmp = 5,
    DDSPhase = 6,
    DDSDetPhase = 7,
    DDSReset = 8,
    DAC = 9,
};

namespace Inst {

/**
 * Command format:
 * TTL all: [#0: 8][val: 32] (5 bytes)
 * TTL1: [#1: 8][t: 2][val: 1][chn: 5] (2 bytes)
 * Wait: [#2: 8][t: 64] (9 bytes)
 * Clock: [#3: 8][period: 8] (2 bytes)
 * DDS Freq: [#4: 8][chn: 8][freq: 32] (6 bytes)
 * DDS Amp: [#5: 8][chn: 8][amp: 16] (4 bytes)
 * DDS Phase: [#6: 8][chn: 8][phase: 16] (4 bytes)
 * DDS Det Phase: [#7: 8][chn: 8][det_phase: 16] (4 bytes)
 * DDS Reset: [#8: 8][chn: 8] (2 bytes)
 * DAC: [#9: 8][chn: 8][val: 16] (4 bytes)
 **/

struct NACS_PACKED TTLAll {
    OpCode op; // 0
    uint32_t val;
};
static_assert(sizeof(TTLAll) == 5, "");

struct NACS_PACKED TTL1 {
    OpCode op; // 1
    uint8_t t: 2; // Total time is `t + 3`
    uint8_t val: 1;
    uint8_t chn: 5;
};
static_assert(sizeof(TTL1) == 2, "");

struct NACS_PACKED Wait {
    OpCode op; // 2
    uint64_t t;
};
static_assert(sizeof(Wait) == 9, "");

struct NACS_PACKED Clock {
    OpCode op; // 3
    uint8_t period;
};
static_assert(sizeof(Clock) == 2, "");

struct NACS_PACKED DDSFreq {
    OpCode op; // 4
    uint8_t chn;
    uint32_t freq;
};
static_assert(sizeof(DDSFreq) == 6, "");

struct NACS_PACKED DDSAmp {
    OpCode op; // 5
    uint8_t chn;
    uint16_t amp;
};
static_assert(sizeof(DDSAmp) == 4, "");

struct NACS_PACKED DDSPhase {
    OpCode op; // 6
    uint8_t chn;
    uint16_t phase;
};
static_assert(sizeof(DDSPhase) == 4, "");

struct NACS_PACKED DDSDetPhase {
    OpCode op; // 7
    uint8_t chn;
    uint16_t det_phase;
};
static_assert(sizeof(DDSDetPhase) == 4, "");

struct NACS_PACKED DDSReset {
    OpCode op; // 8
    uint8_t chn;
};
static_assert(sizeof(DDSReset) == 2, "");

struct NACS_PACKED DAC {
    OpCode op; // 9
    uint8_t chn;
    uint16_t amp;
};
static_assert(sizeof(DAC) == 4, "");

}

static constexpr uint8_t cmd_size[10] = {
    5, 2, // TTL
    9, 2, // wait, clock
    6, 4, 4, 4, 2, // DDS
    4, // DAC
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
 * Similar to most other functions on the cmdlist, this merge wait instructions together
 * and also merge wait into TTL instructions.
 */
void print(std::ostream &stm, const uint8_t *code, size_t code_len, uint32_t ttl_mask=0);
static inline void print(std::ostream &stm, const std::vector<uint8_t> &code,
                         uint32_t ttl_mask=0)
{
    print(stm, &code[0], code.size(), ttl_mask);
}

/**
 * Total time it takes to execute the cmdlist.
 */
uint64_t total_time(const uint8_t *code, size_t code_len);
static inline uint64_t total_time(const std::vector<uint8_t> &code)
{
    return total_time(&code[0], code.size());
}

// Keep track of user state during execution of cmdlist
struct ExeState {
    /**
     * Execute the cmdlist and use the user provided `cb` to generate action.
     * Wait will be merged together and into TTL pulses when possible.
     *
     * The `cb` must have the following field/method.
     *
     * * `ttl(uint32_t ttl, uint64_t t)`:
     *
     *    Generate a TTL pulse with the value `ttl` and wait a total of time `t` cycles.
     *
     * * `ttl1(uint8_t chn, bool val, uint64_t t)`:
     *
     *    Generate a TTL pulse that set the channel `chn` to `val`,
     *    wait a total of time `t` cycles.
     *
     * * `dds_freq(uint8_t chn, uint32_t freq)`:
     *
     *    Generate a DDS frequency pulse. Should take `50` cycles.
     *
     * * `dds_amp(uint8_t chn, uint16_t amp)`:
     *
     *    Generate a DDS amplitude pulse. Should take `50` cycles.
     *
     * * `dds_phase(uint8_t chn, uint16_t phase)`:
     *
     *    Generate a DDS phase pulse. Should take `50` cycles.
     *
     * * `dds_detphase(uint8_t chn, uint16_t detphase)`:
     *
     *    Generate a DDS det phase pulse. Should take `50` cycles.
     *
     * * `dds_reset(uint8_t chn)`:
     *
     *    Generate a DDS reset pulse. Should take `50` cycles.
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
};

template<typename T>
void ExeState::run(T &&cb, const uint8_t *code, size_t code_len)
{
    for (size_t i = 0; i < code_len;) {
        auto *p = &code[i];
        uint8_t op = *p;
        auto inst_len = cmd_size[op];
        i += inst_len;
        // Look forward until the end of the code or a non-wait pulse is find.
        // Return the total time consumed.
        auto consumeAllWait = [&] {
            uint64_t t = 0;
            while (i < code_len) {
                auto *p2 = &code[i];
                uint8_t op2 = *p2;
                if (op2 != OpCode::Wait)
                    break;
                i += sizeof(Inst::Wait);
                auto inst = loadInst<Inst::Wait>(p2);
                t += inst.t;
            }
            return t;
        };
        auto runTTL = [&] (uint32_t ttl) {
            cb.ttl(ttl, consumeAllWait() + 3);
        };
        auto runTTL1 = [&] (uint8_t chn, bool val, uint64_t t) {
            t += consumeAllWait();
            cb.ttl1(chn, val, t);
        };
        switch (op) {
        case OpCode::TTLAll: {
            auto inst = loadInst<Inst::TTLAll>(p);
            runTTL(inst.val);
            break;
        }
        case OpCode::TTL1: {
            auto inst = loadInst<Inst::TTL1>(p);
            runTTL1(inst.chn, inst.val, inst.t + 3);
            break;
        }
        case OpCode::Wait: {
            auto inst = loadInst<Inst::Wait>(p);
            uint64_t t = inst.t;
            cb.wait(t + consumeAllWait());
            break;
        }
        case OpCode::Clock: {
            cb.clock(loadInst<Inst::Clock>(p).period);
            break;
        }
        case OpCode::DDSFreq: {
            auto inst = loadInst<Inst::DDSFreq>(p);
            cb.dds_freq(inst.chn, inst.freq);
            break;
        }
        case OpCode::DDSAmp: {
            auto inst = loadInst<Inst::DDSAmp>(p);
            cb.dds_amp(inst.chn, inst.amp);
            break;
        }
        case OpCode::DDSPhase: {
            auto inst = loadInst<Inst::DDSPhase>(p);
            cb.dds_phase(inst.chn, inst.phase);
            break;
        }
        case OpCode::DDSDetPhase: {
            auto inst = loadInst<Inst::DDSDetPhase>(p);
            cb.dds_detphase(inst.chn, inst.det_phase);
            break;
        }
        case OpCode::DDSReset: {
            auto inst = loadInst<Inst::DDSReset>(p);
            cb.dds_reset(inst.chn);
            break;
        }
        case OpCode::DAC: {
            auto inst = loadInst<Inst::DAC>(p);
            cb.dac(inst.chn, inst.amp);
            break;
        }
        default:
            abort();
        }
    }
}

uint32_t parse(buff_ostream &output, std::istream &input);

}

}
}

#endif
