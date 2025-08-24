/*************************************************************************
 *   Copyright (c) 2018 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include <nacs-utils/mem.h>

#ifndef __NACS_SEQ_ZYNQ_CMDLIST_H__
#define __NACS_SEQ_ZYNQ_CMDLIST_H__

/**
 * This file defines the format used for running simple sequences on the ZYNQ board.
 * This is not for experiment sequences, which uses the bytecode.
 *
 * Comparing to bytecode, this is not as compact and is not optimized for
 * small changes. It's simpler for decoding and encoding
 * (though the bytecode isn't much harder) and support all operations instead of
 * only the ones that can be used in an experimental sequence.
 */

#include "pulse_time.h"
#include "seq_shared.h"

#include <istream>
#include <stdexcept>
#include <vector>

namespace NaCs {

class buff_ostream;

namespace Seq::Zynq::CmdList {

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

namespace InstDefs {

/**
 * Command format (ver 3):
 * TTL all: [#0: 8][t: 5][bank: 3][val: 32] (6 bytes)
 * TTL1: [#1: 8][t: 7][val: 1][chn: 8] (3 bytes)
 * Wait: [#2: 8][t: 64] (9 bytes)
 * Clock: [#3: 8][period: 8] (2 bytes)
 * DDS Freq: [#4: 8][chn: 8][freq: 32] (6 bytes)
 * DDS Amp: [#5: 8][chn: 8][amp: 16] (4 bytes)
 * DDS Phase: [#6: 8][chn: 8][phase: 16] (4 bytes)
 * DDS Det Phase: [#7: 8][chn: 8][det_phase: 16] (4 bytes)
 * DDS Reset: [#8: 8][chn: 8] (2 bytes)
 * DAC: [#9: 8][chn: 8][val: 16] (4 bytes)
 *
 * Old command formats:
 * TTL all (ver 1-2): [#0: 8][val: 32] (5 bytes)
 * TTL1 (ver 1-2): [#1: 8][t: 2][val: 1][chn: 5] (2 bytes)
 **/

struct NACS_PACKED TTLAll_v1 {
    OpCode op; // 0
    uint32_t val;
};
static_assert(sizeof(TTLAll_v1) == 5, "");

struct NACS_PACKED TTL1_v1 {
    OpCode op; // 1
    uint8_t t: 2; // Total time is `t + PulseTime::Min`
    uint8_t val: 1;
    uint8_t chn: 5;
};
static_assert(sizeof(TTL1_v1) == 2, "");

struct NACS_PACKED TTLAll_v3 {
    OpCode op; // 0
    uint8_t t: 5;
    uint8_t bank: 3;
    uint32_t val;
};
static_assert(sizeof(TTLAll_v3) == 6, "");

struct NACS_PACKED TTL1_v3 {
    OpCode op; // 1
    uint8_t t: 7; // Total time is `t + PulseTime::Min`
    uint8_t val: 1;
    uint8_t chn;
};
static_assert(sizeof(TTL1_v3) == 3, "");

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

struct Inst_v1 {
    using TTLAll = InstDefs::TTLAll_v1;
    using TTL1 = InstDefs::TTL1_v1;
    using Wait = InstDefs::Wait;
    using Clock = InstDefs::Clock;
    using DDSFreq = InstDefs::DDSFreq;
    using DDSAmp = InstDefs::DDSAmp;
    using DDSPhase = InstDefs::DDSPhase;
    using DDSDetPhase = InstDefs::DDSDetPhase;
    using DDSReset = InstDefs::DDSReset;
    using DAC = InstDefs::DAC;
    static constexpr uint8_t version = 1;
    static constexpr uint8_t cmd_size[10] = {
        sizeof(TTLAll), sizeof(TTL1), // TTL
        sizeof(Wait), sizeof(Clock), // wait, clock
        sizeof(DDSFreq), sizeof(DDSAmp), sizeof(DDSPhase),
        sizeof(DDSDetPhase), sizeof(DDSReset), // DDS
        sizeof(DAC), // DAC
    };
};

struct Inst_v3 {
    using TTLAll = InstDefs::TTLAll_v3;
    using TTL1 = InstDefs::TTL1_v3;
    using Wait = InstDefs::Wait;
    using Clock = InstDefs::Clock;
    using DDSFreq = InstDefs::DDSFreq;
    using DDSAmp = InstDefs::DDSAmp;
    using DDSPhase = InstDefs::DDSPhase;
    using DDSDetPhase = InstDefs::DDSDetPhase;
    using DDSReset = InstDefs::DDSReset;
    using DAC = InstDefs::DAC;
    static constexpr uint8_t version = 3;
    static constexpr uint8_t cmd_size[10] = {
        sizeof(TTLAll), sizeof(TTL1), // TTL
        sizeof(Wait), sizeof(Clock), // wait, clock
        sizeof(DDSFreq), sizeof(DDSAmp), sizeof(DDSPhase),
        sizeof(DDSDetPhase), sizeof(DDSReset), // DDS
        sizeof(DAC), // DAC
    };
};


/**
 * Returns the number of instructions.
 *
 * Note that this does not merge wait instructions.
 */
size_t count(const uint8_t *code, size_t code_len, uint32_t version);
static inline size_t count(const std::vector<uint8_t> &code, uint32_t version)
{
    return count(&code[0], code.size(), version);
}

/**
 * Print all instructions in a human readable format.
 *
 * Similar to most other functions on the cmdlist, this merge wait instructions together
 * and also merge wait into TTL instructions.
 */
void print(std::ostream &stm, const uint8_t *code, size_t code_len,
           uint32_t ttl_mask, uint32_t version);
void print(std::ostream &stm, const uint8_t *code, size_t code_len,
           const std::vector<uint32_t> &ttl_mask, uint32_t version);
template<typename T>
static inline void print(std::ostream &stm, const std::vector<uint8_t> &code,
                         T &&ttl_mask, uint32_t version)
{
    print(stm, &code[0], code.size(), std::forward<T>(ttl_mask), version);
}

/**
 * Total time it takes to execute the cmdlist.
 */
uint64_t total_time(const uint8_t *code, size_t code_len, uint32_t version);
static inline uint64_t total_time(const std::vector<uint8_t> &code, uint32_t version)
{
    return total_time(&code[0], code.size(), version);
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
     *    Generate a DDS frequency pulse. Should take `PulseTime::DDSFreq` cycles.
     *
     * * `dds_amp(uint8_t chn, uint16_t amp)`:
     *
     *    Generate a DDS amplitude pulse. Should take `PulseTime::DDSAmp` cycles.
     *
     * * `dds_phase(uint8_t chn, uint16_t phase)`:
     *
     *    Generate a DDS phase pulse. Should take `PulseTime::DDSPhase` cycles.
     *
     * * `dds_detphase(uint8_t chn, uint16_t detphase)`:
     *
     *    Generate a DDS det phase pulse. Should take `PulseTime::DDSPhase` cycles.
     *
     * * `dds_reset(uint8_t chn)`:
     *
     *    Generate a DDS reset pulse. Should take `PulseTime::DDSReset` cycles.
     *
     * * `dac(uint8_t chn, uint16_t V)`:
     *
     *    Generate a DAC pulse. Should take `PulseTime::DAC` cycles.
     *
     * * `wait(uint64_t t)`:
     *
     *    Wait a total of time `t` cycles.
     *
     * * `clock(uint8_t period)`:
     *
     *    Generate a clock pulse (`255` is off). Should take `PulseTime::Clock` cycles.
     */
    template<typename T>
    void run(T &&cb, const uint8_t *code, size_t len, uint32_t version);

    uint8_t min_time = PulseTime::Min;
private:
    template<typename Inst, typename T>
    void _run(T &&cb, const uint8_t *code, size_t len);
};

template<typename T>
void ExeState::run(T &&cb, const uint8_t *code, size_t code_len, uint32_t version)
{
    if (version == 0 || version > 3)
        throw std::runtime_error("Invalid CmdList version number.");
    if (version >= 2)
        min_time = PulseTime::Min2;
    if (version >= 3) {
        _run<Inst_v3>(std::forward<T>(cb), code, code_len);
    }
    else {
        _run<Inst_v1>(std::forward<T>(cb), code, code_len);
    }
}

template<typename Inst, typename T>
void ExeState::_run(T &&cb, const uint8_t *code, size_t code_len)
{
    for (size_t i = 0; i < code_len;) {
        auto *p = &code[i];
        uint8_t op = *p;
        auto inst_len = Inst::cmd_size[op];
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
                i += sizeof(typename Inst::Wait);
                auto inst = Mem::load_unalign<typename Inst::Wait>(p2);
                t += inst.t;
            }
            return t;
        };
        auto runTTL = [&] (uint32_t ttl, uint64_t t, int bank) {
            cb.ttl(ttl, consumeAllWait() + t + min_time, bank);
        };
        auto runTTL1 = [&] (uint8_t chn, bool val, uint64_t t) {
            cb.ttl1(chn, val, consumeAllWait() + t + min_time);
        };
        switch (op) {
        case OpCode::TTLAll: {
            auto inst = Mem::load_unalign<typename Inst::TTLAll>(p);
            if constexpr (Inst::version >= 3) {
                runTTL(inst.val, inst.t, inst.bank);
            }
            else {
                runTTL(inst.val, 0, 0);
            }
            break;
        }
        case OpCode::TTL1: {
            auto inst = Mem::load_unalign<typename Inst::TTL1>(p);
            runTTL1(inst.chn, inst.val, inst.t);
            break;
        }
        case OpCode::Wait: {
            auto inst = Mem::load_unalign<typename Inst::Wait>(p);
            uint64_t t = inst.t;
            cb.wait(t + consumeAllWait());
            break;
        }
        case OpCode::Clock: {
            cb.clock(Mem::load_unalign<typename Inst::Clock>(p).period);
            break;
        }
        case OpCode::DDSFreq: {
            auto inst = Mem::load_unalign<typename Inst::DDSFreq>(p);
            cb.dds_freq(inst.chn, inst.freq);
            break;
        }
        case OpCode::DDSAmp: {
            auto inst = Mem::load_unalign<typename Inst::DDSAmp>(p);
            cb.dds_amp(inst.chn, inst.amp);
            break;
        }
        case OpCode::DDSPhase: {
            auto inst = Mem::load_unalign<typename Inst::DDSPhase>(p);
            cb.dds_phase(inst.chn, inst.phase);
            break;
        }
        case OpCode::DDSDetPhase: {
            auto inst = Mem::load_unalign<typename Inst::DDSDetPhase>(p);
            cb.dds_detphase(inst.chn, inst.det_phase);
            break;
        }
        case OpCode::DDSReset: {
            auto inst = Mem::load_unalign<typename Inst::DDSReset>(p);
            cb.dds_reset(inst.chn);
            break;
        }
        case OpCode::DAC: {
            auto inst = Mem::load_unalign<typename Inst::DAC>(p);
            cb.dac(inst.chn, inst.amp);
            break;
        }
        default:
            throw std::runtime_error("Invalid opcode.");
        }
    }
}

SeqMetadata parse(buff_ostream &output, std::istream &input, uint32_t version);

}
}

#endif
