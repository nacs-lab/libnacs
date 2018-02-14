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

#include "pulser.h"

namespace NaCs {
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
    uint8_t t: 2;
    uint8_t val: 1;
    uint16_t chn: 5;
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

}

}
}

#endif
