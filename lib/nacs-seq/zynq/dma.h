/*************************************************************************
 *   Copyright (c) 2026 - 2026 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_SEQ_ZYNQ_DMA_H__
#define __NACS_SEQ_ZYNQ_DMA_H__

#include "utils.h"

#include <nacs-utils/mem.h>

#include <optional>
#include <span>
#include <stdexcept>

namespace NaCs::Seq::Zynq::DMA {

enum OpCode : uint8_t {
    WAIT1 = 0,
    TTL_SET4 = 1,
    CLOCKOUT = 3,

    WAIT2 = 0,
    TTL_SET16 = 1,
    DDS_SET16 = 2,

    WAIT_TRIG = 0,
    TTL_SET32 = 1,
    DDS_SET32 = 2,
    DAC = 3,
};

// Instruction format:
//   [len: 2][opcode: 2][data: 12/28/44]
// Instruction length is fully determined by the first two bits to simplify decoding
//
// All instructions
//
// * 2 bytes (len = 0, 4 + 12 bits)
//
//     *     wait1: [len<0>: 2][opcode<0>: 2][cycle: 12]
//     *  ttl_set4: [len<0>: 2][opcode<1>: 2][bank4_1: 6][val1: 4][<0>: 2]
//     *  clockout: [len<0>: 2][opcode<3>: 2][period: 9][<0>: 3]
//
// * 4 bytes (len = 1, 4 + 28 bits)
//
//     *     wait2: [len<1>: 2][opcode<0>: 2][cycle: 28]
//     * ttl_set16: [len<1>: 2][opcode<1>: 2][bank8_1: 5][val1: 8][bank8_2: 5][val2: 8][<0>: 2]
//     * dds_set16: [len<1>: 2][opcode<2>: 2][bus_id: 1][dds_id: 4][fud: 1][addr: 6][data: 16]
//
// * 6 bytes (len = 2, 4 + 44 bits)
//
//     * wait_trig: [len<2>: 2][opcode<0>: 2][chn: 8][edge: 1][cycle: 35]
//     * ttl_set32: [len<2>: 2][opcode<1>: 2][bank16_1: 4][val1: 16][bank16_2: 4][val2: 16][<0>: 4]
//     * dds_set32: [len<2>: 2][opcode<2>: 2][bus_id: 1][dds_id: 4][fud: 1][addr: 6][data: 32]
//     *       dac: [len<2>: 2][opcode<3>: 2][id: 2][cycle: 9][clk_pha: 1][clk_pol: 1][data: 18][<0>: 13]

namespace Inst_v0 {

struct NACS_PACKED _Header {
    uint8_t len: 2;
    OpCode op: 2;
};

/**
 * 16 bit instructions
 */

struct NACS_PACKED Wait1 {
    uint8_t len: 2;
    OpCode op: 2;
    uint16_t cycle: 12;
    Wait1() = default;
    Wait1(uint16_t cycle)
        : len(0),
          op(OpCode::WAIT1),
          cycle(cycle)
    {
    }
};
static_assert(sizeof(Wait1) == 2, "");

struct NACS_PACKED TTLSet4 {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t bank4_1: 6;
    uint8_t val1: 4;
    TTLSet4() = default;
    TTLSet4(uint8_t bank4_1, uint8_t val1)
        : len(0),
          op(OpCode::TTL_SET4),
          bank4_1(bank4_1),
          val1(val1)
    {
    }
};
static_assert(sizeof(TTLSet4) == 2, "");

struct NACS_PACKED ClockOut {
    uint8_t len: 2;
    OpCode op: 2;
    uint16_t period: 9;
    ClockOut() = default;
    ClockOut(uint16_t period)
        : len(0),
          op(OpCode::CLOCKOUT),
          period(period)
    {
    }
};
static_assert(sizeof(ClockOut) == 2, "");

/**
 * 32 bit instructions
 */

struct NACS_PACKED Wait2 {
    uint8_t len: 2;
    OpCode op: 2;
    uint32_t cycle: 28;
    Wait2() = default;
    Wait2(uint32_t cycle) : len(1), op(OpCode::WAIT2), cycle(cycle)
    {
    }
};
static_assert(sizeof(Wait2) == 4, "");

struct NACS_PACKED TTLSet16 {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t bank8_1: 5;
    uint8_t val1: 8;
    uint8_t bank8_2: 5;
    uint8_t val2: 8;
    TTLSet16() = default;
    TTLSet16(uint8_t bank8_1, uint8_t val1,
             uint8_t bank8_2, uint8_t val2)
    : len(1), op(OpCode::TTL_SET16),
    bank8_1(bank8_1), val1(val1), bank8_2(bank8_2), val2(val2)
    {
    }
};
static_assert(sizeof(TTLSet16) == 4, "");

struct NACS_PACKED DDSSet16 {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t bus_id: 1;
    uint8_t dds_id: 4;
    uint8_t fud: 1;
    uint8_t addr: 6;
    uint16_t data: 16;
    DDSSet16() = default;
    DDSSet16(uint8_t bus_id, uint8_t dds_id, uint8_t fud, uint8_t addr, uint16_t data)
    : len(1), op(OpCode::DDS_SET16),
    bus_id(bus_id), dds_id(dds_id), fud(fud), addr(addr), data(data)
    {
    }
};
static_assert(sizeof(DDSSet16) == 4, "");

/**
 * 48 bit instructions
 */

struct NACS_PACKED WaitTrig {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t chn: 8;
    uint8_t edge: 1;
    uint64_t cycle: 35;
    WaitTrig() = default;
    WaitTrig(uint8_t chn, uint8_t edge, uint64_t cycle)
    : len(2), op(OpCode::WAIT_TRIG), chn(chn), edge(edge), cycle(cycle)
    {
    }
};
static_assert(sizeof(WaitTrig) == 6, "");

struct NACS_PACKED TTLSet32 {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t bank16_1: 4;
    uint16_t val1: 16;
    uint8_t bank16_2: 4;
    uint16_t val2: 16;
    TTLSet32() = default;
    TTLSet32(uint8_t bank16_1, uint16_t val1,
             uint8_t bank16_2, uint16_t val2)
    : len(2), op(OpCode::TTL_SET32),
    bank16_1(bank16_1), val1(val1), bank16_2(bank16_2), val2(val2)
    {
    }
};
static_assert(sizeof(TTLSet32) == 6, "");

struct NACS_PACKED DDSSet32 {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t bus_id: 1;
    uint8_t dds_id: 4;
    uint8_t fud: 1;
    uint8_t addr: 6;
    uint32_t data: 32;
    DDSSet32() = default;
    DDSSet32(uint8_t bus_id, uint8_t dds_id, uint8_t fud, uint8_t addr, uint32_t data)
    : len(2), op(OpCode::DDS_SET32),
    bus_id(bus_id), dds_id(dds_id), fud(fud), addr(addr), data(data)
    {
    }
};
static_assert(sizeof(DDSSet32) == 6, "");

struct NACS_PACKED DAC {
    uint8_t len: 2;
    OpCode op: 2;
    uint8_t id: 2;
    uint16_t cycle: 9;
    uint8_t clk_pha: 1;
    uint8_t clk_pol: 1;
    uint32_t data: 18;
    DAC() = default;
    DAC(uint8_t id, uint16_t cycle, uint8_t clk_pha, uint8_t clk_pol, uint32_t data)
    : len(2), op(OpCode::DAC),
    id(id), cycle(cycle), clk_pha(clk_pha), clk_pol(clk_pol), data(data)
    {
    }
};
static_assert(sizeof(DDSSet32) == 6, "");

}

// No state to keep track of currently
struct ExeState {
    template<typename CB>
    void run(CB &&cb, std::span<const uint8_t> code, uint32_t version);
};

namespace {

template<typename T> struct _tag_t {};
template<typename T> static constexpr _tag_t<T> _tag_v;

}

template<typename CB>
inline void ExeState::run(CB &&cb, std::span<const uint8_t> code, uint32_t version)
{
    if (version != 0)
        throw std::runtime_error("Invalid DMA instructions version number.");
    auto code_len = code.size();
    if (code_len % 2 != 0)
        throw std::runtime_error("Invalid DMA instructions length.");
    for (size_t i = 0; i < code_len;) {
        auto load = [&] <typename T> (_tag_t<T>) -> std::optional<T> {
            if (i + sizeof(T) > code_len)
                return std::optional<T>();
            return Mem::load_unalign<T>(&code[i]);
        };
        auto header = load(_tag_v<Inst_v0::_Header>);
        auto invalid = [&] {
            auto len = (header->len + 1) * 2;
            if (i + len > code_len)
                len = code_len - i;
            cb.invalid(code.subspan(i, len));
            i += len;
        };
        auto handle = [&] <typename T> (_tag_t<T>) {
            auto inst = load(_tag_v<T>);
            if (inst) {
                cb.inst(*inst);
                i += sizeof(T);
            }
            else {
                invalid();
            }
        };
        switch (header->len) {
        case 0:
            switch (header->op) {
            case OpCode::WAIT1:
                handle(_tag_v<Inst_v0::Wait1>);
                break;
            case OpCode::TTL_SET4:
                handle(_tag_v<Inst_v0::TTLSet4>);
                break;
            case OpCode::CLOCKOUT:
                handle(_tag_v<Inst_v0::ClockOut>);
                break;
            default:
                invalid();
            }
            break;
        case 1:
            switch (header->op) {
            case OpCode::WAIT2:
                handle(_tag_v<Inst_v0::Wait2>);
                break;
            case OpCode::TTL_SET16:
                handle(_tag_v<Inst_v0::TTLSet16>);
                break;
            case OpCode::DDS_SET16:
                handle(_tag_v<Inst_v0::DDSSet16>);
                break;
            default:
                invalid();
            }
            break;
        case 2:
            switch (header->op) {
            case OpCode::WAIT_TRIG:
                handle(_tag_v<Inst_v0::WaitTrig>);
                break;
            case OpCode::TTL_SET32:
                handle(_tag_v<Inst_v0::TTLSet32>);
                break;
            case OpCode::DDS_SET32:
                handle(_tag_v<Inst_v0::DDSSet32>);
                break;
            case OpCode::DAC:
                handle(_tag_v<Inst_v0::DAC>);
                break;
            }
            break;
        default:
            invalid();
        }
    }
}

}

#endif
