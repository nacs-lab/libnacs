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

#include "exehelper_p.h"

#include "bytecode.h"

#include <assert.h>
#include <math.h>

namespace NaCs::Seq::Zynq::ByteCode {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len, uint32_t version)
{
    if (version == 0 || version > 3)
        throw std::runtime_error("Invalid ByteCode version number.");
    // Do not use the `ExeState` helper to avoid fusing of instructions.
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b & 0xf;
        assert(op < 14);
        auto len = version == 3 ? Inst_v3::inst_size[op] : Inst_v1::inst_size[op];
        i += len;
        count++;
    }
    return count;
}

NACS_EXPORT() void print_raw(std::ostream &stm, const uint8_t *code, size_t code_len,
                             uint32_t version)
{
    int min_time = version >= 2 ? PulseTime::Min2 : PulseTime::Min;
    for (size_t i = 0; i < code_len;) {
        auto *p = &code[i];
        uint8_t b = *p;
        uint8_t op = b & 0xf;
        auto inst_len = version >= 3 ? Inst_v3::inst_size[op] : Inst_v1::inst_size[op];
        i += inst_len;
        switch (op) {
        case OpCode::TTLAll: {
            if (version >= 3) {
                auto inst = Mem::load_unalign<Inst_v3::TTLAll>(p);
                stm << "ttl[" << int(inst.bank) << "]=0x"
                    << std::hex << inst.val << " t=0x" << int(inst.t)
                    << "(+0x" << min_time << ")" << std::dec << std::endl;
            }
            else {
                auto inst = Mem::load_unalign<Inst_v1::TTLAll>(p);
                stm << "ttl=0x" << std::hex << inst.val << " t=0x" << int(inst.t)
                    << "(+0x" << min_time << ")" << std::dec << std::endl;
            }
            break;
        }
        case OpCode::_TTL_flip_2b: {
            if (version >= 3) {
                auto inst = Mem::load_unalign<Inst_v3::TTL1>(p);
                stm << "ttl1(" << int(inst.val1) << ", bank=" << int(inst.bank)
                    << ") t=0x" << int(inst.t) << "(+0x" << min_time << ")"
                    << std::dec << std::endl;
            }
            else {
                auto inst = Mem::load_unalign<Inst_v1::TTL2>(p);
                stm << "ttl2(" << int(inst.val1) << ", " << int(inst.val2)
                    << ") t=0x" << int(inst.t) << "(+0x" << min_time << ")"
                    << std::dec << std::endl;
            }
            break;
        }
        case OpCode::_TTL_flip_3b: {
            if (version >= 3) {
                auto inst = Mem::load_unalign<Inst_v3::TTL3>(p);
                stm << "ttl3(" << int(inst.val1) << ", " << int(inst.val2) << ", "
                    << int(inst.val3) << ", bank=" << int(inst.bank)
                    << ") t=0x" << int(inst.t) << "(+0x" << min_time << ")"
                    << std::dec << std::endl;
            }
            else {
                auto inst = Mem::load_unalign<Inst_v1::TTL4>(p);
                stm << "ttl4(" << int(inst.val1) << ", " << int(inst.val2) << ", "
                    << int(inst.val3) << ", " << int(inst.val4) << ") t=0x0(+0x"
                    << min_time << ")" << std::dec << std::endl;
            }
            break;
        }
        case OpCode::TTL5: {
            if (version >= 3) {
                auto inst = Mem::load_unalign<Inst_v3::TTL5>(p);
                stm << "ttl5(" << int(inst.val1) << ", " << int(inst.val2) << ", "
                    << int(inst.val3) << ", " << int(inst.val4) << ", " << int(inst.val5)
                    << ", bank=" << int(inst.bank) << ") t=0x0(+0x" << min_time << ")"
                    << std::dec << std::endl;
            }
            else {
                auto inst = Mem::load_unalign<Inst_v1::TTL5>(p);
                stm << "ttl5(" << int(inst.val1) << ", " << int(inst.val2) << ", "
                    << int(inst.val3) << ", " << int(inst.val4) << ", " << int(inst.val5)
                    << ") t=0x" << int(inst.t) << "(+0x" << min_time << ")"
                    << std::dec << std::endl;
            }
            break;
        }
        case OpCode::Wait: {
            auto inst = Mem::load_unalign<InstDefs::Wait>(p);
            stm << "wait(0x" << std::hex << int(inst.t) << std::dec
                << " * 2^(" << int(inst.exp) << " * 3))" << std::endl;
            break;
        }
        case OpCode::Clock: {
            if (b & 0x10) {
                auto inst = Mem::load_unalign<InstDefs::Wait2>(p);
                stm << "wait2(0x" << std::hex << int(inst.t) << std::dec << ")" << std::endl;
            }
            else {
                stm << "clock(" << int(Mem::load_unalign<InstDefs::Clock>(p).period)
                    << ")" << std::endl;
            }
            break;
        }
        case OpCode::DDSFreq: {
            auto inst = Mem::load_unalign<InstDefs::DDSFreq>(p);
            stm << "freq(" << int(inst.chn) << ")=0x"
                << std::hex << inst.freq << std::dec << std::endl;
            break;
        }
        case OpCode::DDSDetFreq2: {
            auto inst = Mem::load_unalign<InstDefs::DDSDetFreq2>(p);
            uint32_t freq = inst.freq;
            char sign = '+';
            if (freq & 0x40) {
                sign = '-';
                freq = uint32_t(-(freq | 0xffffff80));
            }
            stm << "freq2(" << int(inst.chn) << ")" << sign << "=0x"
                << std::hex << freq << std::dec << std::endl;
            break;
        }
        case OpCode::DDSDetFreq3: {
            auto inst = Mem::load_unalign<InstDefs::DDSDetFreq3>(p);
            uint32_t freq = inst.freq;
            char sign = '+';
            if (freq & 0x4000) {
                sign = '-';
                freq = uint32_t(-(freq | 0xffff8000));
            }
            stm << "freq3(" << int(inst.chn) << ")" << sign << "=0x"
                << std::hex << freq << std::dec << std::endl;
            break;
        }
        case OpCode::DDSDetFreq4: {
            auto inst = Mem::load_unalign<InstDefs::DDSDetFreq4>(p);
            uint32_t freq = inst.freq;
            char sign = '+';
            if (freq & 0x400000) {
                sign = '-';
                freq = uint32_t(-(freq | 0xff800000));
            }
            stm << "freq4(" << int(inst.chn) << ")" << sign << "=0x"
                << std::hex << freq << std::dec << std::endl;
            break;
        }
        case OpCode::DDSAmp: {
            auto inst = Mem::load_unalign<InstDefs::DDSAmp>(p);
            stm << "amp(" << int(inst.chn) << ")=0x"
                << std::hex << inst.amp << std::dec << std::endl;
            break;
        }
        case OpCode::DDSDetAmp: {
            auto inst = Mem::load_unalign<InstDefs::DDSDetAmp>(p);
            uint16_t amp = inst.amp;
            char sign = '+';
            if (amp & 0x40) {
                sign = '-';
                amp = uint16_t(-(amp | 0xff80));
            }
            stm << "amp(" << int(inst.chn) << ")" << sign << "=0x"
                << std::hex << amp << std::dec << std::endl;
            break;
        }
        case OpCode::DDSPhase: {
            auto inst = Mem::load_unalign<InstDefs::DDSPhase>(p);
            stm << "phase(" << int(inst.chn) << ")=0x"
                << std::hex << inst.phase << std::dec << std::endl;
            break;
        }
        case OpCode::DAC: {
            auto inst = Mem::load_unalign<InstDefs::DAC>(p);
            stm << "dac(" << int(inst.chn) << ")=0x"
                << std::hex << inst.amp << std::dec << std::endl;
            break;
        }
        case OpCode::DACDet: {
            auto inst = Mem::load_unalign<InstDefs::DACDet>(p);
            uint16_t amp = inst.amp;
            char sign = '+';
            if (amp & 0x200) {
                sign = '-';
                amp = uint16_t(-(amp | 0xfc00));
            }
            stm << "dac(" << int(inst.chn) << ")" << sign << "=0x"
                << std::hex << amp << std::dec << std::endl;
            break;
        }
        default:
            throw std::runtime_error("Invalid opcode.");
        }
    }
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask, uint32_t version)
{
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
    Printer printer{stm};
    ExeState state;
    state.run(printer, code, code_len, version);
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         const std::vector<uint32_t> &ttl_masks, uint32_t version)
{
    stm << "ttl_mask=" << std::hex;
    bool isfirst = true;
    for (auto ttl_mask: ttl_masks) {
        stm << (isfirst ? "0x" : " 0x") << ttl_mask;
        isfirst = false;
    }
    stm << std::dec << std::endl;
    Printer printer{stm};
    ExeState state;
    state.run(printer, code, code_len, version);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len, uint32_t version)
{
    TimeKeeper keeper;
    ExeState state;
    state.run(keeper, code, code_len, version);
    return keeper.total_t;
}

}
