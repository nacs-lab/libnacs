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

#include "dma.h"

#include <concepts>

namespace NaCs::Seq::Zynq::DMA {

namespace {

struct Counter {
    void inst(auto)
    {
        count += 1;
    }
    void invalid(auto)
    {
    }

    size_t count{0};
};

struct Printer {
    void write_hex(std::integral auto v)
    {
        char c = (v < 10) ? ('0' + v) : ('a' + v - 10);
        stm << c;
    }
    void write_nhex(std::integral auto v, int width)
    {
        for (int shift = width - 1; shift >= 0; shift--) {
            write_hex((v >> (shift * 4)) & 0xf);
        }
    }

    void inst(Inst_v0::Wait1 inst)
    {
        stm << "wait1(" << inst.cycle << ")" << std::endl;
    }
    void inst(Inst_v0::TTLSet4 inst)
    {
        stm << "ttl_set4([" << inst.bank4_1 * 4 << ":"
            << (inst.bank4_1 + 1) * 4 - 1 << "]=0x";
        write_hex(inst.val1);
        stm << ")" << std::endl;
    }
    void inst(Inst_v0::ClockOut inst)
    {
        if (inst.period == 0x1ff) {
            stm << "clockout(off)" << std::endl;
        }
        else {
            stm << "clockout(" << inst.period << ")" << std::endl;
        }
    }
    void inst(Inst_v0::Wait2 inst)
    {
        stm << "wait2(" << inst.cycle << ")" << std::endl;
    }
    void inst(Inst_v0::TTLSet16 inst)
    {
        stm << "ttl_set16([" << inst.bank8_1 * 8 << ":"
            << (inst.bank8_1 + 1) * 8 - 1 << "]=0x";
        write_nhex(inst.val1, 2);
        stm << ", [" << inst.bank8_2 * 8 << ":"
            << (inst.bank8_2 + 1) * 8 - 1 << "]=0x";
        write_nhex(inst.val2, 2);
        stm << ")" << std::endl;
    }
    void inst(Inst_v0::DDSSet16 inst)
    {
        stm << "dds_set16(bus=" << int(inst.bus_id) << ", dds=" << int(inst.dds_id)
            << ", fud=" << int(inst.fud) << ", [0x";
        write_nhex(inst.addr << 1, 2);
        stm << ":0x";
        write_nhex((inst.addr << 1) + 1, 2);
        stm << "]=0x";
        write_nhex(inst.data, 4);
        stm << ")" << std::endl;
    }
    void inst(Inst_v0::WaitTrig inst)
    {
        stm << "wait_trig(chn=" << int(inst.chn)
            << ", edge=" << int(inst.edge)
            << ", timeout=" << inst.cycle << ")" << std::endl;
    }
    void inst(Inst_v0::TTLSet32 inst)
    {
        stm << "ttl_set32([" << inst.bank16_1 * 16 << ":"
            << (inst.bank16_1 + 1) * 16 << "]=0x";
        write_nhex(inst.val1, 4);
        stm << ", [" << inst.bank16_2 * 16 << ":"
            << (inst.bank16_2 + 1) * 16 << "]=0x";
        write_nhex(inst.val2, 4);
        stm << ")" << std::endl;
    }
    void inst(Inst_v0::DDSSet32 inst)
    {
        stm << "dds_set32(bus=" << int(inst.bus_id) << ", dds=" << int(inst.dds_id)
            << ", fud=" << int(inst.fud) << ", [0x";
        write_nhex(inst.addr << 1, 2);
        stm << ":0x";
        write_nhex((inst.addr << 1) + 3, 2);
        stm << "]=0x";
        write_nhex(inst.data, 8);
        stm << ")" << std::endl;
    }
    void inst(Inst_v0::DAC inst)
    {
        stm << "dac(id=" << int(inst.id) << ", cycle=" << int(inst.cycle)
            << ", pha=" << int(inst.clk_pha)
            << ", pol=" << int(inst.clk_pol) << ", data=0x";
        write_nhex(inst.data << 1, 5);
        stm << ")" << std::endl;
    }
    void invalid(std::span<const uint8_t> inst)
    {
        stm << "invalid(";
        bool isfirst = true;
        for (auto c: inst) {
            if (isfirst) {
                isfirst = false;
                stm << "0x";
            }
            else {
                stm << ", 0x";
            }
            write_nhex(c, 2);
        }
        stm << ")" << std::endl;
    }

    std::ostream &stm;
};

}

NACS_EXPORT() size_t count(std::span<const uint8_t> code, uint32_t version)
{
    ExeState state;
    Counter counter;
    state.run(counter, code, version);
    return counter.count;
}

NACS_EXPORT() void print(std::ostream &stm, std::span<const uint8_t> code,
                         uint32_t version)
{
    ExeState state;
    Printer printer(stm);
    state.run(printer, code, version);
}

}
