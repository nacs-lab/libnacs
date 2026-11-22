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
#include <type_traits>

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

// Collect the values of a channel.
// The time here is the absolute time (in DMA cycles) since the start of the instructions,
// without the time of the trigger waits.
struct ValueCollector {
    ValueCollector(ChnType type, uint8_t num, int start_ttl, const DDSTiming &timing)
        : type(type),
          num(num),
          start_ttl(start_ttl),
          timing(timing),
          started(start_ttl < 0)
    {
    }

    void record(uint64_t t, double val)
    {
        ts.push_back(t);
        vals.push_back(val);
    }

    void inst(Inst_v0::Wait1 inst)
    {
        t += inst.cycle + 1;
    }
    void inst(Inst_v0::Wait2 inst)
    {
        t += inst.cycle + 1;
    }
    void inst(Inst_v0::WaitTrig)
    {
        // The time of the trigger is unknown.
    }
    void inst(Inst_v0::ClockOut)
    {
    }

    // Set `width` TTL channels starting from `start`.
    void ttl(unsigned start, unsigned width, uint32_t val)
    {
        if (start_ttl >= int(start) && start_ttl < int(start + width)) {
            bool new_val = (val >> (start_ttl - start)) & 1;
            // The sequence starts at the end of the start trigger.
            if (!started && start_val && !new_val) {
                started = true;
                start_t = t;
            }
            start_val = new_val;
        }
        if (type == ChnType::TTL && num >= start && num < start + width) {
            record(t, double((val >> (num - start)) & 1));
        }
    }
    void inst(Inst_v0::TTLSet4 inst)
    {
        ttl(inst.bank4_1 * 4, 4, inst.val1);
    }
    void inst(Inst_v0::TTLSet16 inst)
    {
        ttl(inst.bank8_1 * 8, 8, inst.val1);
        if (inst.bank8_2 != inst.bank8_1) {
            ttl(inst.bank8_2 * 8, 8, inst.val2);
        }
    }
    void inst(Inst_v0::TTLSet32 inst)
    {
        ttl(inst.bank16_1 * 16, 16, inst.val1);
        if (inst.bank16_2 != inst.bank16_1) {
            ttl(inst.bank16_2 * 16, 16, inst.val2);
        }
    }

    // The registers for the frequency (low and high 16 bits), phase and amplitude
    // (byte address 0x2c - 0x33).
    static constexpr uint8_t reg_base = 0x2c >> 1;
    void set_reg(unsigned chn, uint8_t addr, uint16_t val)
    {
        if (addr < reg_base || addr >= reg_base + 4)
            return;
        auto reg = addr - reg_base;
        regs[chn][reg] = val;
        // Bit 0: frequency, bit 1: phase, bit 2: amplitude
        written[chn] = uint8_t(written[chn] | (1 << (reg == 0 ? 0 : reg - 1)));
    }
    void dds(uint8_t bus, uint8_t id, uint8_t addr, bool fud, bool is32, uint32_t data)
    {
        if (id >= 11)
            return;
        auto chn = bus * 11 + id;
        set_reg(chn, addr, uint16_t(data));
        if (is32)
            set_reg(chn, uint8_t(addr | 1), uint16_t(data >> 16));
        if (!fud)
            return;
        // All the DDS's on the bus written since the last FUD are updated at the FUD.
        auto fud_t = t + (is32 ? timing.set32_fud : timing.set16_fud) - timing.fud_to_end;
        for (unsigned c = bus * 11; c < unsigned(bus * 11 + 11); c++) {
            auto w = written[c];
            written[c] = 0;
            if (c != num)
                continue;
            auto &r = regs[c];
            if (type == ChnType::Freq && (w & 1)) {
                record(fud_t, dds_freq_from_mu(uint32_t(r[0]) | (uint32_t(r[1]) << 16)));
            }
            else if (type == ChnType::Phase && (w & 2)) {
                record(fud_t, dds_phase_from_mu(r[2]) * 360); // in degrees
            }
            else if (type == ChnType::Amp && (w & 4)) {
                record(fud_t, dds_amp_from_mu(r[3] & 0xfff));
            }
        }
    }
    void inst(Inst_v0::DDSSet16 inst)
    {
        dds(inst.bus_id, inst.dds_id, inst.addr, inst.fud, false, inst.data);
    }
    void inst(Inst_v0::DDSSet32 inst)
    {
        dds(inst.bus_id, inst.dds_id, inst.addr, inst.fud, true, inst.data);
    }
    void inst(Inst_v0::DAC inst)
    {
        uint32_t data = inst.data;
        if (type == ChnType::DAC && (data >> 16) == num) {
            record(t, dac_from_mu(data & 0xffff));
        }
    }
    void invalid(std::span<const uint8_t>)
    {
        throw std::runtime_error("Invalid DMA instruction.");
    }

    // Convert the times to be relative to the start of the sequence.
    void finalize()
    {
        for (auto &v: ts) {
            if (!started) {
                v = 0;
            }
            else {
                v = v > start_t ? v - start_t : 0;
            }
        }
    }

    const ChnType type;
    const uint8_t num;
    const int start_ttl;
    const DDSTiming timing;
    uint64_t t = 0;
    bool started;
    bool start_val = false;
    uint64_t start_t = 0;
    uint16_t regs[22][4] = {};
    uint8_t written[22] = {};
    std::vector<uint64_t> ts;
    std::vector<double> vals;
};

struct TimeKeeper {
    template<typename T>
    void inst(T inst)
    {
        using T2 = std::remove_cv_t<T>;
        if constexpr (std::same_as<T2,Inst_v0::Wait1> ||
                      std::same_as<T2,Inst_v0::Wait2>) {
            time += inst.cycle + 1;
        }
    }
    void invalid(auto)
    {
    }

    uint64_t time{0};
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

NACS_EXPORT() uint64_t total_time(std::span<const uint8_t> code, uint32_t version)
{
    ExeState state;
    TimeKeeper keeper;
    state.run(keeper, code, version);
    return keeper.time;
}

NACS_EXPORT() std::pair<std::vector<uint64_t>,std::vector<double>>
extract_values(std::span<const uint8_t> code, uint32_t version, ChnType type,
               uint8_t num, int start_ttl, const DDSTiming &timing)
{
    ExeState state;
    ValueCollector collector(type, num, start_ttl, timing);
    state.run(collector, code, version);
    collector.finalize();
    return {std::move(collector.ts), std::move(collector.vals)};
}

}
