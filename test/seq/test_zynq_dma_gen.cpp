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

#define CATCH_CONFIG_MAIN

#include "zynq_helper.h"

#include "../../lib/nacs-seq/zynq/dma.h"

#include <algorithm>
#include <random>
#include <set>
#include <type_traits>

// Test the DMA instruction generation (`BCGen::dma`).
// The generated instructions are decoded back into the groups of outputs
// executed in the same cycle, which are checked against the hardware constraints
// (`check_hw_rules`) and the expected timing of each test.

struct DummyHostSeq : HostSeq {
    DummyHostSeq()
        : HostSeq()
    {
        nconsts = 0;
        nglobals = 0;
        nglobal_vals = 0;
        nchannels = 0;
        nshared = 0;
        npublic_globals = 0;

        seqs.resize(1);
        {
            auto &host_bseq = seqs[0];
            host_bseq.id = 1;

            host_bseq.nmeasure = 0;
            host_bseq.ndirect = 0;
            host_bseq.nneed_order = 0;
            host_bseq.default_branch = -1;
            host_bseq.ndirect_assumes = 0;
        }
    }

    template<typename T>
    void set_value(uint32_t idx, T v)
    {
        reserve(idx + 1);
        if constexpr (std::is_same_v<T,bool>) {
            types[idx] = HostSeq::Type::Bool;
            values[idx].b = v;
        }
        else if constexpr (std::is_same_v<T,int64_t>) {
            types[idx] = HostSeq::Type::Int64;
            values[idx].i64 = v;
        }
        else {
            static_assert(std::is_same_v<T,double>);
            types[idx] = HostSeq::Type::Float64;
            values[idx].f64 = v;
        }
    }
private:
    void reserve(uint32_t len)
    {
        if (len > nconsts) {
            nconsts = len;
            nshared = len;
            values.resize(len);
            types.resize(len);
        }
    }
};

namespace {

// Sequence time unit is 1ps and the DMA cycle is 5ns.
constexpr int64_t ps_per_cycle = 5000;
constexpr int64_t cycles_per_us = 200;
constexpr int64_t us = 1000000; // in sequence time unit (ps)

struct DDSWrite {
    int64_t t;
    uint8_t bus;
    uint8_t id;
    uint8_t addr;
    bool is32;
    bool fud;
    uint32_t data;
};

// DDS parameters applied by a FUD
struct DDSUpdate {
    int64_t t; // (rising) edge of the FUD
    uint8_t chn;
    uint32_t freq;
    uint16_t amp;
    uint16_t phase;
};

struct TTLChange {
    int64_t t;
    int bank;
    uint32_t val; // The value of the bank after the change
};

struct Decoded {
    std::vector<uint32_t> masks;
    std::vector<TTLChange> ttl;
    // Time and width (4, 16 or 32) of the TTL instructions
    std::vector<std::pair<int64_t,int>> ttl_insts;
    std::vector<DDSWrite> dds;
    std::vector<DDSUpdate> dds_updates;
    std::vector<std::pair<int64_t,uint32_t>> dac;
    std::vector<std::pair<int64_t,uint16_t>> clock;
    std::vector<std::pair<int64_t,Zynq::DMA::Inst_v0::WaitTrig>> trig;
    std::vector<std::pair<int64_t,uint64_t>> waits; // time of the wait and its length
    int64_t total = 0;
    bool ends_with_wait = false;

    // First time the bank has the value after `t0`
    int64_t ttl_time(int bank, uint32_t mask, uint32_t val, int64_t t0=0) const
    {
        for (auto &change: ttl) {
            if (change.t >= t0 && change.bank == bank && (change.val & mask) == val) {
                return change.t;
            }
        }
        return -1;
    }
};

struct DMADecoder {
    DMADecoder(Decoded &res, const Zynq::BCGen::DDSTiming &timing)
        : res(res),
          timing(timing)
    {
    }
    // Current group
    uint32_t group_ttl_mask[8] = {};
    bool group_dds[2] = {false, false};
    bool group_dac = false;
    bool group_clock = false;
    uint32_t ttl[8] = {};
    uint16_t dds_regs[22][64] = {};
    // The DDS's written since the last FUD on each bus.
    std::vector<uint8_t> dds_written[2];
    int64_t t = 0;
    Decoded &res;
    const Zynq::BCGen::DDSTiming &timing;

    void set_ttl(unsigned bit, unsigned width, uint32_t val)
    {
        auto bank = bit / 32;
        REQUIRE(bank < res.masks.size());
        auto shift = bit % 32;
        auto mask = uint32_t(((uint64_t(1) << width) - 1) << shift);
        val = (val << shift) & mask;
        // Each sub-bank is written only once in a group
        REQUIRE((group_ttl_mask[bank] & mask) == 0);
        // The value must not have any bit set outside of the DMA TTL mask.
        REQUIRE((val & ~res.masks[bank]) == 0);
        group_ttl_mask[bank] |= mask;
        ttl[bank] = (ttl[bank] & ~mask) | val;
    }
    void end_group(bool is_wait=true)
    {
        for (int bank = 0; bank < 8; bank++) {
            if (group_ttl_mask[bank]) {
                res.ttl.push_back({t, bank, ttl[bank]});
            }
            group_ttl_mask[bank] = 0;
        }
        group_dds[0] = group_dds[1] = false;
        group_dac = false;
        group_clock = false;
        res.ends_with_wait = is_wait;
    }
    void inst(Zynq::DMA::Inst_v0::Wait1 inst)
    {
        end_group();
        res.waits.push_back({t, inst.cycle + 1});
        t += inst.cycle + 1;
    }
    void inst(Zynq::DMA::Inst_v0::Wait2 inst)
    {
        end_group();
        res.waits.push_back({t, inst.cycle + 1});
        t += inst.cycle + 1;
    }
    void inst(Zynq::DMA::Inst_v0::WaitTrig inst)
    {
        end_group();
        // The time the trigger takes is unknown, treat it as no time.
        res.trig.push_back({t, inst});
    }
    void inst(Zynq::DMA::Inst_v0::TTLSet4 inst)
    {
        res.ttl_insts.push_back({t, 4});
        set_ttl(inst.bank4_1 * 4, 4, inst.val1);
        res.ends_with_wait = false;
    }
    void inst(Zynq::DMA::Inst_v0::TTLSet16 inst)
    {
        res.ttl_insts.push_back({t, 16});
        set_ttl(inst.bank8_1 * 8, 8, inst.val1);
        if (inst.bank8_2 != inst.bank8_1) {
            set_ttl(inst.bank8_2 * 8, 8, inst.val2);
        }
        else {
            REQUIRE(inst.val1 == inst.val2);
        }
        res.ends_with_wait = false;
    }
    void inst(Zynq::DMA::Inst_v0::TTLSet32 inst)
    {
        res.ttl_insts.push_back({t, 32});
        set_ttl(inst.bank16_1 * 16, 16, inst.val1);
        if (inst.bank16_2 != inst.bank16_1) {
            set_ttl(inst.bank16_2 * 16, 16, inst.val2);
        }
        else {
            REQUIRE(inst.val1 == inst.val2);
        }
        res.ends_with_wait = false;
    }
    void dds_write(uint8_t bus, uint8_t id, uint8_t addr, bool fud, bool is32,
                   uint32_t data)
    {
        REQUIRE(bus < 2);
        REQUIRE(id < 11);
        // One output per bus per group
        REQUIRE(!group_dds[bus]);
        group_dds[bus] = true;
        res.dds.push_back({t, bus, id, addr, is32, fud, data});
        auto chn = bus * 11 + id;
        dds_regs[chn][addr] = uint16_t(data);
        if (is32) {
            // The address of the 32 bit write must be even
            REQUIRE(addr % 2 == 0);
            dds_regs[chn][addr | 1] = uint16_t(data >> 16);
        }
        auto &written = dds_written[bus];
        if (std::find(written.begin(), written.end(), uint8_t(chn)) == written.end())
            written.push_back(uint8_t(chn));
        if (fud) {
            // The FUD is shared by all the DDS's on the bus.
            auto fud_t = t + (is32 ? timing.set32_fud : timing.set16_fud) - timing.fud_to_end;
            for (auto c: written) {
                auto &regs = dds_regs[c];
                res.dds_updates.push_back({fud_t, c,
                        uint32_t(regs[0x2c >> 1]) | (uint32_t(regs[0x2e >> 1]) << 16),
                        regs[0x32 >> 1], regs[0x30 >> 1]});
            }
            written.clear();
        }
        res.ends_with_wait = false;
    }
    void inst(Zynq::DMA::Inst_v0::DDSSet16 inst)
    {
        dds_write(inst.bus_id, inst.dds_id, inst.addr, inst.fud, false, inst.data);
    }
    void inst(Zynq::DMA::Inst_v0::DDSSet32 inst)
    {
        dds_write(inst.bus_id, inst.dds_id, inst.addr, inst.fud, true, inst.data);
    }
    void inst(Zynq::DMA::Inst_v0::DAC inst)
    {
        REQUIRE(!group_dac);
        group_dac = true;
        REQUIRE(inst.id == 0);
        REQUIRE(inst.cycle == 1);
        REQUIRE(inst.clk_pha == 0);
        REQUIRE(inst.clk_pol == 0);
        res.dac.push_back({t, uint32_t(inst.data)});
        res.ends_with_wait = false;
    }
    void inst(Zynq::DMA::Inst_v0::ClockOut inst)
    {
        REQUIRE(!group_clock);
        group_clock = true;
        res.clock.push_back({t, uint16_t(inst.period)});
        res.ends_with_wait = false;
    }
    void invalid(std::span<const uint8_t>)
    {
        FAIL("Invalid DMA instruction");
    }
};

static std::string dump_dma(const std::vector<uint8_t> &code)
{
    std::stringstream stm;
    uint32_t nbanks;
    memcpy(&nbanks, code.data(), 4);
    for (uint32_t i = 0; i < nbanks; i++) {
        uint32_t mask;
        memcpy(&mask, code.data() + 4 + i * 4, 4);
        stm << "# ttl_mask=0x" << std::hex << mask << std::dec << std::endl;
    }
    Zynq::DMA::print(stm, std::span(code).subspan(4 + nbanks * 4), 0);
    return stm.str();
}

static int64_t dds_write_time(const Zynq::BCGen::DDSTiming &timing, const DDSWrite &write)
{
    if (write.is32)
        return write.fud ? timing.set32_fud : timing.set32;
    return write.fud ? timing.set16_fud : timing.set16;
}

// Check the constraints of the hardware on the decoded output.
static void check_hw_rules(const Zynq::BCGen &bc_gen, const Decoded &res)
{
    REQUIRE(res.ends_with_wait);
    for (int bus = 0; bus < 2; bus++) {
        const DDSWrite *prev = nullptr;
        for (auto &write: res.dds) {
            if (write.bus != bus)
                continue;
            if (prev) {
                INFO("DDS write at " << write.t << " after " << prev->t);
                // The previous write on the bus must have finished.
                REQUIRE(write.t - prev->t >= dds_write_time(bc_gen.dds_timing, *prev));
                if (!prev->fud) {
                    // The write without FUD must be followed by more writes
                    // on the same bus without any long wait
                    // (during which the CPU commands may run).
                    REQUIRE(write.t - prev->t < 128);
                }
            }
            prev = &write;
        }
        if (prev) {
            REQUIRE(prev->fud);
            REQUIRE(res.total >= prev->t + dds_write_time(bc_gen.dds_timing, *prev));
        }
    }
    for (size_t i = 1; i < res.dac.size(); i++)
        REQUIRE(res.dac[i].first - res.dac[i - 1].first >= 76);
    if (!res.dac.empty()) {
        REQUIRE(res.total >= res.dac.back().first + 76);
    }
}

static Decoded generate_dma(Zynq::BCGen &bc_gen, HostSeq &host_seq)
{
    bc_gen.dma = true;
    REQUIRE(bc_gen.version() == 0);
    bc_gen.generate(host_seq);
    auto code = bc_gen.bytecode;
    INFO(dump_dma(code));
    // Generation is repeatable.
    for (int i = 0; i < 3; i++) {
        bc_gen.generate(host_seq);
        REQUIRE(code == bc_gen.bytecode);
    }

    Decoded res;
    REQUIRE(code.size() >= 8);
    uint32_t nbanks;
    memcpy(&nbanks, code.data(), 4);
    REQUIRE(nbanks >= 1);
    REQUIRE(nbanks <= 8);
    REQUIRE(code.size() >= 4 + nbanks * 4);
    res.masks.resize(nbanks);
    memcpy(res.masks.data(), code.data() + 4, nbanks * 4);
    auto insts = std::span(code).subspan(4 + nbanks * 4);
    REQUIRE(insts.size() % 2 == 0);

    DMADecoder decoder(res, bc_gen.dds_timing);
    Zynq::DMA::ExeState state;
    state.run(decoder, insts, 0);
    res.total = decoder.t;
    REQUIRE(uint64_t(res.total) == Zynq::DMA::total_time(insts, 0));
    check_hw_rules(bc_gen, res);
    return res;
}

static Zynq::BCGen::SeqPulse value_pulse(Zynq::ChnType chn_type, uint8_t chn, uint32_t id,
                                         uint32_t time_id, uint32_t endvalue_id)
{
    return {
        .chn_type = chn_type,
        .pulse_type = Zynq::BCGen::PulseType::Value,
        .chn = chn,

        .id = id,
        .time_id = time_id,
        .len_id = uint32_t(-1),
        .endvalue_id = endvalue_id,
        .cond_id = uint32_t(-1),
        .ramp_func = nullptr
    };
}

static void init_bc_gen(Zynq::BCGen &bc_gen)
{
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask[0] = 0;
    bc_gen.fpga_clock_div = ps_per_cycle;
    bc_gen.len_ns = 0;
}

// The start of the sequence, i.e. when the start trigger turns off
// (with `seq_delay == 0`)
static int64_t seq_start(const Decoded &res, int start_chn=0)
{
    auto bank = start_chn / 32;
    auto bit = uint32_t(1) << (start_chn % 32);
    auto t_hi = res.ttl_time(bank, bit, bit);
    REQUIRE(t_hi >= 0);
    auto t_lo = res.ttl_time(bank, bit, 0, t_hi);
    // High for 1us
    REQUIRE(t_lo == t_hi + cycles_per_us);
    return t_lo;
}

static uint32_t mu(Zynq::ChnType type, double val)
{
    return Zynq::BCGen::convert_value(type, val);
}

// The DDS update of `chn` at time `t`
static const DDSUpdate *find_update(const Decoded &res, int chn, int64_t t)
{
    for (auto &update: res.dds_updates) {
        if (update.chn == chn && update.t == t) {
            return &update;
        }
    }
    return nullptr;
}

static std::vector<DDSUpdate> chn_updates(const Decoded &res, int chn)
{
    std::vector<DDSUpdate> updates;
    for (auto &update: res.dds_updates) {
        if (update.chn == chn) {
            updates.push_back(update);
        }
    }
    return updates;
}

// Values of a channel extracted from the instructions (time in DMA cycles
// relative to the start of the sequence).
static std::pair<std::vector<uint64_t>,std::vector<double>>
extract(const Zynq::BCGen &bc_gen, Zynq::ChnType type, uint8_t num)
{
    auto &code = bc_gen.bytecode;
    uint32_t nbanks;
    memcpy(&nbanks, code.data(), 4);
    return Zynq::DMA::extract_values(std::span(code).subspan(4 + nbanks * 4), 0, type, num,
                                     bc_gen.start_ttl_chn, bc_gen.dds_timing);
}

// Check the extracted DDS values against the decoded FUD's.
// Every extracted value must be from a FUD updating the DDS and the parameter
// must not change at the FUD's without an extracted value.
static void check_extract_dds(const Zynq::BCGen &bc_gen, const Decoded &res,
                              int64_t t0, uint8_t chn)
{
    using Zynq::ChnType;
    for (auto type: {ChnType::Freq, ChnType::Amp, ChnType::Phase}) {
        INFO("DDS " << int(chn) << " type " << int(type));
        auto [ts, vals] = extract(bc_gen, type, chn);
        auto param = [&] (const DDSUpdate &update) -> double {
            if (type == ChnType::Freq)
                return Zynq::dds_freq_from_mu(update.freq);
            if (type == ChnType::Amp)
                return Zynq::dds_amp_from_mu(update.amp);
            return Zynq::dds_phase_from_mu(update.phase) * 360;
        };
        size_t idx = 0;
        // The registers that are never written are 0 in the decoder.
        double prev = 0;
        for (auto &update: chn_updates(res, chn)) {
            auto t = uint64_t(std::max(update.t - t0, int64_t(0)));
            auto val = param(update);
            if (idx < ts.size() && ts[idx] == t && vals[idx] == val) {
                idx++;
            }
            else {
                // Not written, no change
                REQUIRE(val == prev);
            }
            prev = val;
        }
        REQUIRE(idx == ts.size());
    }
}

static std::vector<DDSWrite> bus_writes(const Decoded &res, int bus)
{
    std::vector<DDSWrite> writes;
    for (auto &write: res.dds) {
        if (write.bus == bus) {
            writes.push_back(write);
        }
    }
    return writes;
}

}

TEST_CASE("empty") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    auto res = generate_dma(bc_gen, host_seq);
    REQUIRE(res.masks == std::vector<uint32_t>{1});

    bc_gen.generate(host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_dma(bc_gen.bytecode));
    REQUIRE(checker.cmp<uint32_t>(1)); // nbanks
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    using namespace Zynq::DMA::Inst_v0;
    // Initial wait for buffering
    REQUIRE(checker.cmp(Wait1(199)));
    // Start trigger
    REQUIRE(checker.cmp(TTLSet4(0, 1)));
    REQUIRE(checker.cmp(Wait1(199)));
    REQUIRE(checker.cmp(TTLSet4(0, 0)));
    // End
    REQUIRE(checker.cmp(Wait1(0)));
    REQUIRE(checker.is_end());
}

TEST_CASE("concurrent") {
    // TTL, both DDS buses and the DAC at the same time.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask[0] = 0x8;
    bc_gen.len_ns = 10000;
    bc_gen.start_vals[{Zynq::ChnType::TTL, 3}] = 0;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<double>(2, 100e6);
    host_seq.set_value<double>(3, 0.5);
    host_seq.set_value<double>(4, 0.25);
    host_seq.set_value<double>(5, 50e6);
    host_seq.set_value<double>(6, 2.0);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 3, 1, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Freq, 3, 2, 0, 2));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Amp, 3, 3, 0, 3));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Phase, 3, 4, 0, 4));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Freq, 14, 5, 0, 5));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::DAC, 1, 6, 0, 6));
    auto res = generate_dma(bc_gen, host_seq);
    REQUIRE(res.masks == std::vector<uint32_t>{0x9});
    auto &timing = bc_gen.dds_timing;

    auto t0 = seq_start(res);
    auto t = t0 + cycles_per_us;
    // All outputs happen at the same time.
    REQUIRE(res.ttl_time(0, 0x8, 0x8) == t);
    REQUIRE(res.dds.size() == 3);
    // DDS 3 (bus 0): frequency without FUD, then phase and amplitude with FUD,
    // started early so that the FUD happens at `t`.
    auto freq3 = mu(Zynq::ChnType::Freq, 100e6);
    auto amp3 = mu(Zynq::ChnType::Amp, 0.5);
    auto phase3 = mu(Zynq::ChnType::Phase, 0.25);
    auto bus0 = bus_writes(res, 0);
    REQUIRE(bus0.size() == 2);
    REQUIRE(bus0[0].id == 3);
    REQUIRE(bus0[0].addr == 0x2c >> 1);
    REQUIRE(bus0[0].is32);
    REQUIRE(!bus0[0].fud);
    REQUIRE(bus0[0].data == freq3);
    REQUIRE(bus0[0].t == t - timing.set32 - (timing.set32_fud - timing.fud_to_end));
    // Phase and amplitude in one write right after the frequency
    REQUIRE(bus0[1].t == bus0[0].t + timing.set32);
    REQUIRE(bus0[1].addr == 0x30 >> 1);
    REQUIRE(bus0[1].is32);
    REQUIRE(bus0[1].fud);
    // DDS 14 (bus 1)
    auto bus1 = bus_writes(res, 1);
    REQUIRE(bus1.size() == 1);
    REQUIRE(bus1[0].id == 3);
    REQUIRE(bus1[0].fud);
    REQUIRE(bus1[0].data == mu(Zynq::ChnType::Freq, 50e6));
    // All the parameters of DDS 3 are updated with a single FUD at `t`
    // and DDS 14 at the same time.
    REQUIRE(res.dds_updates.size() == 2);
    auto update3 = find_update(res, 3, t);
    REQUIRE(update3);
    REQUIRE(update3->freq == freq3);
    REQUIRE(update3->amp == amp3);
    REQUIRE(update3->phase == phase3);
    auto update14 = find_update(res, 14, t);
    REQUIRE(update14);
    REQUIRE(update14->freq == mu(Zynq::ChnType::Freq, 50e6));
    // DAC at the same time
    REQUIRE(res.dac.size() == 1);
    REQUIRE(res.dac[0].first == t);
    REQUIRE(res.dac[0].second == ((1u << 16) | mu(Zynq::ChnType::DAC, 2.0)));
    REQUIRE(res.total >= t0 + 10 * cycles_per_us);

    // Extracted values
    using Zynq::ChnType;
    auto single = [] (uint64_t t, double v) {
        return std::make_pair(std::vector<uint64_t>{t}, std::vector<double>{v});
    };
    REQUIRE(extract(bc_gen, ChnType::Freq, 3) ==
            single(cycles_per_us, Zynq::dds_freq_from_mu(freq3)));
    REQUIRE(extract(bc_gen, ChnType::Amp, 3) ==
            single(cycles_per_us, Zynq::dds_amp_from_mu(amp3)));
    REQUIRE(extract(bc_gen, ChnType::Phase, 3) ==
            single(cycles_per_us, Zynq::dds_phase_from_mu(phase3) * 360));
    REQUIRE(extract(bc_gen, ChnType::Freq, 14) ==
            single(cycles_per_us, Zynq::dds_freq_from_mu(mu(ChnType::Freq, 50e6))));
    REQUIRE(extract(bc_gen, ChnType::Freq, 4).first.empty());
    REQUIRE(extract(bc_gen, ChnType::DAC, 1) ==
            single(cycles_per_us, Zynq::dac_from_mu(mu(ChnType::DAC, 2.0))));
    REQUIRE(extract(bc_gen, ChnType::DAC, 0).first.empty());
    // TTL 3 is in the same 4-bit sub-bank as the start trigger.
    REQUIRE(extract(bc_gen, ChnType::TTL, 3) ==
            std::make_pair(std::vector<uint64_t>{0, 0, uint64_t(cycles_per_us)},
                           std::vector<double>{0, 0, 1}));
}

TEST_CASE("dds_shared_fud") {
    // Two DDS on the same bus are updated with a single FUD.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 0.5);
    host_seq.set_value<double>(2, 0.7);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Amp, 2, 1, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Amp, 7, 2, 0, 2));
    auto res = generate_dma(bc_gen, host_seq);
    auto &timing = bc_gen.dds_timing;
    auto t = seq_start(res) + cycles_per_us;
    REQUIRE(res.dds.size() == 2);
    REQUIRE(res.dds[0].id == 2);
    REQUIRE(!res.dds[0].is32);
    REQUIRE(!res.dds[0].fud);
    REQUIRE(res.dds[0].data == mu(Zynq::ChnType::Amp, 0.5));
    REQUIRE(res.dds[1].t == res.dds[0].t + timing.set16);
    REQUIRE(res.dds[1].id == 7);
    REQUIRE(!res.dds[1].is32);
    REQUIRE(res.dds[1].fud);
    REQUIRE(res.dds[1].data == mu(Zynq::ChnType::Amp, 0.7));
    REQUIRE(res.dds_updates.size() == 2);
    REQUIRE(find_update(res, 2, t));
    REQUIRE(find_update(res, 2, t)->amp == mu(Zynq::ChnType::Amp, 0.5));
    REQUIRE(find_update(res, 7, t));
    REQUIRE(find_update(res, 7, t)->amp == mu(Zynq::ChnType::Amp, 0.7));
}

TEST_CASE("dds_timing") {
    // The FUD happens at the requested time with a different DDS timing.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.dds_timing = {20, 30, 40, 50, 5};
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<double>(2, 0.7);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Freq, 2, 1, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Amp, 13, 2, 0, 2));
    auto res = generate_dma(bc_gen, host_seq);
    auto t = seq_start(res) + cycles_per_us;
    REQUIRE(res.dds.size() == 2);
    REQUIRE(bus_writes(res, 0)[0].t == t - 45);
    REQUIRE(bus_writes(res, 1)[0].t == t - 25);
    REQUIRE(find_update(res, 2, t));
    REQUIRE(find_update(res, 13, t));

    // Invalid timing
    bc_gen.dds_timing = {20, 30, 40, 50, 0};
    REQUIRE_THROWS_WITH(bc_gen.generate(host_seq), Catch::Contains("Invalid DDS timing"));
    bc_gen.dds_timing = {20, 30, 40, 50, 30};
    REQUIRE_THROWS_WITH(bc_gen.generate(host_seq), Catch::Contains("Invalid DDS timing"));
}

namespace {

static Zynq::BCGen::SeqPulse ramp_pulse(Zynq::ChnType chn_type, uint8_t chn, uint32_t id,
                                        uint32_t time_id, uint32_t len_id,
                                        uint32_t endvalue_id,
                                        double (*func)(double, const void*))
{
    return {
        .chn_type = chn_type,
        .pulse_type = Zynq::BCGen::PulseType::Scalar,
        .chn = chn,

        .id = id,
        .time_id = time_id,
        .len_id = len_id,
        .endvalue_id = endvalue_id,
        .cond_id = uint32_t(-1),
        .ramp_func = (void (*)(void))func
    };
}

}

TEST_CASE("dds_ramp") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 20 * us);
    bc_gen.seq_pulses.push_back(ramp_pulse(Zynq::ChnType::Freq, 3, 1, 0, 2, 1,
                                           [] (double in, const void*) {
                                               return in * 5;
                                           }));
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    auto start = t0 + cycles_per_us;
    auto end = t0 + 21 * cycles_per_us;
    // One sample every `set32_fud` cycles as long as it finishes before
    // the end value needs to be written, plus the end value.
    auto write_t = int64_t(bc_gen.dds_timing.set32_fud);
    REQUIRE(res.dds_updates.size() == size_t((20 * cycles_per_us - write_t) / write_t + 2));
    // The start of the ramp at the requested time and as fast as the bus allows after that.
    REQUIRE(res.dds_updates[0].t == start);
    for (size_t i = 0; i < res.dds_updates.size(); i++) {
        auto &update = res.dds_updates[i];
        REQUIRE(update.chn == 3);
        if (i + 1 == res.dds_updates.size())
            break;
        REQUIRE(update.t == start + int64_t(i) * write_t);
        // Value at the time of the update
        REQUIRE(update.freq == mu(Zynq::ChnType::Freq,
                                  double(update.t - start) * ps_per_cycle * 5));
    }
    for (auto &write: res.dds) {
        REQUIRE(write.fud);
        REQUIRE(write.is32);
        REQUIRE(write.addr == 0x2c >> 1);
    }
    // End value at the end of the ramp
    REQUIRE(res.dds_updates.back().t == end);
    REQUIRE(res.dds_updates.back().freq == mu(Zynq::ChnType::Freq, 100e6));

    auto [ts, vals] = extract(bc_gen, Zynq::ChnType::Freq, 3);
    REQUIRE(ts.size() == res.dds_updates.size());
    for (size_t i = 0; i < ts.size(); i++) {
        REQUIRE(int64_t(ts[i]) == res.dds_updates[i].t - t0);
        REQUIRE(vals[i] == Zynq::dds_freq_from_mu(res.dds_updates[i].freq));
    }
    check_extract_dds(bc_gen, res, t0, 3);
}

TEST_CASE("dds_ramp_freq16") {
    // Slow frequency ramp that only changes the lower 16 bits most of the time
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 20 * us);
    auto func = [] (double in, const void*) { return 10e6 + in * 1e-3; };
    bc_gen.seq_pulses.push_back(ramp_pulse(Zynq::ChnType::Freq, 3, 1, 0, 2, 1, func));
    auto res = generate_dma(bc_gen, host_seq);
    auto start = seq_start(res) + cycles_per_us;
    auto end = start + 20 * cycles_per_us;
    auto &timing = bc_gen.dds_timing;
    // The first write sets the whole frequency
    REQUIRE(res.dds[0].is32);
    size_t n32 = 0;
    for (auto &write: res.dds) {
        REQUIRE(write.fud);
        if (write.is32) {
            REQUIRE(write.addr == 0x2c >> 1);
            n32++;
        }
        else {
            REQUIRE((write.addr == 0x2c >> 1 || write.addr == 0x2e >> 1));
        }
    }
    // Only the start (and maybe when the higher 16 bits changes)
    REQUIRE(n32 <= 2);
    REQUIRE(res.dds_updates.size() > 100);
    REQUIRE(res.dds_updates[0].t == start);
    for (size_t i = 0; i < res.dds_updates.size(); i++) {
        auto &update = res.dds_updates[i];
        if (i + 1 == res.dds_updates.size()) {
            REQUIRE(update.t == end);
            REQUIRE(update.freq == mu(Zynq::ChnType::Freq, 10e6));
            break;
        }
        REQUIRE(update.freq == mu(Zynq::ChnType::Freq,
                                  func(double(update.t - start) * ps_per_cycle, nullptr)));
        if (i > 0) {
            // As fast as the 16-bit write allows.
            REQUIRE(update.t - res.dds_updates[i - 1].t >= timing.set16_fud);
            if (!res.dds[i].is32 && !res.dds[i - 1].is32) {
                REQUIRE(update.t - res.dds_updates[i - 1].t == timing.set16_fud);
            }
        }
    }
}

TEST_CASE("dds_ramp_shared_fud") {
    // Ramps on two DDS's on the same bus are updated together.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 20 * us);
    host_seq.set_value<double>(3, 0.5);
    auto freq_ramp = [] (double in, const void*) { return in * 5; };
    auto amp_ramp = [] (double in, const void*) { return in * 2e-8; };
    bc_gen.seq_pulses.push_back(ramp_pulse(Zynq::ChnType::Freq, 1, 1, 0, 2, 1, freq_ramp));
    bc_gen.seq_pulses.push_back(ramp_pulse(Zynq::ChnType::Amp, 2, 2, 0, 2, 3, amp_ramp));
    auto res = generate_dma(bc_gen, host_seq);
    auto start = seq_start(res) + cycles_per_us;
    auto end = start + 20 * cycles_per_us;
    auto updates1 = chn_updates(res, 1);
    auto updates2 = chn_updates(res, 2);
    REQUIRE(updates1.size() > 10);
    REQUIRE(updates1.front().t == start);
    REQUIRE(updates2.front().t == start);
    REQUIRE(updates1.back().t == end);
    REQUIRE(updates2.back().t == end);
    REQUIRE(updates1.back().freq == mu(Zynq::ChnType::Freq, 100e6));
    REQUIRE(updates2.back().amp == mu(Zynq::ChnType::Amp, 0.5));
    // Every FUD updates both DDS's
    REQUIRE(updates1.size() == updates2.size());
    for (size_t i = 0; i < updates1.size(); i++) {
        REQUIRE(updates1[i].t == updates2[i].t);
        if (i + 1 < updates1.size()) {
            auto t = double(updates1[i].t - start) * ps_per_cycle;
            REQUIRE(updates1[i].freq == mu(Zynq::ChnType::Freq, freq_ramp(t, nullptr)));
            REQUIRE(updates2[i].amp == mu(Zynq::ChnType::Amp, amp_ramp(t, nullptr)));
        }
    }
    // One FUD for both
    size_t nfud = 0;
    for (auto &write: res.dds)
        nfud += write.fud;
    REQUIRE(nfud == updates1.size());
}

TEST_CASE("dds_ramp+fixed") {
    // A DDS update during a ramp of another DDS on the same bus
    // happens at the requested time.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 20 * us);
    host_seq.set_value<int64_t>(3, 10 * us);
    host_seq.set_value<double>(4, 0.3);
    bc_gen.seq_pulses.push_back(ramp_pulse(Zynq::ChnType::Freq, 3, 1, 0, 2, 1,
                                           [] (double in, const void*) {
                                               return in * 5;
                                           }));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::Amp, 5, 2, 3, 4));
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    auto updates5 = chn_updates(res, 5);
    REQUIRE(updates5.size() == 1);
    REQUIRE(updates5[0].t == t0 + 10 * cycles_per_us);
    REQUIRE(updates5[0].amp == mu(Zynq::ChnType::Amp, 0.3));
    auto updates3 = chn_updates(res, 3);
    REQUIRE(updates3.front().t == t0);
    REQUIRE(updates3.back().t == t0 + 20 * cycles_per_us);
    REQUIRE(updates3.back().freq == mu(Zynq::ChnType::Freq, 100e6));
}

TEST_CASE("ttl_multi_bank") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask = {0x8, 0, 0x101};
    bc_gen.start_ttl_chn = 33;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 2 * us);
    host_seq.set_value<bool>(3, false);
    host_seq.set_value<int64_t>(4, 2 * us + 7 * ps_per_cycle);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 3, 1, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 64, 2, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 72, 3, 2, 1));
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 3, 4, 4, 3));
    auto res = generate_dma(bc_gen, host_seq);
    REQUIRE(res.masks == std::vector<uint32_t>{0x8, 0x2, 0x101});
    auto t0 = seq_start(res, 33);
    REQUIRE(res.ttl_time(0, 0x8, 0x8) == t0 + cycles_per_us);
    REQUIRE(res.ttl_time(2, 0x1, 0x1) == t0 + cycles_per_us);
    REQUIRE(res.ttl_time(2, 0x100, 0x100) == t0 + 2 * cycles_per_us);
    // TTL pulses are not delayed even if they are close to each other
    REQUIRE(res.ttl_time(0, 0x8, 0, t0) == t0 + 2 * cycles_per_us + 7);
}

TEST_CASE("first_bseq") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.first_bseq = true;
    bc_gen.ttl_mask[0] = 0x6;
    bc_gen.start_vals[{Zynq::ChnType::TTL, 1}] = 1;
    bc_gen.start_vals[{Zynq::ChnType::TTL, 2}] = 0;
    bc_gen.start_vals[{Zynq::ChnType::Freq, 3}] = 0x12345678;
    bc_gen.start_vals[{Zynq::ChnType::Amp, 3}] = 0x123;
    bc_gen.start_vals[{Zynq::ChnType::Phase, 3}] = 0x4567;
    bc_gen.start_vals[{Zynq::ChnType::Freq, 14}] = 0x23456789;
    bc_gen.start_vals[{Zynq::ChnType::DAC, 2}] = 0x8000;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<bool>(1, false);
    // t=0 pulse merged into the start value
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 2, 1, 0, 1));
    auto res = generate_dma(bc_gen, host_seq);
    // Initialization after the initial wait
    constexpr int64_t init_t = 200;
    REQUIRE(res.ttl_time(0, 0x7, 0x2) == init_t);
    REQUIRE(res.dds.size() == 3);
    REQUIRE(res.dds[0].t == init_t);
    REQUIRE(res.dds[0].bus == 0);
    REQUIRE(!res.dds[0].fud);
    REQUIRE(res.dds[1].t == init_t);
    REQUIRE(res.dds[1].bus == 1);
    REQUIRE(res.dds[1].fud);
    REQUIRE(res.dds[2].t == init_t + bc_gen.dds_timing.set32);
    REQUIRE(res.dds[2].fud);
    REQUIRE(res.dds_updates.size() == 2);
    REQUIRE(res.dds_updates[0].chn == 14);
    REQUIRE(res.dds_updates[0].freq == 0x23456789);
    REQUIRE(res.dds_updates[1].chn == 3);
    REQUIRE(res.dds_updates[1].freq == 0x12345678);
    REQUIRE(res.dds_updates[1].amp == 0x123);
    REQUIRE(res.dds_updates[1].phase == 0x4567);
    REQUIRE(res.dac.size() == 1);
    REQUIRE(res.dac[0].first == init_t);
    REQUIRE(res.dac[0].second == ((2u << 16) | 0x8000));
    // The sequence starts after the initialization
    auto t0 = seq_start(res);
    REQUIRE(t0 >= init_t + 2500 * 2);
}

TEST_CASE("trigger") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask[0] = 0x10;
    bc_gen.seq_delay = 500;
    bc_gen.len_ns = 50000;
    bc_gen.trig_type = Zynq::BCGen::TrigType::Raise;
    bc_gen.trig_channel = 5;
    bc_gen.trig_timeout_cycle = 0x123456;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<bool>(1, true);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 4, 1, 0, 1));
    auto res = generate_dma(bc_gen, host_seq);
    REQUIRE(res.trig.size() == 1);
    auto [trig_t, trig] = res.trig[0];
    REQUIRE(trig.chn == 5);
    REQUIRE(trig.edge == 1);
    REQUIRE(trig.cycle == 0x123456 * 2);
    // Start trigger after the trigger
    auto t_hi = res.ttl_time(0, 1, 1);
    REQUIRE(t_hi == trig_t);
    auto t_lo = res.ttl_time(0, 1, 0, t_hi);
    REQUIRE(t_lo == t_hi + 200);
    // Sequence starts `seq_delay` after the start trigger
    auto t0 = t_lo + 500 * 2;
    REQUIRE(res.ttl_time(0, 0x10, 0x10) == t0 + cycles_per_us);
    REQUIRE(res.total >= t0 + 50 * cycles_per_us);

    // The extracted time is relative to the end of the start trigger,
    // i.e. including `seq_delay`.
    auto [ts, vals] = extract(bc_gen, Zynq::ChnType::TTL, 4);
    REQUIRE(ts.back() == uint64_t(500 * 2 + cycles_per_us));
    REQUIRE(vals.back() == 1);
}

TEST_CASE("long_wait") {
    // Waits longer than a single wait instruction.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask[0] = 0x10;
    bc_gen.len_ns = 5000000000;
    host_seq.set_value<int64_t>(0, 3000000 * us);
    host_seq.set_value<bool>(1, true);
    bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, 4, 1, 0, 1));
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    REQUIRE(res.ttl_time(0, 0x10, 0x10) == t0 + 3000000 * cycles_per_us);
    REQUIRE(res.total == t0 + 5000000 * cycles_per_us);
}

TEST_CASE("clock") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.clocks.push_back({-1 * us, 100});
    bc_gen.clocks.push_back({5 * us, 0});
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    REQUIRE(res.clock.size() == 4);
    // Toggle every 100 x 10ns
    REQUIRE(res.clock[0] == std::make_pair(t0 - cycles_per_us, uint16_t(199)));
    REQUIRE(res.clock[1] == std::make_pair(t0 + 5 * cycles_per_us, uint16_t(0x1ff)));
    // Clock output for 10ms at the end
    REQUIRE(res.clock[2].first == t0 + 5 * cycles_per_us + 1);
    REQUIRE(res.clock[2].second == 19);
    REQUIRE(res.clock[3].first == res.clock[2].first + 10000 * cycles_per_us);
    REQUIRE(res.clock[3].second == 0x1ff);

    // Different length
    bc_gen.clock_tail_ns = 1000000;
    auto res2 = generate_dma(bc_gen, host_seq);
    REQUIRE(res2.clock.size() == 4);
    REQUIRE(res2.clock[3].first == res2.clock[2].first + 1000 * cycles_per_us);
    REQUIRE(res2.clock[3].second == 0x1ff);

    // Disabled
    bc_gen.clock_tail_ns = 0;
    auto res3 = generate_dma(bc_gen, host_seq);
    REQUIRE(res3.clock.size() == 2);
    REQUIRE(res3.total == seq_start(res3) + 5 * cycles_per_us + 1);
}

namespace {

// Construct the instruction on top of garbage to check that the reserved bits are cleared.
template<typename T, typename... Args>
static uint64_t inst_bits(Args... args)
{
    alignas(T) uint8_t buff[sizeof(T)];
    volatile uint8_t *vbuff = buff;
    for (size_t i = 0; i < sizeof(T); i++)
        vbuff[i] = 0xff;
    auto inst = new (buff) T(args...);
    uint64_t res = 0;
    memcpy(&res, inst, sizeof(T));
    return res;
}

}

TEST_CASE("encoding") {
    using namespace Zynq::DMA::Inst_v0;
    // The bit layout from the instruction format, all the reserved bits must be 0.
    REQUIRE(inst_bits<Wait1>(199) == ((199u << 4) | 0));
    REQUIRE(inst_bits<TTLSet4>(3, 0xa) == ((0xau << 10) | (3u << 4) | (1u << 2)));
    REQUIRE(inst_bits<ClockOut>(0x1ff) == ((0x1ffu << 4) | (3u << 2)));
    REQUIRE(inst_bits<TTLSet16>(1, 0xab, 2, 0xcd) ==
            ((0xcdull << 22) | (2ull << 17) | (0xabull << 9) | (1ull << 4) |
             (1ull << 2) | 1ull));
    REQUIRE(inst_bits<TTLSet32>(1, 0x1234, 2, 0x5678) ==
            ((0x5678ull << 28) | (2ull << 24) | (0x1234ull << 8) | (1ull << 4) |
             (1ull << 2) | 2ull));
    REQUIRE(inst_bits<DAC>(0, 1, 0, 0, 0x2abcd) ==
            ((0x2abcdull << 17) | (1ull << 6) | (3ull << 2) | 2ull));

    // Raw bytes of the empty sequence
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.dma = true;
    bc_gen.generate(host_seq);
    REQUIRE(bc_gen.bytecode == std::vector<uint8_t>{
            1, 0, 0, 0, // nbanks
            1, 0, 0, 0, // ttl_mask
            0x70, 0x0c, // wait1(199)
            0x04, 0x04, // ttl_set4([0:3]=1)
            0x70, 0x0c, // wait1(199)
            0x04, 0x00, // ttl_set4([0:3]=0)
            0x00, 0x00, // wait1(0)
        });
}

TEST_CASE("dds_ramp_vector") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 20 * us);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Vector,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 2,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(void (*)(double*, const double*, const void*))
            [] (double *out, const double *in, const void*) {
                for (unsigned i = 0; i < vector_size; i++) {
                    out[i] = in[i] * 5;
                }
            }
        });
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    REQUIRE(res.dds_updates.size() > 10);
    for (auto &update: res.dds_updates) {
        auto expected = mu(Zynq::ChnType::Freq,
                           double(update.t - t0 - cycles_per_us) * ps_per_cycle * 5);
        if (update.t < t0 + 21 * cycles_per_us) {
            REQUIRE(update.freq == expected);
        }
    }
    REQUIRE(res.dds_updates.back().t == t0 + 21 * cycles_per_us);
    REQUIRE(res.dds_updates.back().freq == mu(Zynq::ChnType::Freq, 100e6));
}

TEST_CASE("stress") {
    // Overlapping ramps and fixed pulses on all the buses.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask[0] = 0xff00;
    bc_gen.len_ns = 100000;
    // Values
    host_seq.set_value<double>(0, 100e6);
    host_seq.set_value<double>(1, 0.9);
    host_seq.set_value<double>(2, 0.6);
    host_seq.set_value<double>(3, 5.0);
    host_seq.set_value<double>(4, 0.3);
    host_seq.set_value<double>(5, 30e6);
    host_seq.set_value<double>(6, -2.0);
    host_seq.set_value<bool>(7, true);
    host_seq.set_value<bool>(8, false);
    // Times
    for (int i = 0; i < 40; i++)
        host_seq.set_value<int64_t>(20 + i, int64_t(i) * us + i * 37 * ps_per_cycle);
    // Lengths
    host_seq.set_value<int64_t>(70, 20 * us);
    host_seq.set_value<int64_t>(71, 7 * us);
    host_seq.set_value<int64_t>(72, 30 * us);

    uint32_t id = 1;
    auto freq_ramp = [] (double in, const void*) { return in * 5; };
    auto amp_ramp = [] (double in, const void*) { return in * 3e-8; };
    auto dac_ramp = [] (double in, const void*) { return in * 1e-7 - 1; };
    using Zynq::ChnType;
    // Bus 0: frequency ramps on DDS 1, 2 and amplitude ramp on DDS 4,
    // overlapping with each other and the fixed pulses on DDS 1 and 6.
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 1, id++, 21, 70, 0, freq_ramp));
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 2, id++, 23, 72, 0, freq_ramp));
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Amp, 4, id++, 22, 71, 1, amp_ramp));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Phase, 1, id++, 25, 4));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Amp, 6, id++, 26, 2));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Amp, 1, id++, 26, 1));
    // Bus 1: ramp on DDS 12 and fixed pulses on DDS 15 and 12
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 12, id++, 21, 72, 0, freq_ramp));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Freq, 15, id++, 24, 5));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Phase, 12, id++, 27, 4));
    // DAC ramp on 0 and fixed value on 3
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::DAC, 0, id++, 22, 70, 3, dac_ramp));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::DAC, 3, id++, 28, 6));
    // TTL pulses at the same time as other outputs
    for (int i = 0; i < 30; i++) {
        bc_gen.seq_pulses.push_back(value_pulse(ChnType::TTL, uint8_t(8 + (i / 2) % 8), id++,
                                                uint32_t(20 + i), i % 2 ? 8 : 7));
    }
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);

    // Final values
    std::map<int,DDSUpdate> last_dds;
    for (auto &update: res.dds_updates)
        last_dds[update.chn] = update;
    REQUIRE(last_dds[1].freq == mu(ChnType::Freq, 100e6));
    REQUIRE(last_dds[1].amp == mu(ChnType::Amp, 0.9));
    REQUIRE(last_dds[1].phase == mu(ChnType::Phase, 0.3));
    REQUIRE(last_dds[2].freq == mu(ChnType::Freq, 100e6));
    REQUIRE(last_dds[4].amp == mu(ChnType::Amp, 0.9));
    REQUIRE(last_dds[6].amp == mu(ChnType::Amp, 0.6));
    REQUIRE(last_dds[12].freq == mu(ChnType::Freq, 100e6));
    REQUIRE(last_dds[12].phase == mu(ChnType::Phase, 0.3));
    REQUIRE(last_dds[15].freq == mu(ChnType::Freq, 30e6));
    std::map<int,uint32_t> last_dac;
    for (auto [t, data]: res.dac)
        last_dac[data >> 16] = data & 0xffff;
    REQUIRE(last_dac[0] == mu(ChnType::DAC, 5.0));
    REQUIRE(last_dac[3] == mu(ChnType::DAC, -2.0));

    // The fixed values (including the start and end of the ramps) take effect
    // at the requested time without being delayed by the ramps.
    auto pulse_t = [&] (int i) { return t0 + i * cycles_per_us + i * 37; };
    auto check_update = [&] (int chn, int64_t t) {
        INFO("DDS " << chn << " at " << t);
        auto update = find_update(res, chn, t);
        REQUIRE(update);
        return *update;
    };
    // Start of the ramps
    REQUIRE(check_update(1, pulse_t(1)).freq == mu(ChnType::Freq, 0));
    REQUIRE(check_update(12, pulse_t(1)).freq == mu(ChnType::Freq, 0));
    REQUIRE(check_update(4, pulse_t(2)).amp == mu(ChnType::Amp, 0));
    REQUIRE(check_update(2, pulse_t(3)).freq == mu(ChnType::Freq, 0));
    // Fixed values
    REQUIRE(check_update(15, pulse_t(4)).freq == mu(ChnType::Freq, 30e6));
    REQUIRE(check_update(1, pulse_t(5)).phase == mu(ChnType::Phase, 0.3));
    // Two DDS's on the same bus at the same time with one FUD
    REQUIRE(check_update(6, pulse_t(6)).amp == mu(ChnType::Amp, 0.6));
    REQUIRE(check_update(1, pulse_t(6)).amp == mu(ChnType::Amp, 0.9));
    REQUIRE(check_update(12, pulse_t(7)).phase == mu(ChnType::Phase, 0.3));
    // End of the ramps
    REQUIRE(check_update(4, pulse_t(2) + 7 * cycles_per_us).amp == mu(ChnType::Amp, 0.9));
    REQUIRE(check_update(1, pulse_t(1) + 20 * cycles_per_us).freq ==
            mu(ChnType::Freq, 100e6));
    REQUIRE(check_update(12, pulse_t(1) + 30 * cycles_per_us).freq ==
            mu(ChnType::Freq, 100e6));
    REQUIRE(check_update(2, pulse_t(3) + 30 * cycles_per_us).freq ==
            mu(ChnType::Freq, 100e6));
    // DAC
    auto find_dac = [&] (int64_t t) -> uint32_t {
        for (auto [dac_t, data]: res.dac) {
            if (dac_t == t) {
                return data;
            }
        }
        return 0;
    };
    REQUIRE(find_dac(pulse_t(2)) == mu(ChnType::DAC, -1));
    REQUIRE(find_dac(pulse_t(8)) == ((3u << 16) | mu(ChnType::DAC, -2.0)));
    REQUIRE(find_dac(pulse_t(2) + 20 * cycles_per_us) == mu(ChnType::DAC, 5.0));

    // Extracted values
    for (uint8_t chn: {1, 2, 4, 6, 12, 15})
        check_extract_dds(bc_gen, res, t0, chn);
    {
        auto [ts, vals] = extract(bc_gen, ChnType::DAC, 3);
        REQUIRE(ts == std::vector<uint64_t>{uint64_t(pulse_t(8) - t0)});
        REQUIRE(vals == std::vector<double>{Zynq::dac_from_mu(mu(ChnType::DAC, -2.0))});
    }
    {
        auto [ts, vals] = extract(bc_gen, ChnType::DAC, 0);
        REQUIRE(ts.size() == size_t(std::count_if(res.dac.begin(), res.dac.end(),
                                                  [] (auto &v) { return (v.second >> 16) == 0; })));
        REQUIRE(ts.back() == uint64_t(pulse_t(2) + 20 * cycles_per_us - t0));
        REQUIRE(vals.back() == Zynq::dac_from_mu(mu(ChnType::DAC, 5.0)));
    }
    // TTL at the exact time
    for (int i = 0; i < 30; i++) {
        INFO("TTL pulse " << i);
        // Each channel is turned on and then off.
        auto bit = uint32_t(1) << (8 + (i / 2) % 8);
        auto t = res.ttl_time(0, bit, i % 2 ? 0 : bit, pulse_t(i));
        REQUIRE(t == pulse_t(i));
    }
}

TEST_CASE("dds_late") {
    // The bus is still busy at the time the outputs needs to start.
    // The outputs happen as soon as possible but never before the requested time.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 1 * us + 80 * ps_per_cycle);
    host_seq.set_value<double>(3, 20e6);
    host_seq.set_value<int64_t>(4, 5 * us);
    host_seq.set_value<int64_t>(5, 1 * us + 10 * ps_per_cycle);
    using Zynq::ChnType;
    auto ramp = [] (double in, const void*) { return 1e6 + in * 5; };
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Freq, 1, 1, 0, 1));
    // A value and the start of a ramp on the same bus 80 cycles later,
    // which can't be done on time with the ramp.
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Freq, 2, 2, 2, 3));
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 3, 3, 2, 4, 1, ramp));
    // Two DAC outputs at the same time and another one shortly after
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::DAC, 0, 4, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::DAC, 1, 5, 0, 1));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::DAC, 2, 6, 5, 1));
    auto res = generate_dma(bc_gen, host_seq);
    auto &timing = bc_gen.dds_timing;
    auto t1 = seq_start(res) + cycles_per_us;
    auto t2 = t1 + 80;
    REQUIRE(find_update(res, 1, t1));
    auto updates2 = chn_updates(res, 2);
    auto updates3 = chn_updates(res, 3);
    REQUIRE(updates2.size() == 1);
    // Right after the previous write, with one FUD for both.
    auto t_late = t1 + timing.fud_to_end + timing.set32 +
        (timing.set32_fud - timing.fud_to_end);
    REQUIRE(t_late > t2);
    REQUIRE(updates2[0].t == t_late);
    REQUIRE(updates2[0].freq == mu(ChnType::Freq, 20e6));
    REQUIRE(updates3[0].t == t_late);
    REQUIRE(updates3[0].freq == mu(ChnType::Freq,
                                   ramp(double(t_late - t2) * ps_per_cycle, nullptr)));
    for (size_t i = 0; i + 1 < updates3.size(); i++) {
        REQUIRE(updates3[i].freq == mu(ChnType::Freq,
                                       ramp(double(updates3[i].t - t2) * ps_per_cycle,
                                            nullptr)));
    }

    // DAC outputs one after another
    REQUIRE(res.dac.size() == 3);
    REQUIRE(res.dac[0].first == t1);
    REQUIRE(res.dac[1].first == t1 + 76 + 1);
    REQUIRE(res.dac[2].first == t1 + 2 * (76 + 1));
    std::set<uint32_t> dac_chns;
    for (auto [t, data]: res.dac) {
        REQUIRE((data & 0xffff) == mu(ChnType::DAC, 10e6));
        dac_chns.insert(data >> 16);
    }
    REQUIRE(dac_chns == std::set<uint32_t>{0, 1, 2});
}

TEST_CASE("dds_ramp_truncated") {
    // Ramps interrupted by another pulse on the same channel.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 30000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<int64_t>(1, 5 * us);
    host_seq.set_value<int64_t>(2, 20 * us);
    host_seq.set_value<int64_t>(3, 2 * us);
    host_seq.set_value<double>(4, 100e6);
    host_seq.set_value<double>(5, 0.3);
    using Zynq::ChnType;
    auto ramp1 = [] (double in, const void*) { return in * 5; };
    auto ramp2 = [] (double in, const void*) { return 50e6 - in * 5; };
    auto ramp3 = [] (double in, const void*) { return in * 1e-7; };
    // DDS 1: ramp cut by a value
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 1, 1, 0, 2, 4, ramp1));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Freq, 1, 2, 1, 4));
    // DDS 2: ramp cut by another ramp
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 2, 3, 0, 2, 4, ramp1));
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 2, 4, 1, 3, 4, ramp2));
    // DDS 13: ramp cut by a value, alone on bus 1
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Amp, 13, 5, 0, 2, 5, ramp3));
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Amp, 13, 6, 1, 5));
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    auto start = t0 + cycles_per_us;
    auto cut = t0 + 5 * cycles_per_us;
    auto check_ramp = [&] (int chn, int64_t from, int64_t to, auto func, bool amp) {
        for (auto &update: chn_updates(res, chn)) {
            if (update.t < from || update.t >= to)
                continue;
            auto v = func(double(update.t - from) * ps_per_cycle, nullptr);
            if (amp) {
                REQUIRE(update.amp == mu(ChnType::Amp, v));
            }
            else {
                REQUIRE(update.freq == mu(ChnType::Freq, v));
            }
        }
    };
    auto updates1 = chn_updates(res, 1);
    REQUIRE(updates1.front().t == start);
    REQUIRE(updates1.back().t == cut);
    REQUIRE(updates1.back().freq == mu(ChnType::Freq, 100e6));
    check_ramp(1, start, cut, ramp1, false);
    // No update for the end of the first ramp of DDS 2
    REQUIRE(find_update(res, 2, cut));
    REQUIRE(find_update(res, 2, cut)->freq == mu(ChnType::Freq, 50e6));
    check_ramp(2, start, cut, ramp1, false);
    check_ramp(2, cut, cut + 2 * cycles_per_us, ramp2, false);
    REQUIRE(chn_updates(res, 2).back().t == cut + 2 * cycles_per_us);
    REQUIRE(chn_updates(res, 2).back().freq == mu(ChnType::Freq, 100e6));
    auto updates13 = chn_updates(res, 13);
    REQUIRE(updates13.size() > 10);
    REQUIRE(updates13.back().t == cut);
    REQUIRE(updates13.back().amp == mu(ChnType::Amp, 0.3));
    check_ramp(13, start, cut, ramp3, true);
}

TEST_CASE("dds_ramp_idle_start") {
    // A ramp starting from the current value long after the last output.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.len_ns = 2000000000;
    host_seq.set_value<int64_t>(0, 1 * us);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 1000000 * us);
    host_seq.set_value<int64_t>(3, 2 * us);
    using Zynq::ChnType;
    auto ramp = [] (double in, const void*) { return 10e6 + in * 5; };
    bc_gen.seq_pulses.push_back(value_pulse(ChnType::Freq, 1, 1, 0, 1));
    bc_gen.seq_pulses.push_back(ramp_pulse(ChnType::Freq, 1, 2, 2, 3, 1, ramp));
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res);
    auto start = t0 + 1000000 * cycles_per_us;
    auto updates = chn_updates(res, 1);
    REQUIRE(updates.size() > 3);
    REQUIRE(updates[0].t == t0 + cycles_per_us);
    // No write needed at the start of the ramp, the next sample is done
    // after one retry (no change at the start).
    auto &timing = bc_gen.dds_timing;
    REQUIRE(updates[1].t > start);
    REQUIRE(updates[1].t <= start + timing.set16_fud +
            (timing.set32_fud - timing.fud_to_end));
    for (size_t i = 1; i + 1 < updates.size(); i++) {
        REQUIRE(updates[i].freq == mu(ChnType::Freq,
                                      ramp(double(updates[i].t - start) * ps_per_cycle,
                                           nullptr)));
    }
    REQUIRE(updates.back().t == start + 2 * cycles_per_us);
    // Few waits (i.e. no busy looping in the generator)
    REQUIRE(res.waits.size() < 200);
}

namespace {

// TTL instruction widths in the group at `t`
static std::multiset<int> ttl_insts_at(const Decoded &res, int64_t t)
{
    std::multiset<int> insts;
    for (auto [inst_t, width]: res.ttl_insts) {
        if (inst_t == t) {
            insts.insert(width);
        }
    }
    return insts;
}

}

TEST_CASE("ttl_cross_bank") {
    // TTL changes in multiple banks at the same time share the same instructions.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    init_bc_gen(bc_gen);
    bc_gen.ttl_mask = {0xffffffff, 0xffffffff, 0xffffffff, 0x7fffffff};
    // Start trigger on a channel not used otherwise.
    bc_gen.start_ttl_chn = 127;
    // All the banks are initialized (so that the later changes only need
    // to write the changed sub-banks).
    bc_gen.first_bseq = true;
    for (int bank = 0; bank < 4; bank++) {
        for (int i = 0; i < (bank == 3 ? 31 : 32); i++) {
            bc_gen.start_vals[{Zynq::ChnType::TTL, uint8_t(bank * 32 + i)}] = 0;
        }
    }
    host_seq.set_value<bool>(0, true);
    for (int i = 0; i < 9; i++)
        host_seq.set_value<int64_t>(1 + i, (i + 1) * us);
    uint32_t id = 1;
    auto ttl = [&] (int chn, uint32_t time_id) {
        bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, uint8_t(chn), id++,
                                                time_id, 0));
    };
    // 1us: one 8-bit sub-bank in each of the 4 banks -> 2 x TTLSet16
    for (int bank = 0; bank < 4; bank++) {
        ttl(bank * 32 + 1, 1);
        ttl(bank * 32 + 6, 1);
    }
    // 2us: one 4-bit sub-bank in two banks -> 1 x TTLSet16
    ttl(8, 2);
    ttl(32 + 24, 2);
    // 3us: one 4-bit sub-bank in three banks -> 1 x TTLSet16 + 1 x TTLSet4
    ttl(12, 3);
    ttl(32 + 12, 3);
    ttl(64 + 31, 3);
    // 4us: three 16-bit sub-banks (both 8-bit sub-banks changed in each)
    // -> 1 x TTLSet32 + 1 x TTLSet16 (smaller than 2 x TTLSet32)
    ttl(16, 4);
    ttl(31, 4);
    ttl(32 + 16, 4);
    ttl(32 + 31, 4);
    ttl(64 + 17, 4);
    ttl(64 + 30, 4);
    // 5us: 16-bit sub-bank in two banks and a 4-bit sub-bank in another one
    // -> 1 x TTLSet32 + 1 x TTLSet4
    ttl(0, 5);
    ttl(15, 5);
    ttl(64 + 0, 5);
    ttl(64 + 9, 5);
    ttl(96 + 4, 5);
    // 6us: 8 bits in one bank and 4 bits in another one -> 1 x TTLSet16
    ttl(96 + 16, 6);
    ttl(96 + 20, 6);
    ttl(0 + 28, 6);
    // 7us: 4 bits in one 16-bit sub-bank and 4 bits in each of the 8-bit sub-banks
    // of another one -> 1 x TTLSet32 (fewer instructions than TTLSet4 + TTLSet16)
    ttl(32 + 0, 7);
    ttl(64 + 20, 7);
    ttl(64 + 24, 7);
    // 8us: 4 bits in two 16-bit sub-banks and 8 bits in a third one -> 1 x TTLSet16
    // (the 8-bit sub-bank and one of the 4-bit sub-banks) + 1 x TTLSet4
    ttl(32 + 4, 8);
    ttl(32 + 20, 8);
    ttl(96 + 8, 8);
    ttl(96 + 12, 8);
    // 9us: four 16-bit sub-banks (both 8-bit sub-banks changed in each)
    // -> 2 x TTLSet32
    for (int bank = 0; bank < 3; bank++) {
        ttl(bank * 32 + 2, 9);
        ttl(bank * 32 + 13, 9);
    }
    ttl(96 + 2, 9);
    ttl(96 + 14, 9);
    auto res = generate_dma(bc_gen, host_seq);
    auto t0 = seq_start(res, 127);
    auto at = [&] (int i) { return ttl_insts_at(res, t0 + i * cycles_per_us); };
    REQUIRE(at(1) == std::multiset<int>{16, 16});
    REQUIRE(at(2) == std::multiset<int>{16});
    REQUIRE(at(3) == std::multiset<int>{4, 16});
    REQUIRE(at(4) == std::multiset<int>{16, 32});
    REQUIRE(at(5) == std::multiset<int>{4, 32});
    REQUIRE(at(6) == std::multiset<int>{16});
    REQUIRE(at(7) == std::multiset<int>{32});
    REQUIRE(at(8) == std::multiset<int>{4, 16});
    REQUIRE(at(9) == std::multiset<int>{32, 32});

    // The values are correct
    std::map<int,uint32_t> bank_vals;
    for (auto &change: res.ttl) {
        if (change.t > t0) {
            bank_vals[change.bank] = change.val;
        }
    }
    auto bits = [] (std::initializer_list<int> chns) {
        uint32_t v = 0;
        for (auto chn: chns)
            v |= uint32_t(1) << chn;
        return v;
    };
    REQUIRE(bank_vals[0] == bits({1, 6, 8, 12, 16, 31, 0, 15, 28, 2, 13}));
    REQUIRE(bank_vals[1] == bits({1, 6, 24, 12, 16, 31, 0, 4, 20, 2, 13}));
    REQUIRE(bank_vals[2] == bits({1, 6, 31, 17, 30, 0, 9, 20, 24, 2, 13}));
    REQUIRE(bank_vals[3] == bits({1, 6, 4, 16, 20, 8, 12, 2, 14}));

    // Extracted values of TTL 64 + 31 (set at 3us). The other times are when
    // the sub-bank of the channel is written.
    auto [ts, vals] = extract(bc_gen, Zynq::ChnType::TTL, 64 + 31);
    REQUIRE(!ts.empty());
    for (size_t i = 0; i < ts.size(); i++)
        REQUIRE(vals[i] == (ts[i] >= uint64_t(3 * cycles_per_us) ? 1 : 0));
    REQUIRE(std::count(ts.begin(), ts.end(), uint64_t(3 * cycles_per_us)) == 1);
}

TEST_CASE("ttl_random") {
    // Random TTL changes on multiple banks, checked against the expected values
    // and the minimum size of the instructions.
    std::mt19937 rng(1234);
    for (int iter = 0; iter < 20; iter++) {
        INFO("iteration " << iter);
        DummyHostSeq host_seq;
        Zynq::BCGen bc_gen;
        init_bc_gen(bc_gen);
        int nbanks = 1 + int(rng() % 8);
        bc_gen.ttl_mask.resize(nbanks);
        for (auto &mask: bc_gen.ttl_mask)
            mask = rng() % 4 == 0 ? 0 : uint32_t(rng() | rng());
        bc_gen.ttl_mask[0] |= 1;
        bc_gen.start_ttl_chn = 0;
        host_seq.set_value<bool>(0, false);
        host_seq.set_value<bool>(1, true);
        constexpr int ntimes = 30;
        for (int i = 0; i < ntimes; i++)
            host_seq.set_value<int64_t>(2 + i, (i + 1) * us);
        std::vector<std::map<int,bool>> expected(ntimes);
        uint32_t id = 1;
        for (int i = 0; i < ntimes; i++) {
            auto nchanges = rng() % 20;
            for (uint32_t j = 0; j < nchanges; j++) {
                auto chn = 1 + int(rng() % (nbanks * 32 - 1));
                bool val = rng() % 2;
                expected[i][chn] = val;
                bc_gen.seq_pulses.push_back(value_pulse(Zynq::ChnType::TTL, uint8_t(chn),
                                                        id++, uint32_t(2 + i), val));
            }
        }
        auto res = generate_dma(bc_gen, host_seq);
        auto t0 = seq_start(res);
        // Value of the masked channels after each time.
        std::map<int,bool> state;
        for (int i = 0; i < ntimes; i++) {
            auto t = t0 + (i + 1) * cycles_per_us;
            for (auto [chn, val]: expected[i])
                state[chn] = val;
            for (auto [chn, val]: state) {
                auto bank = chn / 32;
                auto bit = uint32_t(1) << (chn % 32);
                if (!(res.masks[bank] & bit))
                    continue;
                uint32_t bank_val = 0;
                for (auto &change: res.ttl) {
                    if (change.t <= t && change.bank == bank) {
                        bank_val = change.val;
                    }
                }
                INFO("time " << i << " chn " << chn);
                REQUIRE(bool(bank_val & bit) == val);
            }
        }
    }
}
