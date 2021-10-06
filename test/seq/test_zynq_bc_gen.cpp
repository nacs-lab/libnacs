/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

        // Create one sequence since a real sequence always has one
        // even though we don't really need it...
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
    void set_value(uint32_t idx, T v);
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

// The BCGen only uses the value and type array
// and that'll be the only thing we'll fill in the HostSeq.
template<>
void DummyHostSeq::set_value<bool>(uint32_t idx, bool v)
{
    reserve(idx + 1);
    types[idx] = HostSeq::Type::Bool;
    values[idx].b = v;
}

template<>
void DummyHostSeq::set_value<int32_t>(uint32_t idx, int32_t v)
{
    reserve(idx + 1);
    types[idx] = HostSeq::Type::Int32;
    values[idx].i32 = v;
}

template<>
void DummyHostSeq::set_value<double>(uint32_t idx, double v)
{
    reserve(idx + 1);
    types[idx] = HostSeq::Type::Float64;
    values[idx].f64 = v;
}

template<>
void DummyHostSeq::set_value<int64_t>(uint32_t idx, int64_t v)
{
    reserve(idx + 1);
    types[idx] = HostSeq::Type::Int64;
    values[idx].i64 = v;
}

TEST_CASE("static") {
    REQUIRE(Zynq::BCGen::version() == 2);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::TTL, 0) == 0);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::TTL, 1) == 1);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Freq, 0) == 0);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Freq, 1.75e9 / 2) == (1u << 30u));
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Freq, 1.75e9) == (1u << 31u) - 1);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Amp, 0) == 0);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Amp, 0.29) == 1188);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::Amp, 1.0) == 4095);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::DAC, -10) == 0);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::DAC, 0.1) == 33095);
    REQUIRE(Zynq::BCGen::convert_value(Zynq::BCGen::ChnType::DAC, 10) == 65535);
}

static void check_generate(Zynq::BCGen &bc_gen, HostSeq &host_seq)
{
    bc_gen.generate(host_seq);
    auto bc = bc_gen.bytecode;
    REQUIRE(&bc != &bc_gen.bytecode);
    INFO(dump_bytecode(bc));
    for (int i = 0; i < 10; i++) {
        bc_gen.generate(host_seq);
        INFO(dump_bytecode(bc_gen.bytecode));
        REQUIRE(bc == bc_gen.bytecode);
    }
}

TEST_CASE("empty") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 0;
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(0)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("offset") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 100;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 0;
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(0)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 197}));
    REQUIRE(checker.is_end());
}

TEST_CASE("clock") {
    for (int i = 0; i < 2; i++) {
        DummyHostSeq host_seq;
        Zynq::BCGen bc_gen;
        bc_gen.seq_idx = 0;
        // The clock offset will swallow a small sequence delay
        bc_gen.seq_delay = i * 40;
        bc_gen.first_bseq = false;
        bc_gen.start_ttl_chn = 0;
        bc_gen.ttl_mask = 0;
        bc_gen.fpga_clock_div = 10000;
        bc_gen.len_ns = 300;
        bc_gen.clocks.push_back({-100 * 10000, 100});
        bc_gen.clocks.push_back({1000 * 10000, 0});
        check_generate(bc_gen, host_seq);
        Checker checker(bc_gen.bytecode);
        INFO(dump_bytecode(bc_gen.bytecode));
        REQUIRE(checker.cmp<int64_t>(300)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 1095}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
}

TEST_CASE("clock+offset") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 3000;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 300;
    bc_gen.clocks.push_back({-100 * 10000, 100});
    bc_gen.clocks.push_back({1000 * 10000, 0});
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(300)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2900}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 1095}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_set") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x208;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 0;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100800);
    host_seq.set_value<bool>(3, false);
    host_seq.set_value<int64_t>(4, 100700);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 9,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x209)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 1, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 5}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 9}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_set_off") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x208;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 0;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100800);
    host_seq.set_value<bool>(3, false);
    host_seq.set_value<int64_t>(4, 100700);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = 3,
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 9,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x209)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 1, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 5}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 9, 9}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl0") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x4;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<bool>(1, true);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x5)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 2, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl0_first_bseq") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x4;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<bool>(1, true);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x5)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl0_first_bseq_off") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x4;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<bool>(2, false);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = 2,
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x5)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl0+_first_bseq") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x4;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    // Automatic update of start value should only happen if the sequence time is 0
    // not always when the FPGA time step is 0.
    host_seq.set_value<int64_t>(0, 10);
    host_seq.set_value<bool>(1, true);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x5)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 2, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_set+offset") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 40;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x208;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100800);
    host_seq.set_value<bool>(3, false);
    host_seq.set_value<int64_t>(4, 100700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 9,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x209)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 38}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 5}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 9}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_set+clock") {
    for (int i = 0; i < 2; i++) {
        DummyHostSeq host_seq;
        Zynq::BCGen bc_gen;
        bc_gen.seq_idx = 0;
        // The clock offset will swallow a small sequence delay
        bc_gen.seq_delay = i * 40;
        bc_gen.first_bseq = false;
        bc_gen.start_ttl_chn = 0;
        bc_gen.ttl_mask = 0x208;
        bc_gen.fpga_clock_div = 10000;
        bc_gen.len_ns = 987;
        bc_gen.clocks.push_back({-100 * 10000, 100});
        bc_gen.clocks.push_back({1000 * 10000, 0});
        host_seq.set_value<int64_t>(0, 10300);
        host_seq.set_value<bool>(1, true);
        host_seq.set_value<int64_t>(2, 100800);
        host_seq.set_value<bool>(3, false);
        host_seq.set_value<int64_t>(4, 100700);
        bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
        bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 0;
        bc_gen.seq_pulses.push_back({
                .chn_type = Zynq::BCGen::ChnType::TTL,
                .pulse_type = Zynq::BCGen::PulseType::Value,
                .chn = 3,

                .id = 1,
                .time_id = 0,
                .len_id = uint32_t(-1),
                .endvalue_id = 1,
                .cond_id = uint32_t(-1),
                .ramp_func = nullptr
            });
        bc_gen.seq_pulses.push_back({
                .chn_type = Zynq::BCGen::ChnType::TTL,
                .pulse_type = Zynq::BCGen::PulseType::Value,
                .chn = 3,

                .id = 2,
                .time_id = 2,
                .len_id = uint32_t(-1),
                .endvalue_id = 3,
                .cond_id = uint32_t(-1),
                .ramp_func = nullptr
            });
        bc_gen.seq_pulses.push_back({
                .chn_type = Zynq::BCGen::ChnType::TTL,
                .pulse_type = Zynq::BCGen::PulseType::Value,
                .chn = 9,

                .id = 3,
                .time_id = 4,
                .len_id = uint32_t(-1),
                .endvalue_id = 1,
                .cond_id = uint32_t(-1),
                .ramp_func = nullptr
            });
        check_generate(bc_gen, host_seq);
        Checker checker(bc_gen.bytecode);
        INFO(dump_bytecode(bc_gen.bytecode));
        REQUIRE(checker.cmp<int64_t>(987)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x209)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 5}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 3, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 986}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
}

TEST_CASE("dds_set") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<int64_t>(4, 10000700);
    host_seq.set_value<double>(5, 110e6);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 5,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 1, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x7507507}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 449}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 450}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x80bb3ee}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_set_skip") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 100800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<int64_t>(4, 100700);
    host_seq.set_value<double>(5, 110e6);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    // This pulse is only 10 cycles from the next one and should be skipped.
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 5,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x80bb3ee}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_set_skip+offset") {
    // Skipping the first set pulse with a clock with a negative start time.
    // This check if the intial time offset can be handled correctly.
    // This also requires `first_bseq` to be set to have all the channels
    // set to their start value before the sequence start.
    // And a pulse time `> 0` to make sure it didn't get folded into the start sequence.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 1034;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<int64_t>(0, 1);
    host_seq.set_value<double>(1, 0);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence+Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 1131}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_set+offset") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 1034;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<int64_t>(4, 10000700);
    host_seq.set_value<double>(5, 110e6);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 5,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 1032}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x7507507}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 449}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 450}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x80bb3ee}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0xbec17}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp_vector") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Vector,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(void (*)(double*, const double*, const void*))
            [] (double *out, const double *in, const void*) {
                REQUIRE(in[0] >= 0);
                REQUIRE(in[0] <= 12000700);
                for (unsigned i = 0; i < vector_size; i++) {
                    out[i] = in[i] * 0.8149072527885437;
                }
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0xbec17}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp_interrupt") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<bool>(3, true);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = 3,
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = 3,
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0xbb3ee7}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp_interrupt_off") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<bool>(3, false);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = 3,
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0xbec17}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp+offset") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 4875;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 4875}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0xbec17}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp+clock") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    bc_gen.clocks.push_back({-100 * 10000, 100});
    bc_gen.clocks.push_back({500 * 10000, 0});
    bc_gen.clocks.push_back({1000 * 10000, 100});
    bc_gen.clocks.push_back({2500 * 10000, 0});
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 95}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1050000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 45}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0x12c9e7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 1295}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_multiramp") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = 4,
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return 1.0 - in * 2.9304029304029303e-8;
            }
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xfff}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf87}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf0f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe97}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe1f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xda7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xd2f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xcb7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0xbec17}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xc3f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44})); // -60
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_multiramp+clock") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.

    // NOTE: In additional to the ramp, overlapping part of the DDS outputs,
    // should be interrupted by the clock once after each channel
    // (i.e. once after Freq3 and once after Amp7) and the ramp output should restart
    // on the other channel
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    bc_gen.clocks.push_back({600 * 10000, 100});
    bc_gen.clocks.push_back({955 * 10000, 0});
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 12000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = 4,
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 12000700);
                return 1.0 - in * 2.9304029304029303e-8;
            }
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xfff}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf81}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1050000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf09}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe91}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe19}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Clock{5, 0, 255}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1050000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xd9b}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xd23}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0x19a7b7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xcab}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44})); // -60
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0x44}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_slowramp") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 10300);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 5000800);
    host_seq.set_value<double>(3, 0.1);
    host_seq.set_value<double>(4, 15000700);
    host_seq.set_value<double>(5, 35000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 15000700);
                return in * 0.8149072527885437;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = 5,
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 35000700);
                return in * 1.9536019536019536e-10;
            }
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 3, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 490000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 500000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0x65c677})); // - 0x1a3989
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 50}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x19a}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds0") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 300);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(5, 110e6);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 5}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 5,

            .id = 3,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 5,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x7507507}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 5, 0x80bb3ee}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds0_start") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<double>(1, 100e6);
    host_seq.set_value<int64_t>(2, 300);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(5, 110e6);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 5}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 7,

            .id = 2,
            .time_id = 0,
            .len_id = uint32_t(-1),
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 5,

            .id = 3,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 5,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x7507507}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq2{7, 5, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 5, 0x80bb3ee}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp0") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 200);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 8000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 8000700);
                return in * 0.8149072527885437 + 300;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = 4,
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 8000700);
                return 1.0 - in * 2.9304029304029303e-8;
            }
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq3{8, 3, 0x170}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xfc3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf4b}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xed3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe5b}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xde3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xd6b}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xcf3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0xbb3ee7}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("dds_ramp0_start") {
    // Note that the exact lowering of the ramp may depend on the scheduling logic
    // Test failure after changing the scheduling logic does not necessarily indicate a bug.
    // However, the resulting new bytecode should be reviewed manually
    // to make sure it still produce the correct result
    // and the test should be updated accordingly.
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 1230;
    host_seq.set_value<int64_t>(0, 0);
    host_seq.set_value<double>(1, 10e6);
    host_seq.set_value<int64_t>(2, 200);
    host_seq.set_value<double>(3, 0.4);
    host_seq.set_value<double>(4, 8000700);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0;
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Freq,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 3,

            .id = 1,
            .time_id = 0,
            .len_id = 4,
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 8000700);
                return in * 0.8149072527885437 + 300;
            }
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::Amp,
            .pulse_type = Zynq::BCGen::PulseType::Scalar,
            .chn = 7,

            .id = 2,
            .time_id = 2,
            .len_id = 4,
            .endvalue_id = 3,
            .cond_id = uint32_t(-1),
            .ramp_func = (void (*)(void))(double (*)(double, const void*))
            [] (double in, const void*) {
                REQUIRE(in >= 0);
                REQUIRE(in <= 8000700);
                return 1.0 - in * 2.9304029304029303e-8;
            }
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(1230)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq3{8, 3, 0x170}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetAmp{11, 7, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xfff}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 0x7a120}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf87}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xf0f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe97}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xe1f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xda7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xd2f}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSDetFreq4{9, 3, 1000000}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0xcb7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0xbb3ee7}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x666}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 100}));
    REQUIRE(checker.is_end());
}

TEST_CASE("start_vals") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x228;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 5}] = 1;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 1;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0x94ab4e3;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0x348;
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(0)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x229)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL4{2, 0, 5, 9, 9}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("start_vals_start") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x228;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 3}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 5}] = 1;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 9}] = 1;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Freq, 3}] = 0x94ab4e3;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::Amp, 7}] = 0x348;
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(0)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x229)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSFreq{6, 3, 0x94ab4e3}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::DDSAmp{10, 0, 7, 0x348}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL4{2, 0, 5, 9, 9}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 99}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_mgr") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x6;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<bool>(0, false);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100000);
    host_seq.set_value<int64_t>(3, 3000000);
    host_seq.set_value<int64_t>(4, 8000000);
    host_seq.set_value<int64_t>(5, 8100000);
    host_seq.set_value<int64_t>(6, 12000000);
    host_seq.set_value<int64_t>(7, 15000000);
    host_seq.set_value<int64_t>(8, 15500000);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 1}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 1;
    bc_gen.add_ttl_manager(1, 300000, 200000, 1000000, 1300000, false);
    bc_gen.add_ttl_manager(2, 300000, 200000, 1000000, 1300000, true);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 1,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 2,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 4,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 5,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 6,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 7,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });

    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 8,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 9,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 10,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 11,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 12,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 13,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 14,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 266}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 506}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 126}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 266}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_mgr_inverse") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = false;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x6;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<bool>(0, false);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100000);
    host_seq.set_value<int64_t>(3, 3000000);
    host_seq.set_value<int64_t>(4, 8000000);
    host_seq.set_value<int64_t>(5, 8100000);
    host_seq.set_value<int64_t>(6, 12000000);
    host_seq.set_value<int64_t>(7, 15000000);
    host_seq.set_value<int64_t>(8, 15500000);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 1}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 1;
    bc_gen.add_ttl_manager(1, 300000, 200000, 1000000, 1300000, true);
    bc_gen.add_ttl_manager(2, 300000, 200000, 1000000, 1300000, false);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 1,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 2,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 4,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 5,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 6,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 7,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });

    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 8,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 9,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 10,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 11,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 12,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 13,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 14,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 276}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 886}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 306}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 126}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_mgr_start") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x6;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<bool>(0, false);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100000);
    host_seq.set_value<int64_t>(3, 3000000);
    host_seq.set_value<int64_t>(4, 8000000);
    host_seq.set_value<int64_t>(5, 8100000);
    host_seq.set_value<int64_t>(6, 12000000);
    host_seq.set_value<int64_t>(7, 15000000);
    host_seq.set_value<int64_t>(8, 15500000);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 1}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 1;
    bc_gen.add_ttl_manager(1, 300000, 200000, 1000000, 1300000, false);
    bc_gen.add_ttl_manager(2, 300000, 200000, 1000000, 1300000, true);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 1,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 2,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 4,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 5,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 6,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 7,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });

    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 8,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 9,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 10,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 11,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 12,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 13,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 14,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 1}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 267}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 506}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 126}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 266}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_mgr_inverse_start") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x6;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<bool>(0, false);
    host_seq.set_value<bool>(1, true);
    host_seq.set_value<int64_t>(2, 100000);
    host_seq.set_value<int64_t>(3, 3000000);
    host_seq.set_value<int64_t>(4, 8000000);
    host_seq.set_value<int64_t>(5, 8100000);
    host_seq.set_value<int64_t>(6, 12000000);
    host_seq.set_value<int64_t>(7, 15000000);
    host_seq.set_value<int64_t>(8, 15500000);
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 1}] = 0;
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 2}] = 1;
    bc_gen.add_ttl_manager(1, 300000, 200000, 1000000, 1300000, true);
    bc_gen.add_ttl_manager(2, 300000, 200000, 1000000, 1300000, false);
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 1,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 2,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 3,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 4,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 5,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 6,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 1,

            .id = 7,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });

    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 8,
            .time_id = 2,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 9,
            .time_id = 3,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 10,
            .time_id = 4,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 11,
            .time_id = 5,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 12,
            .time_id = 6,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 13,
            .time_id = 7,
            .len_id = uint32_t(-1),
            .endvalue_id = 1,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    bc_gen.seq_pulses.push_back({
            .chn_type = Zynq::BCGen::ChnType::TTL,
            .pulse_type = Zynq::BCGen::PulseType::Value,
            .chn = 2,

            .id = 14,
            .time_id = 8,
            .len_id = uint32_t(-1),
            .endvalue_id = 0,
            .cond_id = uint32_t(-1),
            .ramp_func = nullptr
        });
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 0, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 276}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 886}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 306}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 126}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 2}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}

TEST_CASE("ttl_rate_limit") {
    DummyHostSeq host_seq;
    Zynq::BCGen bc_gen;
    bc_gen.seq_idx = 0;
    bc_gen.seq_delay = 0;
    bc_gen.first_bseq = true;
    bc_gen.start_ttl_chn = 0;
    bc_gen.ttl_mask = 0x6;
    bc_gen.fpga_clock_div = 10000;
    bc_gen.len_ns = 123;
    host_seq.set_value<bool>(0, false);
    host_seq.set_value<bool>(1, true);
    for (int i = 0; i < 100000; i++) {
        host_seq.set_value<int64_t>(i + 2, int64_t(100000) * i);
        bc_gen.seq_pulses.push_back({
                .chn_type = Zynq::BCGen::ChnType::TTL,
                .pulse_type = Zynq::BCGen::PulseType::Value,
                .chn = 1,

                .id = uint32_t(i + 1),
                .time_id = uint32_t(i + 2),
                .len_id = uint32_t(-1),
                .endvalue_id = uint32_t(i % 2), // 0->false, 1->true
                .cond_id = uint32_t(-1),
                .ramp_func = nullptr
            });
    }
    bc_gen.start_vals[{Zynq::BCGen::ChnType::TTL, 1}] = 0;
    check_generate(bc_gen, host_seq);
    Checker checker(bc_gen.bytecode);
    INFO(dump_bytecode(bc_gen.bytecode));
    REQUIRE(checker.cmp<int64_t>(123)); // len_ns
    REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
    // Startup
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait{4, 0, 2500}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 96}));
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 0, 0}));
    // Sequence
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 7}));
    for (int i = 0; i < 1023; i++) {
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 1}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 6}));
    }
    // We should probably make this transition smoother...
    for (int i = 0; i < 385; i++) {
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 1}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 246}));
    }
    for (int i = 0; i < 5956; i++) {
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 1}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 146}));
    }
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::TTL2{1, 3, 1, 1}));
    // Finish
    REQUIRE(checker.cmp(Zynq::ByteCode::Inst::Wait2{5, 1, 97}));
    REQUIRE(checker.is_end());
}
