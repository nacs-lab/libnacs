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

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1
#define CATCH_CONFIG_MAIN

#include "../error_helper.h"

#include "../../lib/nacs-seq/builder.h"
#include "../../lib/nacs-seq/error.h"
#include "../../lib/nacs-seq/seq.h"

#include "../../lib/nacs-utils/llvm/codegen.h"
#include "../../lib/nacs-utils/llvm/utils.h"
#include "../../lib/nacs-utils/log.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

using namespace NaCs;
using namespace std::literals::string_literals;

static auto llvm_ctx = LLVM::new_context();

TEST_CASE("dup_channel_name") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.chnnames.push_back("test_chn1");

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Duplicated channel name"s);
}

TEST_CASE("default_invalid_chn") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[3] = IR::TagVal(2.3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("default_invalid_chn0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[0] = IR::TagVal(2.3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("assign_invalid_slot") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.slots.push_back(IR::Type::Float64);
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.assignments.push_back({41, 3, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid global ID"s);
}

TEST_CASE("assign_invalid_val") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.slots.push_back(IR::Type::Float64);
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.assignments.push_back({41, 0, 2});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("assign_invalid_val0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.slots.push_back(IR::Type::Float64);
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.assignments.push_back({41, 0, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("assign_valarg") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Arg;
    builder.nodes[0].args[0].id = 0;
    builder.slots.push_back(IR::Type::Float64);
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.assignments.push_back({41, 0, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid assignment value"s);
}

TEST_CASE("node_loop") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[0].args[0].id = 2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[1].args[0].id = 1;

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Dependency loop between values"s);
}

TEST_CASE("node_invalid_slot") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[1].args[0].id = 2;
    builder.slots.push_back(IR::Type::Float64);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid global ID"s);
}

TEST_CASE("node_func_invalid_slot") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.nodes[1].op = Seq::Builder::OpCode::Mul;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[1].args[0].id = 2;
    builder.nodes[1].argtypes[1] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[1].f64 = 2.3;

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid global ID"s);
}

TEST_CASE("node_disabled_measure") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Measure;
    builder.nodes[0].args[0].id = 1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 101, 2, 0});
    basicseq.measures[2] = {1 /* time */, 1 /* channel */};

    uint32_t warning_count = 0;
    Log::pushLogger([&] (Log::Level level, const char *func, const char *msg) {
        REQUIRE(warning_count++ == 0);
        REQUIRE(level == Log::Warn);
        REQUIRE(!func);
        REQUIRE(msg == "Measurement 1 taken on disabled channel.\n"s);
    });
    builder.buildseq();
    Log::popLogger();
}

TEST_CASE("node_interp_invalid_data") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Interp;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[0].id = 0;
    builder.nodes[0].argtypes[1] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[1].id = 1;
    builder.nodes[0].argtypes[2] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[2].id = 2;
    builder.nodes[0].data_id = 1;
    builder.slots.push_back(IR::Type::Bool);
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);
    builder.datas.push_back({1, 2, 5, 6});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid interp data ID"s);
}

TEST_CASE("time_loop") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 1, 1, 2});
    basicseq.times.push_back({Seq::Sign::Pos, 2, 2, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Dependency loop between times"s);
}

TEST_CASE("time_invalid_prev") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 1, 1, 3});
    basicseq.times.push_back({Seq::Sign::Pos, 2, 2, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("time_invalid_val") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 1, 2, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("time_invalid_val0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 1, 0, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("time_prev_arg") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Arg;
    builder.nodes[0].args[0].id = 0;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 1, 1, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid time node"s);
}

TEST_CASE("endtime_invalid_time") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 220, 2, 0});
    basicseq.endtimes.push_back(2);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("endtime_invalid_time0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 220, 2, 0});
    basicseq.endtimes.push_back(0);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("measure_dup_id") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 2000;
    builder.seqs.resize(2);
    auto &basicseq1 = builder.seqs[0];
    auto &basicseq2 = builder.seqs[1];
    basicseq1.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq1.measures[2] = {1 /* time */, 1 /* channel */};
    basicseq2.times.push_back({Seq::Sign::Pos, 102, 2, 0});
    basicseq2.measures[2] = {1 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Duplicated measure ID"s);
}

TEST_CASE("measure_invalid_chn") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.measures[2] = {1 /* time */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("measure_invalid_chn0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.measures[2] = {1 /* time */, 0 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("measure_invalid_time") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.measures[2] = {2 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("measure_invalid_time0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.measures[2] = {0 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("output_invalid_chn") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 1 /* value */,
        0 /* cond */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("output_invalid_chn0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 1 /* value */,
        0 /* cond */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("output_invalid_len") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 2 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("output_invalid_val") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 1 /* length */, 2 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("output_invalid_val0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 1 /* length */, 0 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("output_invalid_ramparg") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(3);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 9123.45;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[2].args[0].f64 = 5;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 2;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 3 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Too many arguments for ramp"s);
}

TEST_CASE("output_invalid_ramparg_indirect") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(4);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 9123.45;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[2].args[0].f64 = 5;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 2;
    builder.nodes[3].op = Seq::Builder::OpCode::Mul;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 3;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[3].args[1].id = 1;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 4 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Too many arguments for ramp"s);
}

TEST_CASE("output_invalid_ramparg_len0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(3);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 9123.45;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[2].args[0].f64 = 5;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 0;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 0 /* length */, 3 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid output value"s);
}

TEST_CASE("output_noramp_noold") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(3);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 9123.45;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[0].id = 2;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 0;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 3 /* value */,
        0 /* cond */, 1 /* channel */};
    builder.noramp_chns.push_back(1);

    auto err = expect_error<Seq::Error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NoRamp));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 8);
    REQUIRE(err.type2 == Seq::Error::Type::Channel);
    REQUIRE(err.id2 == 1);
    REQUIRE(err.what() == "Ramp not supported on channel"s);
}

TEST_CASE("output_noramp_old") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.nodes.resize(4);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 9123.45;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[0].id = 2;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 0;
    builder.nodes[3].op = Seq::Builder::OpCode::Add;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 3;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[3].args[1].id = 1;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[9] = {1 /* time */, 2 /* length */, 4 /* value */,
        0 /* cond */, 2 /* channel */};
    builder.noramp_chns.push_back(2);

    auto err = expect_error<Seq::Error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NoRamp));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 9);
    REQUIRE(err.type2 == Seq::Error::Type::Channel);
    REQUIRE(err.id2 == 2);
    REQUIRE(err.what() == "Ramp not supported on channel"s);
}

TEST_CASE("output_invalid_time") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {2 /* time */, 1 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("output_invalid_time0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {0 /* time */, 1 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound time ID"s);
}

TEST_CASE("output_lenarg") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.nodes[1].op = Seq::Builder::OpCode::Add;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.nodes[1].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[1].args[1].id = 1;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 2 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid pulse length value"s);
}

TEST_CASE("branch_invalid") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.resize(2);
    auto &basicseq1 = builder.seqs[0];
    auto &basicseq2 = builder.seqs[1];
    basicseq1.default_target = 2;
    basicseq1.branches.push_back({102 /* branch ID */, 3 /* target */, 1 /* condition */});
    basicseq2.default_target = 0;
    basicseq2.branches.push_back({103 /* branch ID */, 2 /* target */, 2 /* condition */});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid branch target"s);
}

TEST_CASE("branch_invalid_default") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1000;
    builder.seqs.resize(2);
    auto &basicseq1 = builder.seqs[0];
    auto &basicseq2 = builder.seqs[1];
    basicseq1.default_target = 3;
    basicseq1.branches.push_back({102 /* branch ID */, 1 /* target */, 1 /* condition */});
    basicseq2.default_target = 0;
    basicseq2.branches.push_back({103 /* branch ID */, 2 /* target */, 2 /* condition */});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid branch target"s);
}

TEST_CASE("branch_invalid_cond") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.2;
    builder.seqs.resize(2);
    auto &basicseq1 = builder.seqs[0];
    basicseq1.branches.push_back({102 /* branch ID */, 1 /* target */, 2 /* condition */});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Out of bound node ID"s);
}

TEST_CASE("branch_condarg") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Arg;
    builder.nodes[0].args[0].id = 1;
    builder.seqs.resize(2);
    auto &basicseq1 = builder.seqs[0];
    basicseq1.branches.push_back({102 /* branch ID */, 1 /* target */, 1 /* condition */});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid branch condition"s);
}

TEST_CASE("noramp_invalid_chn") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.noramp_chns.push_back(3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}

TEST_CASE("noramp_invalid_chn0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.noramp_chns.push_back(0);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    REQUIRE(err.what() == "Invalid channel ID"s);
}
