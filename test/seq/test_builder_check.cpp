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

#include "../error_helper.h"

#include "../../lib/seq/builder.h"
#include "../../lib/seq/error.h"

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

static void test_dup_channel_name(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.chnnames.push_back("test_chn1");

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Duplicated channel name") == 0);
}

static void test_default_invalid_chn(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[3] = IR::TagVal(2.3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_default_invalid_chn0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[0] = IR::TagVal(2.3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_assign_invalid_slot(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid global ID") == 0);
}

static void test_assign_invalid_val(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_assign_invalid_val0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_assign_valarg(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid assignment value") == 0);
}

static void test_node_loop(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Dependency loop between values") == 0);
}

static void test_node_invalid_slot(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid global ID") == 0);
}

static void test_node_invalid_measure(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 2, 0});
    basicseq.measures[2] = {1 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid measure ID") == 0);
}

static void test_node_interp_invalid_data(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid interp data ID") == 0);
}

static void test_time_loop(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 1, 1, 2});
    basicseq.times.push_back({Seq::EventTime::Pos, 2, 2, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Dependency loop between times") == 0);
}

static void test_time_invalid_prev(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 1, 1, 3});
    basicseq.times.push_back({Seq::EventTime::Pos, 2, 2, 1});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_time_invalid_val(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 1, 2, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_time_invalid_val0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 4.5;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 1, 0, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_time_prev_arg(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Arg;
    builder.nodes[0].args[0].id = 0;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 1, 1, 0});

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid time node") == 0);
}

static void test_endtime_invalid_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 220, 2, 0});
    basicseq.endtimes.push_back(2);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_endtime_invalid_time0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 220, 2, 0});
    basicseq.endtimes.push_back(0);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_measure_dup_id(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq1.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq1.measures[2] = {1 /* time */, 1 /* channel */};
    basicseq2.times.push_back({Seq::EventTime::Pos, 102, 2, 0});
    basicseq2.measures[2] = {1 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Duplicated measure ID") == 0);
}

static void test_measure_invalid_chn(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.measures[2] = {1 /* time */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_measure_invalid_chn0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.measures[2] = {1 /* time */, 0 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_measure_invalid_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.measures[2] = {2 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_measure_invalid_time0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.measures[2] = {0 /* time */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_output_invalid_chn(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 1 /* value */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_output_invalid_chn0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 1 /* value */, 2 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_output_invalid_len(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 2 /* length */, 1 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_output_invalid_val(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 1 /* length */, 2 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_output_invalid_val0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 1 /* length */, 0 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_output_invalid_ramparg(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 3 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Too many arguments for ramp") == 0);
}

static void test_output_invalid_ramparg_indirect(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 4 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Too many arguments for ramp") == 0);
}

static void test_output_invalid_ramparg_len0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 0 /* length */, 3 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid output value") == 0);
}

static void test_output_noramp_noold(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 3 /* value */, 1 /* channel */};
    builder.noramp_chns.push_back(1);

    auto err = expect_error<Seq::Error>([&] {
        builder.buildseq();
    });

    assert(err.code >> 16 == Seq::Error::Pulse);
    assert((err.code & ((1u << 16) - 1)) == Seq::Pulse::NoRamp);
    assert(err.type1 == Seq::Error::Pulse);
    assert(err.id1 == 8);
    assert(err.type2 == Seq::Error::Channel);
    assert(err.id2 == 1);
    assert(strcmp(err.what(), "Ramp not supported on channel") == 0);
}

static void test_output_noramp_old(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[9] = {1 /* time */, 2 /* length */, 4 /* value */, 2 /* channel */};
    builder.noramp_chns.push_back(2);

    auto err = expect_error<Seq::Error>([&] {
        builder.buildseq();
    });

    assert(err.code >> 16 == Seq::Error::Pulse);
    assert((err.code & ((1u << 16) - 1)) == Seq::Pulse::NoRamp);
    assert(err.type1 == Seq::Error::Pulse);
    assert(err.id1 == 9);
    assert(err.type2 == Seq::Error::Channel);
    assert(err.id2 == 2);
    assert(strcmp(err.what(), "Ramp not supported on channel") == 0);
}

static void test_output_invalid_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {2 /* time */, 1 /* length */, 1 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_output_invalid_time0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.nodes.resize(1);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 1000;
    builder.seqs.resize(1);
    auto &basicseq = builder.seqs[0];
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {0 /* time */, 1 /* length */, 1 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Out of bound time ID") == 0);
}

static void test_output_lenarg(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 0});
    basicseq.outputs[29] = {1 /* time */, 2 /* length */, 1 /* value */, 1 /* channel */};

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid pulse length value") == 0);
}

static void test_branch_invalid(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid branch target") == 0);
}

static void test_branch_invalid_default(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid branch target") == 0);
}

static void test_branch_invalid_cond(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Out of bound node ID") == 0);
}

static void test_branch_condarg(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    assert(strcmp(err.what(), "Invalid branch condition") == 0);
}

static void test_noramp_invalid_chn(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.noramp_chns.push_back(3);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

static void test_noramp_invalid_chn0(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.noramp_chns.push_back(0);

    auto err = expect_error<std::runtime_error>([&] {
        builder.buildseq();
    });

    assert(strcmp(err.what(), "Invalid channel ID") == 0);
}

int main()
{
    auto llvm_ctx = LLVM::new_context();

    test_dup_channel_name(*llvm_ctx);
    test_default_invalid_chn(*llvm_ctx);
    test_default_invalid_chn0(*llvm_ctx);
    test_assign_invalid_slot(*llvm_ctx);
    test_assign_invalid_val(*llvm_ctx);
    test_assign_invalid_val0(*llvm_ctx);
    test_assign_valarg(*llvm_ctx);
    test_node_loop(*llvm_ctx);
    test_node_invalid_slot(*llvm_ctx);
    test_node_invalid_measure(*llvm_ctx);
    test_node_interp_invalid_data(*llvm_ctx);
    test_time_loop(*llvm_ctx);
    test_time_invalid_prev(*llvm_ctx);
    test_time_invalid_val(*llvm_ctx);
    test_time_invalid_val0(*llvm_ctx);
    test_time_prev_arg(*llvm_ctx);
    test_endtime_invalid_time(*llvm_ctx);
    test_endtime_invalid_time0(*llvm_ctx);
    test_measure_dup_id(*llvm_ctx);
    test_measure_invalid_chn(*llvm_ctx);
    test_measure_invalid_chn0(*llvm_ctx);
    test_measure_invalid_time(*llvm_ctx);
    test_measure_invalid_time0(*llvm_ctx);
    test_output_invalid_chn(*llvm_ctx);
    test_output_invalid_chn0(*llvm_ctx);
    test_output_invalid_len(*llvm_ctx);
    test_output_invalid_val(*llvm_ctx);
    test_output_invalid_val0(*llvm_ctx);
    test_output_invalid_ramparg(*llvm_ctx);
    test_output_invalid_ramparg_indirect(*llvm_ctx);
    test_output_invalid_ramparg_len0(*llvm_ctx);
    test_output_noramp_noold(*llvm_ctx);
    test_output_noramp_old(*llvm_ctx);
    test_output_invalid_time(*llvm_ctx);
    test_output_invalid_time0(*llvm_ctx);
    test_output_lenarg(*llvm_ctx);
    test_branch_invalid(*llvm_ctx);
    test_branch_invalid_default(*llvm_ctx);
    test_branch_invalid_cond(*llvm_ctx);
    test_branch_condarg(*llvm_ctx);
    test_noramp_invalid_chn(*llvm_ctx);
    test_noramp_invalid_chn0(*llvm_ctx);

    return 0;
}
