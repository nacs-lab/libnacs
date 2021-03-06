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

#include "../../lib/seq/builder.h"
#include "../../lib/seq/error.h"

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

static void test_channels(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");

    builder.buildseq();

    assert(seq.get_chn_id("test_chn2") == 2);
    assert(seq.get_chn_id("test_chn1") == 1);
    assert(seq.get_chn_id("test_chn") == 0); // Not existing channel
    assert(seq.get_chn_names().size() == 2);
    assert(seq.get_chn_names()[0] == "test_chn1");
    assert(seq.get_chn_names()[1] == "test_chn2");
}

static void test_default(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[1] = IR::TagVal(2.3);

    builder.buildseq();

    assert(seq.get_chn_id("test_chn1") == 1);
    auto defval = seq.defval(1);
    assert(defval == 2.3);
}

static void test_global(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Bool);
    builder.slots.push_back(IR::Type::Int32);

    builder.buildseq();

    assert(seq.get_slots().size() == 3);
    auto slot1 = seq.get_slots()[0].get();
    auto slot2 = seq.get_slots()[1].get();
    auto slot3 = seq.get_slots()[2].get();
    assert(slot1 == seq.get_slot(IR::Type::Bool, 0));
    assert(slot2 == seq.get_slot(IR::Type::Int32, 1));
    assert(slot3 == seq.get_slot(IR::Type::Int32, 2));
    assert(slot1->type() == IR::Type::Float64);
    assert(slot2->type() == IR::Type::Bool);
    assert(slot3->type() == IR::Type::Int32);
}

static void test_assign(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[0].id = 0;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[1].args[0].id = 1;
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.assignments.push_back({41, 0, 2});
    basicseq.assignments.push_back({52, 1, 1});

    builder.buildseq();

    assert(seq.get_slots().size() == 2);
    auto slot1 = seq.get_slots()[0].get();
    auto slot2 = seq.get_slots()[1].get();
    assert(slot1 == seq.get_slot(IR::Type::Bool, 0));
    assert(slot2 == seq.get_slot(IR::Type::Int32, 1));
    assert(slot1->type() == IR::Type::Float64);
    assert(slot2->type() == IR::Type::Int32);

    assert(seq.get_basicseqs().size() == 1);

    assert(seq.get_basicseqs().front().id() == 1);

    auto &assigns = seq.get_basicseqs().front().get_assigns();
    assert(assigns.size() == 2);
    assert(assigns.find(0)->second.val.get() == slot2);
    assert(assigns.find(0)->second.id == 41);
    assert(assigns.find(1)->second.val.get() == slot1);
    assert(assigns.find(1)->second.id == 52);
}

static void test_node_cse(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(8);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[0].id = 0;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[1].args[0].id = 1;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[0].id = 1;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[1].id = 2;
    builder.nodes[3].op = Seq::Builder::OpCode::Add;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 1;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[1].id = 2;
    builder.nodes[4].op = Seq::Builder::OpCode::Mul;
    builder.nodes[4].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[4].args[0].id = 3;
    builder.nodes[4].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[4].args[1].id = 4;
    builder.nodes[5].op = Seq::Builder::OpCode::Mul;
    builder.nodes[5].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[0].id = 4;
    builder.nodes[5].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[1].id = 3;
    builder.nodes[6].op = Seq::Builder::OpCode::Mul;
    builder.nodes[6].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[6].args[0].id = 5;
    builder.nodes[6].argtypes[1] = Seq::Builder::ArgType::ConstInt32;
    builder.nodes[6].args[1].i32 = 8;
    builder.nodes[7].op = Seq::Builder::OpCode::Mul;
    builder.nodes[7].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[7].args[0].id = 6;
    builder.nodes[7].argtypes[1] = Seq::Builder::ArgType::ConstInt32;
    builder.nodes[7].args[1].i32 = 8;
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);

    builder.buildseq();

    std::vector<Seq::Var*> vars(seq.env().begin(), seq.env().end());

    assert(vars.size() == 5);

    assert(vars[4]->varid() == 0);
    assert(vars[4]->is_extern());
    assert(vars[4]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0)));
    assert(vars[3]->varid() == 1);
    assert(vars[3]->is_extern());
    assert(vars[3]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(1)));
    assert(vars[2]->varid() == 2);
    assert(vars[2]->is_call());
    assert(vars[2]->get_callee().is_llvm);
    assert(vars[2]->nfreeargs() == 0);
    assert(vars[2]->args().size() == 2);
    assert(vars[2]->args()[0].compare(Seq::Arg::create_var(vars[4])) == 0);
    assert(vars[2]->args()[1].compare(Seq::Arg::create_var(vars[3])) == 0);
    assert(vars[1]->varid() == 3);
    assert(vars[1]->is_call());
    assert(vars[1]->get_callee().is_llvm);
    assert(vars[1]->nfreeargs() == 0);
    assert(vars[1]->args().size() == 2);
    assert(vars[1]->args()[0].compare(Seq::Arg::create_var(vars[2])) == 0);
    assert(vars[1]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    assert(vars[0]->varid() == 4);
    assert(vars[0]->is_call());
    assert(vars[0]->get_callee().is_llvm);
    assert(vars[0]->nfreeargs() == 0);
    assert(vars[0]->args().size() == 2);
    assert(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[1])) == 0);
    assert(vars[0]->args()[1].compare(Seq::Arg::create_const(IR::TagVal(8))) == 0);
}

static void test_node_commute(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.nodes.resize(6);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[0].args[0].id = 0;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[1].args[0].id = 1;
    builder.nodes[2].op = Seq::Builder::OpCode::Add;
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[0].id = 1;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[2].args[1].id = 2;
    builder.nodes[3].op = Seq::Builder::OpCode::Add;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 2;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[1].id = 1;
    builder.nodes[4].op = Seq::Builder::OpCode::CmpLT;
    builder.nodes[4].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[4].args[0].id = 1;
    builder.nodes[4].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[4].args[1].id = 2;
    builder.nodes[5].op = Seq::Builder::OpCode::CmpGT;
    builder.nodes[5].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[0].id = 2;
    builder.nodes[5].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[1].id = 1;
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);

    builder.buildseq();

    std::vector<Seq::Var*> vars(seq.env().begin(), seq.env().end());

    assert(vars.size() == 4);

    assert(vars[3]->varid() == 0);
    assert(vars[3]->is_extern());
    assert(vars[3]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0)));
    assert(vars[2]->varid() == 1);
    assert(vars[2]->is_extern());
    assert(vars[2]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(1)));
    assert(vars[1]->varid() == 2);
    assert(vars[1]->is_call());
    assert(vars[1]->get_callee().is_llvm);
    assert(vars[1]->nfreeargs() == 0);
    assert(vars[1]->args().size() == 2);
    assert(vars[1]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    assert(vars[1]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    assert(vars[0]->varid() == 3);
    assert(vars[0]->is_call());
    assert(vars[0]->get_callee().is_llvm);
    assert(vars[0]->nfreeargs() == 0);
    assert(vars[0]->args().size() == 2);
    assert(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    assert(vars[0]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
}

static void test_node_interp(llvm::LLVMContext &llvm_ctx)
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
    builder.nodes[0].data_id = 0;
    builder.slots.push_back(IR::Type::Bool);
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);
    builder.datas.push_back({1, 2, 5, 6});

    builder.buildseq();

    std::vector<Seq::Var*> vars(seq.env().begin(), seq.env().end());

    assert(vars.size() == 4);

    assert(vars[3]->varid() == 0);
    assert(vars[3]->is_extern());
    assert(vars[3]->get_extern() == std::make_pair(IR::Type::Bool, uint64_t(0)));
    assert(vars[2]->varid() == 1);
    assert(vars[2]->is_extern());
    assert(vars[2]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(1)));
    assert(vars[1]->varid() == 2);
    assert(vars[1]->is_extern());
    assert(vars[1]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(2)));
    assert(vars[0]->varid() == 3);
    assert(vars[0]->is_call());
    assert(vars[0]->get_callee().is_llvm);
    assert(vars[0]->nfreeargs() == 0);
    assert(vars[0]->args().size() == 3);
    assert(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    assert(vars[0]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    assert(vars[0]->args()[2].compare(Seq::Arg::create_var(vars[1])) == 0);
}

static void test_time(llvm::LLVMContext &llvm_ctx)
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
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 2});
    basicseq.times.push_back({Seq::EventTime::Pos, 220, 2, 0});
    basicseq.times.push_back({Seq::EventTime::Unknown, 299, 1, 1});
    basicseq.endtimes.push_back(1);
    basicseq.endtimes.push_back(3);

    builder.buildseq();

    assert(seq.get_basicseqs().size() == 1);
    assert(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    std::vector<decltype(&*bseq.get_assumes().begin())> assumes;
    for (auto &assume: bseq.get_assumes())
        assumes.push_back(&assume);

    assert(assumes.size() == 2);
    assert(assumes[0]->second.sign == Seq::EventTime::Pos);
    assert(assumes[0]->second.id == 101);
    assert(assumes[0]->first->is_const());
    assert(assumes[0]->first.get() == seq.get_const(IR::TagVal(4.2)));
    assert(assumes[1]->second.sign == Seq::EventTime::Pos);
    assert(assumes[1]->second.id == 220);
    assert(assumes[1]->first->is_const());
    assert(assumes[1]->first.get() == seq.get_const(IR::TagVal(1000.0)));

    std::vector<const Seq::EventTime*> endtimes;
    for (auto &endtime: bseq.get_endtimes())
        endtimes.push_back(endtime.get());
    assert(endtimes.size() == 2);
    assert(endtimes[0]->tconst == 0);
    assert(endtimes[0]->terms.size() == 2);
    assert(endtimes[0]->terms[0].sign == Seq::EventTime::Pos);
    assert(endtimes[0]->terms[0].id == 220);
    assert(endtimes[0]->terms[0].var->is_const());
    assert(endtimes[0]->terms[0].var.get() == seq.get_const(IR::TagVal(1000.0)));
    assert(endtimes[0]->terms[1].sign == Seq::EventTime::Pos);
    assert(endtimes[0]->terms[1].id == 101);
    assert(endtimes[0]->terms[1].var->is_const());
    assert(endtimes[0]->terms[1].var.get() == seq.get_const(IR::TagVal(4.2)));
    assert(endtimes[1]->tconst == 0);
    assert(endtimes[1]->terms.size() == 3);
    assert(endtimes[1]->terms[0].sign == Seq::EventTime::Pos);
    assert(endtimes[1]->terms[0].id == 220);
    assert(endtimes[1]->terms[0].var->is_const());
    assert(endtimes[1]->terms[0].var.get() == seq.get_const(IR::TagVal(1000.0)));
    assert(endtimes[1]->terms[1].sign == Seq::EventTime::Pos);
    assert(endtimes[1]->terms[1].id == 101);
    assert(endtimes[1]->terms[1].var->is_const());
    assert(endtimes[1]->terms[1].var.get() == seq.get_const(IR::TagVal(4.2)));
    assert(endtimes[1]->terms[2].sign == Seq::EventTime::Unknown);
    assert(endtimes[1]->terms[2].id == 299);
    assert(endtimes[1]->terms[2].var->is_const());
    assert(endtimes[1]->terms[2].var.get() == seq.get_const(IR::TagVal(4.2)));
}

static void test_measure(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.chnnames.push_back("test_chn3");
    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 9.2;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[1].args[0].f64 = 1020;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 101, 1, 2});
    basicseq.times.push_back({Seq::EventTime::NonNeg, 220, 2, 0});
    basicseq.measures[12] = {1 /* time */, 3 /* channel */};

    builder.buildseq();

    assert(seq.get_basicseqs().size() == 1);
    assert(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    assert(!bseq.get_pulses(1));
    assert(!bseq.get_pulses(2));
    assert(bseq.get_pulses(3));
    assert(bseq.get_pulses(3)->size() == 1);

    auto &measure = bseq.get_pulses(3)->front();
    assert(measure.is_measure());
    assert(measure.id() == 12);
    assert(!measure.len());
    assert(measure.val()->is_extern());
    assert(measure.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 12));
    auto &t = measure.start();
    assert(t.tconst == 0);
    assert(t.terms.size() == 2);
    assert(t.terms[0].sign == Seq::EventTime::NonNeg);
    assert(t.terms[0].id == 220);
    assert(t.terms[0].var->is_const());
    assert(t.terms[0].var.get() == seq.get_const(IR::TagVal(1020.0)));
    assert(t.terms[1].sign == Seq::EventTime::Pos);
    assert(t.terms[1].id == 101);
    assert(t.terms[1].var->is_const());
    assert(t.terms[1].var.get() == seq.get_const(IR::TagVal(9.2)));
}

static void test_measure_val(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.nodes.resize(2);
    builder.nodes[0].op = Seq::Builder::OpCode::Identity;
    builder.nodes[0].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[0].args[0].f64 = 134.1;
    builder.nodes[1].op = Seq::Builder::OpCode::Identity;
    builder.nodes[1].argtypes[0] = Seq::Builder::ArgType::Measure;
    builder.nodes[1].args[0].id = 19;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.measures[19] = {1 /* time */, 2 /* channel */};
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 2 /* value */, 1 /* channel */};

    builder.buildseq();

    assert(seq.get_basicseqs().size() == 1);
    assert(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    assert(bseq.get_pulses(1));
    assert(bseq.get_pulses(1)->size() == 1);
    assert(bseq.get_pulses(2));
    assert(bseq.get_pulses(2)->size() == 1);

    auto &measure = bseq.get_pulses(2)->front();
    assert(measure.is_measure());
    assert(measure.id() == 19);
    assert(!measure.len());
    assert(measure.val()->is_extern());
    assert(measure.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 19));
    auto &t = measure.start();
    assert(t.tconst == 0);
    assert(t.terms.size() == 1);
    assert(t.terms[0].sign == Seq::EventTime::Pos);
    assert(t.terms[0].id == 105);
    assert(t.terms[0].var->is_const());
    assert(t.terms[0].var.get() == seq.get_const(IR::TagVal(134.1)));

    auto &pulse = bseq.get_pulses(1)->front();
    assert(!pulse.is_measure());
    assert(pulse.id() == 29);
    assert(!pulse.len());
    assert(pulse.val()->is_extern());
    assert(pulse.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 19));
    assert(pulse.val() == measure.val());
    assert(!pulse.len());
    auto &t2 = pulse.start();
    assert(t2.tconst == 0);
    assert(t2.terms.size() == 1);
    assert(t2.terms[0].sign == Seq::EventTime::Pos);
    assert(t2.terms[0].id == 105);
    assert(t2.terms[0].var->is_const());
    assert(t2.terms[0].var.get() == seq.get_const(IR::TagVal(134.1)));
    assert(&pulse.start() == &measure.start());
}

static void test_output(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.chnnames.push_back("test_chn3");
    builder.chnnames.push_back("test_chn4");
    builder.nodes.resize(6);
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
    builder.nodes[2].args[1].id = 1;
    builder.nodes[3].op = Seq::Builder::OpCode::Mul;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 2;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[1].id = 3;
    builder.nodes[4].op = Seq::Builder::OpCode::Mul;
    builder.nodes[4].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[4].args[0].id = 1;
    builder.nodes[4].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[4].args[1].id = 0;
    builder.nodes[5].op = Seq::Builder::OpCode::Add;
    builder.nodes[5].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[0].id = 4;
    builder.nodes[5].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[5].args[1].id = 5;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 1 /* value */, 1 /* channel */};
    basicseq.outputs[65] = {1 /* time */, 2 /* length */, 4 /* value */, 2 /* channel */};
    basicseq.outputs[32] = {1 /* time */, 2 /* length */, 5 /* value */, 3 /* channel */};
    basicseq.outputs[49] = {1 /* time */, 2 /* length */, 6 /* value */, 4 /* channel */};

    builder.buildseq();

    assert(seq.get_basicseqs().size() == 1);
    assert(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    assert(bseq.get_pulses(1));
    assert(bseq.get_pulses(1)->size() == 1);
    assert(bseq.get_pulses(2));
    assert(bseq.get_pulses(2)->size() == 1);
    assert(bseq.get_pulses(3));
    assert(bseq.get_pulses(3)->size() == 1);
    assert(bseq.get_pulses(4));
    assert(bseq.get_pulses(4)->size() == 1);

    auto &p1 = bseq.get_pulses(1)->front();
    assert(!p1.is_measure());
    assert(p1.id() == 8);
    assert(!p1.len());
    assert(!p1.needs_oldval());
    assert(p1.val()->is_const());
    assert(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    assert(!p2.is_measure());
    assert(p2.id() == 65);
    assert(!p2.len());
    assert(p2.needs_oldval());
    assert(p2.val()->is_call());

    auto &p3 = bseq.get_pulses(3)->front();
    assert(!p3.is_measure());
    assert(p3.id() == 32);
    assert(p3.len());
    assert(!p3.needs_oldval());
    assert(p3.val()->is_call());

    auto &p4 = bseq.get_pulses(4)->front();
    assert(!p4.is_measure());
    assert(p4.id() == 49);
    assert(p4.len());
    assert(p4.needs_oldval());
    assert(p4.val()->is_call());
}

static void test_output_len0(llvm::LLVMContext &llvm_ctx)
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
    builder.nodes[2].argtypes[0] = Seq::Builder::ArgType::ConstFloat64;
    builder.nodes[2].args[0].f64 = 5;
    builder.nodes[2].argtypes[1] = Seq::Builder::ArgType::Arg;
    builder.nodes[2].args[1].id = 1;
    builder.nodes[3].op = Seq::Builder::OpCode::Mul;
    builder.nodes[3].argtypes[0] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[0].id = 2;
    builder.nodes[3].argtypes[1] = Seq::Builder::ArgType::Node;
    builder.nodes[3].args[1].id = 3;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::EventTime::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 0 /* length */, 1 /* value */, 1 /* channel */};
    basicseq.outputs[65] = {1 /* time */, 0 /* length */, 4 /* value */, 2 /* channel */};

    builder.buildseq();

    assert(seq.get_basicseqs().size() == 1);
    assert(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    assert(bseq.get_pulses(1));
    assert(bseq.get_pulses(1)->size() == 1);
    assert(bseq.get_pulses(2));
    assert(bseq.get_pulses(2)->size() == 1);

    auto &p1 = bseq.get_pulses(1)->front();
    assert(!p1.is_measure());
    assert(p1.id() == 8);
    assert(!p1.len());
    assert(!p1.needs_oldval());
    assert(p1.val()->is_const());
    assert(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    assert(!p2.is_measure());
    assert(p2.id() == 65);
    assert(!p2.len());
    assert(p2.needs_oldval());
    assert(p2.val()->is_call());
}

static void test_branch(llvm::LLVMContext &llvm_ctx)
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
    basicseq1.branches.push_back({102 /* branch ID */, 0 /* target */, 1 /* condition */});
    basicseq2.default_target = 0;
    basicseq2.branches.push_back({103 /* branch ID */, 2 /* target */, 2 /* condition */});

    builder.buildseq();

    std::vector<const Seq::BasicSeq*> seqs;
    for (auto &bs: seq.get_basicseqs())
        seqs.push_back(&bs);
    assert(seqs.size() == 2);

    assert(seqs[0]->id() == 1);
    assert(seqs[0]->get_default_branch() == seqs[1]);
    assert(seqs[0]->get_branches().size() == 1);
    assert(seqs[0]->get_branches()[0].target == nullptr);
    assert(seqs[0]->get_branches()[0].id == 102);
    assert(seqs[0]->get_branches()[0].cond.get() == seq.get_const(IR::TagVal(4.2)));

    assert(seqs[1]->id() == 2);
    assert(seqs[1]->get_default_branch() == nullptr);
    assert(seqs[1]->get_branches().size() == 1);
    assert(seqs[1]->get_branches()[0].target == seqs[1]);
    assert(seqs[1]->get_branches()[0].id == 103);
    assert(seqs[1]->get_branches()[0].cond.get() == seq.get_const(IR::TagVal(1000.0)));
}

int main()
{
    auto llvm_ctx = LLVM::new_context();

    test_channels(*llvm_ctx);
    test_default(*llvm_ctx);
    test_global(*llvm_ctx);
    test_assign(*llvm_ctx);
    test_node_cse(*llvm_ctx);
    test_node_commute(*llvm_ctx);
    test_node_interp(*llvm_ctx);

    test_time(*llvm_ctx);
    test_measure(*llvm_ctx);
    test_measure_val(*llvm_ctx);
    test_output(*llvm_ctx);
    test_output_len0(*llvm_ctx);

    test_branch(*llvm_ctx);

    return 0;
}
