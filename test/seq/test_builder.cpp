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

#include "compiler_helper.h"

#include "../../lib/nacs-seq/builder.h"

#include "../../lib/nacs-utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

using namespace NaCs;

static auto llvm_ctx = LLVM::new_context();

TEST_CASE("channels") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");

    builder.buildseq();

    REQUIRE(seq.get_chn_id("test_chn2") == 2);
    REQUIRE(seq.get_chn_id("test_chn1") == 1);
    REQUIRE(seq.get_chn_id("test_chn") == 0); // Not existing channel
    REQUIRE(seq.get_chn_names().size() == 2);
    REQUIRE(seq.get_chn_names()[0] == "test_chn1");
    REQUIRE(seq.get_chn_names()[1] == "test_chn2");

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("default") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.defvals[1] = IR::TagVal(2.3);

    builder.buildseq();

    REQUIRE(seq.get_chn_id("test_chn1") == 1);
    auto defval = seq.defval(1);
    REQUIRE(defval == 2.3);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("global") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Bool);
    builder.slots.push_back(IR::Type::Int32);

    builder.buildseq();

    REQUIRE(seq.get_slots().size() == 3);
    auto slot1 = seq.get_slots()[0].get();
    auto slot2 = seq.get_slots()[1].get();
    auto slot3 = seq.get_slots()[2].get();
    REQUIRE(slot1 == seq.get_slot(IR::Type::Bool, 0));
    REQUIRE(slot2 == seq.get_slot(IR::Type::Int32, 1));
    REQUIRE(slot3 == seq.get_slot(IR::Type::Int32, 2));
    REQUIRE(slot1->type() == IR::Type::Float64);
    REQUIRE(slot2->type() == IR::Type::Bool);
    REQUIRE(slot3->type() == IR::Type::Int32);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("assign") {
    Seq::Seq seq(*llvm_ctx);
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

    REQUIRE(seq.get_slots().size() == 2);
    auto slot1 = seq.get_slots()[0].get();
    auto slot2 = seq.get_slots()[1].get();
    REQUIRE(slot1 == seq.get_slot(IR::Type::Bool, 0));
    REQUIRE(slot2 == seq.get_slot(IR::Type::Int32, 1));
    REQUIRE(slot1->type() == IR::Type::Float64);
    REQUIRE(slot2->type() == IR::Type::Int32);

    REQUIRE(seq.get_basicseqs().size() == 1);

    REQUIRE(seq.get_basicseqs().front().id() == 1);

    auto &assigns = seq.get_basicseqs().front().get_assigns();
    REQUIRE(assigns.size() == 2);
    REQUIRE(assigns.find(0)->second.val.get() == slot2);
    REQUIRE(assigns.find(0)->second.id == 41);
    REQUIRE(assigns.find(1)->second.val.get() == slot1);
    REQUIRE(assigns.find(1)->second.id == 52);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("node_cse") {
    Seq::Seq seq(*llvm_ctx);
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

    REQUIRE(vars.size() == 5);

    REQUIRE(vars[4]->varid() == 0);
    REQUIRE(vars[4]->is_extern());
    REQUIRE(vars[4]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0)));
    REQUIRE(vars[3]->varid() == 1);
    REQUIRE(vars[3]->is_extern());
    REQUIRE(vars[3]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(1)));
    REQUIRE(vars[2]->varid() == 2);
    REQUIRE(vars[2]->is_call());
    REQUIRE(vars[2]->get_callee().is_llvm);
    REQUIRE(vars[2]->nfreeargs() == 0);
    REQUIRE(vars[2]->args().size() == 2);
    REQUIRE(vars[2]->args()[0].compare(Seq::Arg::create_var(vars[4])) == 0);
    REQUIRE(vars[2]->args()[1].compare(Seq::Arg::create_var(vars[3])) == 0);
    REQUIRE(vars[1]->varid() == 3);
    REQUIRE(vars[1]->is_call());
    REQUIRE(vars[1]->get_callee().is_llvm);
    REQUIRE(vars[1]->nfreeargs() == 0);
    REQUIRE(vars[1]->args().size() == 2);
    REQUIRE(vars[1]->args()[0].compare(Seq::Arg::create_var(vars[2])) == 0);
    REQUIRE(vars[1]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    REQUIRE(vars[0]->varid() == 4);
    REQUIRE(vars[0]->is_call());
    REQUIRE(vars[0]->get_callee().is_llvm);
    REQUIRE(vars[0]->nfreeargs() == 0);
    REQUIRE(vars[0]->args().size() == 2);
    REQUIRE(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[1])) == 0);
    REQUIRE(vars[0]->args()[1].compare(Seq::Arg::create_const(IR::TagVal(8))) == 0);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("node_commute") {
    Seq::Seq seq(*llvm_ctx);
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

    REQUIRE(vars.size() == 4);

    REQUIRE(vars[3]->varid() == 0);
    REQUIRE(vars[3]->is_extern());
    REQUIRE(vars[3]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0)));
    REQUIRE(vars[2]->varid() == 1);
    REQUIRE(vars[2]->is_extern());
    REQUIRE(vars[2]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(1)));
    REQUIRE(vars[1]->varid() == 2);
    REQUIRE(vars[1]->is_call());
    REQUIRE(vars[1]->get_callee().is_llvm);
    REQUIRE(vars[1]->nfreeargs() == 0);
    REQUIRE(vars[1]->args().size() == 2);
    REQUIRE(vars[1]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    REQUIRE(vars[1]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    REQUIRE(vars[0]->varid() == 3);
    REQUIRE(vars[0]->is_call());
    REQUIRE(vars[0]->get_callee().is_llvm);
    REQUIRE(vars[0]->nfreeargs() == 0);
    REQUIRE(vars[0]->args().size() == 2);
    REQUIRE(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    REQUIRE(vars[0]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("node_interp") {
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
    builder.nodes[0].data_id = 0;
    builder.slots.push_back(IR::Type::Bool);
    builder.slots.push_back(IR::Type::Float64);
    builder.slots.push_back(IR::Type::Int32);
    builder.datas.push_back({1, 2, 5, 6});

    builder.buildseq();

    std::vector<Seq::Var*> vars(seq.env().begin(), seq.env().end());

    REQUIRE(vars.size() == 4);

    REQUIRE(vars[3]->varid() == 0);
    REQUIRE(vars[3]->is_extern());
    REQUIRE(vars[3]->get_extern() == std::make_pair(IR::Type::Bool, uint64_t(0)));
    REQUIRE(vars[2]->varid() == 1);
    REQUIRE(vars[2]->is_extern());
    REQUIRE(vars[2]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(1)));
    REQUIRE(vars[1]->varid() == 2);
    REQUIRE(vars[1]->is_extern());
    REQUIRE(vars[1]->get_extern() == std::make_pair(IR::Type::Int32, uint64_t(2)));
    REQUIRE(vars[0]->varid() == 3);
    REQUIRE(vars[0]->is_call());
    REQUIRE(vars[0]->get_callee().is_llvm);
    REQUIRE(vars[0]->nfreeargs() == 0);
    REQUIRE(vars[0]->args().size() == 3);
    REQUIRE(vars[0]->args()[0].compare(Seq::Arg::create_var(vars[3])) == 0);
    REQUIRE(vars[0]->args()[1].compare(Seq::Arg::create_var(vars[2])) == 0);
    REQUIRE(vars[0]->args()[2].compare(Seq::Arg::create_var(vars[1])) == 0);

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("time") {
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
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 2});
    basicseq.times.push_back({Seq::Sign::Pos, 220, 2, 0});
    basicseq.times.push_back({Seq::Sign::Unknown, 299, 1, 1});
    basicseq.endtimes.push_back(1);
    basicseq.endtimes.push_back(3);

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    std::vector<decltype(&*bseq.get_assumes().begin())> assumes;
    for (auto &assume: bseq.get_assumes())
        assumes.push_back(&assume);
    // Assumptions are sorted by variables, sort them by ID so that we can check more easily
    std::sort(assumes.begin(), assumes.end(),
              [&] (auto &a1, auto &a2) { return a1->second.id < a2->second.id; });

    REQUIRE(assumes.size() == 2);
    REQUIRE(assumes[0]->second.sign == Seq::Sign::Pos);
    REQUIRE(assumes[0]->second.id == 101);
    REQUIRE(assumes[0]->first->is_const());
    REQUIRE(assumes[0]->first.get() == seq.get_const(IR::TagVal(4.2)));
    REQUIRE(assumes[1]->second.sign == Seq::Sign::Pos);
    REQUIRE(assumes[1]->second.id == 220);
    REQUIRE(assumes[1]->first->is_const());
    REQUIRE(assumes[1]->first.get() == seq.get_const(IR::TagVal(1000.0)));

    std::vector<const Seq::EventTime*> endtimes;
    for (auto &endtime: bseq.get_endtimes())
        endtimes.push_back(endtime.get());
    REQUIRE(endtimes.size() == 2);
    REQUIRE(endtimes[0]->tconst == 0);
    REQUIRE(endtimes[0]->terms.size() == 2);
    REQUIRE(endtimes[0]->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(endtimes[0]->terms[0].id == 220);
    REQUIRE(endtimes[0]->terms[0].var->is_const());
    REQUIRE(endtimes[0]->terms[0].var.get() == seq.get_const(IR::TagVal(1000.0)));
    REQUIRE(endtimes[0]->terms[1].sign == Seq::Sign::Pos);
    REQUIRE(endtimes[0]->terms[1].id == 101);
    REQUIRE(endtimes[0]->terms[1].var->is_const());
    REQUIRE(endtimes[0]->terms[1].var.get() == seq.get_const(IR::TagVal(4.2)));
    REQUIRE(endtimes[1]->tconst == 0);
    REQUIRE(endtimes[1]->terms.size() == 3);
    REQUIRE(endtimes[1]->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(endtimes[1]->terms[0].id == 220);
    REQUIRE(endtimes[1]->terms[0].var->is_const());
    REQUIRE(endtimes[1]->terms[0].var.get() == seq.get_const(IR::TagVal(1000.0)));
    REQUIRE(endtimes[1]->terms[1].sign == Seq::Sign::Pos);
    REQUIRE(endtimes[1]->terms[1].id == 101);
    REQUIRE(endtimes[1]->terms[1].var->is_const());
    REQUIRE(endtimes[1]->terms[1].var.get() == seq.get_const(IR::TagVal(4.2)));
    REQUIRE(endtimes[1]->terms[2].sign == Seq::Sign::Unknown);
    REQUIRE(endtimes[1]->terms[2].id == 299);
    REQUIRE(endtimes[1]->terms[2].var->is_const());
    REQUIRE(endtimes[1]->terms[2].var.get() == seq.get_const(IR::TagVal(4.2)));

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("measure") {
    Seq::Seq seq(*llvm_ctx);
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
    basicseq.times.push_back({Seq::Sign::Pos, 101, 1, 2});
    basicseq.times.push_back({Seq::Sign::NonNeg, 220, 2, 0});
    basicseq.measures[12] = {1 /* time */, 3 /* channel */};

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    REQUIRE(!bseq.get_pulses(1));
    REQUIRE(!bseq.get_pulses(2));
    REQUIRE(bseq.get_pulses(3));
    REQUIRE(bseq.get_pulses(3)->size() == 1);

    auto &measure = bseq.get_pulses(3)->front();
    REQUIRE(measure.is_measure());
    REQUIRE(measure.id() == 12);
    REQUIRE(!measure.len());
    REQUIRE(measure.val()->is_extern());
    REQUIRE(measure.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 12));
    auto &t = measure.start();
    REQUIRE(t.tconst == 0);
    REQUIRE(t.terms.size() == 2);
    REQUIRE(t.terms[0].sign == Seq::Sign::NonNeg);
    REQUIRE(t.terms[0].id == 220);
    REQUIRE(t.terms[0].var->is_const());
    REQUIRE(t.terms[0].var.get() == seq.get_const(IR::TagVal(1020.0)));
    REQUIRE(t.terms[1].sign == Seq::Sign::Pos);
    REQUIRE(t.terms[1].id == 101);
    REQUIRE(t.terms[1].var->is_const());
    REQUIRE(t.terms[1].var.get() == seq.get_const(IR::TagVal(9.2)));

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("measure_val") {
    Seq::Seq seq(*llvm_ctx);
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
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.measures[19] = {1 /* time */, 2 /* channel */};
    basicseq.outputs[29] = {1 /* time */, 0 /* length */, 2 /* value */,
        0 /* cond */, 1 /* channel */};

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    REQUIRE(bseq.get_pulses(1));
    REQUIRE(bseq.get_pulses(1)->size() == 1);
    REQUIRE(bseq.get_pulses(2));
    REQUIRE(bseq.get_pulses(2)->size() == 1);

    auto &measure = bseq.get_pulses(2)->front();
    REQUIRE(measure.is_measure());
    REQUIRE(measure.id() == 19);
    REQUIRE(!measure.len());
    REQUIRE(measure.val()->is_extern());
    REQUIRE(measure.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 19));
    auto &t = measure.start();
    REQUIRE(t.tconst == 0);
    REQUIRE(t.terms.size() == 1);
    REQUIRE(t.terms[0].sign == Seq::Sign::Pos);
    REQUIRE(t.terms[0].id == 105);
    REQUIRE(t.terms[0].var->is_const());
    REQUIRE(t.terms[0].var.get() == seq.get_const(IR::TagVal(134.1)));

    auto &pulse = bseq.get_pulses(1)->front();
    REQUIRE(!pulse.is_measure());
    REQUIRE(pulse.id() == 29);
    REQUIRE(!pulse.len());
    REQUIRE(pulse.val()->is_extern());
    REQUIRE(pulse.val()->get_extern() == std::make_pair(IR::Type::Float64,
                                                         (uint64_t(1) << 32) | 19));
    REQUIRE(pulse.val() == measure.val());
    REQUIRE(!pulse.len());
    auto &t2 = pulse.start();
    REQUIRE(t2.tconst == 0);
    REQUIRE(t2.terms.size() == 1);
    REQUIRE(t2.terms[0].sign == Seq::Sign::Pos);
    REQUIRE(t2.terms[0].id == 105);
    REQUIRE(t2.terms[0].var->is_const());
    REQUIRE(t2.terms[0].var.get() == seq.get_const(IR::TagVal(134.1)));
    REQUIRE(&pulse.start() == &measure.start());

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("output") {
    Seq::Seq seq(*llvm_ctx);
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
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 2 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};
    basicseq.outputs[65] = {1 /* time */, 2 /* length */, 4 /* value */,
        0 /* cond */,2 /* channel */};
    basicseq.outputs[32] = {1 /* time */, 2 /* length */, 5 /* value */,
        0 /* cond */,3 /* channel */};
    basicseq.outputs[49] = {1 /* time */, 2 /* length */, 6 /* value */,
        0 /* cond */,4 /* channel */};

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    REQUIRE(bseq.get_pulses(1));
    REQUIRE(bseq.get_pulses(1)->size() == 1);
    REQUIRE(bseq.get_pulses(2));
    REQUIRE(bseq.get_pulses(2)->size() == 1);
    REQUIRE(bseq.get_pulses(3));
    REQUIRE(bseq.get_pulses(3)->size() == 1);
    REQUIRE(bseq.get_pulses(4));
    REQUIRE(bseq.get_pulses(4)->size() == 1);

    auto &p1 = bseq.get_pulses(1)->front();
    REQUIRE(!p1.is_measure());
    REQUIRE(p1.id() == 8);
    REQUIRE(!p1.len());
    REQUIRE(!p1.cond());
    REQUIRE(!p1.needs_oldval());
    REQUIRE(p1.val()->is_const());
    REQUIRE(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    REQUIRE(!p2.is_measure());
    REQUIRE(p2.id() == 65);
    REQUIRE(!p2.len());
    REQUIRE(!p2.cond());
    REQUIRE(p2.needs_oldval());
    REQUIRE(p2.val()->is_call());

    auto &p3 = bseq.get_pulses(3)->front();
    REQUIRE(!p3.is_measure());
    REQUIRE(p3.id() == 32);
    REQUIRE(p3.len());
    REQUIRE(!p3.cond());
    REQUIRE(!p3.needs_oldval());
    REQUIRE(p3.val()->is_call());

    auto &p4 = bseq.get_pulses(4)->front();
    REQUIRE(!p4.is_measure());
    REQUIRE(p4.id() == 49);
    REQUIRE(p4.len());
    REQUIRE(!p4.cond());
    REQUIRE(p4.needs_oldval());
    REQUIRE(p4.val()->is_call());

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("output_len0") {
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
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[8] = {1 /* time */, 0 /* length */, 1 /* value */,
        0 /* cond */, 1 /* channel */};
    basicseq.outputs[65] = {1 /* time */, 0 /* length */, 4 /* value */,
        0 /* cond */, 2 /* channel */};

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    REQUIRE(bseq.get_pulses(1));
    REQUIRE(bseq.get_pulses(1)->size() == 1);
    REQUIRE(bseq.get_pulses(2));
    REQUIRE(bseq.get_pulses(2)->size() == 1);

    auto &p1 = bseq.get_pulses(1)->front();
    REQUIRE(!p1.is_measure());
    REQUIRE(p1.id() == 8);
    REQUIRE(!p1.len());
    REQUIRE(!p1.cond());
    REQUIRE(!p1.needs_oldval());
    REQUIRE(p1.val()->is_const());
    REQUIRE(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    REQUIRE(!p2.is_measure());
    REQUIRE(p2.id() == 65);
    REQUIRE(!p2.len());
    REQUIRE(!p2.cond());
    REQUIRE(p2.needs_oldval());
    REQUIRE(p2.val()->is_call());

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("output_cond") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.chnnames.push_back("test_chn1");
    builder.chnnames.push_back("test_chn2");
    builder.slots.push_back(IR::Type::Bool);
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
    builder.nodes[4].op = Seq::Builder::OpCode::Identity;
    builder.nodes[4].argtypes[0] = Seq::Builder::ArgType::Global;
    builder.nodes[4].args[0].id = 0;
    builder.nodes[5].op = Seq::Builder::OpCode::Identity;
    builder.nodes[5].argtypes[0] = Seq::Builder::ArgType::ConstBool;
    builder.nodes[5].args[0].b = false;
    builder.seqs.emplace_back();
    auto &basicseq = builder.seqs.back();
    basicseq.times.push_back({Seq::Sign::Pos, 105, 1, 0});
    basicseq.outputs[3] = {1 /* time */, 0 /* length */, 1 /* value */,
        5 /* cond */, 1 /* channel */};
    basicseq.outputs[57] = {1 /* time */, 0 /* length */, 4 /* value */,
        6 /* cond */, 2 /* channel */};

    builder.buildseq();

    REQUIRE(seq.get_basicseqs().size() == 1);
    REQUIRE(seq.get_basicseqs().front().id() == 1);
    auto &bseq = seq.get_basicseqs().front();

    REQUIRE(bseq.get_pulses(1));
    REQUIRE(bseq.get_pulses(1)->size() == 1);
    REQUIRE(bseq.get_pulses(2));
    REQUIRE(bseq.get_pulses(2)->size() == 1);

    auto &p1 = bseq.get_pulses(1)->front();
    REQUIRE(!p1.is_measure());
    REQUIRE(p1.id() == 3);
    REQUIRE(!p1.len());
    REQUIRE(p1.cond());
    REQUIRE(p1.cond()->is_extern());
    REQUIRE(p1.cond() == seq.get_slots()[0].get());
    REQUIRE(!p1.needs_oldval());
    REQUIRE(p1.val()->is_const());
    REQUIRE(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    REQUIRE(!p2.is_measure());
    REQUIRE(p2.id() == 57);
    REQUIRE(!p2.len());
    REQUIRE(p2.cond());
    REQUIRE(p2.cond()->is_const());
    REQUIRE(p1.cond()->get_const().is(IR::TagVal(false)));
    REQUIRE(p2.needs_oldval());
    REQUIRE(p2.val()->is_call());

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}

TEST_CASE("branch") {
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
    basicseq1.branches.push_back({102 /* branch ID */, 0 /* target */, 1 /* condition */});
    basicseq2.default_target = 0;
    basicseq2.branches.push_back({103 /* branch ID */, 2 /* target */, 2 /* condition */});

    builder.buildseq();

    std::vector<const Seq::BasicSeq*> seqs;
    for (auto &bs: seq.get_basicseqs())
        seqs.push_back(&bs);
    REQUIRE(seqs.size() == 2);

    REQUIRE(seqs[0]->id() == 1);
    REQUIRE(seqs[0]->get_default_branch() == seqs[1]);
    REQUIRE(seqs[0]->get_branches().size() == 1);
    REQUIRE(seqs[0]->get_branches()[0].target == nullptr);
    REQUIRE(seqs[0]->get_branches()[0].id == 102);
    REQUIRE(seqs[0]->get_branches()[0].cond.get() == seq.get_const(IR::TagVal(4.2)));

    REQUIRE(seqs[1]->id() == 2);
    REQUIRE(seqs[1]->get_default_branch() == nullptr);
    REQUIRE(seqs[1]->get_branches().size() == 1);
    REQUIRE(seqs[1]->get_branches()[0].target == seqs[1]);
    REQUIRE(seqs[1]->get_branches()[0].id == 103);
    REQUIRE(seqs[1]->get_branches()[0].cond.get() == seq.get_const(IR::TagVal(1000.0)));

    seq.prepare();
    seq.optimize();
    CompileSeq cs(seq);
}
