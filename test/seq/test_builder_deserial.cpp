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
#include "../../lib/nacs-utils/streams.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

using namespace NaCs;
using namespace std::literals::string_literals;

struct ostream : uvector_ostream {
    using uvector_ostream::write;
    template<typename T>
    void write(T v)
    {
        STATIC_REQUIRE(std::is_trivial_v<T>);
        write((const char*)&v, sizeof(T));
    }
    void write_string(const char *s)
    {
        write(s, strlen(s) + 1);
    }
    void parse(Seq::Builder &builder)
    {
        auto data = get_buf();
        builder.deserialize(data.data(), data.size());
    }
};

static auto llvm_ctx = LLVM::new_context();

TEST_CASE("channels") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 2);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");

    builder.buildseq();

    REQUIRE(seq.get_chn_id("test_chn2") == 2);
    REQUIRE(seq.get_chn_id("test_chn1") == 1);
    REQUIRE(seq.get_chn_id("test_chn") == 0); // Not existing channel
    REQUIRE(seq.get_chn_names().size() == 2);
    REQUIRE(seq.get_chn_names()[0] == "test_chn1");
    REQUIRE(seq.get_chn_names()[1] == "test_chn2");

    seq.prepare();
    seq.optimize();
}

TEST_CASE("default") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(1);
    stm.write_string("test_chn1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.3);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 1);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.defvals.size() == 1);
    REQUIRE(builder.defvals[1].is(IR::TagVal(2.3)));

    builder.buildseq();

    REQUIRE(seq.get_chn_id("test_chn1") == 1);
    auto defval = seq.defval(1);
    REQUIRE(defval == 2.3);

    seq.prepare();
    seq.optimize();
}

TEST_CASE("global") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Bool);
    stm.write(IR::Type::Int32);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.slots.size() == 3);
    REQUIRE(builder.slots[0] == IR::Type::Float64);
    REQUIRE(builder.slots[1] == IR::Type::Bool);
    REQUIRE(builder.slots[2] == IR::Type::Int32);

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
}

TEST_CASE("assign") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Int32);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(41);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(52);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 2);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[0].id == 0);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[1].args[0].id == 1);
    REQUIRE(builder.slots.size() == 2);
    REQUIRE(builder.slots[0] == IR::Type::Float64);
    REQUIRE(builder.slots[1] == IR::Type::Int32);
    REQUIRE(builder.seqs.size() == 1);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.assignments.size() == 2);
    REQUIRE(basicseq.assignments[0].assign_id == 41);
    REQUIRE(basicseq.assignments[0].global_id == 0);
    REQUIRE(basicseq.assignments[0].val_node == 2);
    REQUIRE(basicseq.assignments[1].assign_id == 52);
    REQUIRE(basicseq.assignments[1].global_id == 1);
    REQUIRE(basicseq.assignments[1].val_node == 1);

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
}

TEST_CASE("node_cse") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(8);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 5
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(4);
    // 6
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(4);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    // 7
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(5);
    stm.write(Seq::Builder::ArgType::ConstInt32);
    stm.write<int32_t>(8);
    // 8
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(6);
    stm.write(Seq::Builder::ArgType::ConstInt32);
    stm.write<int32_t>(8);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Int32);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 8);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[0].id == 0);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[1].args[0].id == 1);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[0].id == 1);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[1].id == 2);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 1);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[1].id == 2);
    REQUIRE(builder.nodes[4].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[4].args[0].id == 3);
    REQUIRE(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[4].args[1].id == 4);
    REQUIRE(builder.nodes[5].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[0].id == 4);
    REQUIRE(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[1].id == 3);
    REQUIRE(builder.nodes[6].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[6].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[6].args[0].id == 5);
    REQUIRE(builder.nodes[6].argtypes[1] == Seq::Builder::ArgType::ConstInt32);
    REQUIRE(builder.nodes[6].args[1].i32 == 8);
    REQUIRE(builder.nodes[7].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[7].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[7].args[0].id == 6);
    REQUIRE(builder.nodes[7].argtypes[1] == Seq::Builder::ArgType::ConstInt32);
    REQUIRE(builder.nodes[7].args[1].i32 == 8);
    REQUIRE(builder.slots.size() == 2);
    REQUIRE(builder.slots[0] == IR::Type::Float64);
    REQUIRE(builder.slots[1] == IR::Type::Int32);

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
}

TEST_CASE("node_commute") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    // 5
    stm.write(Seq::Builder::OpCode::CmpLT);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 6
    stm.write(Seq::Builder::OpCode::CmpGT);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Int32);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 6);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[0].id == 0);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[1].args[0].id == 1);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[0].id == 1);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[1].id == 2);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 2);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[1].id == 1);
    REQUIRE(builder.nodes[4].op == Seq::Builder::OpCode::CmpLT);
    REQUIRE(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[4].args[0].id == 1);
    REQUIRE(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[4].args[1].id == 2);
    REQUIRE(builder.nodes[5].op == Seq::Builder::OpCode::CmpGT);
    REQUIRE(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[0].id == 2);
    REQUIRE(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[1].id == 1);
    REQUIRE(builder.slots.size() == 2);
    REQUIRE(builder.slots[0] == IR::Type::Float64);
    REQUIRE(builder.slots[1] == IR::Type::Int32);

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
}

TEST_CASE("node_interp") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(1);
    // 1
    stm.write(Seq::Builder::OpCode::Interp);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Bool);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Int32);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(0);
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(4);
    stm.write<double>(1);
    stm.write<double>(2);
    stm.write<double>(5);
    stm.write<double>(6);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 1);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Interp);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[0].id == 0);
    REQUIRE(builder.nodes[0].argtypes[1] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[1].id == 1);
    REQUIRE(builder.nodes[0].argtypes[2] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[0].args[2].id == 2);
    REQUIRE(builder.nodes[0].data_id == 0);
    REQUIRE(builder.slots.size() == 3);
    REQUIRE(builder.slots[0] == IR::Type::Bool);
    REQUIRE(builder.slots[1] == IR::Type::Float64);
    REQUIRE(builder.slots[2] == IR::Type::Int32);
    REQUIRE(builder.datas.size() == 1);
    double data[] = {1, 2, 5, 6};
    REQUIRE(builder.datas[0].size() == 4);
    REQUIRE(memcmp(builder.datas[0].data(), data, sizeof(data)) == 0);

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
}

TEST_CASE("time") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(4.2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(3);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(101);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(220);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(299);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(3);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 2);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 4.2);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 1000);
    REQUIRE(builder.seqs.size() == 1);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 3);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 101);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 2);
    REQUIRE(basicseq.times[1].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[1].id == 220);
    REQUIRE(basicseq.times[1].delta_node == 2);
    REQUIRE(basicseq.times[1].prev_id == 0);
    REQUIRE(basicseq.times[2].sign == Seq::Sign::Unknown);
    REQUIRE(basicseq.times[2].id == 299);
    REQUIRE(basicseq.times[2].delta_node == 1);
    REQUIRE(basicseq.times[2].prev_id == 1);
    REQUIRE(basicseq.endtimes.size() == 2);
    REQUIRE(basicseq.endtimes[0] == 1);
    REQUIRE(basicseq.endtimes[1] == 3);

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
}

TEST_CASE("measure") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9.2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1020);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(3);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    stm.write_string("test_chn3");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(101);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Seq::Sign::NonNeg);
        stm.write<uint32_t>(220);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(12);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(3);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 3);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.chnnames[2] == "test_chn3");
    REQUIRE(builder.nodes.size() == 2);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 9.2);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 1020);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 2);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 101);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 2);
    REQUIRE(basicseq.times[1].sign == Seq::Sign::NonNeg);
    REQUIRE(basicseq.times[1].id == 220);
    REQUIRE(basicseq.times[1].delta_node == 2);
    REQUIRE(basicseq.times[1].prev_id == 0);
    REQUIRE(basicseq.measures.size() == 1);
    REQUIRE(basicseq.measures[12].time_id == 1);
    REQUIRE(basicseq.measures[12].chn == 3);

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
}

TEST_CASE("measure_val") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Measure);
    stm.write<uint32_t>(19);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(29);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(19);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 2);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.nodes.size() == 2);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Measure);
    REQUIRE(builder.nodes[1].args[0].id == 19);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.measures.size() == 1);
    REQUIRE(basicseq.measures[19].time_id == 1);
    REQUIRE(basicseq.measures[19].chn == 2);
    REQUIRE(basicseq.outputs.size() == 1);
    REQUIRE(basicseq.outputs[29].time_id == 1);
    REQUIRE(basicseq.outputs[29].len_node == 0);
    REQUIRE(basicseq.outputs[29].val_node == 2);
    REQUIRE(basicseq.outputs[29].cond_node == 0);
    REQUIRE(basicseq.outputs[29].chn == 1);

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
}

TEST_CASE("output") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9123.45);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(5);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(1);
    // 4
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    // 5
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 6
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(4);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(5);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(4);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    stm.write_string("test_chn3");
    stm.write_string("test_chn4");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(4);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        // 2
        stm.write<uint32_t>(65);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        // 3
        stm.write<uint32_t>(32);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(3);
        // 4
        stm.write<uint32_t>(49);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(4);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 4);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.chnnames[2] == "test_chn3");
    REQUIRE(builder.chnnames[3] == "test_chn4");
    REQUIRE(builder.nodes.size() == 6);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 9123.45);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[2].args[0].f64 == 5);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[2].args[1].id == 1);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 2);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[1].id == 3);
    REQUIRE(builder.nodes[4].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[4].args[0].id == 1);
    REQUIRE(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[4].args[1].id == 0);
    REQUIRE(builder.nodes[5].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[0].id == 4);
    REQUIRE(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[5].args[1].id == 5);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.outputs.size() == 4);
    REQUIRE(basicseq.outputs[8].time_id == 1);
    REQUIRE(basicseq.outputs[8].len_node == 2);
    REQUIRE(basicseq.outputs[8].val_node == 1);
    REQUIRE(basicseq.outputs[8].cond_node == 0);
    REQUIRE(basicseq.outputs[8].chn == 1);
    REQUIRE(basicseq.outputs[65].time_id == 1);
    REQUIRE(basicseq.outputs[65].len_node == 2);
    REQUIRE(basicseq.outputs[65].val_node == 4);
    REQUIRE(basicseq.outputs[65].cond_node == 0);
    REQUIRE(basicseq.outputs[65].chn == 2);
    REQUIRE(basicseq.outputs[32].time_id == 1);
    REQUIRE(basicseq.outputs[32].len_node == 2);
    REQUIRE(basicseq.outputs[32].val_node == 5);
    REQUIRE(basicseq.outputs[32].cond_node == 0);
    REQUIRE(basicseq.outputs[32].chn == 3);
    REQUIRE(basicseq.outputs[49].time_id == 1);
    REQUIRE(basicseq.outputs[49].len_node == 2);
    REQUIRE(basicseq.outputs[49].val_node == 6);
    REQUIRE(basicseq.outputs[49].cond_node == 0);
    REQUIRE(basicseq.outputs[49].chn == 4);

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
    REQUIRE(!p1.needs_oldval());
    REQUIRE(p1.val()->is_const());
    REQUIRE(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    REQUIRE(!p2.is_measure());
    REQUIRE(p2.id() == 65);
    REQUIRE(!p2.len());
    REQUIRE(p2.needs_oldval());
    REQUIRE(p2.val()->is_call());

    auto &p3 = bseq.get_pulses(3)->front();
    REQUIRE(!p3.is_measure());
    REQUIRE(p3.id() == 32);
    REQUIRE(p3.len());
    REQUIRE(!p3.needs_oldval());
    REQUIRE(p3.val()->is_call());

    auto &p4 = bseq.get_pulses(4)->front();
    REQUIRE(!p4.is_measure());
    REQUIRE(p4.id() == 49);
    REQUIRE(p4.len());
    REQUIRE(p4.needs_oldval());
    REQUIRE(p4.val()->is_call());

    seq.prepare();
    seq.optimize();
}

TEST_CASE("output_len0") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(4);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9123.45);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(5);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(1);
    // 4
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        // 2
        stm.write<uint32_t>(65);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 2);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.nodes.size() == 4);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 9123.45);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[2].args[0].f64 == 5);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[2].args[1].id == 1);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 2);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[1].id == 3);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.outputs.size() == 2);
    REQUIRE(basicseq.outputs[8].time_id == 1);
    REQUIRE(basicseq.outputs[8].len_node == 0);
    REQUIRE(basicseq.outputs[8].val_node == 1);
    REQUIRE(basicseq.outputs[8].cond_node == 0);
    REQUIRE(basicseq.outputs[8].chn == 1);
    REQUIRE(basicseq.outputs[65].time_id == 1);
    REQUIRE(basicseq.outputs[65].len_node == 0);
    REQUIRE(basicseq.outputs[65].val_node == 4);
    REQUIRE(basicseq.outputs[65].cond_node == 0);
    REQUIRE(basicseq.outputs[65].chn == 2);

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
    REQUIRE(!p1.needs_oldval());
    REQUIRE(p1.val()->is_const());
    REQUIRE(p1.val()->get_const().is(IR::TagVal(134.1)));

    auto &p2 = bseq.get_pulses(2)->front();
    REQUIRE(!p2.is_measure());
    REQUIRE(p2.id() == 65);
    REQUIRE(!p2.len());
    REQUIRE(p2.needs_oldval());
    REQUIRE(p2.val()->is_call());

    seq.prepare();
    seq.optimize();
}

TEST_CASE("output_cond") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9123.45);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(5);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(1);
    // 4
    stm.write(Seq::Builder::OpCode::Mul);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    // 5
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 6
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstBool);
    stm.write<bool>(false);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Bool);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(1);
        // 2
        stm.write<uint32_t>(57);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(2);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 2);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.slots.size() == 1);
    REQUIRE(builder.slots[0] == IR::Type::Bool);
    REQUIRE(builder.nodes.size() == 6);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 9123.45);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[2].args[0].f64 == 5);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[2].args[1].id == 1);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Mul);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 2);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[1].id == 3);
    REQUIRE(builder.nodes[4].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Global);
    REQUIRE(builder.nodes[4].args[0].id == 0);
    REQUIRE(builder.nodes[5].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::ConstBool);
    REQUIRE(builder.nodes[5].args[0].b == false);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.outputs.size() == 2);
    REQUIRE(basicseq.outputs[3].time_id == 1);
    REQUIRE(basicseq.outputs[3].len_node == 0);
    REQUIRE(basicseq.outputs[3].val_node == 1);
    REQUIRE(basicseq.outputs[3].cond_node == 5);
    REQUIRE(basicseq.outputs[3].chn == 1);
    REQUIRE(basicseq.outputs[57].time_id == 1);
    REQUIRE(basicseq.outputs[57].len_node == 0);
    REQUIRE(basicseq.outputs[57].val_node == 4);
    REQUIRE(basicseq.outputs[57].cond_node == 6);
    REQUIRE(basicseq.outputs[57].chn == 2);

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
}

TEST_CASE("branch") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(2);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(4.2);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(2);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(102);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        // [default_target: 4B]
        stm.write<uint32_t>(2);
    }
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(103);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(2);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.nodes.size() == 2);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 4.2);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 1000);
    auto &basicseq1 = builder.seqs[0];
    auto &basicseq2 = builder.seqs[1];
    REQUIRE(basicseq1.default_target == 2);
    REQUIRE(basicseq1.branches.size() == 1);
    REQUIRE(basicseq1.branches[0].branch_id == 102);
    REQUIRE(basicseq1.branches[0].target_id == 0);
    REQUIRE(basicseq1.branches[0].cond_node == 1);
    REQUIRE(basicseq2.default_target == 0);
    REQUIRE(basicseq2.branches.size() == 1);
    REQUIRE(basicseq2.branches[0].branch_id == 103);
    REQUIRE(basicseq2.branches[0].target_id == 2);
    REQUIRE(basicseq2.branches[0].cond_node == 2);

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
}

TEST_CASE("output_noramp_noold") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(3);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9123.45);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(1);
    stm.write_string("test_chn1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(1);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 1);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.nodes.size() == 3);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 9123.45);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[0].id == 2);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[2].args[1].id == 0);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.outputs.size() == 1);
    REQUIRE(basicseq.outputs[8].time_id == 1);
    REQUIRE(basicseq.outputs[8].len_node == 2);
    REQUIRE(basicseq.outputs[8].val_node == 3);
    REQUIRE(basicseq.outputs[8].cond_node == 0);
    REQUIRE(basicseq.outputs[8].chn == 1);
    REQUIRE(builder.noramp_chns.size() == 1);
    REQUIRE(builder.noramp_chns[0] == 1);

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

    seq.prepare();
    seq.optimize();
}

TEST_CASE("output_noramp_old") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(4);
    // 1
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(134.1);
    // 2
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(9123.45);
    // 3
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 4
    stm.write(Seq::Builder::OpCode::Add);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(3);
    stm.write(Seq::Builder::ArgType::Arg);
    stm.write<uint32_t>(1);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("test_chn1");
    stm.write_string("test_chn2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(2);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Seq::Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        // 1
        stm.write<uint32_t>(9);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    stm.parse(builder);

    REQUIRE(builder.chnnames.size() == 2);
    REQUIRE(builder.chnnames[0] == "test_chn1");
    REQUIRE(builder.chnnames[1] == "test_chn2");
    REQUIRE(builder.nodes.size() == 4);
    REQUIRE(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[0].args[0].f64 == 134.1);
    REQUIRE(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    REQUIRE(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    REQUIRE(builder.nodes[1].args[0].f64 == 9123.45);
    REQUIRE(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[2].args[0].id == 2);
    REQUIRE(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[2].args[1].id == 0);
    REQUIRE(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    REQUIRE(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    REQUIRE(builder.nodes[3].args[0].id == 3);
    REQUIRE(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Arg);
    REQUIRE(builder.nodes[3].args[1].id == 1);
    auto &basicseq = builder.seqs.back();
    REQUIRE(basicseq.times.size() == 1);
    REQUIRE(basicseq.times[0].sign == Seq::Sign::Pos);
    REQUIRE(basicseq.times[0].id == 105);
    REQUIRE(basicseq.times[0].delta_node == 1);
    REQUIRE(basicseq.times[0].prev_id == 0);
    REQUIRE(basicseq.outputs.size() == 1);
    REQUIRE(basicseq.outputs[9].time_id == 1);
    REQUIRE(basicseq.outputs[9].len_node == 2);
    REQUIRE(basicseq.outputs[9].val_node == 4);
    REQUIRE(basicseq.outputs[9].cond_node == 0);
    REQUIRE(basicseq.outputs[9].chn == 2);
    REQUIRE(builder.noramp_chns.size() == 1);
    REQUIRE(builder.noramp_chns[0] == 2);

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

    seq.prepare();
    seq.optimize();
}
