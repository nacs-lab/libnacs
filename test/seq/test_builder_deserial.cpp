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
#include "../../lib/utils/streams.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

struct ostream : uvector_ostream {
    using uvector_ostream::write;
    template<typename T>
    void write(T v)
    {
        static_assert(std::is_trivial_v<T>);
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

static void test_channels(llvm::LLVMContext &llvm_ctx){

    Seq::Seq seq(llvm_ctx);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 2);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");

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

    stm.parse(builder);

    assert(builder.chnnames.size() == 1);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.defvals.size() == 1);
    assert(builder.defvals[1].is(IR::TagVal(2.3)));

    builder.buildseq();

    assert(seq.get_chn_id("test_chn1") == 1);
    auto defval = seq.defval(1);
    assert(defval == 2.3);
}

static void test_global(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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

    stm.parse(builder);

    assert(builder.slots.size() == 3);
    assert(builder.slots[0] == IR::Type::Float64);
    assert(builder.slots[1] == IR::Type::Bool);
    assert(builder.slots[2] == IR::Type::Int32);

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
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
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

    stm.parse(builder);

    assert(builder.nodes.size() == 2);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[0].id == 0);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[1].args[0].id == 1);
    assert(builder.slots.size() == 2);
    assert(builder.slots[0] == IR::Type::Float64);
    assert(builder.slots[1] == IR::Type::Int32);
    assert(builder.seqs.size() == 1);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.assignments.size() == 2);
    assert(basicseq.assignments[0].assign_id == 41);
    assert(basicseq.assignments[0].global_id == 0);
    assert(basicseq.assignments[0].val_node == 2);
    assert(basicseq.assignments[1].assign_id == 52);
    assert(basicseq.assignments[1].global_id == 1);
    assert(basicseq.assignments[1].val_node == 1);

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

    stm.parse(builder);

    assert(builder.nodes.size() == 8);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[0].id == 0);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[1].args[0].id == 1);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[0].id == 1);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[1].id == 2);
    assert(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[0].id == 1);
    assert(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[1].id == 2);
    assert(builder.nodes[4].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[4].args[0].id == 3);
    assert(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[4].args[1].id == 4);
    assert(builder.nodes[5].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[0].id == 4);
    assert(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[1].id == 3);
    assert(builder.nodes[6].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[6].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[6].args[0].id == 5);
    assert(builder.nodes[6].argtypes[1] == Seq::Builder::ArgType::ConstInt32);
    assert(builder.nodes[6].args[1].i32 == 8);
    assert(builder.nodes[7].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[7].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[7].args[0].id == 6);
    assert(builder.nodes[7].argtypes[1] == Seq::Builder::ArgType::ConstInt32);
    assert(builder.nodes[7].args[1].i32 == 8);
    assert(builder.slots.size() == 2);
    assert(builder.slots[0] == IR::Type::Float64);
    assert(builder.slots[1] == IR::Type::Int32);

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

    stm.parse(builder);

    assert(builder.nodes.size() == 6);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[0].id == 0);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[1].args[0].id == 1);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[0].id == 1);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[1].id == 2);
    assert(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[0].id == 2);
    assert(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[1].id == 1);
    assert(builder.nodes[4].op == Seq::Builder::OpCode::CmpLT);
    assert(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[4].args[0].id == 1);
    assert(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[4].args[1].id == 2);
    assert(builder.nodes[5].op == Seq::Builder::OpCode::CmpGT);
    assert(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[0].id == 2);
    assert(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[1].id == 1);
    assert(builder.slots.size() == 2);
    assert(builder.slots[0] == IR::Type::Float64);
    assert(builder.slots[1] == IR::Type::Int32);

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

    stm.parse(builder);

    assert(builder.nodes.size() == 1);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Interp);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[0].id == 0);
    assert(builder.nodes[0].argtypes[1] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[1].id == 1);
    assert(builder.nodes[0].argtypes[2] == Seq::Builder::ArgType::Global);
    assert(builder.nodes[0].args[2].id == 2);
    assert(builder.nodes[0].data_id == 0);
    assert(builder.slots.size() == 3);
    assert(builder.slots[0] == IR::Type::Bool);
    assert(builder.slots[1] == IR::Type::Float64);
    assert(builder.slots[2] == IR::Type::Int32);
    assert(builder.datas.size() == 1);
    double data[] = {1, 2, 5, 6};
    assert(builder.datas[0].size() == 4);
    assert(memcmp(builder.datas[0].data(), data, sizeof(data)) == 0);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(101);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(220);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>(Seq::EventTime::Unknown);
        stm.write<uint32_t>(299);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(3);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
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

    stm.parse(builder);

    assert(builder.nodes.size() == 2);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 4.2);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 1000);
    assert(builder.seqs.size() == 1);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 3);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 101);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 2);
    assert(basicseq.times[1].sign == Seq::EventTime::Pos);
    assert(basicseq.times[1].id == 220);
    assert(basicseq.times[1].delta_node == 2);
    assert(basicseq.times[1].prev_id == 0);
    assert(basicseq.times[2].sign == Seq::EventTime::Unknown);
    assert(basicseq.times[2].id == 299);
    assert(basicseq.times[2].delta_node == 1);
    assert(basicseq.times[2].prev_id == 1);
    assert(basicseq.endtimes.size() == 2);
    assert(basicseq.endtimes[0] == 1);
    assert(basicseq.endtimes[1] == 3);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(101);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint8_t>(Seq::EventTime::NonNeg);
        stm.write<uint32_t>(220);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 3);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");
    assert(builder.chnnames[2] == "test_chn3");
    assert(builder.nodes.size() == 2);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 9.2);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 1020);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 2);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 101);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 2);
    assert(basicseq.times[1].sign == Seq::EventTime::NonNeg);
    assert(basicseq.times[1].id == 220);
    assert(basicseq.times[1].delta_node == 2);
    assert(basicseq.times[1].prev_id == 0);
    assert(basicseq.measures.size() == 1);
    assert(basicseq.measures[12].time_id == 1);
    assert(basicseq.measures[12].chn == 3);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(29);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 2);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");
    assert(builder.nodes.size() == 2);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 134.1);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::Measure);
    assert(builder.nodes[1].args[0].id == 19);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 1);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 105);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 0);
    assert(basicseq.measures.size() == 1);
    assert(basicseq.measures[19].time_id == 1);
    assert(basicseq.measures[19].chn == 2);
    assert(basicseq.outputs.size() == 1);
    assert(basicseq.outputs[29].time_id == 1);
    assert(basicseq.outputs[29].len_node == 0);
    assert(basicseq.outputs[29].val_node == 2);
    assert(basicseq.outputs[29].chn == 1);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(4);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // 2
        stm.write<uint32_t>(65);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(2);
        // 3
        stm.write<uint32_t>(32);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(3);
        // 4
        stm.write<uint32_t>(49);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(6);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 4);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");
    assert(builder.chnnames[2] == "test_chn3");
    assert(builder.chnnames[3] == "test_chn4");
    assert(builder.nodes.size() == 6);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 134.1);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 9123.45);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[2].args[0].f64 == 5);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[2].args[1].id == 1);
    assert(builder.nodes[3].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[0].id == 2);
    assert(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[1].id == 3);
    assert(builder.nodes[4].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[4].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[4].args[0].id == 1);
    assert(builder.nodes[4].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[4].args[1].id == 0);
    assert(builder.nodes[5].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[5].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[0].id == 4);
    assert(builder.nodes[5].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[5].args[1].id == 5);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 1);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 105);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 0);
    assert(basicseq.outputs.size() == 4);
    assert(basicseq.outputs[8].time_id == 1);
    assert(basicseq.outputs[8].len_node == 2);
    assert(basicseq.outputs[8].val_node == 1);
    assert(basicseq.outputs[8].chn == 1);
    assert(basicseq.outputs[65].time_id == 1);
    assert(basicseq.outputs[65].len_node == 2);
    assert(basicseq.outputs[65].val_node == 4);
    assert(basicseq.outputs[65].chn == 2);
    assert(basicseq.outputs[32].time_id == 1);
    assert(basicseq.outputs[32].len_node == 2);
    assert(basicseq.outputs[32].val_node == 5);
    assert(basicseq.outputs[32].chn == 3);
    assert(basicseq.outputs[49].time_id == 1);
    assert(basicseq.outputs[49].len_node == 2);
    assert(basicseq.outputs[49].val_node == 6);
    assert(basicseq.outputs[49].chn == 4);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // 2
        stm.write<uint32_t>(65);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(4);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 2);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");
    assert(builder.nodes.size() == 4);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 134.1);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 9123.45);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[2].args[0].f64 == 5);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[2].args[1].id == 1);
    assert(builder.nodes[3].op == Seq::Builder::OpCode::Mul);
    assert(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[0].id == 2);
    assert(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[1].id == 3);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 1);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 105);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 0);
    assert(basicseq.outputs.size() == 2);
    assert(basicseq.outputs[8].time_id == 1);
    assert(basicseq.outputs[8].len_node == 0);
    assert(basicseq.outputs[8].val_node == 1);
    assert(basicseq.outputs[8].chn == 1);
    assert(basicseq.outputs[65].time_id == 1);
    assert(basicseq.outputs[65].len_node == 0);
    assert(basicseq.outputs[65].val_node == 4);
    assert(basicseq.outputs[65].chn == 2);

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
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
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
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
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

    stm.parse(builder);

    assert(builder.nodes.size() == 2);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 4.2);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 1000);
    auto &basicseq1 = builder.seqs[0];
    auto &basicseq2 = builder.seqs[1];
    assert(basicseq1.default_target == 2);
    assert(basicseq1.branches.size() == 1);
    assert(basicseq1.branches[0].branch_id == 102);
    assert(basicseq1.branches[0].target_id == 0);
    assert(basicseq1.branches[0].cond_node == 1);
    assert(basicseq2.default_target == 0);
    assert(basicseq2.branches.size() == 1);
    assert(basicseq2.branches[0].branch_id == 103);
    assert(basicseq2.branches[0].target_id == 2);
    assert(basicseq2.branches[0].cond_node == 2);

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

static void test_output_noramp_noold(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        // 1
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(3);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 1);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.nodes.size() == 3);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 134.1);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 9123.45);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[0].id == 2);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[2].args[1].id == 0);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 1);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 105);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 0);
    assert(basicseq.outputs.size() == 1);
    assert(basicseq.outputs[8].time_id == 1);
    assert(basicseq.outputs[8].len_node == 2);
    assert(basicseq.outputs[8].val_node == 3);
    assert(basicseq.outputs[8].chn == 1);
    assert(builder.noramp_chns.size() == 1);
    assert(builder.noramp_chns[0] == 1);

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
        stm.write<uint8_t>(Seq::EventTime::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        // 1
        stm.write<uint32_t>(9);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(4);
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

    stm.parse(builder);

    assert(builder.chnnames.size() == 2);
    assert(builder.chnnames[0] == "test_chn1");
    assert(builder.chnnames[1] == "test_chn2");
    assert(builder.nodes.size() == 4);
    assert(builder.nodes[0].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[0].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[0].args[0].f64 == 134.1);
    assert(builder.nodes[1].op == Seq::Builder::OpCode::Identity);
    assert(builder.nodes[1].argtypes[0] == Seq::Builder::ArgType::ConstFloat64);
    assert(builder.nodes[1].args[0].f64 == 9123.45);
    assert(builder.nodes[2].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[2].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[2].args[0].id == 2);
    assert(builder.nodes[2].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[2].args[1].id == 0);
    assert(builder.nodes[3].op == Seq::Builder::OpCode::Add);
    assert(builder.nodes[3].argtypes[0] == Seq::Builder::ArgType::Node);
    assert(builder.nodes[3].args[0].id == 3);
    assert(builder.nodes[3].argtypes[1] == Seq::Builder::ArgType::Arg);
    assert(builder.nodes[3].args[1].id == 1);
    auto &basicseq = builder.seqs.back();
    assert(basicseq.times.size() == 1);
    assert(basicseq.times[0].sign == Seq::EventTime::Pos);
    assert(basicseq.times[0].id == 105);
    assert(basicseq.times[0].delta_node == 1);
    assert(basicseq.times[0].prev_id == 0);
    assert(basicseq.outputs.size() == 1);
    assert(basicseq.outputs[9].time_id == 1);
    assert(basicseq.outputs[9].len_node == 2);
    assert(basicseq.outputs[9].val_node == 4);
    assert(basicseq.outputs[9].chn == 2);
    assert(builder.noramp_chns.size() == 1);
    assert(builder.noramp_chns[0] == 2);

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
    test_output_noramp_noold(*llvm_ctx);
    test_output_noramp_old(*llvm_ctx);

    return 0;
}
