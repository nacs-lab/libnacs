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
#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/streams.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

using namespace NaCs;

template<typename T>
static bool compare_vector_sorted(const std::vector<T> &v1, const std::vector<T> &v2)
{
    if (v1.size() != v2.size())
        return false;
    auto s1 = v1;
    auto s2 = v2;
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());
    return s1 == s2;
}

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

namespace NaCs::Seq {

static bool operator==(const HostSeq::Assignment &a1, const HostSeq::Assignment &a2)
{
    return std::make_pair(a1.global_id, a1.value) == std::make_pair(a2.global_id, a2.value);
}

static bool operator<(const HostSeq::Assignment &a1, const HostSeq::Assignment &a2)
{
    return std::make_pair(a1.global_id, a1.value) < std::make_pair(a2.global_id, a2.value);
}

static bool operator==(const HostSeq::Assumption &a1, const HostSeq::Assumption &a2)
{
    return std::make_tuple(a1.sign, a1.value, a1.id) ==
        std::make_tuple(a2.sign, a2.value, a2.id);
}

}

static auto llvm_ctx = LLVM::new_context();

TEST_CASE("const_merge") {
    // Make sure the constants are allocated correctly.
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(0)),
                             seq.get_const(IR::TagVal(1.0)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
                                 builder.createRet(builder.createMul(0, builder.getConst(3.0)));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0)}, 1);
                             }());
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto p2 = bs2->add_pulse(1, 2, bs2->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_const(IR::TagVal(3.0)));
    REQUIRE(p2->id() == 2);
    REQUIRE(!p2->needs_oldval());
    bs1->set_default_branch(bs2);

    seq.prepare();
    seq.optimize();

    // The event time tracking is basic sequence specific
    // so the two times that belongs to two basic sequences can never have the same addresses.
    REQUIRE(&p1->start() != &p2->start());
    REQUIRE(p1->start() == p2->start());
    REQUIRE(p1->start().terms.empty());

    REQUIRE(p1->endval()->is_const());
    REQUIRE(p1->endval()->get_const().is(IR::TagVal(3.0)));
    REQUIRE(p2->val()->is_const());
    REQUIRE(p2->val()->get_const().is(IR::TagVal(3.0)));

    // The endval is produced by the optimizer and right now we do not merge
    // constants at this level.
    REQUIRE(p1->endval() != p2->val());
    REQUIRE(p2->endval() == p2->val());

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 5);
    REQUIRE(cs.host_seq.nglobals == 0);
    REQUIRE(cs.host_seq.npublic_globals == 0);
    REQUIRE(cs.host_seq.nglobal_vals == 0);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 6);

    REQUIRE(cs.host_seq.values.size() == 6);
    REQUIRE(cs.host_seq.types.size() == 5);
    REQUIRE(cs.host_seq.global_evals.size() == 0);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    auto slot_f1 = cs.compiler.var_slot(p1->len());
    REQUIRE(slot_f1 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1].f64 == 1.0);

    REQUIRE(p1->val() != p1->endval());
    REQUIRE(p1->val()->is_call());
    auto slot_p1val = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_p1val < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_p1val] == Seq::HostSeq::Type::Float64);

    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    auto slot_f3 = cs.compiler.var_slot(p1->endval());
    REQUIRE(slot_f3 < cs.host_seq.nconsts);
    REQUIRE(cs.compiler.var_slot(p2->val()) == slot_f3);
    REQUIRE(cs.compiler.var_slot(bs2->startval(1)) == slot_f3);
    REQUIRE(cs.host_seq.types[slot_f3] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f3].f64 == 3.0);

    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.compiler.time_slot(&p2->start()) == slot_i0);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    REQUIRE(cs.host_seq.seqs[0].id == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses.size() == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].id == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].time == slot_i0);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].measure == uint32_t(-1));
    REQUIRE(cs.host_seq.seqs[0].pulses[0].len == slot_f1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].value == slot_p1val);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].chn == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].ramp_func);
    for (int i = 0; i < 10; i++) {
        auto t = i / 10.0;
        REQUIRE(cs.host_seq.seqs[0].pulses[0].ramp_func(t, nullptr) == t * 3);
    }
    REQUIRE(cs.host_seq.seqs[0].pulses[0].endvalue == slot_f3);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].cond == uint32_t(-1));

    REQUIRE(cs.host_seq.seqs[1].id == 2);
    REQUIRE(cs.host_seq.seqs[1].pulses.size() == 1);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].id == 2);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].time == slot_i0);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].measure == uint32_t(-1));
    REQUIRE(cs.host_seq.seqs[1].pulses[0].len == uint32_t(-1));
    REQUIRE(cs.host_seq.seqs[1].pulses[0].value == slot_f3);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].chn == 1);
    REQUIRE(!cs.host_seq.seqs[1].pulses[0].ramp_func);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].endvalue == slot_f3);
    REQUIRE(cs.host_seq.seqs[1].pulses[0].cond == uint32_t(-1));
}

TEST_CASE("global_value") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    Seq::EventTime t;
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto p1 = bs1->add_pulse(1, 4, bs1->track_time(t), seq.get_slot(IR::Type::Float64, 1),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto v = builder.createMul(0, builder.getConst(2.0));
                                 builder.createRet(builder.createAdd(v, 1));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0),
                                                      Seq::Arg::create_var(g0)}, 1);
                             }());
    t.terms.clear();

    seq.prepare();
    seq.optimize();

    REQUIRE(p1->endval()->is_call());
    REQUIRE(p1->val()->is_call());

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 1);
    REQUIRE(cs.host_seq.nglobals == 2);
    REQUIRE(cs.host_seq.npublic_globals == 2);
    REQUIRE(cs.host_seq.nglobal_vals == 3);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 7);

    REQUIRE(cs.host_seq.values.size() == 7);
    REQUIRE(cs.host_seq.depends.size() == 3);
    REQUIRE(cs.host_seq.types.size() == 6);
    REQUIRE(cs.host_seq.global_evals.size() == 3);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == 1);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(slot_g1 == 2);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_p1val = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_p1val >= gv_offset);
    REQUIRE(slot_p1val < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1val] == Seq::HostSeq::Type::Float64);
    // This is `dummy_eval_func` and it can take a null pointer
    cs.host_seq.global_evals[slot_p1val - gv_offset](nullptr);
    REQUIRE(cs.host_seq.depends[slot_p1val - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1val - gv_offset][0] == slot_g0 - gl_offset);

    auto slot_p1endval = cs.compiler.var_slot(p1->endval());
    REQUIRE(slot_p1endval >= gv_offset);
    REQUIRE(slot_p1endval < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1endval] == Seq::HostSeq::Type::Float64);
    for (int i = 0; i < 11; i++) {
        double g0 = (i - 5) * 0.5;
        cs.host_seq.values[slot_g0].f64 = g0;
        for (int j = 0; j < 11; j++) {
            double g1 = (j - 5) * 0.5;
            cs.host_seq.values[slot_g1].f64 = g1;
            cs.host_seq.global_evals[slot_p1endval - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p1endval].f64 == g0 + g1 * 2);
        }
    }
    REQUIRE(cs.host_seq.depends[slot_p1endval - gv_offset].size() == 2);
    REQUIRE(cs.host_seq.depends[slot_p1endval - gv_offset][0] == slot_g0 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_p1endval - gv_offset][1] == slot_g1 - gl_offset);

    auto slot_p1time = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_p1time >= gv_offset);
    REQUIRE(slot_p1time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g1 = i * 0.4;
        cs.host_seq.values[slot_g1].f64 = g1;
        cs.host_seq.global_evals[slot_p1time - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p1time].i64 == round<int64_t>(g1));
    }
    REQUIRE(cs.host_seq.depends[slot_p1time - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1time - gv_offset][0] == slot_g1 - gl_offset);

    REQUIRE(cs.host_seq.seqs[0].id == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses.size() == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].id == 4);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].time == slot_p1time);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].measure == uint32_t(-1));
    REQUIRE(cs.host_seq.seqs[0].pulses[0].len == slot_g1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].value == slot_p1val);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].chn == 1);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].ramp_func);
    for (int i = 0; i < 11; i++) {
        double g0 = (i - 5) * 0.5;
        cs.host_seq.values[slot_g0].f64 = g0;
        for (int j = 0; j < 10; j++) {
            auto t = j / 10.0;
            auto ramp_func = cs.host_seq.seqs[0].pulses[0].ramp_func;
            REQUIRE(ramp_func(t, cs.host_seq.values.data()) == g0 + t * 2);
        }
    }
    REQUIRE(cs.host_seq.seqs[0].pulses[0].endvalue == slot_p1endval);
    REQUIRE(cs.host_seq.seqs[0].pulses[0].cond == uint32_t(-1));
}

TEST_CASE("bseq_value") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);
    auto g3 = seq.get_slot(IR::Type::Float64, 3);

    auto p1 = bs1->add_pulse(1, 2, bs1->track_time(Seq::EventTime(0)), g1,
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64,
                                                      IR::Type::Float64});
                                 auto v = builder.createMul(0, 2);
                                 builder.createRet(builder.createAdd(v, 1));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0),
                                                      Seq::Arg::create_arg(1),
                                                      Seq::Arg::create_var(g0)}, 2);
                             }());
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, g2);
    auto p2 = bs1->add_pulse(1, 3, bs1->track_time(t), seq.get_const(IR::TagVal(1.0)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 builder.createRet(builder.createSub(1, 0));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0),
                                                      Seq::Arg::create_arg(1)}, 2);
                             }());
    Seq::EventTime t2;
    t2.add_term(Seq::Sign::Pos, g3);
    auto p3 = bs1->add_pulse(1, 4, bs1->track_time(t2), seq.get_const(IR::TagVal(1.5)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 builder.createRet(builder.createAdd(1, 0));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0),
                                                      Seq::Arg::create_arg(1)}, 2);
                             }());
    bs1->add_branch([&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, builder.getConstFloat(1.0)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g3)}, 0);
    }(), bs1, 8);
    t.terms.clear();
    t2.terms.clear();

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].id == 8);
    REQUIRE(bs1->get_branches()[0].target == bs1);

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 3);
    REQUIRE(cs.host_seq.nglobals == 4);
    REQUIRE(cs.host_seq.npublic_globals == 4);
    REQUIRE(cs.host_seq.nglobal_vals == 3);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 11);

    REQUIRE(cs.host_seq.values.size() == 19);
    REQUIRE(cs.host_seq.depends.size() == 3);
    REQUIRE(cs.host_seq.types.size() == 10);
    REQUIRE(cs.host_seq.global_evals.size() == 3);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_f1 = cs.compiler.var_slot(p2->len());
    REQUIRE(slot_f1 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1].f64 == 1);

    auto slot_f1_5 = cs.compiler.var_slot(p3->len());
    REQUIRE(slot_f1_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_5].f64 == 1.5);

    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);
    auto slot_g3 = cs.compiler.var_slot(g3);
    REQUIRE(slot_g3 == gl_offset + 3);
    REQUIRE(cs.host_seq.types[slot_g3] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_p2time = cs.compiler.time_slot(&p2->start());
    REQUIRE(slot_p2time >= gv_offset);
    REQUIRE(slot_p2time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p2time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g2 = i * 0.4;
        cs.host_seq.values[slot_g2].f64 = g2;
        cs.host_seq.global_evals[slot_p2time - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p2time].i64 == round<int64_t>(g2));
    }
    REQUIRE(cs.host_seq.depends[slot_p2time - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p2time - gv_offset][0] == slot_g2 - gl_offset);

    auto slot_p3time = cs.compiler.time_slot(&p3->start());
    REQUIRE(slot_p3time >= gv_offset);
    REQUIRE(slot_p3time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p3time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g3 = i * 0.4;
        cs.host_seq.values[slot_g3].f64 = g3;
        cs.host_seq.global_evals[slot_p3time - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p3time].i64 == round<int64_t>(g3));
    }
    REQUIRE(cs.host_seq.depends[slot_p3time - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p3time - gv_offset][0] == slot_g3 - gl_offset);

    auto slot_brcond = cs.compiler.var_slot(bs1->get_branches()[0].cond.get());
    REQUIRE(slot_brcond >= gv_offset);
    REQUIRE(slot_brcond < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_brcond] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g3 = i * 0.4 - 1.1;
        cs.host_seq.values[slot_g3].f64 = g3;
        cs.host_seq.global_evals[slot_brcond - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_brcond].b == (g3 < 1.0));
    }
    REQUIRE(cs.host_seq.depends[slot_brcond - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_brcond - gv_offset][0] == slot_g3 - gl_offset);

    // Channel:
    auto chn_offset = cs.host_seq.nconsts + cs.host_seq.nglobals + cs.host_seq.nglobal_vals;
    auto slot_ch1 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_ch1 >= chn_offset);
    REQUIRE(slot_ch1 < chn_offset + cs.host_seq.nchannels);

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);
        auto &info1 = bs1->get_channel_infos().find(1)->second;

        REQUIRE(host_bseq->nmeasure == 2);
        REQUIRE(host_bseq->ndirect == 2);
        REQUIRE(host_bseq->nneed_order == 4);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_p2time - gv_offset, slot_p3time - gv_offset}));

        REQUIRE(std::is_sorted(host_bseq->cond_global_refs.begin(),
                               host_bseq->cond_global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->cond_global_refs, {slot_brcond - gv_offset}));

        REQUIRE(host_bseq->types.size() == 6);
        REQUIRE(host_bseq->reverse_depends.size() == 2);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        REQUIRE(info1.pulse_start_vars.find(p2) != info1.pulse_start_vars.end());
        auto p2oldval = info1.pulse_start_vars.find(p2)->second.get();
        auto slot_p2oldval = cs.compiler.var_slot(p2oldval);
        REQUIRE(slot_p2oldval >= mea_offset);
        REQUIRE(slot_p2oldval < mea_offset + host_bseq->nmeasure);

        REQUIRE(info1.pulse_start_vars.find(p3) != info1.pulse_start_vars.end());
        auto p3oldval = info1.pulse_start_vars.find(p3)->second.get();
        auto slot_p3oldval = cs.compiler.var_slot(p3oldval);
        REQUIRE(slot_p3oldval >= mea_offset);
        REQUIRE(slot_p3oldval < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;
        auto slot_p1val = cs.compiler.var_slot(p1->val());
        REQUIRE(slot_p1val >= dir_offset);
        REQUIRE(slot_p1val < dir_offset + host_bseq->ndirect);
        REQUIRE(host_bseq->types[slot_p1val - dir_offset] == Seq::HostSeq::Type::Float64);
        // This is `dummy_eval_func` and it can take a null pointer
        host_bseq->evals[slot_p1val - dir_offset](nullptr);

        auto slot_p1endval = cs.compiler.var_slot(p1->endval());
        REQUIRE(slot_p1endval >= dir_offset);
        REQUIRE(slot_p1endval < dir_offset + host_bseq->ndirect);
        REQUIRE(host_bseq->types[slot_p1endval - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double g0 = i * 0.3 - 1.2;
            cs.host_seq.values[slot_g0].f64 = g0;
            for (int j = 0; j < 11; j++) {
                double g1 = j * 0.3 + 0.1;
                cs.host_seq.values[slot_g1].f64 = g1;
                for (int k = 0; k < 11; k++) {
                    double ch1 = k * 0.5 - 0.9;
                    cs.host_seq.values[slot_ch1].f64 = ch1;
                    host_bseq->evals[slot_p1endval - dir_offset](cs.host_seq.values.data());
                    REQUIRE(approx(cs.host_seq.values[slot_p1endval].f64, g1 * g0 + ch1));
                }
            }
        }

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        auto slot_p2val = cs.compiler.var_slot(p2->val());
        REQUIRE(slot_p2val >= no_offset);
        REQUIRE(slot_p2val < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p2val - dir_offset] == Seq::HostSeq::Type::Float64);
        // This is `dummy_eval_func` and it can take a null pointer
        host_bseq->evals[slot_p2val - dir_offset](nullptr);
        REQUIRE(host_bseq->deps_count[slot_p2val - no_offset] == 1);

        auto slot_p2endval = cs.compiler.var_slot(p2->endval());
        REQUIRE(slot_p2endval >= no_offset);
        REQUIRE(slot_p2endval < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p2endval - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double p2old = i * 0.4 - 1.9;
            cs.host_seq.values[slot_p2oldval].f64 = p2old;
            host_bseq->evals[slot_p2endval - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p2endval].f64 == p2old - 1);
        }
        REQUIRE(host_bseq->deps_count[slot_p2endval - no_offset] == 1);

        auto slot_p3val = cs.compiler.var_slot(p3->val());
        REQUIRE(slot_p3val >= no_offset);
        REQUIRE(slot_p3val < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p3val - dir_offset] == Seq::HostSeq::Type::Float64);
        // This is `dummy_eval_func` and it can take a null pointer
        host_bseq->evals[slot_p3val - dir_offset](nullptr);
        REQUIRE(host_bseq->deps_count[slot_p3val - no_offset] == 1);

        auto slot_p3endval = cs.compiler.var_slot(p3->endval());
        REQUIRE(slot_p3endval >= no_offset);
        REQUIRE(slot_p3endval < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p3endval - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double p3old = i * 0.4 - 1.9;
            cs.host_seq.values[slot_p3oldval].f64 = p3old;
            host_bseq->evals[slot_p3endval - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p3endval].f64 == p3old + 1.5);
        }
        REQUIRE(host_bseq->deps_count[slot_p3endval - no_offset] == 1);

        REQUIRE(host_bseq->assumptions_idx ==
                (std::vector<uint32_t>{uint32_t(-1), uint32_t(-1), uint32_t(-1), uint32_t(-1)}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_p2oldval - mea_offset],
                                      {slot_p2val - no_offset, slot_p2endval - no_offset}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_p3oldval - mea_offset],
                                      {slot_p3val - no_offset, slot_p3endval - no_offset}));

        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.empty());

        REQUIRE(host_bseq->branches.size() == 1);
        REQUIRE(host_bseq->branches[0].cond == slot_brcond);
        REQUIRE(host_bseq->branches[0].target == 0);
        REQUIRE(host_bseq->default_branch == uint32_t(-1));

        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 3);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 2);
        REQUIRE(host_p1->time == slot_i0);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == slot_g1);
        REQUIRE(host_p1->value == slot_p1val);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(host_p1->ramp_func);
        for (int i = 0; i < 11; i++) {
            double g0 = i * 0.3 - 1.2;
            cs.host_seq.values[slot_g0].f64 = g0;
            for (int j = 0; j < 11; j++) {
                double t = j * 0.3 + 0.1;
                for (int k = 0; k < 11; k++) {
                    double ch1 = k * 0.5 - 0.9;
                    cs.host_seq.values[slot_ch1].f64 = ch1;
                    REQUIRE(approx(host_p1->ramp_func(t, cs.host_seq.values.data()),
                                   t * g0 + ch1));
                }
            }
        }
        REQUIRE(host_p1->endvalue == slot_p1endval);
        REQUIRE(host_p1->cond == uint32_t(-1));

        auto host_p2 = pulses[1];
        REQUIRE(host_p2->id == 3);
        REQUIRE(host_p2->time == slot_p2time);
        REQUIRE(host_p2->measure == slot_p2oldval - mea_offset);
        REQUIRE(host_p2->len == slot_f1);
        REQUIRE(host_p2->value == slot_p2val);
        REQUIRE(host_p2->chn == 1);
        REQUIRE(host_p2->ramp_func);
        for (int i = 0; i < 11; i++) {
            double p2old = i * 0.3 - 1.2;
            cs.host_seq.values[slot_p2oldval].f64 = p2old;
            for (int j = 0; j < 11; j++) {
                double t = j * 0.3 + 0.1;
                REQUIRE(host_p2->ramp_func(t, cs.host_seq.values.data()) == p2old - t);
            }
        }
        REQUIRE(host_p2->endvalue == slot_p2endval);
        REQUIRE(host_p2->cond == uint32_t(-1));

        auto host_p3 = pulses[2];
        REQUIRE(host_p3->id == 4);
        REQUIRE(host_p3->time == slot_p3time);
        REQUIRE(host_p3->measure == slot_p3oldval - mea_offset);
        REQUIRE(host_p3->len == slot_f1_5);
        REQUIRE(host_p3->value == slot_p3val);
        REQUIRE(host_p3->chn == 1);
        REQUIRE(host_p3->ramp_func);
        for (int i = 0; i < 11; i++) {
            double p3old = i * 0.3 - 1.2;
            cs.host_seq.values[slot_p3oldval].f64 = p3old;
            for (int j = 0; j < 11; j++) {
                double t = j * 0.3 + 0.1;
                REQUIRE(host_p3->ramp_func(t, cs.host_seq.values.data()) == p3old + t);
            }
        }
        REQUIRE(host_p3->endvalue == slot_p3endval);
        REQUIRE(host_p3->cond == uint32_t(-1));
    }
}

TEST_CASE("cond_pulse_and_measure") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);

    auto p1cond = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, builder.getConstFloat(0)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g0)}, 0);
    }();

    auto p1 = bs1->add_pulse(1, 7, bs1->track_time(Seq::EventTime(0)), g1,
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64,
                                                      IR::Type::Float64});
                                 auto v = builder.createMul(0, 2);
                                 builder.createRet(builder.createSub(v, 1));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_arg(0),
                                                      Seq::Arg::create_arg(1),
                                                      Seq::Arg::create_var(g0)}, 2);
                             }(), p1cond);
    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(5)),
                               bs1->new_measure(seq.env(), 1));
    auto p2 = bs1->add_pulse(1, 12, bs1->track_time(Seq::EventTime(10)), nullptr,
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 builder.createRet(builder.createSub(1, 0));
                                 return seq.get_call(builder.get(),
                                                     {Seq::Arg::create_var(m1->val()),
                                                      Seq::Arg::create_arg(1)}, 2);
                             }());
    bs1->add_branch([&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, builder.getConstFloat(-1.0)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g2)}, 0);
    }(), nullptr, 87);
    bs1->set_default_branch(bs1);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].id == 87);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_default_branch() == bs1);

    REQUIRE(p2->val() == p2->endval());

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 3);
    REQUIRE(cs.host_seq.nglobals == 3);
    REQUIRE(cs.host_seq.npublic_globals == 3);
    REQUIRE(cs.host_seq.nglobal_vals == 2);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 9);

    REQUIRE(cs.host_seq.values.size() == 14);
    REQUIRE(cs.host_seq.depends.size() == 2);
    REQUIRE(cs.host_seq.types.size() == 8);
    REQUIRE(cs.host_seq.global_evals.size() == 2);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    auto slot_i5 = cs.compiler.time_slot(&m1->start());
    REQUIRE(slot_i5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i5] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i5].i64 == 5);

    auto slot_i10 = cs.compiler.time_slot(&p2->start());
    REQUIRE(slot_i10 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i10] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i10].i64 == 10);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_p1cond = cs.compiler.var_slot(p1cond);
    REQUIRE(slot_p1cond >= gv_offset);
    REQUIRE(slot_p1cond < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1cond] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4 - 1.1;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_p1cond - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p1cond].b == (g0 < 0));
    }
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset][0] == slot_g0 - gl_offset);

    auto slot_brcond = cs.compiler.var_slot(bs1->get_branches()[0].cond.get());
    REQUIRE(slot_brcond >= gv_offset);
    REQUIRE(slot_brcond < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_brcond] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g2 = i * 0.4 - 1.1;
        cs.host_seq.values[slot_g2].f64 = g2;
        cs.host_seq.global_evals[slot_brcond - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_brcond].b == (g2 < -1.0));
    }
    REQUIRE(cs.host_seq.depends[slot_brcond - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_brcond - gv_offset][0] == slot_g2 - gl_offset);

    // Channel:
    auto chn_offset = cs.host_seq.nconsts + cs.host_seq.nglobals + cs.host_seq.nglobal_vals;
    auto slot_ch1 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_ch1 >= chn_offset);
    REQUIRE(slot_ch1 < chn_offset + cs.host_seq.nchannels);

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);
        auto &info1 = bs1->get_channel_infos().find(1)->second;

        REQUIRE(host_bseq->nmeasure == 2);
        REQUIRE(host_bseq->ndirect == 2);
        REQUIRE(host_bseq->nneed_order == 1);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs, {slot_p1cond - gv_offset}));

        REQUIRE(std::is_sorted(host_bseq->cond_global_refs.begin(),
                               host_bseq->cond_global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->cond_global_refs, {slot_brcond - gv_offset}));

        REQUIRE(host_bseq->types.size() == 3);
        REQUIRE(host_bseq->reverse_depends.size() == 2);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        REQUIRE(info1.pulse_start_vars.find(p2) != info1.pulse_start_vars.end());
        auto p2oldval = info1.pulse_start_vars.find(p2)->second.get();
        auto slot_p2oldval = cs.compiler.var_slot(p2oldval);
        REQUIRE(slot_p2oldval >= mea_offset);
        REQUIRE(slot_p2oldval < mea_offset + host_bseq->nmeasure);

        auto slot_m1val = cs.compiler.var_slot(m1->val());
        REQUIRE(slot_m1val >= mea_offset);
        REQUIRE(slot_m1val < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;
        auto slot_p1val = cs.compiler.var_slot(p1->val());
        REQUIRE(slot_p1val >= dir_offset);
        REQUIRE(slot_p1val < dir_offset + host_bseq->ndirect);
        REQUIRE(host_bseq->types[slot_p1val - dir_offset] == Seq::HostSeq::Type::Float64);
        // This is `dummy_eval_func` and it can take a null pointer
        host_bseq->evals[slot_p1val - dir_offset](nullptr);

        auto slot_p1endval = cs.compiler.var_slot(p1->endval());
        REQUIRE(slot_p1endval >= dir_offset);
        REQUIRE(slot_p1endval < dir_offset + host_bseq->ndirect);
        REQUIRE(host_bseq->types[slot_p1endval - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double g0 = i * 0.3 - 1.2;
            cs.host_seq.values[slot_g0].f64 = g0;
            for (int j = 0; j < 11; j++) {
                double g1 = j * 0.3 + 0.1;
                cs.host_seq.values[slot_g1].f64 = g1;
                for (int k = 0; k < 11; k++) {
                    double ch1 = k * 0.5 - 0.9;
                    cs.host_seq.values[slot_ch1].f64 = ch1;
                    host_bseq->evals[slot_p1endval - dir_offset](cs.host_seq.values.data());
                    REQUIRE(approx(cs.host_seq.values[slot_p1endval].f64, g1 * g0 - ch1));
                }
            }
        }

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        auto slot_p2val = cs.compiler.var_slot(p2->val());
        REQUIRE(slot_p2val >= no_offset);
        REQUIRE(slot_p2val < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p2val - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double p2old = i * 0.4 - 1.9;
            cs.host_seq.values[slot_p2oldval].f64 = p2old;
            for (int j = 0; j < 11; j++) {
                double m1 = j * 0.6 - 2.3;
                cs.host_seq.values[slot_m1val].f64 = m1;
                host_bseq->evals[slot_p2val - dir_offset](cs.host_seq.values.data());
                REQUIRE(cs.host_seq.values[slot_p2val].f64 == p2old - m1);
            }
        }
        REQUIRE(host_bseq->deps_count[slot_p2val - no_offset] == 2);

        REQUIRE(host_bseq->assumptions_idx == (std::vector<uint32_t>{uint32_t(-1)}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_m1val - mea_offset],
                                      {slot_p2val - no_offset}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_p2oldval - mea_offset],
                                      {slot_p2val - no_offset}));

        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.empty());

        REQUIRE(host_bseq->branches.size() == 1);
        REQUIRE(host_bseq->branches[0].cond == slot_brcond);
        REQUIRE(host_bseq->branches[0].target == uint32_t(-1));
        REQUIRE(host_bseq->default_branch == 0);

        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 3);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 7);
        REQUIRE(host_p1->time == slot_i0);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == slot_g1);
        REQUIRE(host_p1->value == slot_p1val);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(host_p1->ramp_func);
        for (int i = 0; i < 11; i++) {
            double g0 = i * 0.3 - 1.2;
            cs.host_seq.values[slot_g0].f64 = g0;
            for (int j = 0; j < 11; j++) {
                double t = j * 0.3 + 0.1;
                for (int k = 0; k < 11; k++) {
                    double ch1 = k * 0.5 - 0.9;
                    cs.host_seq.values[slot_ch1].f64 = ch1;
                    REQUIRE(approx(host_p1->ramp_func(t, cs.host_seq.values.data()),
                                   t * g0 - ch1));
                }
            }
        }
        REQUIRE(host_p1->endvalue == slot_p1endval);
        REQUIRE(host_p1->cond == slot_p1cond);

        auto host_m1 = pulses[1];
        REQUIRE(host_m1->id == 10);
        REQUIRE(host_m1->time == slot_i5);
        REQUIRE(host_m1->measure == slot_m1val - mea_offset);
        REQUIRE(host_m1->len == uint32_t(-2));
        REQUIRE(host_m1->value == uint32_t(-1));
        REQUIRE(host_m1->chn == 1);
        REQUIRE(!host_m1->ramp_func);
        REQUIRE(host_m1->endvalue == uint32_t(-1));
        REQUIRE(host_m1->cond == uint32_t(-1));

        auto host_p2 = pulses[2];
        REQUIRE(host_p2->id == 12);
        REQUIRE(host_p2->time == slot_i10);
        REQUIRE(host_p2->measure == slot_p2oldval - mea_offset);
        REQUIRE(host_p2->len == uint32_t(-1));
        REQUIRE(host_p2->value == slot_p2val);
        REQUIRE(host_p2->chn == 1);
        REQUIRE(!host_p2->ramp_func);
        REQUIRE(host_p2->endvalue == slot_p2val);
        REQUIRE(host_p2->cond == uint32_t(-1));
    }
}

TEST_CASE("bseq_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);

    auto m1 = bs1->add_measure(1, 2, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    Seq::EventTime t1(2);
    t1.add_term(Seq::Sign::NonNeg, m1->val());
    auto p1 = bs1->add_pulse(1, 5, bs1->track_time(t1), nullptr,
                             seq.get_const(IR::TagVal(1.5)));
    auto m2 = bs1->add_measure(1, 9, bs1->track_time(Seq::EventTime(5)),
                               bs1->new_measure(seq.env(), 2));
    Seq::EventTime t2(10);
    t2.add_term(Seq::Sign::NonNeg, m2->val());
    auto p2 = bs1->add_pulse(1, 10, bs1->track_time(t2), nullptr,
                             seq.get_const(IR::TagVal(1.0)));
    auto brcond1 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, builder.getConstFloat(3.4)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g0)}, 0);
    }();
    bs1->add_branch(brcond1, bs1, 79);

    // The reference of measure in the times may mess with the compiler.
    t1.terms.clear();
    t2.terms.clear();

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 3);
    REQUIRE(cs.host_seq.nglobals == 1);
    REQUIRE(cs.host_seq.npublic_globals == 1);
    REQUIRE(cs.host_seq.nglobal_vals == 1);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 6);

    REQUIRE(cs.host_seq.values.size() == 9);
    REQUIRE(cs.host_seq.depends.size() == 1);
    REQUIRE(cs.host_seq.types.size() == 5);
    REQUIRE(cs.host_seq.global_evals.size() == 1);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_i5 = cs.compiler.time_slot(&m2->start());
    REQUIRE(slot_i5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i5] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i5].i64 == 5);

    auto slot_f1 = cs.compiler.var_slot(p2->val());
    REQUIRE(slot_f1 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1].f64 == 1);

    auto slot_f1_5 = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_f1_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_5].f64 == 1.5);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_brcond1 = cs.compiler.var_slot(brcond1);
    REQUIRE(slot_brcond1 >= gv_offset);
    REQUIRE(slot_brcond1 < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_brcond1] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4 + 0.7;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_brcond1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_brcond1].b == (g0 < 3.4));
    }
    REQUIRE(cs.host_seq.depends[slot_brcond1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_brcond1 - gv_offset][0] == slot_g0 - gl_offset);

    // Channel:
    auto chn_offset = cs.host_seq.nconsts + cs.host_seq.nglobals + cs.host_seq.nglobal_vals;
    auto slot_ch1 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_ch1 >= chn_offset);
    REQUIRE(slot_ch1 < chn_offset + cs.host_seq.nchannels);

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 1);
        REQUIRE(host_bseq->ndirect == 1);
        REQUIRE(host_bseq->nneed_order == 1);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(host_bseq->global_refs.empty());

        REQUIRE(std::is_sorted(host_bseq->cond_global_refs.begin(),
                               host_bseq->cond_global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->cond_global_refs, {slot_brcond1 - gv_offset}));


        REQUIRE(host_bseq->types.size() == 2);
        REQUIRE(host_bseq->reverse_depends.size() == 1);
        REQUIRE(host_bseq->assumptions_idx.size() == 1);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        auto slot_m2val = cs.compiler.var_slot(m2->val());
        REQUIRE(slot_m2val >= mea_offset);
        REQUIRE(slot_m2val < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;
        auto slot_p1time = cs.compiler.time_slot(&p1->start());
        REQUIRE(slot_p1time >= dir_offset);
        REQUIRE(slot_p1time < dir_offset + host_bseq->ndirect);
        REQUIRE(host_bseq->types[slot_p1time - dir_offset] == Seq::HostSeq::Type::Int64);
        for (int i = 0; i < 11; i++) {
            double ch1 = i * 1.4 + 9.9;
            cs.host_seq.values[slot_ch1].f64 = ch1;
            host_bseq->evals[slot_p1time - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p1time].i64 == round<int64_t>(ch1) + 2);
        }

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        auto slot_p2time = cs.compiler.time_slot(&p2->start());
        REQUIRE(slot_p2time >= no_offset);
        REQUIRE(slot_p2time < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_p2time - dir_offset] == Seq::HostSeq::Type::Int64);
        for (int i = 0; i < 11; i++) {
            double m2 = i * 4.8 + 4.7;
            cs.host_seq.values[slot_m2val].f64 = m2;
            host_bseq->evals[slot_p2time - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p2time].i64 == round<int64_t>(m2) + 10);
        }
        REQUIRE(host_bseq->deps_count[slot_p2time - no_offset] == 1);

        REQUIRE(host_bseq->assumptions_idx == (std::vector<uint32_t>{uint32_t(-1)}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_m2val - mea_offset],
                                      {slot_p2time - no_offset}));

        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.empty());

        REQUIRE(host_bseq->branches.size() == 1);
        REQUIRE(host_bseq->branches[0].cond == slot_brcond1);
        REQUIRE(host_bseq->branches[0].target == 0);
        REQUIRE(host_bseq->default_branch == uint32_t(-1));

        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 3);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 5);
        REQUIRE(host_p1->time == slot_p1time);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == uint32_t(-1));
        REQUIRE(host_p1->value == slot_f1_5);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(!host_p1->ramp_func);
        REQUIRE(host_p1->endvalue == slot_f1_5);
        REQUIRE(host_p1->cond == uint32_t(-1));

        auto host_m2 = pulses[1];
        REQUIRE(host_m2->id == 9);
        REQUIRE(host_m2->time == slot_i5);
        REQUIRE(host_m2->measure == slot_m2val - mea_offset);
        REQUIRE(host_m2->len == uint32_t(-2));
        REQUIRE(host_m2->value == uint32_t(-1));
        REQUIRE(host_m2->chn == 1);
        REQUIRE(!host_m2->ramp_func);
        REQUIRE(host_m2->endvalue == uint32_t(-1));
        REQUIRE(host_m2->cond == uint32_t(-1));

        auto host_p2 = pulses[2];
        REQUIRE(host_p2->id == 10);
        REQUIRE(host_p2->time == slot_p2time);
        REQUIRE(host_p2->measure == uint32_t(-1));
        REQUIRE(host_p2->len == uint32_t(-1));
        REQUIRE(host_p2->value == slot_f1);
        REQUIRE(host_p2->chn == 1);
        REQUIRE(!host_p2->ramp_func);
        REQUIRE(host_p2->endvalue == slot_f1);
        REQUIRE(host_p2->cond == uint32_t(-1));
    }
}

TEST_CASE("default_value") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    seq.set_defval(1, 123.4);

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 1);
    REQUIRE(cs.host_seq.nglobals == 0);
    REQUIRE(cs.host_seq.npublic_globals == 0);
    REQUIRE(cs.host_seq.nglobal_vals == 0);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 2);

    REQUIRE(cs.host_seq.values.size() == 2);
    REQUIRE(cs.host_seq.depends.size() == 0);
    REQUIRE(cs.host_seq.types.size() == 1);
    REQUIRE(cs.host_seq.global_evals.size() == 0);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 123.4);

    // Constants:
    auto slot_def = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_def < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_def] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_def].f64 == 123.4);
}

TEST_CASE("cond_clone") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);

    auto p1cond = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::ge, 0, builder.getConstFloat(-2)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g0)}, 0);
    }();

    // Use the condition to disable forwarding of measure values
    auto p1 = bs1->add_pulse(1, 5, bs1->track_time(Seq::EventTime(0)), nullptr,
                             seq.get_const(IR::TagVal(1.5)), p1cond);
    auto m1 = bs1->add_measure(1, 9, bs1->track_time(Seq::EventTime(5)),
                               bs1->new_measure(seq.env(), 1));
    auto p2 = bs1->add_pulse(1, 19, bs1->track_time(Seq::EventTime(10)), nullptr, g1);
    auto m2 = bs1->add_measure(1, 28, bs1->track_time(Seq::EventTime(20)),
                               bs1->new_measure(seq.env(), 2));
    // Branch that use measures and globals
    auto brcond1 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g2),
                Seq::Arg::create_var(m2->val())}, 0);
    }();
    bs1->add_branch(brcond1, nullptr, 79);
    // Branch that use only measures
    auto brcond2 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, builder.getConstFloat(-1.0)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }();
    bs1->add_branch(brcond2, bs1, 103);

    seq.prepare();
    seq.optimize();

    brcond1 = bs1->get_branches()[0].cond.get();
    REQUIRE(brcond1->is_call());
    REQUIRE(brcond1->args().size() == 2);
    REQUIRE(brcond1->args()[0].is_var());
    REQUIRE(brcond1->args()[0].get_var() == g2);
    REQUIRE(brcond1->args()[1].is_var());
    auto m2_clone = brcond1->args()[1].get_var();
    REQUIRE(m2_clone->is_extern());
    REQUIRE(m2_clone->get_extern().first == IR::Type::Float64);
    REQUIRE((m2_clone->get_extern().second >> 32) == 0);
    REQUIRE(bs1->get_assigns().find(uint32_t(m2_clone->get_extern().second)) !=
            bs1->get_assigns().end());
    auto &assign_m2 =
        bs1->get_assigns().find(uint32_t(m2_clone->get_extern().second))->second;
    // m2 should have been eliminated and the value replaced with the value of p2
    REQUIRE(assign_m2.val.get() == g1);

    auto brcond2_clone = bs1->get_branches()[1].cond.get();
    REQUIRE(brcond2_clone != brcond2);
    REQUIRE(brcond2_clone->is_extern());
    REQUIRE(brcond2_clone->get_extern().first == IR::Type::Bool);
    REQUIRE((brcond2_clone->get_extern().second >> 32) == 0);
    REQUIRE(bs1->get_assigns().find(uint32_t(brcond2_clone->get_extern().second)) !=
            bs1->get_assigns().end());
    auto &assign_brcond2 =
        bs1->get_assigns().find(uint32_t(brcond2_clone->get_extern().second))->second;
    REQUIRE(assign_brcond2.val.get() == brcond2);

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 4);
    REQUIRE(cs.host_seq.nglobals == 5);
    REQUIRE(cs.host_seq.npublic_globals == 3);
    REQUIRE(cs.host_seq.nglobal_vals == 2);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 12);

    REQUIRE(cs.host_seq.values.size() == 14);
    REQUIRE(cs.host_seq.depends.size() == 2);
    REQUIRE(cs.host_seq.types.size() == 11);
    REQUIRE(cs.host_seq.global_evals.size() == 2);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    auto slot_i5 = cs.compiler.time_slot(&m1->start());
    REQUIRE(slot_i5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i5] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i5].i64 == 5);

    auto slot_i10 = cs.compiler.time_slot(&p2->start());
    REQUIRE(slot_i10 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i10] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i10].i64 == 10);

    auto slot_f1_5 = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_f1_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_5].f64 == 1.5);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);
    auto slot_brcond2_clone = cs.compiler.var_slot(brcond2_clone);
    REQUIRE(slot_brcond2_clone >= gl_offset);
    REQUIRE(slot_brcond2_clone < gl_offset + cs.host_seq.nglobals);
    REQUIRE(slot_brcond2_clone == gl_offset + brcond2_clone->get_extern().second);
    REQUIRE(cs.host_seq.types[slot_brcond2_clone] == Seq::HostSeq::Type::Bool);
    auto slot_m2_clone = cs.compiler.var_slot(m2_clone);
    REQUIRE(slot_m2_clone >= gl_offset);
    REQUIRE(slot_m2_clone < gl_offset + cs.host_seq.nglobals);
    REQUIRE(slot_m2_clone == gl_offset + m2_clone->get_extern().second);
    REQUIRE(cs.host_seq.types[slot_m2_clone] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_p1cond = cs.compiler.var_slot(p1cond);
    REQUIRE(slot_p1cond >= gv_offset);
    REQUIRE(slot_p1cond < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1cond] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4 - 3.1;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_p1cond - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p1cond].b == (g0 >= -2));
    }
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset][0] == slot_g0 - gl_offset);

    auto slot_brcond1 = cs.compiler.var_slot(brcond1);
    REQUIRE(slot_brcond1 >= gv_offset);
    REQUIRE(slot_brcond1 < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_brcond1] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g2 = i * 0.4 - 1.1;
        cs.host_seq.values[slot_g2].f64 = g2;
        for (int j = 0; j < 11; j++) {
            double m2_clone = j * 0.4 - 0.9;
            cs.host_seq.values[slot_m2_clone].f64 = m2_clone;
            cs.host_seq.global_evals[slot_brcond1 - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_brcond1].b == (g2 < m2_clone));
        }
    }
    REQUIRE(cs.host_seq.depends[slot_brcond1 - gv_offset].size() == 2);
    REQUIRE(cs.host_seq.depends[slot_brcond1 - gv_offset][0] == slot_g2 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_brcond1 - gv_offset][1] == slot_m2_clone - gl_offset);

    // Channel:
    auto chn_offset = cs.host_seq.nconsts + cs.host_seq.nglobals + cs.host_seq.nglobal_vals;
    auto slot_ch1 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_ch1 >= chn_offset);
    REQUIRE(slot_ch1 < chn_offset + cs.host_seq.nchannels);

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 1);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 1);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs, {slot_p1cond - gv_offset}));

        REQUIRE(std::is_sorted(host_bseq->cond_global_refs.begin(),
                               host_bseq->cond_global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->cond_global_refs, {slot_brcond1 - gv_offset}));

        REQUIRE(host_bseq->types.size() == 1);
        REQUIRE(host_bseq->reverse_depends.size() == 1);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        auto slot_m1val = cs.compiler.var_slot(m1->val());
        REQUIRE(slot_m1val >= mea_offset);
        REQUIRE(slot_m1val < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        auto slot_brcond2 = cs.compiler.var_slot(brcond2);
        REQUIRE(slot_brcond2 >= no_offset);
        REQUIRE(slot_brcond2 < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_brcond2 - dir_offset] == Seq::HostSeq::Type::Bool);
        for (int i = 0; i < 11; i++) {
            double m1 = i * 0.6 - 2.3;
            cs.host_seq.values[slot_m1val].f64 = m1;
            host_bseq->evals[slot_brcond2 - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_brcond2].b == (m1 < -1.0));
        }
        REQUIRE(host_bseq->deps_count[slot_brcond2 - no_offset] == 1);

        REQUIRE(host_bseq->assumptions_idx == (std::vector<uint32_t>{uint32_t(-1)}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_m1val - mea_offset],
                                      {slot_brcond2 - no_offset}));

        REQUIRE(host_bseq->assignments.size() == 2);
        REQUIRE(compare_vector_sorted(host_bseq->assignments,
                                      {{slot_brcond2_clone - gl_offset, slot_brcond2},
                                       {slot_m2_clone - gl_offset, slot_g1}}));

        REQUIRE(host_bseq->assumptions.empty());

        REQUIRE(host_bseq->branches.size() == 2);
        REQUIRE(host_bseq->branches[0].cond == slot_brcond1);
        REQUIRE(host_bseq->branches[0].target == uint32_t(-1));
        REQUIRE(host_bseq->branches[1].cond == slot_brcond2_clone);
        REQUIRE(host_bseq->branches[1].target == 0);
        REQUIRE(host_bseq->default_branch == uint32_t(-1));

        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 3);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 5);
        REQUIRE(host_p1->time == slot_i0);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == uint32_t(-1));
        REQUIRE(host_p1->value == slot_f1_5);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(!host_p1->ramp_func);
        REQUIRE(host_p1->endvalue == slot_f1_5);
        REQUIRE(host_p1->cond == slot_p1cond);

        auto host_m1 = pulses[1];
        REQUIRE(host_m1->id == 9);
        REQUIRE(host_m1->time == slot_i5);
        REQUIRE(host_m1->measure == slot_m1val - mea_offset);
        REQUIRE(host_m1->len == uint32_t(-2));
        REQUIRE(host_m1->value == uint32_t(-1));
        REQUIRE(host_m1->chn == 1);
        REQUIRE(!host_m1->ramp_func);
        REQUIRE(host_m1->endvalue == uint32_t(-1));
        REQUIRE(host_m1->cond == uint32_t(-1));

        auto host_p2 = pulses[2];
        REQUIRE(host_p2->id == 19);
        REQUIRE(host_p2->time == slot_i10);
        REQUIRE(host_p2->measure == uint32_t(-1));
        REQUIRE(host_p2->len == uint32_t(-1));
        REQUIRE(host_p2->value == slot_g1);
        REQUIRE(host_p2->chn == 1);
        REQUIRE(!host_p2->ramp_func);
        REQUIRE(host_p2->endvalue == slot_g1);
        REQUIRE(host_p2->cond == uint32_t(-1));
    }
}

TEST_CASE("assignment") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);
    auto p1 = bs1->add_pulse(1, 11, bs1->track_time(Seq::EventTime(0)), nullptr,
                             seq.get_const(IR::TagVal(1.2)), g1);
    auto m1 = bs1->add_measure(1, 15, bs1->track_time(Seq::EventTime(5)),
                               bs1->new_measure(seq.env(), 1));
    auto g0_assign = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createSub(0, builder.getConstFloat(1.4)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }();
    bs1->assign_global(0, g0_assign, 10);
    bs1->assign_global(1, g0, 13);
    auto v_f2_3 = seq.get_const(IR::TagVal(2.3));
    bs1->assign_global(2, v_f2_3, 19);

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 5);
    REQUIRE(cs.host_seq.nglobals == 3);
    REQUIRE(cs.host_seq.npublic_globals == 3);
    REQUIRE(cs.host_seq.nglobal_vals == 0);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 9);

    REQUIRE(cs.host_seq.values.size() == 11);
    REQUIRE(cs.host_seq.depends.size() == 0);
    REQUIRE(cs.host_seq.types.size() == 8);
    REQUIRE(cs.host_seq.global_evals.size() == 0);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    auto slot_i5 = cs.compiler.time_slot(&m1->start());
    REQUIRE(slot_i5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i5] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i5].i64 == 5);

    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    auto slot_f1_2 = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_f1_2 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_2] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_2].f64 == 1.2);

    auto slot_f2_3 = cs.compiler.var_slot(v_f2_3);
    REQUIRE(slot_f2_3 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f2_3] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f2_3].f64 == 2.3);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);

    // Global values:

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 1);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 1);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(host_bseq->global_refs.empty());
        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.size() == 1);
        REQUIRE(host_bseq->reverse_depends.size() == 1);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        auto slot_m1val = cs.compiler.var_slot(m1->val());
        REQUIRE(slot_m1val >= mea_offset);
        REQUIRE(slot_m1val < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        auto slot_g0_assign = cs.compiler.var_slot(g0_assign);
        REQUIRE(slot_g0_assign >= no_offset);
        REQUIRE(slot_g0_assign < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_g0_assign - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double m1 = i * 0.6 - 2.3;
            cs.host_seq.values[slot_m1val].f64 = m1;
            host_bseq->evals[slot_g0_assign - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_g0_assign].f64 == m1 - 1.4);
        }
        REQUIRE(host_bseq->deps_count[slot_g0_assign - no_offset] == 1);

        REQUIRE(host_bseq->assumptions_idx == (std::vector<uint32_t>{uint32_t(-1)}));
        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_m1val - mea_offset],
                                      {slot_g0_assign - no_offset}));

        REQUIRE(host_bseq->assignments.size() == 3);
        REQUIRE(compare_vector_sorted(host_bseq->assignments,
                                      {{slot_g0 - gl_offset, slot_g0_assign},
                                       {slot_g1 - gl_offset, slot_g0},
                                       {slot_g2 - gl_offset, slot_f2_3}}));

        REQUIRE(host_bseq->assumptions.empty());

        REQUIRE(host_bseq->branches.empty());

        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 2);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 11);
        REQUIRE(host_p1->time == slot_i0);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == uint32_t(-1));
        REQUIRE(host_p1->value == slot_f1_2);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(!host_p1->ramp_func);
        REQUIRE(host_p1->endvalue == slot_f1_2);
        REQUIRE(host_p1->cond == slot_g1);

        auto host_m1 = pulses[1];
        REQUIRE(host_m1->id == 15);
        REQUIRE(host_m1->time == slot_i5);
        REQUIRE(host_m1->measure == slot_m1val - mea_offset);
        REQUIRE(host_m1->len == uint32_t(-2));
        REQUIRE(host_m1->value == uint32_t(-1));
        REQUIRE(host_m1->chn == 1);
        REQUIRE(!host_m1->ramp_func);
        REQUIRE(host_m1->endvalue == uint32_t(-1));
        REQUIRE(host_m1->cond == uint32_t(-1));
    }
}

TEST_CASE("assumption") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);

    auto p1cond = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::ge, 0, builder.getConstFloat(-2)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g0)}, 0);
    }();

    // Use the condition to disable forwarding of measure values
    auto p1 = bs1->add_pulse(1, 5, bs1->track_time(Seq::EventTime(0)), nullptr,
                             seq.get_const(IR::TagVal(1.5)), p1cond);
    auto m1 = bs1->add_measure(1, 9, bs1->track_time(Seq::EventTime(5)),
                               bs1->new_measure(seq.env(), 1));
    auto assume1val = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        auto v = builder.createMul(0, builder.getConstFloat(23.4));
        builder.createRet(builder.createSub(v, builder.getConstFloat(1.2)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(g1)}, 0);
    }();
    bs1->add_assume(Seq::Sign::Pos, assume1val, 10);
    auto assume2val = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        auto v = builder.createMul(0, builder.getConstFloat(2));
        builder.createRet(builder.createAdd(v, builder.getConstFloat(1.2)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }();
    bs1->add_assume(Seq::Sign::NonNeg, assume2val, 19);
    auto assume3val = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        auto v = builder.createMul(0, builder.createAdd(1, builder.getConstFloat(4.9)));
        builder.createRet(builder.createAdd(v, builder.getConstFloat(101.4)));
        return seq.get_call(builder.get(),
                            {Seq::Arg::create_var(m1->val()), Seq::Arg::create_var(g2)}, 0);
    }();
    bs1->add_assume(Seq::Sign::Pos, assume3val, 49);

    seq.prepare();
    seq.optimize();

    std::vector<const std::pair<const Seq::Var::Ref,Seq::BasicSeq::Assumption>*> bs1_assumes;
    for (auto &assume: bs1->get_assumes())
        bs1_assumes.push_back(&assume);
    std::sort(bs1_assumes.begin(), bs1_assumes.end(), [] (auto *a1, auto *a2) {
        return a1->second.id < a2->second.id;
    });
    REQUIRE(bs1_assumes.size() == 3);
    REQUIRE(bs1_assumes[0]->first.get() == assume1val);
    REQUIRE(bs1_assumes[0]->second.sign == Seq::Sign::Pos);
    REQUIRE(bs1_assumes[0]->second.id == 10);
    REQUIRE(bs1_assumes[1]->first.get() == assume2val);
    REQUIRE(bs1_assumes[1]->second.sign == Seq::Sign::NonNeg);
    REQUIRE(bs1_assumes[1]->second.id == 19);
    REQUIRE(bs1_assumes[2]->first.get() == assume3val);
    REQUIRE(bs1_assumes[2]->second.sign == Seq::Sign::Pos);
    REQUIRE(bs1_assumes[2]->second.id == 49);

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 4);
    REQUIRE(cs.host_seq.nglobals == 3);
    REQUIRE(cs.host_seq.npublic_globals == 3);
    REQUIRE(cs.host_seq.nglobal_vals == 3);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 11);

    REQUIRE(cs.host_seq.values.size() == 16);
    REQUIRE(cs.host_seq.depends.size() == 3);
    REQUIRE(cs.host_seq.types.size() == 10);
    REQUIRE(cs.host_seq.global_evals.size() == 3);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_i0 = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_i0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i0] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i0].i64 == 0);

    auto slot_i5 = cs.compiler.time_slot(&m1->start());
    REQUIRE(slot_i5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_i5] == Seq::HostSeq::Type::Int64);
    REQUIRE(cs.host_seq.values[slot_i5].i64 == 5);

    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    auto slot_f1_5 = cs.compiler.var_slot(p1->val());
    REQUIRE(slot_f1_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_5].f64 == 1.5);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_p1cond = cs.compiler.var_slot(p1cond);
    REQUIRE(slot_p1cond >= gv_offset);
    REQUIRE(slot_p1cond < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1cond] == Seq::HostSeq::Type::Bool);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4 - 3.1;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_p1cond - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p1cond].b == (g0 >= -2));
    }
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1cond - gv_offset][0] == slot_g0 - gl_offset);

    auto slot_assume1 = cs.compiler.var_slot(assume1val);
    REQUIRE(slot_assume1 >= gv_offset);
    REQUIRE(slot_assume1 < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_assume1] == Seq::HostSeq::Type::Float64);
    for (int i = 0; i < 11; i++) {
        double g1 = i * 0.4 - 1.1;
        cs.host_seq.values[slot_g1].f64 = g1;
        cs.host_seq.global_evals[slot_assume1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(approx(cs.host_seq.values[slot_assume1].f64,
                       g1 * 23.4 - 1.2));
    }
    REQUIRE(cs.host_seq.depends[slot_assume1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_assume1 - gv_offset][0] == slot_g1 - gl_offset);

    // Currently we cannot merge the computation with the conversion to int64
    auto slot_assume1i = cs.compiler.assume_slot(bs1_assumes[0]);
    REQUIRE(slot_assume1i > slot_assume1);
    REQUIRE(slot_assume1i >= gv_offset);
    REQUIRE(slot_assume1i < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_assume1i] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double assume1 = i * 1.4 - 5.0;
        cs.host_seq.values[slot_assume1].f64 = assume1;
        cs.host_seq.global_evals[slot_assume1i - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_assume1i].i64 == round<int64_t>(assume1));
    }
    REQUIRE(cs.host_seq.depends[slot_assume1i - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_assume1i - gv_offset][0] == slot_g1 - gl_offset);

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 1);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 4);
        REQUIRE(host_bseq->ndirect_assumes == 1);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_p1cond - gv_offset, slot_assume1 - gv_offset,
                                       slot_assume1i - gv_offset}));

        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.size() == 4);
        REQUIRE(host_bseq->reverse_depends.size() == 1);
        REQUIRE(host_bseq->assumptions_idx.size() == 4);

        // Measure
        auto mea_offset = cs.host_seq.nshared;
        auto slot_m1val = cs.compiler.var_slot(m1->val());
        REQUIRE(slot_m1val >= mea_offset);
        REQUIRE(slot_m1val < mea_offset + host_bseq->nmeasure);

        // Direct
        auto dir_offset = cs.host_seq.nshared + host_bseq->nmeasure;

        // Need order
        auto no_offset = dir_offset + host_bseq->ndirect;
        std::vector<uint32_t> assume_idxs;
        auto slot_assume2 = cs.compiler.var_slot(assume2val);
        REQUIRE(slot_assume2 >= no_offset);
        REQUIRE(slot_assume2 < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_assume2 - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double m1 = i * 0.8 - 0.9;
            cs.host_seq.values[slot_m1val].f64 = m1;
            host_bseq->evals[slot_assume2 - dir_offset](cs.host_seq.values.data());
            REQUIRE(approx(cs.host_seq.values[slot_assume2].f64, m1 * 2 + 1.2));
        }
        REQUIRE(host_bseq->deps_count[slot_assume2 - no_offset] == 1);
        REQUIRE(host_bseq->assumptions_idx[slot_assume2 - no_offset] == uint32_t(-1));

        auto slot_assume2i = cs.compiler.assume_slot(bs1_assumes[1]);
        REQUIRE(slot_assume2i > slot_assume2);
        REQUIRE(slot_assume2i >= no_offset);
        REQUIRE(slot_assume2i < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_assume2i - dir_offset] == Seq::HostSeq::Type::Int64);
        for (int i = 0; i < 11; i++) {
            double assume2 = i * 0.8 - 3.8;
            cs.host_seq.values[slot_assume2].f64 = assume2;
            host_bseq->evals[slot_assume2i - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_assume2i].i64 == round<int64_t>(assume2));
        }
        REQUIRE(host_bseq->deps_count[slot_assume2i - no_offset] == 1);
        assume_idxs.push_back(host_bseq->assumptions_idx[slot_assume2i - no_offset]);

        auto slot_assume3 = cs.compiler.var_slot(assume3val);
        REQUIRE(slot_assume3 >= no_offset);
        REQUIRE(slot_assume3 < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_assume3 - dir_offset] == Seq::HostSeq::Type::Float64);
        for (int i = 0; i < 11; i++) {
            double m1 = i * 0.8 - 0.9;
            cs.host_seq.values[slot_m1val].f64 = m1;
            for (int j = 0; j < 11; j++) {
                double g2 = j * 1.4 - 5.9;
                cs.host_seq.values[slot_g2].f64 = g2;
                host_bseq->evals[slot_assume3 - dir_offset](cs.host_seq.values.data());
                REQUIRE(approx(cs.host_seq.values[slot_assume3].f64, m1 * (g2 + 4.9) + 101.4));
            }
        }
        REQUIRE(host_bseq->deps_count[slot_assume3 - no_offset] == 1);
        REQUIRE(host_bseq->assumptions_idx[slot_assume3 - no_offset] == uint32_t(-1));

        auto slot_assume3i = cs.compiler.assume_slot(bs1_assumes[2]);
        REQUIRE(slot_assume3i > slot_assume3);
        REQUIRE(slot_assume3i >= no_offset);
        REQUIRE(slot_assume3i < no_offset + host_bseq->nneed_order);
        REQUIRE(host_bseq->types[slot_assume3i - dir_offset] == Seq::HostSeq::Type::Int64);
        for (int i = 0; i < 11; i++) {
            double assume3 = i * 0.8 - 3.8;
            cs.host_seq.values[slot_assume3].f64 = assume3;
            host_bseq->evals[slot_assume3i - dir_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_assume3i].i64 == round<int64_t>(assume3));
        }
        REQUIRE(host_bseq->deps_count[slot_assume3i - no_offset] == 1);
        assume_idxs.push_back(host_bseq->assumptions_idx[slot_assume3i - no_offset]);

        REQUIRE(compare_vector_sorted(host_bseq->reverse_depends[slot_m1val - mea_offset],
                                      {slot_assume2 - no_offset, slot_assume2i - no_offset,
                                       slot_assume3 - no_offset, slot_assume3i - no_offset}));

        REQUIRE(host_bseq->assignments.empty());

        REQUIRE(compare_vector_sorted(assume_idxs, {1, 2}));
        REQUIRE(host_bseq->assumptions.size() == 3);
        REQUIRE(host_bseq->assumptions[0] ==
                (Seq::HostSeq::Assumption{Seq::Sign::Pos, slot_assume1i, 10}));
        REQUIRE(host_bseq->assumptions[assume_idxs[0]] ==
                (Seq::HostSeq::Assumption{Seq::Sign::NonNeg, slot_assume2i, 19}));
        REQUIRE(host_bseq->assumptions[assume_idxs[1]] ==
                (Seq::HostSeq::Assumption{Seq::Sign::Pos, slot_assume3i, 49}));

        REQUIRE(host_bseq->branches.empty());
        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 2);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 5);
        REQUIRE(host_p1->time == slot_i0);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == uint32_t(-1));
        REQUIRE(host_p1->value == slot_f1_5);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(!host_p1->ramp_func);
        REQUIRE(host_p1->endvalue == slot_f1_5);
        REQUIRE(host_p1->cond == slot_p1cond);

        auto host_m1 = pulses[1];
        REQUIRE(host_m1->id == 9);
        REQUIRE(host_m1->time == slot_i5);
        REQUIRE(host_m1->measure == slot_m1val - mea_offset);
        REQUIRE(host_m1->len == uint32_t(-2));
        REQUIRE(host_m1->value == uint32_t(-1));
        REQUIRE(host_m1->chn == 1);
        REQUIRE(!host_m1->ramp_func);
        REQUIRE(host_m1->endvalue == uint32_t(-1));
        REQUIRE(host_m1->cond == uint32_t(-1));
    }
}

TEST_CASE("endtime") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);
    Seq::EventTime t1(4);
    t1.add_term(Seq::Sign::NonNeg, g0);
    Seq::EventTime t2(9);
    t2.add_term(Seq::Sign::NonNeg, g1);
    t2.add_term(Seq::Sign::NonNeg, g2);
    t2.add_term(Seq::Sign::NonNeg, g2);
    auto &rt1 = bs1->track_time(t1);
    auto &rt2 = bs1->track_time(t2);
    bs1->add_endtime(rt1);
    bs1->add_endtime(rt2);
    t1.terms.clear();
    t2.terms.clear();

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 1);
    REQUIRE(cs.host_seq.nglobals == 3);
    REQUIRE(cs.host_seq.npublic_globals == 3);
    REQUIRE(cs.host_seq.nglobal_vals == 2);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 7);

    REQUIRE(cs.host_seq.values.size() == 7);
    REQUIRE(cs.host_seq.depends.size() == 2);
    REQUIRE(cs.host_seq.types.size() == 6);
    REQUIRE(cs.host_seq.global_evals.size() == 2);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    auto slot_t1 = cs.compiler.time_slot(&rt1);
    REQUIRE(slot_t1 >= gv_offset);
    REQUIRE(slot_t1 < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_t1] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_t1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_t1].i64 == round<int64_t>(g0) + 4);
    }
    REQUIRE(cs.host_seq.depends[slot_t1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_t1 - gv_offset][0] == slot_g0 - gl_offset);

    auto slot_t2 = cs.compiler.time_slot(&rt2);
    REQUIRE(slot_t2 >= gv_offset);
    REQUIRE(slot_t2 < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_t2] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g1 = i * 0.4;
        cs.host_seq.values[slot_g1].f64 = g1;
        for (int j = 0; j < 11; j++) {
            double g2 = j * 0.4;
            cs.host_seq.values[slot_g2].f64 = g2;
            cs.host_seq.global_evals[slot_t2 - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_t2].i64 ==
                    round<int64_t>(g1) + round<int64_t>(g2) * 2 + 9);
        }
    }
    REQUIRE(cs.host_seq.depends[slot_t2 - gv_offset].size() == 2);
    REQUIRE(cs.host_seq.depends[slot_t2 - gv_offset][0] == slot_g1 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_t2 - gv_offset][1] == slot_g2 - gl_offset);

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 0);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 0);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_t1 - gv_offset, slot_t2 - gv_offset}));
        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.empty());
        REQUIRE(host_bseq->reverse_depends.empty());

        // Measure

        // Direct

        // Need order

        REQUIRE(host_bseq->assumptions_idx.empty());
        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.empty());
        REQUIRE(host_bseq->branches.empty());

        REQUIRE(compare_vector_sorted(host_bseq->endtimes, {slot_t1, slot_t2}));

        REQUIRE(host_bseq->pulses.empty());
    }
}

TEST_CASE("time_chain") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto g0 = seq.get_slot(IR::Type::Float64, 0);
    auto g1 = seq.get_slot(IR::Type::Float64, 1);
    auto g2 = seq.get_slot(IR::Type::Float64, 2);
    auto g3 = seq.get_slot(IR::Type::Float64, 3);

    Seq::EventTime t(5);
    t.add_term(Seq::Sign::Pos, g0);
    auto p1 = bs1->add_pulse(1, 5, bs1->track_time(t), nullptr,
                             seq.get_const(IR::TagVal(0.0)));
    t.tconst = 10;
    auto p2 = bs1->add_pulse(1, 6, bs1->track_time(t), nullptr,
                             seq.get_const(IR::TagVal(1.0)));
    Seq::EventTime t2(t);
    t2.tconst = 12;
    t2.add_term(Seq::Sign::Pos, g1);
    auto p3 = bs1->add_pulse(1, 7, bs1->track_time(t2), nullptr,
                             seq.get_const(IR::TagVal(0.0)));
    t.tconst = 15;
    t.add_term(Seq::Sign::Pos, g2);
    auto p4 = bs1->add_pulse(1, 8, bs1->track_time(t), nullptr,
                             seq.get_const(IR::TagVal(0.5)));
    t.tconst = 20;
    t.add_term(Seq::Sign::Pos, g3);
    auto p5 = bs1->add_pulse(1, 9, bs1->track_time(t), nullptr,
                             seq.get_const(IR::TagVal(1.5)));
    t.terms.clear();
    t2.terms.clear();

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 4);
    REQUIRE(cs.host_seq.nglobals == 4);
    REQUIRE(cs.host_seq.npublic_globals == 4);
    REQUIRE(cs.host_seq.nglobal_vals == 5);
    REQUIRE(cs.host_seq.nchannels == 1);
    REQUIRE(cs.host_seq.nshared == 14);

    REQUIRE(cs.host_seq.values.size() == 14);
    REQUIRE(cs.host_seq.depends.size() == 5);
    REQUIRE(cs.host_seq.types.size() == 13);
    REQUIRE(cs.host_seq.global_evals.size() == 5);
    REQUIRE(cs.host_seq.default_values.size() == 1);
    REQUIRE(cs.host_seq.default_values[0].f64 == 0);

    // Constants:
    auto slot_f0 = cs.compiler.var_slot(bs1->startval(1));
    REQUIRE(slot_f0 < cs.host_seq.nconsts);
    REQUIRE(slot_f0 == cs.compiler.var_slot(p1->val()));
    REQUIRE(slot_f0 == cs.compiler.var_slot(p3->val()));
    REQUIRE(cs.host_seq.types[slot_f0] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0].f64 == 0);

    auto slot_f0_5 = cs.compiler.var_slot(p4->val());
    REQUIRE(slot_f0_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f0_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f0_5].f64 == 0.5);

    auto slot_f1 = cs.compiler.var_slot(p2->val());
    REQUIRE(slot_f1 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1].f64 == 1);

    auto slot_f1_5 = cs.compiler.var_slot(p5->val());
    REQUIRE(slot_f1_5 < cs.host_seq.nconsts);
    REQUIRE(cs.host_seq.types[slot_f1_5] == Seq::HostSeq::Type::Float64);
    REQUIRE(cs.host_seq.values[slot_f1_5].f64 == 1.5);

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = cs.compiler.var_slot(g0);
    REQUIRE(slot_g0 == gl_offset + 0);
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);
    auto slot_g1 = cs.compiler.var_slot(g1);
    REQUIRE(slot_g1 == gl_offset + 1);
    REQUIRE(cs.host_seq.types[slot_g1] == Seq::HostSeq::Type::Float64);
    auto slot_g2 = cs.compiler.var_slot(g2);
    REQUIRE(slot_g2 == gl_offset + 2);
    REQUIRE(cs.host_seq.types[slot_g2] == Seq::HostSeq::Type::Float64);
    auto slot_g3 = cs.compiler.var_slot(g3);
    REQUIRE(slot_g3 == gl_offset + 3);
    REQUIRE(cs.host_seq.types[slot_g3] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = cs.host_seq.nconsts + cs.host_seq.nglobals;
    // p1time = round(g0) + 5
    auto slot_p1time = cs.compiler.time_slot(&p1->start());
    REQUIRE(slot_p1time >= gv_offset);
    REQUIRE(slot_p1time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p1time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g0 = i * 0.4;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_p1time - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p1time].i64 == round<int64_t>(g0) + 5);
    }
    REQUIRE(cs.host_seq.depends[slot_p1time - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p1time - gv_offset][0] == slot_g0 - gl_offset);
    // p2time = p1time + 5
    auto slot_p2time = cs.compiler.time_slot(&p2->start());
    REQUIRE(slot_p2time >= gv_offset);
    REQUIRE(slot_p2time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p2time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        int64_t p1time = i * 21;
        cs.host_seq.values[slot_p1time].i64 = p1time;
        cs.host_seq.global_evals[slot_p2time - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_p2time].i64 == p1time + 5);
    }
    REQUIRE(cs.host_seq.depends[slot_p2time - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_p2time - gv_offset][0] == slot_g0 - gl_offset);
    // p3time = p1time + round(g1) + 7
    auto slot_p3time = cs.compiler.time_slot(&p3->start());
    REQUIRE(slot_p3time >= gv_offset);
    REQUIRE(slot_p3time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p3time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g1 = i * 1.4;
        cs.host_seq.values[slot_g1].f64 = g1;
        for (int j = 0; j < 11; j++) {
            int64_t p1time = j * 25;
            cs.host_seq.values[slot_p1time].i64 = p1time;
            cs.host_seq.global_evals[slot_p3time - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p3time].i64 == p1time + round<int64_t>(g1) + 7);
        }
    }
    REQUIRE(cs.host_seq.depends[slot_p3time - gv_offset].size() == 2);
    REQUIRE(cs.host_seq.depends[slot_p3time - gv_offset][0] == slot_g0 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_p3time - gv_offset][1] == slot_g1 - gl_offset);
    // p4time = p1time + round(g2) + 10
    auto slot_p4time = cs.compiler.time_slot(&p4->start());
    REQUIRE(slot_p4time >= gv_offset);
    REQUIRE(slot_p4time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p4time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g2 = i * 2.8;
        cs.host_seq.values[slot_g2].f64 = g2;
        for (int j = 0; j < 11; j++) {
            int64_t p1time = j * 45;
            cs.host_seq.values[slot_p1time].i64 = p1time;
            cs.host_seq.global_evals[slot_p4time - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p4time].i64 == p1time + round<int64_t>(g2) + 10);
        }
    }
    REQUIRE(cs.host_seq.depends[slot_p4time - gv_offset].size() == 2);
    REQUIRE(cs.host_seq.depends[slot_p4time - gv_offset][0] == slot_g0 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_p4time - gv_offset][1] == slot_g2 - gl_offset);
    // p5time = p4time + round(g3) + 5
    auto slot_p5time = cs.compiler.time_slot(&p5->start());
    REQUIRE(slot_p5time >= gv_offset);
    REQUIRE(slot_p5time < gv_offset + cs.host_seq.nglobal_vals);
    REQUIRE(cs.host_seq.types[slot_p5time] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 11; i++) {
        double g3 = i * 28.8;
        cs.host_seq.values[slot_g3].f64 = g3;
        for (int j = 0; j < 11; j++) {
            int64_t p4time = j * 51;
            cs.host_seq.values[slot_p4time].i64 = p4time;
            cs.host_seq.global_evals[slot_p5time - gv_offset](cs.host_seq.values.data());
            REQUIRE(cs.host_seq.values[slot_p5time].i64 == p4time + round<int64_t>(g3) + 5);
        }
    }
    REQUIRE(cs.host_seq.depends[slot_p5time - gv_offset].size() == 3);
    REQUIRE(cs.host_seq.depends[slot_p5time - gv_offset][0] == slot_g0 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_p5time - gv_offset][1] == slot_g2 - gl_offset);
    REQUIRE(cs.host_seq.depends[slot_p5time - gv_offset][2] == slot_g3 - gl_offset);

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 0);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 0);
        REQUIRE(host_bseq->ndirect_assumes == 0);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_p1time - gv_offset, slot_p2time - gv_offset,
                                       slot_p3time - gv_offset, slot_p4time - gv_offset,
                                       slot_p5time - gv_offset}));

        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.empty());
        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assumptions_idx.empty());

        // Measure

        // Direct

        // Need order

        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.empty());
        REQUIRE(host_bseq->branches.empty());
        REQUIRE(host_bseq->endtimes.empty());

        REQUIRE(host_bseq->pulses.size() == 5);

        std::vector<Seq::HostSeq::Pulse*> pulses;
        for (auto &pulse: host_bseq->pulses)
            pulses.push_back(&pulse);
        std::sort(pulses.begin(), pulses.end(), [] (auto *p1, auto *p2) {
            return p1->id < p2->id;
        });
        auto host_p1 = pulses[0];
        REQUIRE(host_p1->id == 5);
        REQUIRE(host_p1->time == slot_p1time);
        REQUIRE(host_p1->measure == uint32_t(-1));
        REQUIRE(host_p1->len == uint32_t(-1));
        REQUIRE(host_p1->value == slot_f0);
        REQUIRE(host_p1->chn == 1);
        REQUIRE(!host_p1->ramp_func);
        REQUIRE(host_p1->endvalue == slot_f0);
        REQUIRE(host_p1->cond == uint32_t(-1));

        auto host_p2 = pulses[1];
        REQUIRE(host_p2->id == 6);
        REQUIRE(host_p2->time == slot_p2time);
        REQUIRE(host_p2->measure == uint32_t(-1));
        REQUIRE(host_p2->len == uint32_t(-1));
        REQUIRE(host_p2->value == slot_f1);
        REQUIRE(host_p2->chn == 1);
        REQUIRE(!host_p2->ramp_func);
        REQUIRE(host_p2->endvalue == slot_f1);
        REQUIRE(host_p2->cond == uint32_t(-1));

        auto host_p3 = pulses[2];
        REQUIRE(host_p3->id == 7);
        REQUIRE(host_p3->time == slot_p3time);
        REQUIRE(host_p3->measure == uint32_t(-1));
        REQUIRE(host_p3->len == uint32_t(-1));
        REQUIRE(host_p3->value == slot_f0);
        REQUIRE(host_p3->chn == 1);
        REQUIRE(!host_p3->ramp_func);
        REQUIRE(host_p3->endvalue == slot_f0);
        REQUIRE(host_p3->cond == uint32_t(-1));

        auto host_p4 = pulses[3];
        REQUIRE(host_p4->id == 8);
        REQUIRE(host_p4->time == slot_p4time);
        REQUIRE(host_p4->measure == uint32_t(-1));
        REQUIRE(host_p4->len == uint32_t(-1));
        REQUIRE(host_p4->value == slot_f0_5);
        REQUIRE(host_p4->chn == 1);
        REQUIRE(!host_p4->ramp_func);
        REQUIRE(host_p4->endvalue == slot_f0_5);
        REQUIRE(host_p4->cond == uint32_t(-1));

        auto host_p5 = pulses[4];
        REQUIRE(host_p5->id == 9);
        REQUIRE(host_p5->time == slot_p5time);
        REQUIRE(host_p5->measure == uint32_t(-1));
        REQUIRE(host_p5->len == uint32_t(-1));
        REQUIRE(host_p5->value == slot_f1_5);
        REQUIRE(host_p5->chn == 1);
        REQUIRE(!host_p5->ramp_func);
        REQUIRE(host_p5->endvalue == slot_f1_5);
        REQUIRE(host_p5->cond == uint32_t(-1));
    }
}

TEST_CASE("Assumption offset") {
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
    stm.write(Seq::Builder::OpCode::Select);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1e9);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 2
    stm.write(Seq::Builder::OpCode::Rint);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    // 3
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1e9);
    // 4
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(4);
        stm.write<uint8_t>((uint8_t)Seq::Sign::NonNeg);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(1);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(2);
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

    REQUIRE(builder.slots.size() == 1);
    REQUIRE(builder.slots[0] == IR::Type::Float64);

    builder.buildseq();

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 0);
    REQUIRE(cs.host_seq.nglobals == 1);
    REQUIRE(cs.host_seq.npublic_globals == 1);
    REQUIRE(cs.host_seq.nglobal_vals == 3);
    REQUIRE(cs.host_seq.nchannels == 0);
    REQUIRE(cs.host_seq.nshared == 4);

    REQUIRE(cs.host_seq.values.size() == 4);
    REQUIRE(cs.host_seq.depends.size() == 3);
    REQUIRE(cs.host_seq.types.size() == 4);
    REQUIRE(cs.host_seq.global_evals.size() == 3);
    REQUIRE(cs.host_seq.default_values.size() == 0);

    // Constants:

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = gl_offset + 0;
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = gl_offset + cs.host_seq.nglobals;
    // len1 = select(g0, 1e9, 0)
    auto slot_len1 = gv_offset + 0;
    REQUIRE(cs.host_seq.types[slot_len1] == Seq::HostSeq::Type::Float64);
    for (int i = -3; i < 4; i++) {
        double g0 = i * 0.4;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_len1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_len1].f64 == (i ? 1e9 : 0));
    }
    REQUIRE(cs.host_seq.depends[slot_len1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_len1 - gv_offset][0] == slot_g0 - gl_offset);
    // ilen1 = round(len1)
    auto slot_ilen1 = gv_offset + 1;
    REQUIRE(cs.host_seq.types[slot_ilen1] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 1000; i++) {
        double len1 = i * 0.4;
        cs.host_seq.values[slot_len1].f64 = len1;
        cs.host_seq.global_evals[slot_ilen1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_ilen1].i64 == round<int64_t>(len1));
    }
    REQUIRE(cs.host_seq.depends[slot_ilen1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_ilen1 - gv_offset][0] == slot_g0 - gl_offset);
    // ilen2 = ilen1 + 1e9
    auto slot_ilen2 = gv_offset + 2;
    REQUIRE(cs.host_seq.types[slot_ilen2] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 1000; i++) {
        int64_t ilen1 = i * 3;
        cs.host_seq.values[slot_ilen1].i64 = ilen1;
        cs.host_seq.global_evals[slot_ilen2 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_ilen2].i64 == ilen1 + 1000000000);
    }
    REQUIRE(cs.host_seq.depends[slot_ilen2 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_ilen2 - gv_offset][0] == slot_g0 - gl_offset);

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 0);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 0);
        REQUIRE(host_bseq->ndirect_assumes == 1);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_len1 - gv_offset, slot_ilen1 - gv_offset,
                                       slot_ilen2 - gv_offset}));

        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.empty());
        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assumptions_idx.empty());

        // Measure

        // Direct

        // Need order

        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.size() == 1);
        REQUIRE(host_bseq->assumptions[0] ==
                (Seq::HostSeq::Assumption{Seq::Sign::NonNeg, slot_ilen1, 0}));
        REQUIRE(host_bseq->branches.empty());
        REQUIRE(host_bseq->endtimes.size() == 1);
        REQUIRE(compare_vector_sorted(host_bseq->endtimes, {slot_ilen2}));

        REQUIRE(host_bseq->pulses.size() == 0);
    }
}

TEST_CASE("Assumption offset2") {
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(7);
    // 1
    stm.write(Seq::Builder::OpCode::CmpEQ);
    stm.write(Seq::Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 2
    stm.write(Seq::Builder::OpCode::Select);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(5e9);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 3
    stm.write(Seq::Builder::OpCode::Rint);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Seq::Builder::OpCode::Select);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(1);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(5e9);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 5
    stm.write(Seq::Builder::OpCode::Rint);
    stm.write(Seq::Builder::ArgType::Node);
    stm.write<uint32_t>(4);
    // 6
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 7
    stm.write(Seq::Builder::OpCode::Identity);
    stm.write(Seq::Builder::ArgType::ConstFloat64);
    stm.write<double>(1e9);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(0);
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(5);
        stm.write<uint8_t>((uint8_t)Seq::Sign::NonNeg);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);

        stm.write<uint8_t>((uint8_t)Seq::Sign::NonNeg);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(1);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(2);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(0);

        stm.write<uint8_t>((uint8_t)Seq::Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(1);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
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

    REQUIRE(builder.slots.size() == 1);
    REQUIRE(builder.slots[0] == IR::Type::Float64);

    builder.buildseq();

    seq.prepare();
    seq.optimize();

    CompileSeq cs(seq);

    REQUIRE(cs.host_seq.nconsts == 0);
    REQUIRE(cs.host_seq.nglobals == 1);
    REQUIRE(cs.host_seq.npublic_globals == 1);
    REQUIRE(cs.host_seq.nglobal_vals == 3);
    REQUIRE(cs.host_seq.nchannels == 0);
    REQUIRE(cs.host_seq.nshared == 4);

    REQUIRE(cs.host_seq.values.size() == 4);
    REQUIRE(cs.host_seq.depends.size() == 3);
    REQUIRE(cs.host_seq.types.size() == 4);
    REQUIRE(cs.host_seq.global_evals.size() == 3);
    REQUIRE(cs.host_seq.default_values.size() == 0);

    // Constants:

    // Globals:
    auto gl_offset = cs.host_seq.nconsts;
    auto slot_g0 = gl_offset + 0;
    REQUIRE(cs.host_seq.types[slot_g0] == Seq::HostSeq::Type::Float64);

    // Global values:
    auto gv_offset = gl_offset + cs.host_seq.nglobals;
    // len1 = select(g0 == 0, 5e9, 0)
    auto slot_len1 = gv_offset + 0;
    REQUIRE(cs.host_seq.types[slot_len1] == Seq::HostSeq::Type::Float64);
    for (int i = -3; i < 4; i++) {
        double g0 = i * 0.4;
        cs.host_seq.values[slot_g0].f64 = g0;
        cs.host_seq.global_evals[slot_len1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_len1].f64 == (i == 0 ? 5e9 : 0));
    }
    REQUIRE(cs.host_seq.depends[slot_len1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_len1 - gv_offset][0] == slot_g0 - gl_offset);
    // ilen1 = round(len1)
    auto slot_ilen1 = gv_offset + 1;
    REQUIRE(cs.host_seq.types[slot_ilen1] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 1000; i++) {
        double len1 = i * 0.4;
        cs.host_seq.values[slot_len1].f64 = len1;
        cs.host_seq.global_evals[slot_ilen1 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_ilen1].i64 == round<int64_t>(len1));
    }
    REQUIRE(cs.host_seq.depends[slot_ilen1 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_ilen1 - gv_offset][0] == slot_g0 - gl_offset);
    // ilen2 = ilen1 + ilen1
    auto slot_ilen2 = gv_offset + 2;
    REQUIRE(cs.host_seq.types[slot_ilen2] == Seq::HostSeq::Type::Int64);
    for (int i = 0; i < 1000; i++) {
        int64_t ilen1 = i * 3;
        cs.host_seq.values[slot_ilen1].i64 = ilen1;
        cs.host_seq.global_evals[slot_ilen2 - gv_offset](cs.host_seq.values.data());
        REQUIRE(cs.host_seq.values[slot_ilen2].i64 == ilen1 * 2);
    }
    REQUIRE(cs.host_seq.depends[slot_ilen2 - gv_offset].size() == 1);
    REQUIRE(cs.host_seq.depends[slot_ilen2 - gv_offset][0] == slot_g0 - gl_offset);

    // Channel:

    {
        auto host_bseq = &cs.host_seq.seqs[0];
        REQUIRE(host_bseq->id == 1);

        REQUIRE(host_bseq->nmeasure == 0);
        REQUIRE(host_bseq->ndirect == 0);
        REQUIRE(host_bseq->nneed_order == 0);
        REQUIRE(host_bseq->ndirect_assumes == 1);

        REQUIRE(std::is_sorted(host_bseq->global_refs.begin(), host_bseq->global_refs.end()));
        REQUIRE(compare_vector_sorted(host_bseq->global_refs,
                                      {slot_len1 - gv_offset, slot_ilen1 - gv_offset,
                                       slot_ilen2 - gv_offset}));

        REQUIRE(host_bseq->cond_global_refs.empty());

        REQUIRE(host_bseq->types.empty());
        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assumptions_idx.empty());

        // Measure

        // Direct

        // Need order

        REQUIRE(host_bseq->reverse_depends.empty());
        REQUIRE(host_bseq->assignments.empty());
        REQUIRE(host_bseq->assumptions.size() == 1);
        REQUIRE(host_bseq->assumptions[0] ==
                (Seq::HostSeq::Assumption{Seq::Sign::NonNeg, slot_ilen1, 0}));
        REQUIRE(host_bseq->branches.empty());
        REQUIRE(host_bseq->endtimes.size() == 1);
        REQUIRE(compare_vector_sorted(host_bseq->endtimes, {slot_ilen2}));

        REQUIRE(host_bseq->pulses.size() == 0);
    }
}
