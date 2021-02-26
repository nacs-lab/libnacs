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

#include "../../lib/seq/seq.h"
#include "../../lib/seq/error.h"

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/utils.h"
#include "../../lib/utils/number.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

static void test_basic_single_channel(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1.45)));
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(0.3)));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p2->id() == 2);
    assert(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    assert(p1->start() == p2->start());
    assert(p1->start().tconst == 6);
    assert(p1->start().terms.size() == 2);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    // The sequence cache the slot variable so we can compare pointer identity.
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p1->start().terms[1].sign == Seq::EventTime::NonNeg);
    assert(p1->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(approx(p1->endval()->get_const().get<double>(), 0.3));
    assert(approx(p2->endval()->get_const().get<double>(), -0.2));
    assert(!p1->needs_oldval());
    assert(!p2->needs_oldval());
    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(!bs->needs_oldval(p1));
    assert(!bs->needs_oldval(p2));
    assert(bs->startval(1)->get_const().is(IR::TagVal(0.0)));
    assert(approx(bs->endval(1)->get_const().get<double>(), -0.2));
    assert(!bs->get_default_branch());
    assert(bs->get_branches().empty());
    auto slots = seq.get_slots();
    assert(slots.size() == 2);
    assert((slots[0]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0))));
    assert(slots[0].get() == seq.get_slot(IR::Type::Float64, 0));
    assert((slots[1]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(1))));
    assert(slots[1].get() == seq.get_slot(IR::Type::Float64, 1));
}

static void test_basic_single_channel_measure_known_time1(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(0)), nullptr,
                            [&] {
                                IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
                                builder.createRet(builder.createSub(0, builder.getConst(2)));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p2->id() == 2);
    assert(p2->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(30)));
    auto m1 = bs->add_measure(1, 3, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 3);
    assert(m1->is_measure());
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 1));
    auto p3 = bs->add_pulse(1, 4, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64,
                                                     IR::Type::Float64});
                                auto trel = builder.createFDiv(0, builder.getConst(1000));
                                auto vold = builder.createMul(1, builder.createSub(
                                                                  builder.getConst(1), trel));
                                auto vnew = builder.createMul(2, trel);
                                builder.createRet(builder.createAdd(vold, vnew));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1),
                                        Seq::Arg::create_var(m1->val())}, 2);
                            }());
    assert(p3->id() == 4);
    assert(p3->needs_oldval());
    assert(bs->get_pulses(1)->size() == 4);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 3);
    assert(p1->start() == Seq::EventTime(0));
    assert(p2->start().tconst == 5);
    assert(p2->start().terms.size() == 1);
    assert(p2->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p2->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p3->start().tconst == 35);
    assert(p3->start().terms.size() == 2);
    assert(p3->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p3->start().terms[1].sign == Seq::EventTime::Pos);
    if (p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0)) {
        assert(p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
        assert(p3->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    }
    else {
        assert(p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 1));
        assert(p3->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 0));
    }
    assert(approx(p1->endval()->get_const().get<double>(), 0.3));
    assert(approx(p2->endval()->get_const().get<double>(), -0.2));
    assert(approx(p3->endval()->get_const().get<double>(), 0.285));
    assert(!p1->needs_oldval());
    assert(!p2->needs_oldval());
    assert(!p3->needs_oldval());
    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(p3->is_ordered);
    assert(!bs->needs_oldval(p1));
    assert(!bs->needs_oldval(p2));
    assert(!bs->needs_oldval(p3));
    assert(!p1->len());
    assert(p2->len());
    assert(p3->len());
    assert(p1->val()->is_const());
    assert(p2->val()->is_call());
    assert(!p2->val()->argument_unused(0));
    assert(!p3->val()->argument_unused(0));
    // The start value is exact since we set it directly.
    // No approx comparison needed.
    assert(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    assert(approx(bs->endval(1)->get_const().get<double>(), 0.285));
}

static void test_basic_single_channel_measure_known_time2(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1000)));
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 2));
    auto p2 = bs->add_pulse(1, 3, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64,
                                                     IR::Type::Float64});
                                auto trel = builder.createFDiv(0, builder.getConst(1000));
                                auto vold = builder.createMul(1, builder.createSub(
                                                                  builder.getConst(1), trel));
                                auto vnew = builder.createMul(2, trel);
                                builder.createRet(builder.createAdd(vold, vnew));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1),
                                        Seq::Arg::create_var(m1->val())}, 2);
                            }());
    assert(p2->id() == 3);
    assert(p2->needs_oldval());
    assert(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 2);

    assert(p1->start().tconst == 5);
    assert(p1->start().terms.size() == 1);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p2->start().tconst == 1005);
    assert(p2->start().terms.size() == 3);
    assert(approx(p1->endval()->get_const().get<double>(), 1.8));
    assert(approx(p2->endval()->get_const().get<double>(), 1.8));
    assert(!p1->needs_oldval());
    assert(!p2->needs_oldval());
    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(!bs->needs_oldval(p1));
    assert(!bs->needs_oldval(p2));
    assert(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    assert(approx(bs->endval(1)->get_const().get<double>(), 1.8));
}

static void test_basic_single_channel_measure_unknown_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 2));
    auto p2 = bs->add_pulse(1, 3, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64,
                                                     IR::Type::Float64});
                                auto trel = builder.createFDiv(0, builder.getConst(1000));
                                auto vold = builder.createMul(1, builder.createSub(
                                                                  builder.getConst(1), trel));
                                auto vnew = builder.createMul(2, trel);
                                builder.createRet(builder.createAdd(vold, vnew));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1),
                                        Seq::Arg::create_var(m1->val())}, 2);
                            }());
    assert(p2->id() == 3);
    assert(p2->needs_oldval());
    assert(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 2);

    assert(p1->start().tconst == 5);
    assert(p1->start().terms.size() == 1);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p2->start().tconst == 5);
    assert(p2->start().terms.size() == 3);
    assert(approx(p1->endval()->get_const().get<double>(), 1.8));
    assert(p2->endval()->is_call());
    assert(!p1->needs_oldval());
    assert(!p2->needs_oldval());
    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(!bs->needs_oldval(p1));
    assert(!bs->needs_oldval(p2));
    assert(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    assert(bs->endval(1) == p2->endval());

    auto endval = p2->endval();
    assert(endval->args().size() == 1);
    assert(endval->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(endval->get_callee().is_llvm);

    auto var1 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(500))})->ref();
    auto var2 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(1000))})->ref();
    auto var3 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(1500))})->ref();

    seq.env().optimize();
    assert(var1->is_const());
    assert(approx(var1->get_const().get<double>(), 2.05));
    assert(var2->is_const());
    assert(approx(var2->get_const().get<double>(), 1.8));
    assert(var3->is_const());
    assert(approx(var3->get_const().get<double>(), 1.8));
}

static void test_basic_single_channel_unused_measure(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Unknown, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    assert(bs->get_pulses(1)->size() == 2);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 1);

    assert(p1->start().tconst == 5);
    assert(p1->start().terms.size() == 1);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p1->is_ordered);
    assert(approx(p1->endval()->get_const().get<double>(), 1.8));
    assert(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    assert(approx(bs->endval(1)->get_const().get<double>(), 1.8));
}

static void test_basic_two_channel_measure_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    assert(seq.get_chn_id("test_chn2", true) == 2);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1000)));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(10)));
    t.add_term(Seq::EventTime::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(1234)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }());
    auto p2 = bs->add_pulse(2, 3, bs->track_time(t), nullptr, m1->val());
    assert(p2->id() == 3);
    assert(!p2->needs_oldval());
    assert(bs->get_pulses(1)->size() == 2);
    assert(bs->get_pulses(2)->size() == 1);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 1);
    assert(bs->get_pulses(2)->size() == 1);

    assert(p1->start().tconst == 5);
    assert(p1->start().terms.size() == 1);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(p2->start().tconst == 3236); // 1000 + 10 + 5 + 1.8 * 1234
    assert(p2->start().terms.size() == 1);
    assert(p2->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p2->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(!p1->needs_oldval());
    assert(!bs->needs_oldval(p1));
    assert(p1->val()->is_call());
    assert(!p2->needs_oldval());
    assert(!bs->needs_oldval(p2));
    assert(p2->val()->is_const());
    assert(approx(p2->val()->get_const().get<double>(), 1.8));
    assert(p1->endval()->is_const());
    assert(approx(p1->endval()->get_const().get<double>(), 1.8));
    assert(approx(p2->endval()->get_const().get<double>(), 1.8));
    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    assert(bs->startval(2)->get_const().is(IR::TagVal(0.0)));
    assert(approx(bs->endval(1)->get_const().get<double>(), 1.8));
    assert(approx(bs->endval(2)->get_const().get<double>(), 1.8));
}

static void test_basic_single_channel_measure_unknown_order(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    seq.set_defval(1, 2.3);

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    Seq::EventTime t2;
    t2.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1200)));

    // Test that we do not optimize measure when we cannot determine the order
    // Note that in this case p1 and p2 commutes
    // so we can in theory determine the measure value.
    // This test needs to be fixed to account for that if we acquires such capability.
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t2), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0006));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p2->id() == 2);
    assert(p2->needs_oldval());

    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1200)));
    auto m1 = bs->add_measure(1, 3, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 3);
    assert(m1->is_measure());

    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 1));
    auto p3 = bs->add_pulse(1, 4, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64,
                                                     IR::Type::Float64});
                                auto trel = builder.createFDiv(0, builder.getConst(1000));
                                auto vold = builder.createMul(1, builder.createSub(
                                                                  builder.getConst(1), trel));
                                auto vnew = builder.createMul(2, trel);
                                builder.createRet(builder.createAdd(vold, vnew));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1),
                                        Seq::Arg::create_var(m1->val())}, 2);
                            }());
    assert(p3->id() == 4);
    assert(p3->needs_oldval());
    assert(bs->get_pulses(1)->size() == 4);

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1)->size() == 4);
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(bs->needs_oldval(p1));
    assert(p1->endval());
    assert(p2->id() == 2);
    assert(!p2->needs_oldval());
    assert(bs->needs_oldval(p2));
    assert(p2->endval());
    assert(m1->id() == 3);
    assert(m1->is_measure());
    assert(p3->id() == 4);
    assert(!p3->needs_oldval());
    assert(bs->needs_oldval(p3));
    assert(p3->endval());
    assert(!p1->is_ordered);
    assert(!p2->is_ordered);
    assert(p3->is_ordered);
    assert(bs->endval(1));
}

static void test_basic_single_channel_measure_time_violation(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t;
    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.05));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    t.add_term(Seq::EventTime::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(1234)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 2);
    auto p2 = bs->add_pulse(1, 3, bs->track_time(t), nullptr, m1->val());
    assert(p2->id() == 3);
    assert(!p2->needs_oldval());
    assert(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    assert(err.code >> 16 == Seq::Error::EventTime);
    assert((err.code & ((1u << 16) - 1)) == Seq::EventTime::NonPosTime);
    assert(err.type1 == Seq::Error::EventTime);
    assert(err.id1 == 2);
    assert(strcmp(err.what(), "Positive time expected.") == 0);
}

static void test_basic_single_channel_output_unknown_order(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    seq.set_defval(1, 2.3);

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0));
    Seq::EventTime t2;
    t2.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1200)));

    // Test that we do not compute end value of basic sequence
    // when we cannot determine which one is the last output pulse.
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr,
                            seq.get_const(IR::TagVal(1.0)));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t2), nullptr,
                            seq.get_const(IR::TagVal(2.0)));
    assert(p2->id() == 2);
    assert(!p2->needs_oldval());

    seq.prepare();
    seq.optimize();
    assert(bs->get_pulses(1)->size() == 2);
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs->needs_oldval(p1));
    assert(p1->endval()->is_const());
    assert(p2->id() == 2);
    assert(!p2->needs_oldval());
    assert(!bs->needs_oldval(p2));
    assert(p2->endval()->is_const());
    assert(!p1->is_ordered);
    assert(!p2->is_ordered);
    assert(!bs->endval(1));
}

static void test_basic_single_channel_neg_time_offset(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t(-10);
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 29);
    t.add_term(Seq::EventTime::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(-1000)));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());

    seq.prepare();
    seq.optimize();

    assert(bs->get_pulses(1));
    assert(bs->get_pulses(1)->size() == 1);
    assert(&*bs->get_pulses(1)->begin() == p1);
    assert(p1->start().tconst == -1010);
    assert(p1->start().terms.size() == 1);
    assert(p1->start().terms[0].sign == Seq::EventTime::Pos);
    assert(p1->start().terms[0].id == 29);
    assert(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
}

static void test_basic_single_channel_neg_time_output_violation(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t(10);
    t.add_term(Seq::EventTime::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto p1 = bs->add_pulse(1, 14, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(-1000)));
    assert(p1->id() == 14);
    assert(!p1->needs_oldval());

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    assert(err.code >> 16 == Seq::Error::Pulse);
    assert((err.code & ((1u << 16) - 1)) == Seq::Pulse::NegTime);
    assert(err.type1 == Seq::Error::Pulse);
    assert(err.id1 == 14);
    assert(strcmp(err.what(), "Pulse time must be positive.") == 0);
}

static void test_basic_single_channel_neg_time_measure_violation(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);
    Seq::EventTime t(10);
    t.add_term(Seq::EventTime::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto m1 = bs->add_measure(1, 19, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 19);
    assert(m1->is_measure());
    // Make sure the value is used.
    // Otherwise we don't guarantee the throwing of the error.
    seq.get_slot(IR::Type::Float64, 0);
    bs->assign_global(0, m1->val(), 20);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    assert(err.code >> 16 == Seq::Error::Pulse);
    assert((err.code & ((1u << 16) - 1)) == Seq::Pulse::NegTime);
    assert(err.type1 == Seq::Error::Pulse);
    assert(err.id1 == 19);
    assert(strcmp(err.what(), "Pulse time must be positive.") == 0);
}

static void test_single_channel_branch(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    Seq::EventTime t;
    seq.set_defval(1, 2.3);
    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(t), seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0008));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());

    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);
    auto p2 = bs2->add_pulse(1, 1, bs2->track_time(Seq::EventTime(0)),
                             seq.get_const(IR::TagVal(2000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.0002));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 1);
    assert(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    assert(bs1->get_default_branch() == bs2);
    assert(bs2->get_default_branch() == nullptr);

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs2->get_pulses(1)->size() == 1);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(approx(p1->endval()->get_const().get<double>(), 1.5));
    assert(p2->id() == 1);
    assert(!p2->needs_oldval());
    assert(!bs2->needs_oldval(p2));
    assert(p2->endval());
    assert(approx(p2->endval()->get_const().get<double>(), 1.9));
    assert(p1->is_ordered);
    assert(p2->is_ordered);

    assert(approx(bs1->startval(1)->get_const().get<double>(), 2.3));
    assert(approx(bs1->endval(1)->get_const().get<double>(), 1.5));
    assert(approx(bs2->startval(1)->get_const().get<double>(), 1.5));
    assert(approx(bs2->endval(1)->get_const().get<double>(), 1.9));
}

static void test_single_channel_branch_merge(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);
    auto bs4 = seq.add_basicseq(4);
    assert(bs4->id() == 4);

    seq.set_defval(1, 2.3);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.003));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs3, 1);

    auto p2 = bs2->add_pulse(1, 1, bs2->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.0008));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 1);
    assert(p2->needs_oldval());
    bs2->set_default_branch(bs4);

    auto p3 = bs3->add_pulse(1, 1, bs3->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.0008));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p3->id() == 1);
    assert(p3->needs_oldval());
    bs3->set_default_branch(bs4);

    auto p4 = bs4->add_pulse(1, 1, bs4->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.001));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p4->id() == 1);
    assert(p4->needs_oldval());
    // no-op but should work
    bs4->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs2->get_pulses(1)->size() == 1);
    assert(bs3->get_pulses(1)->size() == 1);
    assert(bs4->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 4);

    assert(bs1->get_default_branch() == bs2);
    assert(bs2->get_default_branch() == bs4);
    assert(bs3->get_default_branch() == bs4);
    assert(!bs4->get_default_branch());

    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs3);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs2->get_branches().empty());
    assert(bs3->get_branches().empty());
    assert(bs4->get_branches().empty());

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(approx(p1->endval()->get_const().get<double>(), -0.7));

    assert(p2->id() == 1);
    assert(!p2->needs_oldval());
    assert(!bs2->needs_oldval(p2));
    assert(p2->endval());
    assert(approx(p2->endval()->get_const().get<double>(), 0.1));

    assert(p3->id() == 1);
    assert(!p3->needs_oldval());
    assert(!bs3->needs_oldval(p3));
    assert(p3->endval());
    assert(approx(p3->endval()->get_const().get<double>(), 0.1));

    assert(p2->endval()->get_const().get<double>() == p3->endval()->get_const().get<double>());

    assert(p4->id() == 1);
    assert(!p4->needs_oldval());
    assert(!bs4->needs_oldval(p4));
    assert(p4->endval());
    assert(approx(p4->endval()->get_const().get<double>(), 1.1));

    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(p3->is_ordered);
    assert(p4->is_ordered);

    assert(approx(bs1->startval(1)->get_const().get<double>(), 2.3));
    assert(approx(bs1->endval(1)->get_const().get<double>(), -0.7));

    assert(approx(bs2->startval(1)->get_const().get<double>(), -0.7));
    assert(approx(bs2->endval(1)->get_const().get<double>(), 0.1));
    assert(approx(bs3->startval(1)->get_const().get<double>(), -0.7));
    assert(approx(bs3->endval(1)->get_const().get<double>(), 0.1));

    assert(approx(bs4->startval(1)->get_const().get<double>(), 0.1));
    assert(approx(bs4->endval(1)->get_const().get<double>(), 1.1));
}

static void test_single_channel_branch_not_merge(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);
    auto bs4 = seq.add_basicseq(4);
    assert(bs4->id() == 4);

    seq.set_defval(1, 2.3);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.003));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs3, 1);

    auto p2 = bs2->add_pulse(1, 1, bs2->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.0008));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 1);
    assert(p2->needs_oldval());
    bs2->set_default_branch(bs4);

    auto p3 = bs3->add_pulse(1, 1, bs3->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0008));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p3->id() == 1);
    assert(p3->needs_oldval());
    bs3->set_default_branch(bs4);

    auto p4 = bs4->add_pulse(1, 1, bs4->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(0.001));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p4->id() == 1);
    assert(p4->needs_oldval());

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs2->get_pulses(1)->size() == 1);
    assert(bs3->get_pulses(1)->size() == 1);
    assert(bs4->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 4);

    assert(bs1->get_default_branch() == bs2);
    assert(bs2->get_default_branch() == bs4);
    assert(bs3->get_default_branch() == bs4);
    assert(!bs4->get_default_branch());

    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs3);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs2->get_branches().empty());
    assert(bs3->get_branches().empty());
    assert(bs4->get_branches().empty());

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs2->get_pulses(1)->size() == 1);
    assert(bs3->get_pulses(1)->size() == 1);
    assert(bs4->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 4);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(approx(p1->endval()->get_const().get<double>(), -0.7));

    assert(p2->id() == 1);
    assert(!p2->needs_oldval());
    assert(!bs2->needs_oldval(p2));
    assert(p2->endval());
    assert(approx(p2->endval()->get_const().get<double>(), 0.1));

    assert(p3->id() == 1);
    assert(!p3->needs_oldval());
    assert(!bs3->needs_oldval(p3));
    assert(p3->endval());
    assert(approx(p3->endval()->get_const().get<double>(), -1.5));

    assert(p4->id() == 1);
    assert(!p4->needs_oldval());
    assert(!bs4->needs_oldval(p4));
    assert(p4->endval());

    assert(p1->is_ordered);
    assert(p2->is_ordered);
    assert(p3->is_ordered);
    assert(p4->is_ordered);

    assert(approx(bs1->startval(1)->get_const().get<double>(), 2.3));
    assert(approx(bs1->endval(1)->get_const().get<double>(), -0.7));

    assert(approx(bs2->startval(1)->get_const().get<double>(), -0.7));
    assert(approx(bs2->endval(1)->get_const().get<double>(), 0.1));
    assert(approx(bs3->startval(1)->get_const().get<double>(), -0.7));
    assert(approx(bs3->endval(1)->get_const().get<double>(), -1.5));

    assert(bs4->startval(1)->is_extern());
    assert(bs4->endval(1));
}

static void test_single_channel_branch_nonconst(llvm::LLVMContext &llvm_ctx)
{
    // Non-const cannot be forwarded (in general)
    // since it most likely either depend on a measure,
    // which has the wrong scope,
    // or a global variable, which may be reassigned between sequence.
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             seq.get_slot(IR::Type::Float64, 0));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    bs1->set_default_branch(bs2);

    auto p2 = bs2->add_pulse(1, 1, bs2->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 builder.createRet(
                                     builder.createAdd(1, builder.getConst(0.0008)));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 1);
    assert(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs2->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 2);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(p1->endval() == seq.get_slot(IR::Type::Float64, 0));

    assert(p2->id() == 1);
    assert(!p2->needs_oldval());
    assert(!bs2->needs_oldval(p2));
    assert(p2->endval());
    assert(!p2->len());

    assert(p1->is_ordered);
    assert(p2->is_ordered);

    assert(bs1->startval(1)->get_const().is(IR::TagVal(0.0)));
    assert(bs1->endval(1));
    assert(bs1->endval(1) == seq.get_slot(IR::Type::Float64, 0));

    assert(bs2->startval(1)->is_extern());
    assert(bs2->endval(1));
}

static void test_branch_cond_clone_point(llvm::LLVMContext &llvm_ctx)
{
    // Branch condition that
    // * uses only global value: no cloning
    // * uses only measure: clone as a whole
    // * uses both global and measure: clone from interface
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 2, bs1->track_time(Seq::EventTime(1)), m1val);
    assert(m1->id() == 2);
    assert(m1->is_measure());
    assert(m1->val() == m1val);

    auto gv1 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createAdd(0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_const(0.2),
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1))}, 0);
    }();
    auto gv2 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_const(0.2),
                Seq::Arg::create_var(gv1)}, 0);
    }();

    auto lv1 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createAdd(0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_const(0.6),
                Seq::Arg::create_var(m1val)}, 0);
    }();
    auto lv2 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_const(0.4),
                Seq::Arg::create_var(lv1)}, 0);
    }();

    auto glv1 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createAdd(0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(lv2),
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1))}, 0);
    }();
    auto glv2 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {Seq::Arg::create_const(0),
                Seq::Arg::create_var(glv1)}, 0);
    }();
    bs1->add_branch(gv2, bs2, 1);
    bs1->add_branch(lv2, nullptr, 2);
    bs1->add_branch(glv2, bs1, 3);

    assert(bs1->get_branches().size() == 3);
    assert(bs1->get_branches()[0].target == bs2);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == gv2);
    assert(bs1->get_branches()[1].target == nullptr);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() == lv2);
    assert(bs1->get_branches()[2].target == bs1);
    assert(bs1->get_branches()[2].id == 3);
    assert(bs1->get_branches()[2].cond.get() == glv2);

    seq.prepare();

    assert(seq.get_slots().size() == 2);
    assert(bs1->get_branches().size() == 3);
    assert(bs1->get_branches()[0].target == bs2);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == gv2);
    assert(bs1->get_branches()[1].target == nullptr);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() != lv2);
    auto new_lv2 = bs1->get_branches()[1].cond.get();
    assert(new_lv2->is_extern());
    assert(new_lv2->get_extern().first == IR::Type::Float64);
    assert(new_lv2->get_extern().second == 2);
    assert(bs1->get_branches()[2].target == bs1);
    assert(bs1->get_branches()[2].id == 3);
    assert(bs1->get_branches()[2].cond.get() != glv2);
    auto new_glv2 = bs1->get_branches()[2].cond.get();
    assert(new_glv2->is_call());
    assert(new_glv2->get_callee().is_llvm);
    assert(new_glv2->get_callee().llvm == glv2->get_callee().llvm);
    assert(new_glv2->args().size() == 2);
    assert(new_glv2->args()[1].is_var());
    assert(new_glv2->args()[1].get_var() != glv1);
    auto new_glv1 = new_glv2->args()[1].get_var();
    assert(new_glv1->is_call());
    assert(new_glv1->get_callee().is_llvm);
    assert(new_glv1->get_callee().llvm == glv1->get_callee().llvm);
    assert(new_glv1->args().size() == 2);
    assert(new_glv1->args()[0].is_var());
    assert(new_glv1->args()[0].get_var() == new_lv2);
    assert(new_glv1->args()[1].is_var());
    assert(new_glv1->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 1));

    assert(bs1->get_assigns().size() == 1);
    assert(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(2)->second.id == 0);
    assert(bs1->get_assigns().find(2)->second.val.get() == lv2);

    seq.optimize();

    assert(seq.get_slots().size() == 3);
    assert(bs1->get_branches().size() == 3);
    assert(bs1->get_branches()[0].target == bs2);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == gv2);
    assert(bs1->get_branches()[1].target == nullptr);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 2));
    new_lv2 = seq.get_slot(IR::Type::Float64, 2);
    assert(bs1->get_branches()[2].target == bs1);
    assert(bs1->get_branches()[2].id == 3);
    assert(bs1->get_branches()[2].cond.get() == new_glv2);
    assert(new_glv2->is_call());
    assert(new_glv2->get_callee().is_llvm);
    assert(new_glv2->args().size() == 2);
    assert(new_glv2->args()[0].is_var());
    assert(new_glv2->args()[1].is_var());
    auto arg0 = new_glv2->args()[0].get_var();
    auto arg1 = new_glv2->args()[1].get_var();
    assert((arg0 == seq.get_slot(IR::Type::Float64, 1) &&
            arg1 == seq.get_slot(IR::Type::Float64, 2)) ||
           (arg0 == seq.get_slot(IR::Type::Float64, 2) &&
            arg1 == seq.get_slot(IR::Type::Float64, 1)));

    assert(bs1->get_assigns().size() == 1);
    assert(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(2)->second.id == 0);
    assert(bs1->get_assigns().find(2)->second.val.get() == lv2);
}

static void test_branch_cond_clone_const_opt(llvm::LLVMContext &llvm_ctx)
{
    // * If the cloned value used in a condition is a constant it can be forwarded.
    // * The forwarded condition should not affect the one in a different basic sequence
    //   with the same temporary ID.
    // * The cloned slot after it should be compacted down.
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);

    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 1, bs1->track_time(Seq::EventTime(0)), m1val);
    assert(m1->id() == 1);
    assert(m1->is_measure());
    assert(m1->val() == m1val); // optimizable to constant
    auto p1 = bs1->add_pulse(1, 2, bs1->track_time(Seq::EventTime(1)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    assert(p1->id() == 2);
    assert(!p1->needs_oldval());
    auto m2val = bs1->new_measure(seq.env(), 1);
    auto m2 = bs1->add_measure(1, 3, bs1->track_time(Seq::EventTime(2)), m2val);
    assert(m2->id() == 3);
    assert(m2->is_measure());
    assert(m2->val() == m2val); // not optimizable to constant (slot 0)

    auto glv1 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m1val)}, 0);
    }();
    auto glv2 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m2val)}, 0);
    }();

    bs1->add_branch(glv1, nullptr, 1);
    bs1->add_branch(glv2, bs2, 2);
    assert(bs1->get_branches().size() == 2);
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == glv1);
    assert(bs1->get_branches()[1].target == bs2);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() == glv2);

    auto p2 = bs2->add_pulse(1, 4, bs2->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    assert(p2->id() == 4);
    assert(!p2->needs_oldval());
    auto m3val = bs2->new_measure(seq.env(), 1);
    auto m3 = bs2->add_measure(1, 5, bs2->track_time(Seq::EventTime(1)), m3val);
    assert(m3->id() == 5);
    assert(m3->is_measure());
    assert(m3->val() == m3val);

    auto glv3 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m3val)}, 0);
    }();
    bs2->add_branch(glv3, bs2, 3);
    assert(bs2->get_branches().size() == 1);
    assert(bs2->get_branches()[0].target == bs2);
    assert(bs2->get_branches()[0].id == 3);
    assert(bs2->get_branches()[0].cond.get() == glv3);

    seq.prepare();

    assert(seq.get_slots().size() == 2);

    assert(bs1->get_branches().size() == 2);
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() != glv1);
    auto new_glv1 = bs1->get_branches()[0].cond.get();
    assert(new_glv1->is_call());
    assert(new_glv1->get_callee().is_llvm);
    assert(new_glv1->get_callee().llvm == glv1->get_callee().llvm);
    assert(new_glv1->args().size() == 2);
    assert(new_glv1->args()[0].is_var());
    assert(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv1->args()[1].is_var());
    auto new_m1var = new_glv1->args()[1].get_var();
    assert(new_m1var->is_extern());
    assert(new_m1var->get_extern().first == IR::Type::Float64);
    assert(new_m1var->get_extern().second == 2);
    assert(bs1->get_branches()[1].target == bs2);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() != glv2);
    auto new_glv2 = bs1->get_branches()[1].cond.get();
    assert(new_glv2->is_call());
    assert(new_glv2->get_callee().is_llvm);
    assert(new_glv2->get_callee().llvm == glv2->get_callee().llvm);
    assert(new_glv2->args().size() == 2);
    assert(new_glv2->args()[0].is_var());
    assert(new_glv2->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv2->args()[1].is_var());
    auto new_m2var = new_glv2->args()[1].get_var();
    assert(new_m2var->is_extern());
    assert(new_m2var->get_extern().first == IR::Type::Float64);
    assert(new_m2var->get_extern().second == 3);

    assert(bs1->get_assigns().size() == 2);
    assert(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(2)->second.id == 0);
    assert(bs1->get_assigns().find(2)->second.val.get() == m1val);
    assert(bs1->get_assigns().find(3) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(3)->second.id == 0);
    assert(bs1->get_assigns().find(3)->second.val.get() == m2val);

    assert(bs2->get_branches().size() == 1);
    assert(bs2->get_branches()[0].target == bs2);
    assert(bs2->get_branches()[0].id == 3);
    assert(bs2->get_branches()[0].cond.get() != glv3);
    auto new_glv3 = bs2->get_branches()[0].cond.get();
    assert(new_glv3->is_call());
    assert(new_glv3->get_callee().is_llvm);
    assert(new_glv3->get_callee().llvm == glv3->get_callee().llvm);
    assert(new_glv3->args().size() == 2);
    assert(new_glv3->args()[0].is_var());
    assert(new_glv3->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv3->args()[1].is_var());
    auto new_m3var = new_glv3->args()[1].get_var();
    assert(new_m3var->is_extern());
    assert(new_m3var->get_extern().first == IR::Type::Float64);
    assert(new_m3var->get_extern().second == 2);

    assert(bs2->get_assigns().size() == 1);
    assert(bs2->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs2->get_assigns().find(2)->second.id == 0);
    assert(bs2->get_assigns().find(2)->second.val.get() == m3val);

    seq.optimize();

    assert(seq.get_slots().size() == 3);

    assert(bs1->get_branches().size() == 2);
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == new_glv1);
    assert(new_glv1->is_call());
    assert(new_glv1->args().size() == 1);
    assert(new_glv1->get_callee().is_llvm);
    assert(new_glv1->args()[0].is_var());
    assert(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs1->get_branches()[1].target == bs2);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[1].cond.get() == new_glv2);
    assert(new_glv2->is_call());
    assert(new_glv2->get_callee().is_llvm);
    assert(new_glv2->args()[0].is_var());
    assert(new_glv2->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv2->args()[1].is_var());
    assert(new_glv2->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 2));

    assert(bs1->get_assigns().size() == 1);
    assert(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(2)->second.id == 0);
    assert(bs1->get_assigns().find(2)->second.val.get() == seq.get_slot(IR::Type::Float64, 0));

    assert(bs2->get_branches().size() == 1);
    assert(bs2->get_branches()[0].target == bs2);
    assert(bs2->get_branches()[0].id == 3);
    assert(bs2->get_branches()[0].cond.get() == new_glv3);
    assert(new_glv3->is_call());
    assert(new_glv3->get_callee().is_llvm);
    assert(new_glv3->args()[0].is_var());
    assert(new_glv3->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv3->args()[1].is_var());
    assert(new_glv3->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 2));

    assert(bs2->get_assigns().size() == 1);
    assert(bs2->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs2->get_assigns().find(2)->second.id == 0);
    assert(bs2->get_assigns().find(2)->second.val.get() == seq.get_slot(IR::Type::Float64, 0));
}

static void test_branch_cond_clone_unused_opt(llvm::LLVMContext &llvm_ctx)
{
    // Unused slot will be optimized out.
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 2, bs1->track_time(Seq::EventTime(1)), m1val);
    assert(m1->id() == 2);
    assert(m1->is_measure());
    assert(m1->val() == m1val);

    auto glv1 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m1val)}, 0);
    }();

    bs1->add_branch(glv1, nullptr, 1);
    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() == glv1);

    seq.prepare();

    assert(seq.get_slots().size() == 2);

    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[0].cond.get() != glv1);
    auto new_glv1 = bs1->get_branches()[0].cond.get();
    assert(new_glv1->is_call());
    assert(new_glv1->get_callee().is_llvm);
    assert(new_glv1->get_callee().llvm == glv1->get_callee().llvm);
    assert(new_glv1->args().size() == 2);
    assert(new_glv1->args()[0].is_var());
    assert(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    assert(new_glv1->args()[1].is_var());
    auto new_m1var = new_glv1->args()[1].get_var();
    assert(new_m1var->is_extern());
    assert(new_m1var->get_extern().first == IR::Type::Float64);
    assert(new_m1var->get_extern().second == 2);

    assert(bs1->get_assigns().size() == 1);
    assert(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    assert(bs1->get_assigns().find(2)->second.id == 0);
    assert(bs1->get_assigns().find(2)->second.val.get() == m1val);

    seq.optimize();

    assert(seq.get_slots().size() == 2);

    assert(bs1->get_branches().empty());
    assert(bs1->get_assigns().empty());
}

static void test_branch_loop_const(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    assert(seq.get_chn_id("test_chn2", true) == 2);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);

    seq.set_defval(1, 2.3);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.003));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs2);

    auto m1 = bs2->add_measure(1, 1, bs2->track_time(Seq::EventTime(100)),
                               bs2->new_measure(seq.env(), 1));
    assert(m1->id() == 1);
    assert(m1->is_measure());
    auto p2 = bs2->add_pulse(2, 2, bs2->track_time(Seq::EventTime(200)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64,
                                                      IR::Type::Float64});
                                 auto v = builder.createMul(0, builder.getConst(0.0009));
                                 builder.createRet(builder.createAdd(v, 2));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1),
                                         Seq::Arg::create_var(m1->val())}, 2);
                             }());
    assert(p2->id() == 2);
    assert(p2->needs_oldval());
    bs2->set_default_branch(bs2);
    bs2->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    assert(bs1->has_output(1));
    assert(!bs1->has_output(2));
    assert(!bs2->has_output(1));
    assert(bs2->has_output(2));

    seq.prepare();
    seq.optimize();

    assert(bs1->has_output(1));
    assert(!bs1->has_output(2));
    assert(!bs2->has_output(1));
    assert(bs2->has_output(2));

    assert(bs1->get_pulses(1)->size() == 1);
    assert(bs1->get_pulses(2)->empty());
    assert(bs2->get_pulses(1)->empty());
    assert(bs2->get_pulses(2)->size() == 1);
    assert(seq.get_basicseqs().size() == 2);

    assert(bs1->get_default_branch() == bs2);
    assert(bs2->get_default_branch() == bs2);

    assert(bs1->get_branches().empty());
    assert(bs2->get_branches().size() == 1);
    assert(bs2->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs2->get_branches()[0].target == nullptr);
    assert(bs2->get_branches()[0].id == 10);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(approx(p1->endval()->get_const().get<double>(), -0.7));

    assert(p2->id() == 2);
    assert(!p2->needs_oldval());
    assert(!bs2->needs_oldval(p2));
    assert(p2->endval());
    assert(approx(p2->endval()->get_const().get<double>(), 0.2));

    assert(p1->is_ordered);
    assert(p2->is_ordered);

    assert(approx(bs1->startval(1)->get_const().get<double>(), 2.3));
    assert(approx(bs1->endval(1)->get_const().get<double>(), -0.7));
    assert(bs1->startval(2)->get_const().is(IR::TagVal(0.0)));
    assert(bs1->endval(2)->get_const().is(IR::TagVal(0.0)));

    assert(approx(bs2->startval(1)->get_const().get<double>(), -0.7));
    assert(approx(bs2->endval(1)->get_const().get<double>(), -0.7));
    assert(bs2->startval(2)->is_extern());
    assert(approx(bs2->endval(2)->get_const().get<double>(), 0.2));
}

static void test_single_channel_loop_first(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    seq.set_defval(1, 2.3);
    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto v = builder.createSub(1, builder.getConst(2.0));
                                 builder.createRet(builder.createCall(IR::Builtins::abs, {v}));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 1);

    assert(bs1->get_default_branch() == bs1);
    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 10);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(!p1->len());
    assert(p1->is_ordered);

    assert(bs1->startval(1)->is_extern());
    assert(bs1->endval(1));
}

#if 0
// Currently the optimizer does not look into each pulse to do constant propagation
// so this optimization won't happen.
static void test_single_channel_loop_first2(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    seq.set_defval(1, 1.0);
    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto v = builder.createSub(1, builder.getConst(2.0));
                                 builder.createRet(builder.createCall(IR::Builtins::abs, {v}));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 1);

    assert(bs1->get_default_branch() == bs1);
    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 10);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(p1->endval()->is_const());
    assert(p1->endval()->get_const().is(IR::TagVal(1.0)));
    assert(!p1->len());
    assert(p1->is_ordered);

    assert(bs1->startval(1)->is_const());
    assert(bs1->startval(1)->get_const().is(IR::TagVal(1.0)));
    assert(bs1->endval(1)->is_const());
    assert(bs1->endval(1)->get_const().is(IR::TagVal(1.0)));
}
#endif

static void test_single_channel_loop_first3(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    seq.set_defval(1, 1.0);
    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto v = builder.createSub(0, builder.getConst(1000.0));
                                 builder.createRet(builder.createAdd(1, v));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    assert(bs1->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 1);

    assert(bs1->get_default_branch() == bs1);
    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == nullptr);
    assert(bs1->get_branches()[0].id == 10);

    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(p1->endval()->is_const());
    assert(p1->endval()->get_const().is(IR::TagVal(1.0)));
    assert(p1->len());
    assert(p1->is_ordered);

    assert(bs1->startval(1)->is_const());
    assert(bs1->startval(1)->get_const().is(IR::TagVal(1.0)));
    assert(bs1->endval(1)->is_const());
    assert(bs1->endval(1)->get_const().is(IR::TagVal(1.0)));
}

static void test_optimize_cfg_shortcut(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs2, 2);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs3, 3);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    assert(seq.get_basicseqs().size() == 2);
    assert(bs1->id() == 1);
    assert(bs2->id() == 2);

    assert(bs1->get_default_branch() == bs2);
    assert(!bs2->get_default_branch());

    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs1);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs2->get_branches().empty());
}

static void test_optimize_cfg_delete(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_const(IR::TagVal(0)), bs2, 2);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs3, 3);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    assert(seq.get_basicseqs().size() == 2);
    assert(bs1->id() == 1);
    assert(bs3->id() == 3);

    assert(bs1->get_default_branch() == bs3);
    assert(!bs3->get_default_branch());

    assert(bs1->get_branches().size() == 1);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs1);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs3->get_branches().empty());
}

static void test_optimize_cfg_merge(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 1), bs3, 2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 2), bs2, 3);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 3), bs2, 4);
    bs1->set_default_branch(bs2);

    seq.prepare();
    seq.optimize();

    assert(seq.get_basicseqs().size() == 3);
    assert(bs1->id() == 1);
    assert(bs2->id() == 2);
    assert(bs3->id() == 3);

    assert(bs1->get_default_branch() == bs2);
    assert(!bs2->get_default_branch());
    assert(!bs3->get_default_branch());

    assert(bs1->get_branches().size() == 2);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs1);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs1->get_branches()[1].target == bs3);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs2->get_branches().empty());
    assert(bs3->get_branches().empty());
}

static void test_optimize_cfg_merge_null(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    assert(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 1), bs3, 2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 2), bs2, 3);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 3), nullptr, 4);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    assert(seq.get_basicseqs().size() == 3);
    assert(bs1->id() == 1);
    assert(bs2->id() == 2);
    assert(bs3->id() == 3);

    assert(!bs1->get_default_branch());
    assert(!bs2->get_default_branch());
    assert(!bs3->get_default_branch());

    assert(bs1->get_branches().size() == 3);
    assert(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs1->get_branches()[0].target == bs1);
    assert(bs1->get_branches()[0].id == 1);
    assert(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs1->get_branches()[1].target == bs3);
    assert(bs1->get_branches()[1].id == 2);
    assert(bs1->get_branches()[2].cond.get() == seq.get_slot(IR::Type::Float64, 2));
    assert(bs1->get_branches()[2].target == bs2);
    assert(bs1->get_branches()[2].id == 3);
    assert(bs2->get_branches().empty());
    assert(bs3->get_branches().empty());
}

static void test_basic_single_channel_assumption(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t;
    seq.set_defval(1, 23);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.024));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    bs->add_assume(Seq::EventTime::Pos, m1->val(), 43);
    bs->add_assume(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 3);

    seq.prepare();
    seq.optimize();

    assert(bs->get_assumes().size() == 1);
    assert(bs->get_assumes().begin()->second.sign == Seq::EventTime::Pos);
    assert(bs->get_assumes().begin()->first.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_assumes().begin()->second.id == 3);
}

static void test_basic_single_channel_assumption_merge(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t;
    seq.set_defval(1, 23);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr,
                            seq.get_slot(IR::Type::Float64, 0));
    assert(p1->id() == 1);
    assert(!p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    bs->add_assume(Seq::EventTime::Pos, m1->val(), 43);
    bs->add_assume(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 0), 3);

    seq.prepare();
    seq.optimize();

    assert(bs->get_assumes().size() == 1);
    assert(bs->get_assumes().begin()->second.sign == Seq::EventTime::Pos);
    assert(bs->get_assumes().begin()->first.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_assumes().begin()->second.id == 43);
}

static void test_basic_single_channel_assumption_violation(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t;
    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    assert(p1->id() == 1);
    assert(p1->needs_oldval());
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    assert(m1->id() == 2);
    assert(m1->is_measure());
    // `m1->val()` should be `0.14` but will be rounded to `0`
    // to be consistent with the time check.
    bs->add_assume(Seq::EventTime::Pos, m1->val(), 43);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    assert(err.code >> 16 == Seq::Error::EventTime);
    assert((err.code & ((1u << 16) - 1)) == Seq::EventTime::NonPosTime);
    assert(err.type1 == Seq::Error::EventTime);
    assert(err.id1 == 43);
    assert(strcmp(err.what(), "Positive time expected.") == 0);
}

static void test_endtimes_normalize(llvm::LLVMContext &llvm_ctx)
{
    // Make sure we correctly nomalize the time even if we have only one to start with.
    Seq::Seq seq(llvm_ctx);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1), 8);
    bs->add_endtime(bs->track_time(t));

    assert(bs->get_endtimes().size() == 1);

    seq.prepare();
    seq.optimize();

    assert(bs->get_endtimes().size() == 1);
    assert(bs->get_endtimes().front()->tconst == 11);
    assert(bs->get_endtimes().front()->terms.size() == 2);
    assert(bs->get_endtimes().front()->terms[0].sign == Seq::EventTime::Pos);
    assert(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_endtimes().front()->terms[0].id == 4);
    assert(bs->get_endtimes().front()->terms[1].sign == Seq::EventTime::NonNeg);
    assert(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs->get_endtimes().front()->terms[1].id == 8);
}

static void test_endtimes_elim_diamond(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    bs->add_endtime(bs->track_time(t));
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    bs->add_endtime(bs->track_time(t));
    Seq::EventTime t2(t);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    bs->add_endtime(bs->track_time(t));
    t2.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1), 6);
    bs->add_endtime(bs->track_time(t2));
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1), 8);
    bs->add_endtime(bs->track_time(t));

    assert(bs->get_endtimes().size() == 5);

    seq.prepare();
    seq.optimize();

    assert(bs->get_endtimes().size() == 1);
    assert(bs->get_endtimes().front()->tconst == 11);
    assert(bs->get_endtimes().front()->terms.size() == 2);
    assert(bs->get_endtimes().front()->terms[0].sign == Seq::EventTime::Pos);
    assert(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_endtimes().front()->terms[0].id == 4);
    assert(bs->get_endtimes().front()->terms[1].sign == Seq::EventTime::NonNeg);
    assert(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs->get_endtimes().front()->terms[1].id == 8);
}

static void test_endtimes_elim_Y(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    bs->add_endtime(bs->track_time(t));
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    bs->add_endtime(bs->track_time(t));
    Seq::EventTime t2(t);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    bs->add_endtime(bs->track_time(t));
    t2.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1), 6);
    bs->add_endtime(bs->track_time(t2));

    assert(bs->get_endtimes().size() == 4);

    seq.prepare();
    seq.optimize();

    assert(bs->get_endtimes().size() == 2);

    // t2 having smaller tconst will be sorted to the front.
    assert(bs->get_endtimes().front()->tconst == 6);
    assert(bs->get_endtimes().front()->terms.size() == 2);
    assert(bs->get_endtimes().front()->terms[0].sign == Seq::EventTime::Pos);
    assert(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_endtimes().front()->terms[0].id == 4);
    assert(bs->get_endtimes().front()->terms[1].sign == Seq::EventTime::NonNeg);
    assert(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs->get_endtimes().front()->terms[1].id == 6);

    assert(bs->get_endtimes().back()->tconst == 11);
    assert(bs->get_endtimes().back()->terms.size() == 1);
    assert(bs->get_endtimes().back()->terms[0].sign == Seq::EventTime::Pos);
    assert(bs->get_endtimes().back()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_endtimes().back()->terms[0].id == 4);
}

static void test_endtimes_elim_repeat(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    t.add_term(Seq::EventTime::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    t.add_term(Seq::EventTime::NonNeg, seq.get_slot(IR::Type::Float64, 1), 7);
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));

    assert(bs->get_endtimes().size() == 5);

    seq.prepare();
    seq.optimize();

    assert(bs->get_endtimes().size() == 1);

    // t2 having smaller tconst will be sorted to the front.
    assert(bs->get_endtimes().front()->tconst == 11);
    assert(bs->get_endtimes().front()->terms.size() == 2);
    assert(bs->get_endtimes().front()->terms[0].sign == Seq::EventTime::Pos);
    assert(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    assert(bs->get_endtimes().front()->terms[0].id == 4);
    assert(bs->get_endtimes().front()->terms[1].sign == Seq::EventTime::NonNeg);
    assert(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    assert(bs->get_endtimes().front()->terms[1].id == 7);
}

static void test_endtimes_elim_nonpos(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto bs = seq.add_basicseq(1);
    assert(bs->id() == 1);

    Seq::EventTime t(0);
    bs->add_endtime(bs->track_time(t));

    assert(bs->get_endtimes().size() == 1);

    seq.prepare();
    seq.optimize();

    assert(bs->get_endtimes().size() == 0);
}

int main()
{
    auto llvm_ctx = LLVM::new_context();

    test_basic_single_channel(*llvm_ctx);
    test_basic_single_channel_measure_known_time1(*llvm_ctx);
    test_basic_single_channel_measure_known_time2(*llvm_ctx);
    test_basic_single_channel_measure_unknown_time(*llvm_ctx);
    test_basic_single_channel_unused_measure(*llvm_ctx);
    test_basic_two_channel_measure_time(*llvm_ctx);
    test_basic_single_channel_measure_unknown_order(*llvm_ctx);
    test_basic_single_channel_measure_time_violation(*llvm_ctx);
    test_basic_single_channel_output_unknown_order(*llvm_ctx);
    test_basic_single_channel_neg_time_offset(*llvm_ctx);
    test_basic_single_channel_neg_time_output_violation(*llvm_ctx);
    test_basic_single_channel_neg_time_measure_violation(*llvm_ctx);

    test_single_channel_branch(*llvm_ctx);
    test_single_channel_branch_merge(*llvm_ctx);
    test_single_channel_branch_not_merge(*llvm_ctx);
    test_single_channel_branch_nonconst(*llvm_ctx);
    test_branch_cond_clone_point(*llvm_ctx);
    test_branch_cond_clone_const_opt(*llvm_ctx);
    test_branch_cond_clone_unused_opt(*llvm_ctx);
    test_branch_loop_const(*llvm_ctx);
    test_single_channel_loop_first(*llvm_ctx);
    // test_single_channel_loop_first2(*llvm_ctx);
    test_single_channel_loop_first3(*llvm_ctx);

    test_optimize_cfg_shortcut(*llvm_ctx);
    test_optimize_cfg_delete(*llvm_ctx);
    test_optimize_cfg_merge(*llvm_ctx);
    test_optimize_cfg_merge_null(*llvm_ctx);

    test_basic_single_channel_assumption(*llvm_ctx);
    test_basic_single_channel_assumption_merge(*llvm_ctx);
    test_basic_single_channel_assumption_violation(*llvm_ctx);

    test_endtimes_normalize(*llvm_ctx);
    test_endtimes_elim_diamond(*llvm_ctx);
    test_endtimes_elim_Y(*llvm_ctx);
    test_endtimes_elim_repeat(*llvm_ctx);
    test_endtimes_elim_nonpos(*llvm_ctx);

    return 0;
}
