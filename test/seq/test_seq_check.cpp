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

static void test_check_extern_measure_value(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto p2 = bs2->add_pulse(1, 7, bs2->track_time(Seq::EventTime(0)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_var(m1->val())}, 1);
                             }());
    assert(p2->id() == 7);
    assert(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(p2->id() == 7);

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Pulse);
    assert(err.id2 == p2->id());
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_length(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto len = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }();
    auto p2 = bs2->add_pulse(1, 5, bs2->track_time(Seq::EventTime(0)), len,
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 5);
    assert(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(p2->id() == 5);

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasureLength);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Pulse);
    assert(err.id2 == p2->id());
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 1289);
    auto p2 = bs2->add_pulse(1, 58, bs2->track_time(t), seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_arg(1)}, 2);
                             }());
    assert(p2->id() == 58);
    assert(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(p2->id() == 58);

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::EventTime);
    assert(err.id2 == 1289);
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_assign(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    seq.get_slot(IR::Type::Float64, 0);

    bs2->assign_global(0, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 1009);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Assignment);
    assert(err.id2 == 1009);
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_assume(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    bs2->add_assume(Seq::EventTime::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 1123);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::EventTime);
    assert(err.id2 == 1123);
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_branch(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    bs2->add_branch([&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), bs2, 9803);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Branch);
    assert(err.id2 == 9803);
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_extern_measure_endtime(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    assert(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    seq.get_slot(IR::Type::Float64, 0);

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 12189);
    bs2->add_endtime(bs2->track_time(t));

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::ExternMeasure);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::EventTime);
    assert(err.id2 == 12189);
    assert(strcmp(err.what(), "Measurement value can only be used "
                  "in the same basic sequence") == 0);
}

static void test_check_measure_order_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, seq.get_slot(IR::Type::Float64, 0), 998);
    auto p1 = bs1->add_pulse(1, 11, bs1->track_time(t), seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_var(m1->val())}, 1);
                             }());
    assert(p1->id() == 11);
    assert(!p1->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(p1->id() == 11);

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::MeasureOrder);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Pulse);
    assert(err.id2 == 11);
    assert(strcmp(err.what(), "Use of measure must be later than the measure.") == 0);
}

static void test_check_measure_order_id(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto p1 = bs1->add_pulse(1, 9, bs1->track_time(Seq::EventTime(100)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_var(m1->val())}, 1);
                             }());
    assert(p1->id() == 9);
    assert(!p1->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(p1->id() == 9);

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::MeasureOrder);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Pulse);
    assert(err.id2 == 9);
    assert(strcmp(err.what(), "Use of measure must be later than the measure.") == 0);
}

static void test_check_measure_order_id_ok(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto p1 = bs1->add_pulse(1, 11, bs1->track_time(Seq::EventTime(100)),
                             seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_var(m1->val())}, 1);
                             }());
    assert(p1->id() == 11);
    assert(!p1->is_measure());
    assert(!p1->needs_oldval());

    seq.prepare();
    seq.optimize();

    assert(bs1->has_output(1));

    assert(bs1->get_pulses(1)->size() == 1);
    assert(seq.get_basicseqs().size() == 1);

    assert(p1->id() == 11);
    assert(!p1->needs_oldval());
    assert(!bs1->needs_oldval(p1));
    assert(p1->endval());
    assert(approx(p1->endval()->get_const().get<double>(), -2.4));

    assert(bs1->startval(1)->get_const().is(IR::TagVal(0.0)));
    assert(approx(bs1->endval(1)->get_const().get<double>(), -2.4));
}

static void test_check_measure_order_measure_time(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    assert(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    assert(m1->id() == 10);
    assert(m1->is_measure());
    assert(m1->val()->is_extern());
    assert(m1->val()->get_extern().first == IR::Type::Float64);
    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t(100);
    t.add_term(Seq::EventTime::NonNeg, m1->val(), 998);
    auto m2 = bs1->add_measure(1, 9, bs1->track_time(t), bs1->new_measure(seq.env(), 2));
    assert(m2->id() == 9);
    assert(m2->is_measure());
    assert(m2->val()->is_extern());
    assert(m2->val()->get_extern().first == IR::Type::Float64);
    assert(m2->val()->get_extern().second == ((uint64_t(1) << 32) | 2));

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    assert(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    assert(m2->val()->get_extern().second == ((uint64_t(1) << 32) | 2));

    assert(err.code >> 16 == Seq::Error::BasicSeq);
    assert((err.code & ((1u << 16) - 1)) == Seq::BasicSeq::MeasureOrder);
    assert(err.type1 == Seq::Error::Measure);
    assert(err.id1 == m1->val()->get_extern().second);
    assert(err.type2 == Seq::Error::Measure);
    assert(err.id2 == m2->val()->get_extern().second);
    assert(strcmp(err.what(), "Use of measure must be later than the measure.") == 0);
}

int main()
{
    auto llvm_ctx = LLVM::new_context();

    test_check_extern_measure_value(*llvm_ctx);
    test_check_extern_measure_length(*llvm_ctx);
    test_check_extern_measure_time(*llvm_ctx);
    test_check_extern_measure_assign(*llvm_ctx);
    test_check_extern_measure_assume(*llvm_ctx);
    test_check_extern_measure_branch(*llvm_ctx);
    test_check_extern_measure_endtime(*llvm_ctx);

    test_check_measure_order_time(*llvm_ctx);
    test_check_measure_order_id(*llvm_ctx);
    test_check_measure_order_id_ok(*llvm_ctx);
    test_check_measure_order_measure_time(*llvm_ctx);

    return 0;
}
