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

#include "../../lib/nacs-seq/seq.h"
#include "../../lib/nacs-seq/error.h"

#include "../../lib/nacs-utils/llvm/codegen.h"
#include "../../lib/nacs-utils/llvm/utils.h"
#include "../../lib/nacs-utils/number.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

using namespace NaCs;
using namespace std::literals::string_literals;

static auto llvm_ctx = LLVM::new_context();

TEST_CASE("check_extern_measure_value") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

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
    REQUIRE(p2->id() == 7);
    REQUIRE(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p2->id() == 7);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == p2->id());
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_length") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

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
    REQUIRE(p2->id() == 5);
    REQUIRE(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p2->id() == 5);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasureLength));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == p2->id());
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_cond") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto cond = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(1)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }();
    auto p2 = bs2->add_pulse(1, 5, bs2->track_time(Seq::EventTime(0)), nullptr,
                             seq.get_const(IR::TagVal(1000)), cond);
    REQUIRE(p2->id() == 5);
    REQUIRE(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p2->id() == 5);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasureCond));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == p2->id());
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, [&] {
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
    REQUIRE(p2->id() == 58);
    REQUIRE(!p2->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p2->id() == 58);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::EventTime);
    REQUIRE(err.id2 == 1289);
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_assign") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    seq.get_slot(IR::Type::Float64, 0);

    bs2->assign_global(0, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 1009);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Assignment);
    REQUIRE(err.id2 == 1009);
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_assume") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    bs2->add_assume(Seq::Sign::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 1123);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::EventTime);
    REQUIRE(err.id2 == 1123);
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_branch") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    bs2->add_branch([&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), bs2, 9803);

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Branch);
    REQUIRE(err.id2 == 9803);
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_extern_measure_endtime") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    bs1->set_default_branch(bs2);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(0)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    seq.get_slot(IR::Type::Float64, 0);

    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(20000)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 12189);
    bs2->add_endtime(bs2->track_time(t));

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::ExternMeasure));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::EventTime);
    REQUIRE(err.id2 == 12189);
    REQUIRE(err.what() == "Measurement value can only be used in the same basic sequence"s);
}

TEST_CASE("check_measure_order_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 998);
    auto p1 = bs1->add_pulse(1, 11, bs1->track_time(t), seq.get_const(IR::TagVal(1000)),
                             [&] {
                                 IR::Builder builder(IR::Type::Float64,
                                                     {IR::Type::Float64, IR::Type::Float64});
                                 auto delta = builder.createMul(0, builder.getConst(-0.0024));
                                 builder.createRet(builder.createAdd(1, delta));
                                 return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                         Seq::Arg::create_var(m1->val())}, 1);
                             }());
    REQUIRE(p1->id() == 11);
    REQUIRE(!p1->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p1->id() == 11);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::MeasureOrder));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == 11);
    REQUIRE(err.what() == "Use of measure must be later than the measure."s);
}

TEST_CASE("check_measure_order_id") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

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
    REQUIRE(p1->id() == 9);
    REQUIRE(!p1->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p1->id() == 9);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::MeasureOrder));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == 9);
    REQUIRE(err.what() == "Use of measure must be later than the measure."s);
}

TEST_CASE("check_measure_order_id_only") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    auto p1 = bs1->add_pulse(1, 9, bs1->track_time(Seq::EventTime(1000)),
                             seq.get_const(IR::TagVal(1000)),
                             seq.get_const(IR::TagVal(0.3)),
                             m1->val());
    REQUIRE(p1->id() == 9);
    REQUIRE(!p1->is_measure());

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(p1->id() == 9);

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::MeasureOrder));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Pulse);
    REQUIRE(err.id2 == 9);
    REQUIRE(err.what() == "Use of measure must be later than the measure."s);
}

TEST_CASE("check_measure_order_id_ok") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

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
    REQUIRE(p1->id() == 11);
    REQUIRE(!p1->is_measure());
    REQUIRE(!p1->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->has_output(1));

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 1);

    REQUIRE(p1->id() == 11);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(-2.4));

    REQUIRE(bs1->startval(1)->get_const().is(IR::TagVal(0.0)));
    REQUIRE(bs1->endval(1)->get_const().get<double>() == Approx(-2.4));
}

TEST_CASE("check_measure_order_measure_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto m1 = bs1->add_measure(1, 10, bs1->track_time(Seq::EventTime(100)),
                               bs1->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 10);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val()->is_extern());
    REQUIRE(m1->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));

    Seq::EventTime t(100);
    t.add_term(Seq::Sign::NonNeg, m1->val(), 998);
    auto m2 = bs1->add_measure(1, 9, bs1->track_time(t), bs1->new_measure(seq.env(), 2));
    REQUIRE(m2->id() == 9);
    REQUIRE(m2->is_measure());
    REQUIRE(m2->val()->is_extern());
    REQUIRE(m2->val()->get_extern().first == IR::Type::Float64);
    REQUIRE(m2->val()->get_extern().second == ((uint64_t(1) << 32) | 2));

    auto err = expect_error<Seq::Error>([&] {
        seq.prepare();
    });

    REQUIRE(m1->val()->get_extern().second == ((uint64_t(1) << 32) | 1));
    REQUIRE(m2->val()->get_extern().second == ((uint64_t(1) << 32) | 2));

    REQUIRE(err.type == Seq::Error::Type::BasicSeq);
    REQUIRE(err.code == uint16_t(Seq::Error::BasicSeq::MeasureOrder));
    REQUIRE(err.type1 == Seq::Error::Type::Measure);
    REQUIRE(err.id1 == m1->val()->get_extern().second);
    REQUIRE(err.type2 == Seq::Error::Type::Measure);
    REQUIRE(err.id2 == m2->val()->get_extern().second);
    REQUIRE(err.what() == "Use of measure must be later than the measure."s);
}
