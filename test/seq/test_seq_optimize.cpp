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

TEST_CASE("basic_single_channel") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1.45)));
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(0.3)));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    REQUIRE(p2->id() == 2);
    REQUIRE(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(p1->start() == p2->start());
    REQUIRE(p1->start().tconst == 6);
    REQUIRE(p1->start().terms.size() == 2);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    // The sequence cache the slot variable so we can compare pointer identity.
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->start().terms[1].sign == Seq::Sign::NonNeg);
    REQUIRE(p1->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(0.3));
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(-0.2));
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!p2->needs_oldval());
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(0.0)));
    REQUIRE(bs->endval(1)->get_const().get<double>() == Approx(-0.2));
    REQUIRE(!bs->get_default_branch());
    REQUIRE(bs->get_branches().empty());
    auto slots = seq.get_slots();
    REQUIRE(slots.size() == 2);
    REQUIRE((slots[0]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(0))));
    REQUIRE(slots[0].get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE((slots[1]->get_extern() == std::make_pair(IR::Type::Float64, uint64_t(1))));
    REQUIRE(slots[1].get() == seq.get_slot(IR::Type::Float64, 1));
}

TEST_CASE("basic_single_channel_measure_known_time1") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));

    seq.set_defval(1, 2.3);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(0)), nullptr,
                            [&] {
                                IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
                                builder.createRet(builder.createSub(0, builder.getConst(2)));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(1)}, 2);
                            }());
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0005));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    REQUIRE(p2->id() == 2);
    REQUIRE(p2->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(30)));
    auto m1 = bs->add_measure(1, 3, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 3);
    REQUIRE(m1->is_measure());
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 1));
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
    REQUIRE(p3->id() == 4);
    REQUIRE(p3->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 4);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 3);
    REQUIRE(p1->start() == Seq::EventTime(0));
    REQUIRE(p2->start().tconst == 5);
    REQUIRE(p2->start().terms.size() == 1);
    REQUIRE(p2->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p2->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p3->start().tconst == 35);
    REQUIRE(p3->start().terms.size() == 2);
    REQUIRE(p3->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p3->start().terms[1].sign == Seq::Sign::Pos);
    if (p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0)) {
        REQUIRE(p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
        REQUIRE(p3->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    }
    else {
        REQUIRE(p3->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 1));
        REQUIRE(p3->start().terms[1].var.get() == seq.get_slot(IR::Type::Float64, 0));
    }
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(0.3));
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(-0.2));
    REQUIRE(p3->endval()->get_const().get<double>() == Approx(0.285));
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!p3->needs_oldval());
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(p3->is_ordered);
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(!bs->needs_oldval(p3));
    REQUIRE(!p1->len());
    REQUIRE(p2->len());
    REQUIRE(p3->len());
    REQUIRE(p1->val()->is_const());
    REQUIRE(p2->val()->is_call());
    REQUIRE(!p2->val()->argument_unused(0));
    REQUIRE(!p3->val()->argument_unused(0));
    // The start value is exact since we set it directly.
    // No approx comparison needed.
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    REQUIRE(bs->endval(1)->get_const().get<double>() == Approx(0.285));
}

TEST_CASE("basic_single_channel_measure_known_time2") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1000)));
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 2));
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
    REQUIRE(p2->id() == 3);
    REQUIRE(p2->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 2);

    REQUIRE(p1->start().tconst == 5);
    REQUIRE(p1->start().terms.size() == 1);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p2->start().tconst == 1005);
    REQUIRE(p2->start().terms.size() == 3);
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!p2->needs_oldval());
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    REQUIRE(bs->endval(1)->get_const().get<double>() == Approx(1.8));
}

TEST_CASE("basic_single_channel_measure_unknown_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 2));
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
    REQUIRE(p2->id() == 3);
    REQUIRE(p2->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 2);

    REQUIRE(p1->start().tconst == 5);
    REQUIRE(p1->start().terms.size() == 1);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p2->start().tconst == 5);
    REQUIRE(p2->start().terms.size() == 3);
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(p2->endval()->is_call());
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!p2->needs_oldval());
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    REQUIRE(bs->endval(1) == p2->endval());

    auto endval = p2->endval();
    REQUIRE(endval->args().size() == 1);
    REQUIRE(endval->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(endval->get_callee().is_llvm);

    auto var1 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(500))})->ref();
    auto var2 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(1000))})->ref();
    auto var3 = seq.env().new_call(endval->get_callee().llvm,
                                   {Seq::Arg::create_const(IR::TagVal(1500))})->ref();

    seq.env().optimize();
    REQUIRE(var1->is_const());
    REQUIRE(var1->get_const().get<double>() == Approx(2.05));
    REQUIRE(var2->is_const());
    REQUIRE(var2->get_const().get<double>() == Approx(1.8));
    REQUIRE(var3->is_const());
    REQUIRE(var3->get_const().get<double>() == Approx(1.8));
}

TEST_CASE("basic_single_channel_unused_measure") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Unknown, seq.get_slot(IR::Type::Float64, 1));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    REQUIRE(bs->get_pulses(1)->size() == 2);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 1);

    REQUIRE(p1->start().tconst == 5);
    REQUIRE(p1->start().terms.size() == 1);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->is_ordered);
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    REQUIRE(bs->endval(1)->get_const().get<double>() == Approx(1.8));
}

TEST_CASE("basic_two_channel_measure_time") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    REQUIRE(seq.get_chn_id("test_chn2", true) == 2);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)));

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1000)));
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(10)));
    t.add_term(Seq::Sign::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(1234)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }());
    auto p2 = bs->add_pulse(2, 3, bs->track_time(t), nullptr, m1->val());
    REQUIRE(p2->id() == 3);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 2);
    REQUIRE(bs->get_pulses(2)->size() == 1);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 1);
    REQUIRE(bs->get_pulses(2)->size() == 1);

    REQUIRE(p1->start().tconst == 5);
    REQUIRE(p1->start().terms.size() == 1);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p2->start().tconst == 3236); // 1000 + 10 + 5 + 1.8 * 1234
    REQUIRE(p2->start().terms.size() == 1);
    REQUIRE(p2->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p2->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(p1->val()->is_call());
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(p2->val()->is_const());
    REQUIRE(p2->val()->get_const().get<double>() == Approx(1.8));
    REQUIRE(p1->endval()->is_const());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(1.8));
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(bs->startval(1)->get_const().is(IR::TagVal(2.3)));
    REQUIRE(bs->startval(2)->get_const().is(IR::TagVal(0.0)));
    REQUIRE(bs->endval(1)->get_const().get<double>() == Approx(1.8));
    REQUIRE(bs->endval(2)->get_const().get<double>() == Approx(1.8));
}

TEST_CASE("basic_single_channel_measure_unknown_order") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    seq.set_defval(1, 2.3);

    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    Seq::EventTime t2;
    t2.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1200)));

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t2), seq.get_const(IR::TagVal(1000)),
                            [&] {
                                IR::Builder builder(IR::Type::Float64,
                                                    {IR::Type::Float64, IR::Type::Float64});
                                auto delta = builder.createMul(0, builder.getConst(-0.0006));
                                builder.createRet(builder.createAdd(1, delta));
                                return seq.get_call(builder.get(), {Seq::Arg::create_arg(0),
                                        Seq::Arg::create_arg(1)}, 2);
                            }());
    REQUIRE(p2->id() == 2);
    REQUIRE(p2->needs_oldval());

    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1200)));
    auto m1 = bs->add_measure(1, 3, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 3);
    REQUIRE(m1->is_measure());

    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 1));
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
    REQUIRE(p3->id() == 4);
    REQUIRE(p3->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 4);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 4);
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(bs->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p2->id() == 2);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(bs->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(m1->id() == 3);
    REQUIRE(m1->is_measure());
    REQUIRE(p3->id() == 4);
    REQUIRE(!p3->needs_oldval());
    REQUIRE(bs->needs_oldval(p3));
    REQUIRE(p3->endval());
    REQUIRE(!p1->is_ordered);
    REQUIRE(!p2->is_ordered);
    REQUIRE(p3->is_ordered);
    REQUIRE(bs->endval(1));
}

TEST_CASE("basic_single_channel_measure_time_violation") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    t.add_term(Seq::Sign::Pos, [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createMul(0, builder.getConst(1234)));
        return seq.get_call(builder.get(), {Seq::Arg::create_var(m1->val())}, 0);
    }(), 2);
    auto p2 = bs->add_pulse(1, 3, bs->track_time(t), nullptr, m1->val());
    REQUIRE(p2->id() == 3);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(bs->get_pulses(1)->size() == 3);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NonPosTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 2);
    REQUIRE(err.what() == "Positive time expected."s);
}

TEST_CASE("basic_single_channel_output_unknown_order") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    seq.set_defval(1, 2.3);

    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0));
    Seq::EventTime t2;
    t2.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1200)));

    // Test that we do not compute end value of basic sequence
    // when we cannot determine which one is the last output pulse.
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr,
                            seq.get_const(IR::TagVal(1.0)));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto p2 = bs->add_pulse(1, 2, bs->track_time(t2), nullptr,
                            seq.get_const(IR::TagVal(2.0)));
    REQUIRE(p2->id() == 2);
    REQUIRE(!p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->size() == 2);
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs->needs_oldval(p1));
    REQUIRE(p1->endval()->is_const());
    REQUIRE(p2->id() == 2);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs->needs_oldval(p2));
    REQUIRE(p2->endval()->is_const());
    REQUIRE(!p1->is_ordered);
    REQUIRE(!p2->is_ordered);
    REQUIRE(!bs->endval(1));
}

TEST_CASE("basic_single_channel_neg_time_offset") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t(-10);
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 29);
    t.add_term(Seq::Sign::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(-1000)));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1));
    REQUIRE(bs->get_pulses(1)->size() == 1);
    REQUIRE(&*bs->get_pulses(1)->begin() == p1);
    REQUIRE(p1->start().tconst == -1010);
    REQUIRE(p1->start().terms.size() == 1);
    REQUIRE(p1->start().terms[0].sign == Seq::Sign::Pos);
    REQUIRE(p1->start().terms[0].id == 29);
    REQUIRE(p1->start().terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
}

TEST_CASE("basic_single_channel_neg_time_output_violation") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t(10);
    t.add_term(Seq::Sign::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto p1 = bs->add_pulse(1, 14, bs->track_time(t), nullptr, seq.get_const(IR::TagVal(-1000)));
    REQUIRE(p1->id() == 14);
    REQUIRE(!p1->needs_oldval());

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 14);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("basic_single_channel_neg_time_measure_violation") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    Seq::EventTime t(10);
    t.add_term(Seq::Sign::Unknown, seq.get_const(IR::TagVal(-1000)), 33);
    auto m1 = bs->add_measure(1, 19, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 19);
    REQUIRE(m1->is_measure());
    // Make sure the value is used.
    // Otherwise we don't guarantee the throwing of the error.
    seq.get_slot(IR::Type::Float64, 0);
    bs->assign_global(0, m1->val(), 20);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 19);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("single_channel_branch") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());

    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
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
    REQUIRE(p2->id() == 1);
    REQUIRE(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(bs2->get_default_branch() == nullptr);

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs2->get_pulses(1)->size() == 1);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(1.5));
    REQUIRE(p2->id() == 1);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs2->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(1.9));
    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);

    REQUIRE(bs1->startval(1)->get_const().get<double>() == Approx(2.3));
    REQUIRE(bs1->endval(1)->get_const().get<double>() == Approx(1.5));
    REQUIRE(bs2->startval(1)->get_const().get<double>() == Approx(1.5));
    REQUIRE(bs2->endval(1)->get_const().get<double>() == Approx(1.9));
}

TEST_CASE("single_channel_branch_merge") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);
    auto bs4 = seq.add_basicseq(4);
    REQUIRE(bs4->id() == 4);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
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
    REQUIRE(p2->id() == 1);
    REQUIRE(p2->needs_oldval());
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
    REQUIRE(p3->id() == 1);
    REQUIRE(p3->needs_oldval());
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
    REQUIRE(p4->id() == 1);
    REQUIRE(p4->needs_oldval());
    // no-op but should work
    bs4->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs2->get_pulses(1)->size() == 1);
    REQUIRE(bs3->get_pulses(1)->size() == 1);
    REQUIRE(bs4->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 4);

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(bs2->get_default_branch() == bs4);
    REQUIRE(bs3->get_default_branch() == bs4);
    REQUIRE(!bs4->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs3);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs2->get_branches().empty());
    REQUIRE(bs3->get_branches().empty());
    REQUIRE(bs4->get_branches().empty());

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(-0.7));

    REQUIRE(p2->id() == 1);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs2->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(0.1));

    REQUIRE(p3->id() == 1);
    REQUIRE(!p3->needs_oldval());
    REQUIRE(!bs3->needs_oldval(p3));
    REQUIRE(p3->endval());
    REQUIRE(p3->endval()->get_const().get<double>() == Approx(0.1));

    REQUIRE(p2->endval()->get_const().get<double>() == p3->endval()->get_const().get<double>());

    REQUIRE(p4->id() == 1);
    REQUIRE(!p4->needs_oldval());
    REQUIRE(!bs4->needs_oldval(p4));
    REQUIRE(p4->endval());
    REQUIRE(p4->endval()->get_const().get<double>() == Approx(1.1));

    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(p3->is_ordered);
    REQUIRE(p4->is_ordered);

    REQUIRE(bs1->startval(1)->get_const().get<double>() == Approx(2.3));
    REQUIRE(bs1->endval(1)->get_const().get<double>() == Approx(-0.7));

    REQUIRE(bs2->startval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs2->endval(1)->get_const().get<double>() == Approx(0.1));
    REQUIRE(bs3->startval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs3->endval(1)->get_const().get<double>() == Approx(0.1));

    REQUIRE(bs4->startval(1)->get_const().get<double>() == Approx(0.1));
    REQUIRE(bs4->endval(1)->get_const().get<double>() == Approx(1.1));
}

TEST_CASE("single_channel_branch_not_merge") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);
    auto bs4 = seq.add_basicseq(4);
    REQUIRE(bs4->id() == 4);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
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
    REQUIRE(p2->id() == 1);
    REQUIRE(p2->needs_oldval());
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
    REQUIRE(p3->id() == 1);
    REQUIRE(p3->needs_oldval());
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
    REQUIRE(p4->id() == 1);
    REQUIRE(p4->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs2->get_pulses(1)->size() == 1);
    REQUIRE(bs3->get_pulses(1)->size() == 1);
    REQUIRE(bs4->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 4);

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(bs2->get_default_branch() == bs4);
    REQUIRE(bs3->get_default_branch() == bs4);
    REQUIRE(!bs4->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs3);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs2->get_branches().empty());
    REQUIRE(bs3->get_branches().empty());
    REQUIRE(bs4->get_branches().empty());

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs2->get_pulses(1)->size() == 1);
    REQUIRE(bs3->get_pulses(1)->size() == 1);
    REQUIRE(bs4->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 4);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(-0.7));

    REQUIRE(p2->id() == 1);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs2->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(0.1));

    REQUIRE(p3->id() == 1);
    REQUIRE(!p3->needs_oldval());
    REQUIRE(!bs3->needs_oldval(p3));
    REQUIRE(p3->endval());
    REQUIRE(p3->endval()->get_const().get<double>() == Approx(-1.5));

    REQUIRE(p4->id() == 1);
    REQUIRE(!p4->needs_oldval());
    REQUIRE(!bs4->needs_oldval(p4));
    REQUIRE(p4->endval());

    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);
    REQUIRE(p3->is_ordered);
    REQUIRE(p4->is_ordered);

    REQUIRE(bs1->startval(1)->get_const().get<double>() == Approx(2.3));
    REQUIRE(bs1->endval(1)->get_const().get<double>() == Approx(-0.7));

    REQUIRE(bs2->startval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs2->endval(1)->get_const().get<double>() == Approx(0.1));
    REQUIRE(bs3->startval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs3->endval(1)->get_const().get<double>() == Approx(-1.5));

    REQUIRE(bs4->startval(1)->is_extern());
    REQUIRE(bs4->endval(1));
}

TEST_CASE("single_channel_branch_nonconst") {
    // Non-const cannot be forwarded (in general)
    // since it most likely either depend on a measure,
    // which has the wrong scope,
    // or a global variable, which may be reassigned between sequence.
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(10)),
                             seq.get_const(IR::TagVal(1000)),
                             seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
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
    REQUIRE(p2->id() == 1);
    REQUIRE(p2->needs_oldval());

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs2->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 2);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval() == seq.get_slot(IR::Type::Float64, 0));

    REQUIRE(p2->id() == 1);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs2->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(!p2->len());

    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);

    REQUIRE(bs1->startval(1)->get_const().is(IR::TagVal(0.0)));
    REQUIRE(bs1->endval(1));
    REQUIRE(bs1->endval(1) == seq.get_slot(IR::Type::Float64, 0));

    REQUIRE(bs2->startval(1)->is_extern());
    REQUIRE(bs2->endval(1));
}

TEST_CASE("branch_cond_clone_point") {
    // Branch condition that
    // * uses only global value: no cloning
    // * uses only measure: clone as a whole
    // * uses both global and measure: clone from interface
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 2, bs1->track_time(Seq::EventTime(1)), m1val);
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val() == m1val);

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

    REQUIRE(bs1->get_branches().size() == 3);
    REQUIRE(bs1->get_branches()[0].target == bs2);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == gv2);
    REQUIRE(bs1->get_branches()[1].target == nullptr);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() == lv2);
    REQUIRE(bs1->get_branches()[2].target == bs1);
    REQUIRE(bs1->get_branches()[2].id == 3);
    REQUIRE(bs1->get_branches()[2].cond.get() == glv2);

    seq.prepare();

    REQUIRE(seq.get_slots().size() == 2);
    REQUIRE(bs1->get_branches().size() == 3);
    REQUIRE(bs1->get_branches()[0].target == bs2);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == gv2);
    REQUIRE(bs1->get_branches()[1].target == nullptr);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() != lv2);
    auto new_lv2 = bs1->get_branches()[1].cond.get();
    REQUIRE(new_lv2->is_extern());
    REQUIRE(new_lv2->get_extern().first == IR::Type::Bool);
    REQUIRE(new_lv2->get_extern().second == 2);
    REQUIRE(bs1->get_branches()[2].target == bs1);
    REQUIRE(bs1->get_branches()[2].id == 3);
    REQUIRE(bs1->get_branches()[2].cond.get() != glv2);
    // glv1 is merged into glv2
    auto new_glv2 = bs1->get_branches()[2].cond.get();
    REQUIRE(new_glv2->is_call());
    REQUIRE(new_glv2->get_callee().is_llvm);
    REQUIRE(new_glv2->get_callee().llvm == glv2->get_callee().llvm);
    REQUIRE(new_glv2->args().size() == 2);
    REQUIRE(new_glv2->args()[0].is_var());
    REQUIRE(new_glv2->args()[0].get_var() == new_lv2);
    REQUIRE(new_glv2->args()[1].is_var());
    REQUIRE(new_glv2->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs1->get_assigns().size() == 1);
    REQUIRE(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(2)->second.val.get() == lv2);

    seq.optimize();

    REQUIRE(seq.get_slots().size() == 3);
    REQUIRE(bs1->get_branches().size() == 3);
    REQUIRE(bs1->get_branches()[0].target == bs2);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == gv2);
    REQUIRE(bs1->get_branches()[1].target == nullptr);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 2));
    new_lv2 = seq.get_slot(IR::Type::Float64, 2);
    REQUIRE(bs1->get_branches()[2].target == bs1);
    REQUIRE(bs1->get_branches()[2].id == 3);
    REQUIRE(bs1->get_branches()[2].cond.get() == new_glv2);
    REQUIRE(new_glv2->is_call());
    REQUIRE(new_glv2->get_callee().is_llvm);
    REQUIRE(new_glv2->args().size() == 2);
    REQUIRE(new_glv2->args()[0].is_var());
    REQUIRE(new_glv2->args()[1].is_var());
    auto arg0 = new_glv2->args()[0].get_var();
    auto arg1 = new_glv2->args()[1].get_var();
    REQUIRE(((arg0 == seq.get_slot(IR::Type::Float64, 1) &&
              arg1 == seq.get_slot(IR::Type::Float64, 2)) ||
             (arg0 == seq.get_slot(IR::Type::Float64, 2) &&
              arg1 == seq.get_slot(IR::Type::Float64, 1))));

    REQUIRE(bs1->get_assigns().size() == 1);
    REQUIRE(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(2)->second.val.get() == lv2);
}

TEST_CASE("branch_cond_clone_const_opt") {
    // * If the cloned value used in a condition is a constant it can be forwarded.
    // * The forwarded condition should not affect the one in a different basic sequence
    //   with the same temporary ID.
    // * The cloned slot after it should be compacted down.
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);

    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 1, bs1->track_time(Seq::EventTime(0)), m1val);
    REQUIRE(m1->id() == 1);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val() == m1val); // optimizable to constant
    auto p1 = bs1->add_pulse(1, 2, bs1->track_time(Seq::EventTime(1)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->id() == 2);
    REQUIRE(!p1->needs_oldval());
    auto m2val = bs1->new_measure(seq.env(), 1);
    auto m2 = bs1->add_measure(1, 3, bs1->track_time(Seq::EventTime(2)), m2val);
    REQUIRE(m2->id() == 3);
    REQUIRE(m2->is_measure());
    REQUIRE(m2->val() == m2val); // not optimizable to constant (slot 0)

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
    REQUIRE(bs1->get_branches().size() == 2);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == glv1);
    REQUIRE(bs1->get_branches()[1].target == bs2);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() == glv2);

    auto p2 = bs2->add_pulse(1, 4, bs2->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p2->id() == 4);
    REQUIRE(!p2->needs_oldval());
    auto m3val = bs2->new_measure(seq.env(), 1);
    auto m3 = bs2->add_measure(1, 5, bs2->track_time(Seq::EventTime(1)), m3val);
    REQUIRE(m3->id() == 5);
    REQUIRE(m3->is_measure());
    REQUIRE(m3->val() == m3val);

    auto glv3 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m3val)}, 0);
    }();
    bs2->add_branch(glv3, bs2, 3);
    REQUIRE(bs2->get_branches().size() == 1);
    REQUIRE(bs2->get_branches()[0].target == bs2);
    REQUIRE(bs2->get_branches()[0].id == 3);
    REQUIRE(bs2->get_branches()[0].cond.get() == glv3);

    seq.prepare();

    REQUIRE(seq.get_slots().size() == 2);

    REQUIRE(bs1->get_branches().size() == 2);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() != glv1);
    auto new_glv1 = bs1->get_branches()[0].cond.get();
    REQUIRE(new_glv1->is_call());
    REQUIRE(new_glv1->get_callee().is_llvm);
    REQUIRE(new_glv1->get_callee().llvm == glv1->get_callee().llvm);
    REQUIRE(new_glv1->args().size() == 2);
    REQUIRE(new_glv1->args()[0].is_var());
    REQUIRE(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv1->args()[1].is_var());
    auto new_m1var = new_glv1->args()[1].get_var();
    REQUIRE(new_m1var->is_extern());
    REQUIRE(new_m1var->get_extern().first == IR::Type::Float64);
    REQUIRE(new_m1var->get_extern().second == 2);
    REQUIRE(bs1->get_branches()[1].target == bs2);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() != glv2);
    auto new_glv2 = bs1->get_branches()[1].cond.get();
    REQUIRE(new_glv2->is_call());
    REQUIRE(new_glv2->get_callee().is_llvm);
    REQUIRE(new_glv2->get_callee().llvm == glv2->get_callee().llvm);
    REQUIRE(new_glv2->args().size() == 2);
    REQUIRE(new_glv2->args()[0].is_var());
    REQUIRE(new_glv2->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv2->args()[1].is_var());
    auto new_m2var = new_glv2->args()[1].get_var();
    REQUIRE(new_m2var->is_extern());
    REQUIRE(new_m2var->get_extern().first == IR::Type::Float64);
    REQUIRE(new_m2var->get_extern().second == 3);

    REQUIRE(bs1->get_assigns().size() == 2);
    REQUIRE(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(2)->second.val.get() == m1val);
    REQUIRE(bs1->get_assigns().find(3) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(3)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(3)->second.val.get() == m2val);

    REQUIRE(bs2->get_branches().size() == 1);
    REQUIRE(bs2->get_branches()[0].target == bs2);
    REQUIRE(bs2->get_branches()[0].id == 3);
    REQUIRE(bs2->get_branches()[0].cond.get() != glv3);
    auto new_glv3 = bs2->get_branches()[0].cond.get();
    REQUIRE(new_glv3->is_call());
    REQUIRE(new_glv3->get_callee().is_llvm);
    REQUIRE(new_glv3->get_callee().llvm == glv3->get_callee().llvm);
    REQUIRE(new_glv3->args().size() == 2);
    REQUIRE(new_glv3->args()[0].is_var());
    REQUIRE(new_glv3->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv3->args()[1].is_var());
    auto new_m3var = new_glv3->args()[1].get_var();
    REQUIRE(new_m3var->is_extern());
    REQUIRE(new_m3var->get_extern().first == IR::Type::Float64);
    REQUIRE(new_m3var->get_extern().second == 2);

    REQUIRE(bs2->get_assigns().size() == 1);
    REQUIRE(bs2->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs2->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs2->get_assigns().find(2)->second.val.get() == m3val);

    seq.optimize();

    REQUIRE(seq.get_slots().size() == 3);

    REQUIRE(bs1->get_branches().size() == 2);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == new_glv1);
    REQUIRE(new_glv1->is_call());
    REQUIRE(new_glv1->args().size() == 1);
    REQUIRE(new_glv1->get_callee().is_llvm);
    REQUIRE(new_glv1->args()[0].is_var());
    REQUIRE(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs1->get_branches()[1].target == bs2);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[1].cond.get() == new_glv2);
    REQUIRE(new_glv2->is_call());
    REQUIRE(new_glv2->get_callee().is_llvm);
    REQUIRE(new_glv2->args()[0].is_var());
    REQUIRE(new_glv2->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv2->args()[1].is_var());
    REQUIRE(new_glv2->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 2));

    REQUIRE(bs1->get_assigns().size() == 1);
    REQUIRE(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(2)->second.val.get() == seq.get_slot(IR::Type::Float64, 0));

    REQUIRE(bs2->get_branches().size() == 1);
    REQUIRE(bs2->get_branches()[0].target == bs2);
    REQUIRE(bs2->get_branches()[0].id == 3);
    REQUIRE(bs2->get_branches()[0].cond.get() == new_glv3);
    REQUIRE(new_glv3->is_call());
    REQUIRE(new_glv3->get_callee().is_llvm);
    REQUIRE(new_glv3->args()[0].is_var());
    REQUIRE(new_glv3->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv3->args()[1].is_var());
    REQUIRE(new_glv3->args()[1].get_var() == seq.get_slot(IR::Type::Float64, 2));

    REQUIRE(bs2->get_assigns().size() == 1);
    REQUIRE(bs2->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs2->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs2->get_assigns().find(2)->second.val.get() == seq.get_slot(IR::Type::Float64, 0));
}

TEST_CASE("branch_cond_clone_unused_opt") {
    // Unused slot will be optimized out.
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

    auto p1 = bs1->add_pulse(1, 1, bs1->track_time(Seq::EventTime(0)),
                             nullptr, seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    auto m1val = bs1->new_measure(seq.env(), 1);
    auto m1 = bs1->add_measure(1, 2, bs1->track_time(Seq::EventTime(1)), m1val);
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    REQUIRE(m1->val() == m1val);

    auto glv1 = [&] {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        return seq.get_call(builder.get(), {
                Seq::Arg::create_var(seq.get_slot(IR::Type::Float64, 1)),
                Seq::Arg::create_var(m1val)}, 0);
    }();

    bs1->add_branch(glv1, nullptr, 1);
    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == glv1);

    seq.prepare();

    REQUIRE(seq.get_slots().size() == 2);

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() != glv1);
    auto new_glv1 = bs1->get_branches()[0].cond.get();
    REQUIRE(new_glv1->is_call());
    REQUIRE(new_glv1->get_callee().is_llvm);
    REQUIRE(new_glv1->get_callee().llvm == glv1->get_callee().llvm);
    REQUIRE(new_glv1->args().size() == 2);
    REQUIRE(new_glv1->args()[0].is_var());
    REQUIRE(new_glv1->args()[0].get_var() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(new_glv1->args()[1].is_var());
    auto new_m1var = new_glv1->args()[1].get_var();
    REQUIRE(new_m1var->is_extern());
    REQUIRE(new_m1var->get_extern().first == IR::Type::Float64);
    REQUIRE(new_m1var->get_extern().second == 2);

    REQUIRE(bs1->get_assigns().size() == 1);
    REQUIRE(bs1->get_assigns().find(2) != bs1->get_assigns().end());
    REQUIRE(bs1->get_assigns().find(2)->second.id == 0);
    REQUIRE(bs1->get_assigns().find(2)->second.val.get() == m1val);

    seq.optimize();

    REQUIRE(seq.get_slots().size() == 2);

    REQUIRE(bs1->get_branches().empty());
    REQUIRE(bs1->get_assigns().empty());
}

TEST_CASE("branch_loop_const") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    REQUIRE(seq.get_chn_id("test_chn2", true) == 2);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    bs1->set_default_branch(bs2);

    auto m1 = bs2->add_measure(1, 1, bs2->track_time(Seq::EventTime(100)),
                               bs2->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 1);
    REQUIRE(m1->is_measure());
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
    REQUIRE(p2->id() == 2);
    REQUIRE(p2->needs_oldval());
    bs2->set_default_branch(bs2);
    bs2->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    REQUIRE(bs1->has_output(1));
    REQUIRE(!bs1->has_output(2));
    REQUIRE(!bs2->has_output(1));
    REQUIRE(bs2->has_output(2));

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->has_output(1));
    REQUIRE(!bs1->has_output(2));
    REQUIRE(!bs2->has_output(1));
    REQUIRE(bs2->has_output(2));

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(bs1->get_pulses(2)->empty());
    REQUIRE(bs2->get_pulses(1)->empty());
    REQUIRE(bs2->get_pulses(2)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 2);

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(bs2->get_default_branch() == bs2);

    REQUIRE(bs1->get_branches().empty());
    REQUIRE(bs2->get_branches().size() == 1);
    REQUIRE(bs2->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs2->get_branches()[0].target == nullptr);
    REQUIRE(bs2->get_branches()[0].id == 10);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->get_const().get<double>() == Approx(-0.7));

    REQUIRE(p2->id() == 2);
    REQUIRE(!p2->needs_oldval());
    REQUIRE(!bs2->needs_oldval(p2));
    REQUIRE(p2->endval());
    REQUIRE(p2->endval()->get_const().get<double>() == Approx(0.2));

    REQUIRE(p1->is_ordered);
    REQUIRE(p2->is_ordered);

    REQUIRE(bs1->startval(1)->get_const().get<double>() == Approx(2.3));
    REQUIRE(bs1->endval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs1->startval(2)->get_const().is(IR::TagVal(0.0)));
    REQUIRE(bs1->endval(2)->get_const().is(IR::TagVal(0.0)));

    REQUIRE(bs2->startval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs2->endval(1)->get_const().get<double>() == Approx(-0.7));
    REQUIRE(bs2->startval(2)->is_extern());
    REQUIRE(bs2->endval(2)->get_const().get<double>() == Approx(0.2));
}

TEST_CASE("single_channel_loop_first") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 1);

    REQUIRE(bs1->get_default_branch() == bs1);
    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 10);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(!p1->len());
    REQUIRE(p1->is_ordered);

    REQUIRE(bs1->startval(1)->is_extern());
    REQUIRE(bs1->endval(1));
}

#if 0
// Currently the optimizer does not look into each pulse to do constant propagation
// so this optimization won't happen.
TEST_CASE("single_channel_loop_first2") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 1);

    REQUIRE(bs1->get_default_branch() == bs1);
    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 10);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->is_const());
    REQUIRE(p1->endval()->get_const().is(IR::TagVal(1.0)));
    REQUIRE(!p1->len());
    REQUIRE(p1->is_ordered);

    REQUIRE(bs1->startval(1)->is_const());
    REQUIRE(bs1->startval(1)->get_const().is(IR::TagVal(1.0)));
    REQUIRE(bs1->endval(1)->is_const());
    REQUIRE(bs1->endval(1)->get_const().is(IR::TagVal(1.0)));
}
#endif

TEST_CASE("single_channel_loop_first3") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    bs1->set_default_branch(bs1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), nullptr, 10);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs1->get_pulses(1)->size() == 1);
    REQUIRE(seq.get_basicseqs().size() == 1);

    REQUIRE(bs1->get_default_branch() == bs1);
    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == nullptr);
    REQUIRE(bs1->get_branches()[0].id == 10);

    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(!bs1->needs_oldval(p1));
    REQUIRE(p1->endval());
    REQUIRE(p1->endval()->is_const());
    REQUIRE(p1->endval()->get_const().is(IR::TagVal(1.0)));
    REQUIRE(p1->len());
    REQUIRE(p1->is_ordered);

    REQUIRE(bs1->startval(1)->is_const());
    REQUIRE(bs1->startval(1)->get_const().is(IR::TagVal(1.0)));
    REQUIRE(bs1->endval(1)->is_const());
    REQUIRE(bs1->endval(1)->get_const().is(IR::TagVal(1.0)));
}

TEST_CASE("optimize_cfg_shortcut") {
    Seq::Seq seq(*llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs2, 2);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs3, 3);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    REQUIRE(seq.get_basicseqs().size() == 2);
    REQUIRE(bs1->id() == 1);
    REQUIRE(bs2->id() == 2);

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(!bs2->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs1);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs2->get_branches().empty());
}

TEST_CASE("optimize_cfg_delete") {
    Seq::Seq seq(*llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_const(IR::TagVal(0)), bs2, 2);
    bs1->add_branch(seq.get_const(IR::TagVal(0.4)), bs3, 3);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    REQUIRE(seq.get_basicseqs().size() == 2);
    REQUIRE(bs1->id() == 1);
    REQUIRE(bs3->id() == 3);

    REQUIRE(bs1->get_default_branch() == bs3);
    REQUIRE(!bs3->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 1);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs1);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs3->get_branches().empty());
}

TEST_CASE("optimize_cfg_merge") {
    Seq::Seq seq(*llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 1), bs3, 2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 2), bs2, 3);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 3), bs2, 4);
    bs1->set_default_branch(bs2);

    seq.prepare();
    seq.optimize();

    REQUIRE(seq.get_basicseqs().size() == 3);
    REQUIRE(bs1->id() == 1);
    REQUIRE(bs2->id() == 2);
    REQUIRE(bs3->id() == 3);

    REQUIRE(bs1->get_default_branch() == bs2);
    REQUIRE(!bs2->get_default_branch());
    REQUIRE(!bs3->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 2);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs1);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs1->get_branches()[1].target == bs3);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs2->get_branches().empty());
    REQUIRE(bs3->get_branches().empty());
}

TEST_CASE("optimize_cfg_merge_null") {
    Seq::Seq seq(*llvm_ctx);
    auto bs1 = seq.add_basicseq(1);
    REQUIRE(bs1->id() == 1);
    auto bs2 = seq.add_basicseq(2);
    REQUIRE(bs2->id() == 2);
    auto bs3 = seq.add_basicseq(3);
    REQUIRE(bs3->id() == 3);

    bs1->add_branch(seq.get_slot(IR::Type::Float64, 0), bs1, 1);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 1), bs3, 2);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 2), bs2, 3);
    bs1->add_branch(seq.get_slot(IR::Type::Float64, 3), nullptr, 4);
    bs1->set_default_branch(nullptr);

    seq.prepare();
    seq.optimize();

    REQUIRE(seq.get_basicseqs().size() == 3);
    REQUIRE(bs1->id() == 1);
    REQUIRE(bs2->id() == 2);
    REQUIRE(bs3->id() == 3);

    REQUIRE(!bs1->get_default_branch());
    REQUIRE(!bs2->get_default_branch());
    REQUIRE(!bs3->get_default_branch());

    REQUIRE(bs1->get_branches().size() == 3);
    REQUIRE(bs1->get_branches()[0].cond.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs1->get_branches()[0].target == bs1);
    REQUIRE(bs1->get_branches()[0].id == 1);
    REQUIRE(bs1->get_branches()[1].cond.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs1->get_branches()[1].target == bs3);
    REQUIRE(bs1->get_branches()[1].id == 2);
    REQUIRE(bs1->get_branches()[2].cond.get() == seq.get_slot(IR::Type::Float64, 2));
    REQUIRE(bs1->get_branches()[2].target == bs2);
    REQUIRE(bs1->get_branches()[2].id == 3);
    REQUIRE(bs2->get_branches().empty());
    REQUIRE(bs3->get_branches().empty());
}

TEST_CASE("basic_single_channel_assumption") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    bs->add_assume(Seq::Sign::Pos, m1->val(), 43);
    bs->add_assume(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 3);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_assumes().size() == 1);
    REQUIRE(bs->get_assumes().begin()->second.sign == Seq::Sign::Pos);
    REQUIRE(bs->get_assumes().begin()->first.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_assumes().begin()->second.id == 3);
}

TEST_CASE("basic_single_channel_assumption_merge") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t;
    seq.set_defval(1, 23);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(t), nullptr,
                            seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    bs->add_assume(Seq::Sign::Pos, m1->val(), 43);
    bs->add_assume(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 0), 3);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_assumes().size() == 1);
    REQUIRE(bs->get_assumes().begin()->second.sign == Seq::Sign::Pos);
    REQUIRE(bs->get_assumes().begin()->first.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_assumes().begin()->second.id == 43);
}

TEST_CASE("basic_single_channel_assumption_violation") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

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
    REQUIRE(p1->id() == 1);
    REQUIRE(p1->needs_oldval());
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(900)), 1);
    auto m1 = bs->add_measure(1, 2, bs->track_time(t), bs->new_measure(seq.env(), 1));
    REQUIRE(m1->id() == 2);
    REQUIRE(m1->is_measure());
    // `m1->val()` should be `0.14` but will be rounded to `0`
    // to be consistent with the time check.
    bs->add_assume(Seq::Sign::Pos, m1->val(), 43);

    seq.prepare();
    auto err = expect_error<Seq::Error>([&] {
        seq.optimize();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NonPosTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 43);
    REQUIRE(err.what() == "Positive time expected."s);
}

TEST_CASE("endtimes_normalize") {
    // Make sure we correctly nomalize the time even if we have only one to start with.
    Seq::Seq seq(*llvm_ctx);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1), 8);
    bs->add_endtime(bs->track_time(t));

    REQUIRE(bs->get_endtimes().size() == 1);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_endtimes().size() == 1);
    REQUIRE(bs->get_endtimes().front()->tconst == 11);
    REQUIRE(bs->get_endtimes().front()->terms.size() == 2);
    REQUIRE(bs->get_endtimes().front()->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_endtimes().front()->terms[0].id == 4);
    REQUIRE(bs->get_endtimes().front()->terms[1].sign == Seq::Sign::NonNeg);
    REQUIRE(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs->get_endtimes().front()->terms[1].id == 8);
}

TEST_CASE("endtimes_elim_diamond") {
    Seq::Seq seq(*llvm_ctx);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    bs->add_endtime(bs->track_time(t));
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    bs->add_endtime(bs->track_time(t));
    Seq::EventTime t2(t);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    bs->add_endtime(bs->track_time(t));
    t2.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1), 6);
    bs->add_endtime(bs->track_time(t2));
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1), 8);
    bs->add_endtime(bs->track_time(t));

    REQUIRE(bs->get_endtimes().size() == 5);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_endtimes().size() == 1);
    REQUIRE(bs->get_endtimes().front()->tconst == 11);
    REQUIRE(bs->get_endtimes().front()->terms.size() == 2);
    REQUIRE(bs->get_endtimes().front()->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_endtimes().front()->terms[0].id == 4);
    REQUIRE(bs->get_endtimes().front()->terms[1].sign == Seq::Sign::NonNeg);
    REQUIRE(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs->get_endtimes().front()->terms[1].id == 8);
}

TEST_CASE("endtimes_elim_Y") {
    Seq::Seq seq(*llvm_ctx);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    bs->add_endtime(bs->track_time(t));
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    bs->add_endtime(bs->track_time(t));
    Seq::EventTime t2(t);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    bs->add_endtime(bs->track_time(t));
    t2.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1), 6);
    bs->add_endtime(bs->track_time(t2));

    REQUIRE(bs->get_endtimes().size() == 4);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_endtimes().size() == 2);

    // t2 having smaller tconst will be sorted to the front.
    REQUIRE(bs->get_endtimes().back()->tconst == 6);
    REQUIRE(bs->get_endtimes().back()->terms.size() == 2);
    REQUIRE(bs->get_endtimes().back()->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(bs->get_endtimes().back()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_endtimes().back()->terms[0].id == 4);
    REQUIRE(bs->get_endtimes().back()->terms[1].sign == Seq::Sign::NonNeg);
    REQUIRE(bs->get_endtimes().back()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs->get_endtimes().back()->terms[1].id == 6);

    REQUIRE(bs->get_endtimes().front()->tconst == 11);
    REQUIRE(bs->get_endtimes().front()->terms.size() == 1);
    REQUIRE(bs->get_endtimes().front()->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_endtimes().front()->terms[0].id == 4);
}

TEST_CASE("endtimes_elim_repeat") {
    Seq::Seq seq(*llvm_ctx);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t(5);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(1.45)), 2);
    t.add_term(Seq::Sign::Pos, seq.get_slot(IR::Type::Float64, 0), 4);
    t.add_term(Seq::Sign::Pos, seq.get_const(IR::TagVal(5.45)), 5);
    t.add_term(Seq::Sign::NonNeg, seq.get_slot(IR::Type::Float64, 1), 7);
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));
    bs->add_endtime(bs->track_time(t));

    REQUIRE(bs->get_endtimes().size() == 5);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_endtimes().size() == 1);

    // t2 having smaller tconst will be sorted to the front.
    REQUIRE(bs->get_endtimes().front()->tconst == 11);
    REQUIRE(bs->get_endtimes().front()->terms.size() == 2);
    REQUIRE(bs->get_endtimes().front()->terms[0].sign == Seq::Sign::Pos);
    REQUIRE(bs->get_endtimes().front()->terms[0].var.get() == seq.get_slot(IR::Type::Float64, 0));
    REQUIRE(bs->get_endtimes().front()->terms[0].id == 4);
    REQUIRE(bs->get_endtimes().front()->terms[1].sign == Seq::Sign::NonNeg);
    REQUIRE(bs->get_endtimes().front()->terms[1].var.get() == seq.get_slot(IR::Type::Float64, 1));
    REQUIRE(bs->get_endtimes().front()->terms[1].id == 7);
}

TEST_CASE("endtimes_elim_nonpos") {
    Seq::Seq seq(*llvm_ctx);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);

    Seq::EventTime t(0);
    bs->add_endtime(bs->track_time(t));

    REQUIRE(bs->get_endtimes().size() == 1);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_endtimes().size() == 0);
}

TEST_CASE("cond_elim_true") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(10)), nullptr,
                            seq.get_const(IR::TagVal(0.3)), seq.get_const(IR::TagVal(true)));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(p1->cond());

    seq.prepare();
    seq.optimize();

    REQUIRE(p1->is_ordered);
    REQUIRE(p1->val()->is_const());
    REQUIRE(p1->val()->get_const().is(IR::TagVal(0.3)));
    REQUIRE(!p1->cond());
}

TEST_CASE("cond_elim_false") {
    Seq::Seq seq(*llvm_ctx);
    REQUIRE(seq.get_chn_id("test_chn", true) == 1);
    auto bs = seq.add_basicseq(1);
    REQUIRE(bs->id() == 1);
    auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(10)), nullptr,
                            seq.get_const(IR::TagVal(0.3)), seq.get_const(IR::TagVal(false)));
    REQUIRE(p1->id() == 1);
    REQUIRE(!p1->needs_oldval());
    REQUIRE(p1->cond());
    REQUIRE(bs->get_pulses(1)->size() == 1);

    seq.prepare();
    seq.optimize();

    REQUIRE(bs->get_pulses(1)->empty());
}

TEST_CASE("cond_no_forward_pulse") {
    auto do_test = [&] (bool has_cond) {
        // Note that it is actually legal to do this forwarding with a condition
        // but we needs to be very careful not to truncate the previous pulse
        // when the condition is faulse.
        Seq::Seq seq(*llvm_ctx);
        REQUIRE(seq.get_chn_id("test_chn", true) == 1);
        auto bs = seq.add_basicseq(1);
        REQUIRE(bs->id() == 1);
        auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(10)), nullptr,
                                seq.get_const(IR::TagVal(0.3)),
                                has_cond ? seq.get_slot(IR::Type::Bool, 0) : nullptr);
        REQUIRE(p1->id() == 1);
        REQUIRE(!p1->needs_oldval());
        if (has_cond)
            REQUIRE(p1->cond());
        else
            REQUIRE(!p1->cond());
        auto p2 = bs->add_pulse(1, 2, bs->track_time(Seq::EventTime(20)), nullptr,
                                [&] {
                                    IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
                                    builder.createRet(builder.createSub(0, builder.getConst(2)));
                                    return seq.get_call(builder.get(),
                                                        {Seq::Arg::create_arg(1)}, 2);
                                }());
        REQUIRE(p2->id() == 2);
        REQUIRE(p2->needs_oldval());
        REQUIRE(!p2->cond());

        seq.prepare();
        seq.optimize();

        REQUIRE(!p1->needs_oldval());
        REQUIRE(!p2->needs_oldval());
        REQUIRE(!bs->needs_oldval(p1));
        if (has_cond)
            REQUIRE(bs->needs_oldval(p2));
        else
            REQUIRE(!bs->needs_oldval(p2));
    };
    do_test(false);
    do_test(true);
}

TEST_CASE("cond_no_forward_measure") {
    auto do_test = [&] (bool has_cond) {
        // Note that it is actually legal to do this forwarding with a condition
        // but we needs to be very careful not to truncate the previous pulse
        // when the condition is faulse.
        Seq::Seq seq(*llvm_ctx);
        REQUIRE(seq.get_chn_id("test_chn", true) == 1);
        auto bs = seq.add_basicseq(1);
        REQUIRE(bs->id() == 1);
        auto p1 = bs->add_pulse(1, 1, bs->track_time(Seq::EventTime(10)), nullptr,
                                seq.get_const(IR::TagVal(0.3)),
                                has_cond ? seq.get_slot(IR::Type::Bool, 0) : nullptr);
        REQUIRE(p1->id() == 1);
        REQUIRE(!p1->needs_oldval());
        if (has_cond)
            REQUIRE(p1->cond());
        else
            REQUIRE(!p1->cond());
        auto m1 = bs->add_measure(1, 5, bs->track_time(Seq::EventTime(20)),
                                  bs->new_measure(seq.env(), 1));
        REQUIRE(m1->id() == 5);
        auto p2 = bs->add_pulse(1, 7, bs->track_time(Seq::EventTime(30)), nullptr,
                                [&] {
                                    IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
                                    builder.createRet(builder.createSub(0, builder.getConst(2)));
                                    return seq.get_call(builder.get(),
                                                        {Seq::Arg::create_var(m1->val())}, 2);
                                }());
        REQUIRE(p2->id() == 7);
        REQUIRE(!p2->needs_oldval());
        REQUIRE(!p2->cond());
        REQUIRE(bs->get_pulses(1)->size() == 3);

        seq.prepare();
        seq.optimize();

        REQUIRE(bs->get_pulses(1)->size() == (has_cond ? 3 : 2));
        if (has_cond) {
            REQUIRE(!p1->needs_oldval());
            REQUIRE(!p2->needs_oldval());
            REQUIRE(!bs->needs_oldval(p1));
            REQUIRE(!bs->needs_oldval(p2));
            REQUIRE(p1->cond());
            REQUIRE(p2->val()->is_call());
            REQUIRE(p2->val()->args().size() == 1);
            REQUIRE(p2->val()->args()[0].is_var());
            REQUIRE(p2->val()->args()[0].get_var() == m1->val());
        }
        else {
            REQUIRE(p2->val()->is_const());
            REQUIRE(p2->val()->get_const().get<double>() == Approx(-1.7));
        }
    };
    do_test(false);
    do_test(true);
}
