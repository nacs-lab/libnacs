/*************************************************************************
 *   Copyright (c) 2016 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/nacs-seq/env.h"

#include "../../lib/nacs-utils/llvm/codegen.h"
#include "../../lib/nacs-utils/llvm/utils.h"
#include "../../lib/nacs-utils/number.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

using namespace NaCs;

static Seq::Env env = [] {
    static auto llvm_ctx = LLVM::new_context();
    return Seq::Env(*llvm_ctx);
}();

TEST_CASE("const_call") {
    auto v1 = env.new_const(true);
    auto v2 = env.new_call(v1, {Seq::Arg::create_var(v1)}, 1)->ref();
    REQUIRE(v2->argument_unused(0));
    env.optimize();
    REQUIRE(env.num_vars() == 1);
    REQUIRE(v2->varid() == 0);
    REQUIRE(v2->type() == IR::Type::Bool);
    REQUIRE(v2->is_const());
    REQUIRE(v2->get_const().is(true));
}

TEST_CASE("return_arg") {
    auto v1 = env.new_extern({IR::Type::Float64, 1234})->ref();
    auto v2 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Bool, IR::Type::Float64});
        auto pass_bb = builder.newBB();
        auto fail_bb = builder.newBB();
        builder.createBr(0, pass_bb, fail_bb);
        builder.curBB() = pass_bb;
        builder.createRet(builder.createCall(IR::Builtins::exp10, {1}));
        builder.curBB() = fail_bb;
        builder.createRet(1);
        return env.new_call(builder.get(), {Seq::Arg::create_const(false),
                Seq::Arg::create_var(v1.get())});
    }()->ref();
    env.optimize();
    REQUIRE(env.num_vars() == 2);
    REQUIRE(v2->varid() == 1);
    REQUIRE(v2->type() == IR::Type::Float64);
    REQUIRE(v2->is_call());
    REQUIRE(v2->get_callee().var == v1.get());
    REQUIRE(v2->get_assigned_var() == v1.get());
}

TEST_CASE("return_arg_indirect") {
    // The argument v1 dependency can only be removed after
    // LLVM optimization of the v3 function
    // Test if the local optimization is correctly handling this rescanning.
    // See the caller of `Var::optimize_call`.
    auto v1 = env.new_extern({IR::Type::Float64, 1234});
    auto v2 = env.new_extern({IR::Type::Float64, 1000});
    auto v3 = [&] {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Float64, IR::Type::Float64,
                             IR::Type::Float64, IR::Type::Float64});
        auto zero = builder.createSub(0, 1);
        zero = builder.createMul(zero, 2);
        auto one = builder.createAdd(zero, builder.getConst(1.0));
        builder.createRet(builder.createAdd(one, 3));
        return env.new_call(builder.get(), {Seq::Arg::create_arg(0),
                Seq::Arg::create_arg(0), Seq::Arg::create_var(v1),
                Seq::Arg::create_var(v2)}, 1);
    }()->ref();
    env.optimize();
    REQUIRE(env.num_vars() == 2);
    REQUIRE(v2->varid() == 0);
    REQUIRE(v2->type() == IR::Type::Float64);
    REQUIRE(v3->varid() == 1);
    REQUIRE(v3->type() == IR::Type::Float64);
    REQUIRE(v3->is_call());
    REQUIRE(v3->get_callee().is_llvm);
    REQUIRE(v3->nfreeargs() == 0);
    auto arg3 = v3->args();
    REQUIRE(arg3.size() == 1);
    REQUIRE(arg3[0].is_var());
    REQUIRE(arg3[0].get_var() == v2);
}

TEST_CASE("return_const") {
    auto v1 = env.new_const(3);
    auto v2 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Bool, IR::Type::Float64});
        auto pass_bb = builder.newBB();
        auto fail_bb = builder.newBB();
        builder.createBr(0, pass_bb, fail_bb);
        builder.curBB() = pass_bb;
        builder.createRet(builder.createCall(IR::Builtins::exp10, {1}));
        builder.curBB() = fail_bb;
        builder.createRet(1);
        return env.new_call(builder.get(), {Seq::Arg::create_const(true),
                Seq::Arg::create_var(v1)});
    }()->ref();
    env.optimize();
    REQUIRE(env.num_vars() == 1);
    REQUIRE(v2->varid() == 0);
    REQUIRE(v2->type() == IR::Type::Float64);
    REQUIRE(v2->is_const());
    REQUIRE(v2->get_const().get<double>() == Approx(1000));
}

TEST_CASE("unused_arg") {
    auto v1 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(builder.createAdd(0, 0));
        return env.new_call(builder.get(), {Seq::Arg::create_arg(1)}, 2);
    }()->ref();
    auto v2 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createAdd(1, 1));
        return env.new_call(builder.get(), {Seq::Arg::create_arg(0), Seq::Arg::create_arg(1)}, 2);
    }()->ref();

    env.optimize();

    REQUIRE(v1->nfreeargs() == 2);
    REQUIRE(v1->argument_unused(0));

    REQUIRE(v2->nfreeargs() == 2);
    REQUIRE(v2->argument_unused(0));
}

TEST_CASE("cse") {
    auto v1 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createFDiv(builder.createCall(IR::Builtins::exp10, {0}),
                                             builder.createCall(IR::Builtins::exp10, {1})));
        return env.new_call(builder.get(), {Seq::Arg::create_arg(0),
                Seq::Arg::create_arg(0)}, 1);
    }()->ref();
    auto v2 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
        auto s = builder.createSub(builder.createCall(IR::Builtins::exp10, {0}),
                                   builder.createCall(IR::Builtins::exp10, {1}));
        builder.createRet(builder.createMul(s, 0));
        return env.new_call(builder.get(), {Seq::Arg::create_arg(0),
                Seq::Arg::create_arg(0)}, 1);
    }()->ref();
    env.optimize();
    REQUIRE(env.num_vars() == 2);
    REQUIRE(v1->varid() == 0);
    REQUIRE(v1->type() == IR::Type::Float64);
    REQUIRE(v1->is_const());
    REQUIRE(v1->get_const().is(1.0));
    REQUIRE(v2->varid() == 1);
    REQUIRE(v2->type() == IR::Type::Float64);
    REQUIRE(v2->is_const());
    REQUIRE(v2->get_const().is(0.0));
}
