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

#include "../../lib/seq/env.h"

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/utils.h"
#include "../../lib/utils/number.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

static void test_const_call(Seq::Env &env)
{
    auto v1 = env.new_const(true);
    auto v2 = env.new_call(v1, {Seq::Arg::create_var(v1)}, 1)->ref();
    assert(v2->argument_unused(0));
    env.optimize();
    assert(env.num_vars() == 1);
    assert(v2->varid() == 0);
    assert(v2->type() == IR::Type::Bool);
    assert(v2->is_const());
    assert(v2->get_const().is(true));
}

static void test_return_arg(Seq::Env &env)
{
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
    assert(env.num_vars() == 2);
    assert(v2->varid() == 1);
    assert(v2->type() == IR::Type::Float64);
    assert(v2->is_call());
    assert(v2->get_callee().var == v1.get());
    assert(v2->get_assigned_var() == v1.get());
}

static void test_return_arg_indirect(Seq::Env &env)
{
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
    assert(env.num_vars() == 2);
    assert(v2->varid() == 0);
    assert(v2->type() == IR::Type::Float64);
    assert(v3->varid() == 1);
    assert(v3->type() == IR::Type::Float64);
    assert(v3->is_call());
    assert(v3->get_callee().is_llvm);
    assert(v3->nfreeargs() == 0);
    auto arg3 = v3->args();
    assert(arg3.size() == 1);
    assert(arg3[0].is_var());
    assert(arg3[0].get_var() == v2);
}

static void test_return_const(Seq::Env &env)
{
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
    assert(env.num_vars() == 1);
    assert(v2->varid() == 0);
    assert(v2->type() == IR::Type::Float64);
    assert(v2->is_const());
    assert(approx(v2->get_const().get<double>(), 1000));
}

static void test_unused_arg(Seq::Env &env)
{
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

    assert(v1->nfreeargs() == 2);
    assert(v1->argument_unused(0));

    assert(v2->nfreeargs() == 2);
    assert(v2->argument_unused(0));
}

static void test_cse(Seq::Env &env)
{
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
    assert(env.num_vars() == 2);
    assert(v1->varid() == 0);
    assert(v1->type() == IR::Type::Float64);
    assert(v1->is_const());
    assert(v1->get_const().is(1.0));
    assert(v2->varid() == 1);
    assert(v2->type() == IR::Type::Float64);
    assert(v2->is_const());
    assert(v2->get_const().is(0.0));
}

int main()
{
    auto llvm_ctx = LLVM::new_context();
    Seq::Env env(*llvm_ctx);

    test_const_call(env);
    test_return_arg(env);
    test_return_arg_indirect(env);
    test_return_const(env);
    test_unused_arg(env);
    test_cse(env);

    return 0;
}
