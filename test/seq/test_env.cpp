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

#include "../utils/ir_helper.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

int main()
{
    auto llvm_ctx = LLVM::new_context();
    Seq::Env env(*llvm_ctx);

    auto v1 = env.new_const(true);
    test_str_eq(*v1, "%0 = Bool true");
    auto v2 = env.new_const(-12);
    test_str_eq(*v2, "%1 = Int32 -12");
    auto v3 = [&] {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Bool,
                IR::Type::Float64, IR::Type::Float64});
        auto pass_bb = builder.newBB();
        auto fail_bb = builder.newBB();
        builder.createBr(0, pass_bb, fail_bb);
        builder.curBB() = pass_bb;
        builder.createRet(builder.createFDiv(builder.createCall(IR::Builtins::exp10, {1}),
                                             builder.createCall(IR::Builtins::exp10, {2})));
        builder.curBB() = fail_bb;
        builder.createRet(builder.getConstFloat(3.4));
        return env.new_call(builder.get(), {Seq::Arg::create_var(v1),
                Seq::Arg::create_arg(0), Seq::Arg::create_arg(0)}, 1);
    }()->ref();
    assert(!v3->argument_unused(0));
    test_str_eq(*v3, "*%2(arg[0]) = Float64 llvm:@f(%0, arg[0], arg[0])");
    auto v4 = env.new_extern({IR::Type::Float64, -1});
    test_str_eq(*v4, "%3 = Float64 extern(0xffffffffffffffff)");

    assert(env.num_vars() == 4);
    std::vector<Seq::Var*> vars = {v4, v3.get(), v2, v1};
    unsigned i = 0;
    for (auto v: env)
        assert(v == vars[i++]);
    assert(i == 4);
    assert(vars == std::vector<Seq::Var*>(env.begin(), env.end()));

    assert(v1->type() == IR::Type::Bool);
    assert(v2->type() == IR::Type::Int32);
    assert(v3->type() == IR::Type::Float64);
    assert(v4->type() == IR::Type::Float64);

    assert(v1->used(true));
    assert(!v2->used(true));
    assert(v3->used(true));
    assert(!v4->used(true));

    assert(v1->used(false));
    assert(!v2->used(false));
    assert(!v3->used(false));
    assert(!v4->used(false));

    assert(v1->varid() == 0);
    assert(v2->varid() == 1);
    assert(v3->varid() == 2);
    assert(v4->varid() == 3);

    env.gc();

    assert(env.num_vars() == 2);

    assert(v1->used(true));
    assert(v3->used(true));

    assert(v1->used(false));
    assert(!v3->used(false));

    assert(v1->varid() == 0);
    assert(v3->varid() == 1);

    // std::cout << env << std::endl;

    return 0;
}