/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#ifdef NDEBUG
#  undef NDEBUG
#endif

#include <nacs-utils/ir.h>
#include <nacs-utils/timer.h>
#include <assert.h>
#include <iostream>

using namespace NaCs;

int
main()
{
    assert(IR::validate(IR::Type::Bool));
    assert(IR::validate(IR::Type::Int32));
    assert(IR::validate(IR::Type::Float64));
    assert(!IR::validate(IR::Type(4)));

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(0);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({1.2});
        ctx.eval().dump();
        ctx.reset({4.2});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Bool, {});
        builder.createRet(IR::Consts::False);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64, {});
        builder.createRet(builder.getConstFloat(1.1));
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Int32, {});
        builder.createRet(builder.getConstInt(42));
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Bool,
                    IR::Type::Float64});
        auto pass_bb = builder.newBB();
        auto fail_bb = builder.newBB();
        builder.createBr(0, pass_bb, fail_bb);
        builder.curBB() = pass_bb;
        builder.createRet(1);
        builder.curBB() = fail_bb;
        builder.createRet(builder.getConstFloat(3.4));
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({true, 1.3});
        ctx.eval().dump();
        ctx.reset({false, 1.3});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Float64, IR::Type::Float64});
        auto val1 = builder.createAdd(builder.getConstFloat(3.4), 0);
        auto val2 = builder.createMul(val1, 1);
        auto val3 = builder.createSub(val1, val2);
        builder.createRet(val3);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({2.3, 1.3});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Int32, IR::Type::Int32});
        auto val1 = builder.createFDiv(0, 1);
        builder.createRet(val1);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({3, 2});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Int32,
                    IR::Type::Float64});
        auto pass_bb = builder.newBB();
        auto fail_bb = builder.newBB();
        builder.createBr(builder.createCmp(IR::CmpType::ge, 0, 1),
                         pass_bb, fail_bb);
        builder.curBB() = pass_bb;
        builder.createRet(1);
        builder.curBB() = fail_bb;
        builder.createRet(0);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({20, 1.3});
        ctx.eval().dump();
        ctx.reset({-10, 1.3});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Int32,
                            {IR::Type::Int32, IR::Type::Int32});
        auto loop_bb = builder.newBB();
        auto ret_bb = builder.newBB();
        builder.createBr(loop_bb);
        builder.curBB() = loop_bb;
        auto i = builder.createPhi(IR::Type::Int32, 2);
        auto s = builder.createPhi(IR::Type::Int32, 2);
        builder.addPhiInput(i.second, 0, 0);
        builder.addPhiInput(s.second, 0, builder.getConstInt(0));
        auto i2 = builder.createAdd(i.first, builder.getConstInt(1));
        auto s2 = builder.createAdd(s.first, i.first);
        auto cond = builder.createCmp(IR::CmpType::gt, i2, 1);
        builder.createBr(cond, ret_bb, loop_bb);
        builder.addPhiInput(i.second, loop_bb, i2);
        builder.addPhiInput(s.second, loop_bb, s2);
        builder.curBB() = ret_bb;
        builder.createRet(s2);
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({1, 3});
        ctx.eval().dump();
        ctx.reset({2, 1000});
        ctx.eval().dump();
    }
    std::cout << std::endl;

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Int32});
        builder.createRet(
            builder.createAdd(builder.createCall(IR::Builtins::sin, {0}),
                              builder.createCall(IR::Builtins::sin, {
                                      builder.createMul(0,
                                                        builder.getConstInt(2))
                                  })));
        builder.get().dump();
        IR::EvalContext ctx(builder.get());
        ctx.reset({1});
        ctx.eval().dump();
        ctx.reset({2});
        ctx.eval().dump();

        auto data = builder.get().serialize();
        IR::Function newfunc(data);
        newfunc.dump();
        IR::EvalContext ctx2(newfunc);
        ctx2.reset({1});
        ctx2.eval().dump();

        tic();
        for (int i = 0;i < 1000000;i++) {
            ctx2.reset({1});
            ctx2.eval();
        }
        printToc();
    }
    std::cout << std::endl;

    return 0;
}
