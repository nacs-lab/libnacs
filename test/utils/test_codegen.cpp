/*************************************************************************
 *   Copyright (c) 2016 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "codegen_helper.h"

#include "../../lib/utils/number.h"
#include "../../lib/utils/streams.h"
#include "../../lib/utils/timer.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <math.h>

template<typename T>
static void test_str_eq(T &&v, std::string val)
{
    auto str = sprint(std::forward<T>(v));
    if (val == str)
        return;
    std::cerr << "Test failed:" << std::endl;
    std::cerr << "  Expect: \"" << val << "\"" << std::endl;
    std::cerr << "  Got: \"" << str << "\"" << std::endl;
    abort();
}

int main()
{
    TestCtx ctx;

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        test_codegen<double(double)>(ctx, builder.get(),
                                     test_set(1.2, 1.2),
                                     test_set(4.2, 4.2));
    }

    {
        IR::Builder builder(IR::Type::Bool, {});
        builder.createRet(IR::Consts::False);
        test_str_eq(builder.get(), "Bool () {\n"
                    "L0:\n"
                    "  ret Bool false\n"
                    "}");
        test_codegen<bool()>(ctx, builder.get(), test_set(false));
    }

    {
        IR::Builder builder(IR::Type::Float64, {});
        builder.createRet(builder.getConstFloat(1.1));
        test_str_eq(builder.get(), "Float64 () {\n"
                    "L0:\n"
                    "  ret Float64 1.1\n"
                    "}");
        test_codegen<double()>(ctx, builder.get(), test_set(1.1));
    }

    {
        IR::Builder builder(IR::Type::Int32, {});
        builder.createRet(builder.getConstInt(42));
        test_str_eq(builder.get(), "Int32 () {\n"
                    "L0:\n"
                    "  ret Int32 42\n"
                    "}");
        test_codegen<int()>(ctx, builder.get(), test_set(42));
    }

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
        test_str_eq(builder.get(), "Float64 (Bool %0, Float64 %1) {\n"
                    "L0:\n"
                    "  br Bool %0, L1, L2\n"
                    "L1:\n"
                    "  ret Float64 %1\n"
                    "L2:\n"
                    "  ret Float64 3.4\n"
                    "}");
        test_codegen<double(bool, double)>(ctx, builder.get(),
                                           test_set(1.3, true, 1.3),
                                           test_set(3.4, false, 1.3));
    }

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Float64, IR::Type::Float64});
        auto val1 = builder.createAdd(builder.getConstFloat(3.4), 0);
        auto val2 = builder.createMul(val1, 1);
        auto val3 = builder.createSub(val1, val2);
        builder.createRet(val3);
        test_str_eq(builder.get(), "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %2 = add Float64 3.4, Float64 %0\n"
                    "  Float64 %3 = mul Float64 %2, Float64 %1\n"
                    "  Float64 %4 = sub Float64 %2, Float64 %3\n"
                    "  ret Float64 %4\n"
                    "}");
        test_codegen<double(double, double)>(ctx, builder.get(),
                                             test_set(-1.71, 2.3, 1.3));
    }

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Int32, IR::Type::Int32});
        auto val1 = builder.createFDiv(0, 1);
        builder.createRet(val1);
        test_str_eq(builder.get(), "Float64 (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  Float64 %2 = fdiv Int32 %0, Int32 %1\n"
                    "  ret Float64 %2\n"
                    "}");
        test_codegen<double(int, int)>(ctx, builder.get(),
                                       test_set(1.5, 3, 2));
    }

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
        test_str_eq(builder.get(), "Float64 (Int32 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp ge Int32 %0, Float64 %1\n"
                    "  br Bool %2, L1, L2\n"
                    "L1:\n"
                    "  ret Float64 %1\n"
                    "L2:\n"
                    "  ret Int32 %0\n"
                    "}");
        test_codegen<double(int, double)>(ctx, builder.get(),
                                          test_set(1.3, 20, 1.3),
                                          test_set(-10, -10, 1.3));
    }

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
        test_str_eq(builder.get(), "Int32 (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  br L1\n"
                    "L1:\n"
                    "  Int32 %2 = phi [ L0: Int32 %0 ], [ L1: Int32 %4 ]\n"
                    "  Int32 %3 = phi [ L0: Int32 0 ], [ L1: Int32 %5 ]\n"
                    "  Int32 %4 = add Int32 %2, Int32 1\n"
                    "  Int32 %5 = add Int32 %3, Int32 %2\n"
                    "  Bool %6 = cmp gt Int32 %4, Int32 %1\n"
                    "  br Bool %6, L2, L1\n"
                    "L2:\n"
                    "  ret Int32 %5\n"
                    "}");
        auto test = test_codegen<int(int, int)>(ctx, builder.get(),
                                                test_set(6, 1, 3),
                                                test_set(500499, 2, 1000));

        auto f = test.get_interp_func();
        Timer timer;
        f(2, 1000);
        timer.print();
        auto f2 = test.get_llvm_func();
        timer.restart();
        f2(2, 1000);
        timer.print();
    }

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Int32});
        auto v2 = builder.createMul(0, builder.getConstInt(2));
        auto s1 = builder.createCall(IR::Builtins::sin, {v2});
        auto s2 = builder.createCall(IR::Builtins::sin, {0});
        builder.createRet(builder.createAdd(s2, s1));
        test_str_eq(builder.get(), "Float64 (Int32 %0) {\n"
                    "L0:\n"
                    "  Int32 %1 = mul Int32 %0, Int32 2\n"
                    "  Float64 %2 = call sin(Int32 %1)\n"
                    "  Float64 %3 = call sin(Int32 %0)\n"
                    "  Float64 %4 = add Float64 %3, Float64 %2\n"
                    "  ret Float64 %4\n"
                    "}");
        auto test = test_codegen<double(int), true>(ctx, builder.get(),
                                                    test_set(sin(1) + sin(2), 1),
                                                    test_set(sin(4) + sin(2), 2));

        auto f1 = test.get_interp_func();
        Timer timer;
        for (int i = 0;i < 1000000;i++)
            f1(1);
        timer.print();
        auto f2 = test.get_llvm_func();
        timer.restart();
        for (int i = 0;i < 1000000;i++)
            f2(1);
        timer.print();
    }

    {
        const int32_t data[] = {
            3, 2, 7, 50529027, 197379, 2, 3, 0, 1072693248, 3, 0, 1073741824,
            1, 22, 4, 5, -3, 0, 5, 4, -3, 5, 5, 6, -4, 0, 3, 3, 4, 6, 6, 2, 3,
            -3, 1, 2
        };
        IR::Function newfunc((const uint32_t*)data, sizeof(data) / 4);
        test_str_eq(newfunc, "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %5 = sub Float64 1, Float64 %0\n"
                    "  Float64 %4 = mul Float64 1, Float64 %5\n"
                    "  Float64 %6 = mul Float64 2, Float64 %0\n"
                    "  Float64 %3 = add Float64 %4, Float64 %6\n"
                    "  Float64 %2 = fdiv Float64 %3, Float64 1\n"
                    "  ret Float64 %2\n"
                    "}");
        test_codegen<double(double, double)>(ctx, newfunc, test_set(3, 2, 3));
    }

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        const double data[] = {0.0, 0.1, 0.2, 0.6};
        auto val1 = builder.createInterp(0, 2, 3, sizeof(data) / sizeof(double), data);
        builder.createRet(val1);
        test_str_eq(builder.get(), "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  Float64 %1 = interp [2, (4) +3] (Float64 %0) {0, 0.1, 0.2, 0.6}\n"
                    "  ret Float64 %1\n"
                    "}");
        test_codegen<double(double), true>(
            ctx, builder.get(),
            test_set(linearInterpolate(2.3, 2, 3, 4, data), 2.3),
            test_set(linearInterpolate(3.5, 2, 3, 4, data), 3.5),
            test_set(linearInterpolate(4.4, 2, 3, 4, data), 4.4));
    }

    {
        const int32_t data[] = {
            3, 2, 4, 50529027, 1, 3, 0, 1072693248, 1, 13, 10, 3, 0, -3, -3, 0, 4, 3, 2, 3, 1,
            1, 2, 4, 0, 1072693248, 0, 1073741824, 0, 1074003968, 0, 1072693248
        };
        IR::Function newfunc((const uint32_t*)data, sizeof(data) / 4);
        test_str_eq(newfunc, "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %3 = interp [1, (4) +1] (Float64 %0) {1, 2, 2.5, 1}\n"
                    "  Float64 %2 = add Float64 %3, Float64 %1\n"
                    "  ret Float64 %2\n"
                    "}");
        test_codegen<double(double, double)>(ctx, newfunc, test_set(4.25, 1.5, 2));
    }

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Float64, IR::Type::Float64});
        auto val1 = builder.createAdd(builder.getConstFloat(3.4), 0);
        auto val2 = builder.createMul(val1, 1);
        auto val3 = builder.createSub(val1, val2);
        builder.createRet(val3);
        test_str_eq(builder.get(), "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %2 = add Float64 3.4, Float64 %0\n"
                    "  Float64 %3 = mul Float64 %2, Float64 %1\n"
                    "  Float64 %4 = sub Float64 %2, Float64 %3\n"
                    "  ret Float64 %4\n"
                    "}");
        auto test = test_codegen<double(double, double)>(ctx, builder.get(),
                                                         test_set(-1.71, 2.3, 1.3),
                                                         test_set(-51.3, 2.3, 10.0),
                                                         test_set(0.0, 1.3, 1.0));

        LLVM::Codegen::Wrapper wrap0{true};
        auto test0 = test.get_llvm_test(wrap0);
        auto f0 = (double(*)(double, double, IR::GenVal*))test0.get_ptr();
        test.test_res(-1.71, "NULL closure", f0(2.3, 1.3, nullptr), 2.3, 1.3);
        test.test_res(-51.3, "NULL closure", f0(2.3, 10.0, nullptr), 2.3, 10.0);
        test.test_res(0.0, "NULL closure", f0(1.3, 1.0, nullptr), 1.3, 1.0);

        IR::GenVal vals[2];

        LLVM::Codegen::Wrapper wrap11{true};
        wrap11.add_closure(0, 0);
        auto test11 = test.get_llvm_test(wrap11);
        auto f11 = (double(*)(double, IR::GenVal*))test11.get_ptr();
        vals[0] = IR::TagVal(2.3).val;
        test.test_res(-1.71, "Closure 11", f11(1.3, vals), 2.3, 1.3);
        test.test_res(-51.3, "Closure 11", f11(10.0, vals), 2.3, 10.0);
        vals[0] = IR::TagVal(1.3).val;
        test.test_res(0.0, "Closure 11", f11(1.0, vals), 1.3, 1.0);

        LLVM::Codegen::Wrapper wrap12{true};
        wrap12.add_closure(1, 0);
        auto test12 = test.get_llvm_test(wrap12);
        auto f12 = (double(*)(double, IR::GenVal*))test12.get_ptr();
        vals[0] = IR::TagVal(1.3).val;
        test.test_res(-1.71, "Closure 12", f12(2.3, vals), 2.3, 1.3);
        vals[0] = IR::TagVal(10.0).val;
        test.test_res(-51.3, "Closure 12", f12(2.3, vals), 2.3, 10.0);
        vals[0] = IR::TagVal(1.0).val;
        test.test_res(0.0, "Closure 12", f12(1.3, vals), 1.3, 1.0);

        LLVM::Codegen::Wrapper wrap2{true};
        wrap2.add_closure(1, 0)
            .add_closure(0, 1);
        auto test2 = test.get_llvm_test(wrap2);
        auto f2 = (double(*)(IR::GenVal*))test2.get_ptr();
        vals[0] = IR::TagVal(1.3).val;
        vals[1] = IR::TagVal(2.3).val;
        test.test_res(-1.71, "Closure 2", f2(vals), 2.3, 1.3);
        vals[0] = IR::TagVal(10.0).val;
        test.test_res(-51.3, "Closure 2", f2(vals), 2.3, 10.0);
        vals[0] = IR::TagVal(1.0).val;
        vals[1] = IR::TagVal(1.3).val;
        test.test_res(0.0, "Closure 2", f2(vals), 1.3, 1.0);
    }

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Bool, IR::Type::Int32, IR::Type::Float64});
        builder.createRet(builder.createSelect(0, 1, 2));
        test_str_eq(builder.get(), "Float64 (Bool %0, Int32 %1, Float64 %2) {\n"
                    "L0:\n"
                    "  Float64 %3 = select Bool %0, Int32 %1, Float64 %2\n"
                    "  ret Float64 %3\n"
                    "}");
        test_codegen<double(bool, int, double)>(ctx, builder.get(),
                                                test_set(1.0, true, 1, 2.3),
                                                test_set(2.3, false, 1, 2.3));
    }

    {
        IR::Builder builder(IR::Type::Int32, {IR::Type::Float64});
        builder.createRet(builder.createConvert(IR::Type::Int32, 0));
        test_str_eq(builder.get(), "Int32 (Float64 %0) {\n"
                    "L0:\n"
                    "  Int32 %1 = convert(Float64 %0)\n"
                    "  ret Int32 %1\n"
                    "}");
        test_codegen<int(double)>(ctx, builder.get(),
                                  test_set(2, 2.3),
                                  test_set(2, 2.9),
                                  test_set(10, 10));
    }

    {
        IR::Builder builder(IR::Type::Float64,
                            {IR::Type::Bool, IR::Type::Int32, IR::Type::Float64});
        auto c1 = builder.createCall(IR::Builtins::cos, {1});
        auto s2 = builder.createCall(IR::Builtins::sin, {2});
        builder.createRet(builder.createSelect(0, c1, s2));
        test_str_eq(builder.get(), "Float64 (Bool %0, Int32 %1, Float64 %2) {\n"
                    "L0:\n"
                    "  Float64 %3 = call cos(Int32 %1)\n"
                    "  Float64 %4 = call sin(Float64 %2)\n"
                    "  Float64 %5 = select Bool %0, Float64 %3, Float64 %4\n"
                    "  ret Float64 %5\n"
                    "}");
        test_codegen<double(bool, int, double), true>(
            ctx, builder.get(),
            test_set(cos(1), true, 1, 2.3),
            test_set(sin(2.3), false, 1, 2.3));
    }

    return 0;
}
