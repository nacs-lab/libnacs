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

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN

#include "codegen_helper.h"

#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/streams.h"
#include <iostream>
#include <sstream>
#include <math.h>

TEST_CASE("CodeGen") {
    TestCtx ctx;

#if IR_TESTSET == 0
    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        TestCodegen<double(double)>([] (double v) { return v; },
                                    {1.2, 4.2}, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Bool, {});
        builder.createRet(IR::Consts::False);
        test_str_eq(builder.get(), "Bool () {\n"
                    "L0:\n"
                    "  ret Bool false\n"
                    "}");
        TestCodegen<bool()>([] { return false; }, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Float64, {});
        builder.createRet(builder.getConstFloat(1.1));
        test_str_eq(builder.get(), "Float64 () {\n"
                    "L0:\n"
                    "  ret Float64 1.1\n"
                    "}");
        TestCodegen<double()>([] { return 1.1; }, ctx, builder.get());
    }
#elif IR_TESTSET == 1
    {
        IR::Builder builder(IR::Type::Int32, {});
        builder.createRet(builder.getConstInt(42));
        test_str_eq(builder.get(), "Int32 () {\n"
                    "L0:\n"
                    "  ret Int32 42\n"
                    "}");
        TestCodegen<int()>([] { return 42; }, ctx, builder.get());
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
        TestCodegen<double(bool, double)>([] (bool b, double v) {
                                              return b ? v : 3.4;
                                          }, {false, true}, {1.3, 4.5}, ctx, builder.get());
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
        TestCodegen<double(double, double), true>(
            [] (double v0, double v1) {
                auto v2 = 3.4 + v0;
                return v2 - v2 * v1;
            }, {-1.71, 2.3, 1.3}, {-1.71, 2.3, 1.3}, ctx, builder.get());
    }
#elif IR_TESTSET == 2
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
        TestCodegen<double(int, int)>(
            [] (int v0, int v1) { return double(v0) / v1; },
            {2, 3, 4}, {1, 2, 3}, ctx, builder.get());
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
        TestCodegen<double(int, double)>([] (int i, double v) { return i > v ? v : i; },
                                         {20, -10}, {1.3, 5.6}, ctx, builder.get());
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
        auto test = TestCodegen<int(int, int)>(
            [] (int a, int b) {
                int s = 0;
                for (; a <= b; a++)
                    s += a;
                return s;
            }, {1, 2}, {3, 1000}, ctx, builder.get());

        auto f = test->get_interp_func();
        BENCHMARK("Interp") {
            return f(2, 1000);
        };
        auto f2 = test->get_llvm_func();
        BENCHMARK("Compile") {
            return f2(2, 1000);
        };
    }
#elif IR_TESTSET == 3
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
        auto test = TestCodegen<double(int), true>(
            [] (int v) { return sin(v) + sin(2 * v); }, {1, 2}, ctx, builder.get());

        auto f1 = test->get_interp_func();
        BENCHMARK("Interp") {
            return f1(1);
        };
        auto f2 = test->get_llvm_func();
        BENCHMARK("Compile") {
            return f2(1);
        };
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
        TestCodegen<double(double, double)>(
            [] (double t, double) {
                return (1 - t) + 2 * t;
            }, {2, 3, 4}, {3, 4, 5, 6}, ctx, newfunc);
    }

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        const double data[] = {0.0, 0.1, 0.2, 0.6};
        auto val1 = builder.createInterp(0, builder.getConstFloat(2), builder.getConstFloat(3),
                                         sizeof(data) / sizeof(double), data);
        builder.createRet(val1);
        test_str_eq(builder.get(), "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  Float64 %1 = interp [Float64 2, (4) +Float64 3] (Float64 %0) {0, 0.1, 0.2, 0.6}\n"
                    "  ret Float64 %1\n"
                    "}");
        TestCodegen<double(double), true>(
            [&] (double x) { return linearInterpolate(x, 2, 3, 4, data); },
            {0.7, 2.3, 3.5, 4.4, 5.5}, ctx, builder.get());
    }
#elif IR_TESTSET == 4
    {
        const int32_t data[] = {
            3, 2, 4, 50529027, 1, 3, 0, 1072693248, 1, 13, 10, 3, 0, -3, -3, 0, 4, 3, 2, 3, 1,
            1, 2, 4, 0, 1072693248, 0, 1073741824, 0, 1074003968, 0, 1072693248
        };
        IR::Function newfunc((const uint32_t*)data, sizeof(data) / 4);
        test_str_eq(newfunc, "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %3 = interp [Float64 1, (4) +Float64 1] (Float64 %0) {1, 2, 2.5, 1}\n"
                    "  Float64 %2 = add Float64 %3, Float64 %1\n"
                    "  ret Float64 %2\n"
                    "}");
        const double interp_data[] = {1, 2, 2.5, 1};
        TestCodegen<double(double, double)>(
            [&] (double x, double y) {
                return linearInterpolate(x, 1, 1, 4, interp_data) + y;
            }, {0.1, 1.5, 2.3}, {2, 3, 5}, ctx, newfunc);
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
        TestCodegen<double(bool, int, double)>(
            [] (bool b, int i, double v) { return b ? i : v; },
            {true, false}, {1, 2}, {1.2, 2.3}, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Int32, {IR::Type::Float64});
        builder.createRet(builder.createConvert(IR::Type::Int32, 0));
        test_str_eq(builder.get(), "Int32 (Float64 %0) {\n"
                    "L0:\n"
                    "  Int32 %1 = convert(Float64 %0)\n"
                    "  ret Int32 %1\n"
                    "}");
        TestCodegen<int(double)>([] (double v) { return (int)v; },
                                 {2.3, 2.9, 10}, ctx, builder.get());
    }
#elif IR_TESTSET == 5
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
        TestCodegen<double(bool, int, double), true>(
            [] (bool b, int i, double v) {
                return b ? cos(i) : sin(v);
            }, {true, false}, {1, 2}, {2.3, 3.8}, ctx, builder.get());
    }
#endif
}
