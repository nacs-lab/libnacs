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

#define CATCH_CONFIG_MAIN

#include "vec_codegen_helper.h"

using namespace NaCs;

TEST_CASE("Compare") {
    TestCtx ctx;

    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::eq, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp eq Float64 %0, Float64 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(double, double)>([] (double v1, double v2) {
            return v1 == v2;
        }, {1.0, 1.2, 2.0}, {1.0, 1.2, 2.0}, ctx, builder.get());
        test->test_vec_allarg();
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Int32, IR::Type::Int32});
        builder.createRet(builder.createCmp(IR::CmpType::ne, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp ne Int32 %0, Int32 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(int, int)>([] (int v1, int v2) {
            return v1 != v2;
        }, {1, 12, 2}, {1, 2, 3}, ctx, builder.get());
        test->test_vec_allarg();
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::gt, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp gt Float64 %0, Float64 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(double, double)>([] (double v1, double v2) {
            return v1 > v2;
        }, {1.0, 1.2, 2.0}, {1.0, 1.2, 2.0}, ctx, builder.get());
        test->test_vec_allarg();
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Int32, IR::Type::Int32});
        builder.createRet(builder.createCmp(IR::CmpType::lt, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp lt Int32 %0, Int32 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(int, int)>([] (int v1, int v2) {
            return v1 < v2;
        }, {1, 12, 2}, {1, 2, 3}, ctx, builder.get());
        test->test_vec_allarg();
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64, IR::Type::Float64});
        builder.createRet(builder.createCmp(IR::CmpType::ge, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp ge Float64 %0, Float64 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(double, double)>([] (double v1, double v2) {
            return v1 >= v2;
        }, {1.0, 1.2, 2.0}, {1.0, 1.2, 2.0}, ctx, builder.get());
        test->test_vec_allarg();
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Int32, IR::Type::Int32});
        builder.createRet(builder.createCmp(IR::CmpType::le, 0, 1));
        test_str_eq(builder.get(),
                    "Bool (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  Bool %2 = cmp le Int32 %0, Int32 %1\n"
                    "  ret Bool %2\n"
                    "}");
        auto test = TestVec<bool(int, int)>([] (int v1, int v2) {
            return v1 <= v2;
        }, {1, 12, 2}, {1, 2, 3}, ctx, builder.get());
        test->test_vec_allarg();
    }
}
