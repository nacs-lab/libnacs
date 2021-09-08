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

#include "ir_helper.h"

using namespace NaCs;

TEST_CASE("Convert") {
    auto ctx = IR::ExeContext::get();
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Int32});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Bool (Int32 %0) {\n"
                    "L0:\n"
                    "  ret Int32 %0\n"
                    "}");
        TestIR<bool(int)>([] (int v) {
            return v != 0;
        }, {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}, *ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Bool (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        TestIR<bool(double)>([] (double v) {
            return v != 0;
        }, {-2.0, -1.1, -1, -0.3, 0, 0.3, 1, 1.1, 2.0}, *ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Int32, {IR::Type::Bool});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Int32 (Bool %0) {\n"
                    "L0:\n"
                    "  ret Bool %0\n"
                    "}");
        TestIR<int(bool)>([] (bool v) {
            return v ? 1 : 0;
        }, {false, true}, *ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Int32, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Int32 (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        TestIR<int(double)>([] (double v) {
            return (int)v;
        }, {-2.0, -1.9, -1.1, -1.0, -0.9, -0.1, 0, 0.1, 0.9, 1.0, 1.1, 1.9, 2.0},
            *ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Bool});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Bool %0) {\n"
                    "L0:\n"
                    "  ret Bool %0\n"
                    "}");
        TestIR<double(bool)>([] (bool v) {
            return v ? 1 : 0;
        }, {false, true}, *ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Int32});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Int32 %0) {\n"
                    "L0:\n"
                    "  ret Int32 %0\n"
                    "}");
        TestIR<double(int)>([] (int v) {
            return (double)v;
        }, {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}, *ctx, builder.get());
    }
}
