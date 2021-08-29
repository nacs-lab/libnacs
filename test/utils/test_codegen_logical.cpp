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

#include "codegen_helper.h"

using namespace NaCs;

TEST_CASE("Logical") {
    TestCtx ctx;

    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Bool, IR::Type::Bool});
        builder.createRet(builder.createAnd(0, 1));
        test_str_eq(builder.get(),
                    "Bool (Bool %0, Bool %1) {\n"
                    "L0:\n"
                    "  Bool %2 = and Bool %0, Bool %1\n"
                    "  ret Bool %2\n"
                    "}");
        TestCodegen<bool(bool, bool)>([] (bool v1, bool v2) {
            return v1 && v2;
        }, {false, true}, {false, true}, ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Bool, IR::Type::Bool});
        builder.createRet(builder.createOr(0, 1));
        test_str_eq(builder.get(),
                    "Bool (Bool %0, Bool %1) {\n"
                    "L0:\n"
                    "  Bool %2 = or Bool %0, Bool %1\n"
                    "  ret Bool %2\n"
                    "}");
        TestCodegen<bool(bool, bool)>([] (bool v1, bool v2) {
            return v1 || v2;
        }, {false, true}, {false, true}, ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Bool, IR::Type::Bool});
        builder.createRet(builder.createXor(0, 1));
        test_str_eq(builder.get(),
                    "Bool (Bool %0, Bool %1) {\n"
                    "L0:\n"
                    "  Bool %2 = xor Bool %0, Bool %1\n"
                    "  ret Bool %2\n"
                    "}");
        TestCodegen<bool(bool, bool)>([] (bool v1, bool v2) {
            return v1 ^ v2;
        }, {false, true}, {false, true}, ctx, builder.get());
    }
    {
        IR::Builder builder(IR::Type::Bool, {IR::Type::Bool});
        builder.createRet(builder.createNot(0));
        test_str_eq(builder.get(),
                    "Bool (Bool %0) {\n"
                    "L0:\n"
                    "  Bool %1 = not Bool %0\n"
                    "  ret Bool %1\n"
                    "}");
        TestCodegen<bool(bool)>([] (bool v) {
            return !v;
        }, {false, true}, ctx, builder.get());
    }
}
