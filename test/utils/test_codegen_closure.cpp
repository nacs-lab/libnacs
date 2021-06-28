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

#define CATCH_CONFIG_MAIN

#include "codegen_helper.h"

#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/streams.h"
#include "../../lib/nacs-utils/timer.h"
#include <iostream>
#include <sstream>
#include <math.h>

NACS_EXPORT() IR::GenVal global_vals[10];

TEST_CASE("Closure") {
    TestCtx ctx;

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
    auto test = TestCodegen<double(double, double), true>(
        [] (double t, double len) {
            double v2 = (3.4 + t);
            return v2 - v2 * len;
        }, {2.3, 1.3}, {1.3, 10, 1.0}, ctx, builder.get());

    LLVM::Codegen::Wrapper wrap0{true};
    auto test0 = test->get_llvm_test(wrap0);
    auto f0 = (double(*)(double, double, IR::GenVal*))test0.get_ptr();
    test->testeach_args([&] (double a, double b) { return f0(a, b, nullptr); },
                        "NULL closure");

    {
        IR::GenVal vals[3];

        LLVM::Codegen::Wrapper wrap11{true};
        wrap11.add_closure(0, 0);
        auto test11 = test->get_llvm_test(wrap11);
        auto f11 = (double(*)(double, IR::GenVal*))test11.get_ptr();
        test->testeach_args([&] (double a, double b) {
            vals[0] = IR::TagVal(a).val;
            return f11(b, vals);
        }, "Closure 11");

        LLVM::Codegen::Wrapper wrap12{true};
        wrap12.add_closure(1, 0);
        auto test12 = test->get_llvm_test(wrap12);
        auto f12 = (double(*)(double, IR::GenVal*))test12.get_ptr();
        test->testeach_args([&] (double a, double b) {
            vals[0] = IR::TagVal(b).val;
            return f12(a, vals);
        }, "Closure 12");

        LLVM::Codegen::Wrapper wrap2{true};
        wrap2.add_closure(1, 0)
            .add_closure(0, 1);
        auto test2 = test->get_llvm_test(wrap2);
        auto f2 = (double(*)(IR::GenVal*))test2.get_ptr();
        test->testeach_args([&] (double a, double b) {
            vals[0] = IR::TagVal(b).val;
            vals[1] = IR::TagVal(a).val;
            return f2(vals);
        }, "Closure 2");

        LLVM::Codegen::Wrapper wrap_ret{true};
        wrap_ret.add_closure(-1, 0)
            .add_closure(0, 1);
        auto test_ret = test->get_llvm_test(wrap_ret);
        auto f_ret = (void(*)(double, IR::GenVal*))test_ret.get_ptr();
        test->testeach_args([&] (double a, double b) {
            vals[1] = IR::TagVal(a).val;
            f_ret(b, vals);
            return vals[0].f64;
        }, "Closure Return");

        LLVM::Codegen::Wrapper wrap_ret2{true};
        wrap_ret2.add_closure(-1, 0)
            .add_closure(0, 2)
            .add_closure(1, 1);
        auto test_ret2 = test->get_llvm_test(wrap_ret2);
        auto f_ret2 = (void(*)(IR::GenVal*))test_ret2.get_ptr();
        test->testeach_args([&] (double a, double b) {
            vals[1] = IR::TagVal(b).val;
            vals[2] = IR::TagVal(a).val;
            f_ret2(vals);
            return vals[0].f64;
        }, "Closure Return 2");
    }

    {
        auto gv_ty = LLVM::get_array_type(test0.T_i8, 10 * sizeof(double));
        auto set_closure_ptr = [&] (auto &wrap, auto &test) {
            wrap.closure_ptr = LLVM::new_global_variable(*test.mod, gv_ty, false,
                                                         llvm::GlobalValue::ExternalLinkage,
                                                         nullptr, "global_vals");
        };

        auto test11 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap11{true};
        set_closure_ptr(wrap11, test11);
        wrap11.add_closure(0, 0);
        test11.set_function(test->func, wrap11);
        auto f11 = (double(*)(double))test11.get_ptr();
        test->testeach_args([&] (double a, double b) {
            global_vals[0] = IR::TagVal(a).val;
            return f11(b);
        }, "Closure 11");

        auto test12 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap12{true};
        set_closure_ptr(wrap12, test12);
        wrap12.add_closure(1, 0);
        test12.set_function(test->func, wrap12);
        auto f12 = (double(*)(double))test12.get_ptr();
        test->testeach_args([&] (double a, double b) {
            global_vals[0] = IR::TagVal(b).val;
            return f12(a);
        }, "Closure 12");

        auto test2 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap2{true};
        set_closure_ptr(wrap2, test2);
        wrap2.add_closure(1, 0)
            .add_closure(0, 1);
        test2.set_function(test->func, wrap2);
        auto f2 = (double(*)())test2.get_ptr();
        test->testeach_args([&] (double a, double b) {
            global_vals[0] = IR::TagVal(b).val;
            global_vals[1] = IR::TagVal(a).val;
            return f2();
        }, "Closure 2");

        auto test_ret = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap_ret{true};
        set_closure_ptr(wrap_ret, test_ret);
        wrap_ret.add_closure(-1, 0)
            .add_closure(0, 1);
        test_ret.set_function(test->func, wrap_ret);
        auto f_ret = (void(*)(double))test_ret.get_ptr();
        test->testeach_args([&] (double a, double b) {
            global_vals[1] = IR::TagVal(a).val;
            f_ret(b);
            return global_vals[0].f64;
        }, "Closure Return");

        auto test_ret2 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap_ret2{true};
        set_closure_ptr(wrap_ret2, test_ret2);
        wrap_ret2.add_closure(-1, 0)
            .add_closure(0, 2)
            .add_closure(1, 1);
        test_ret2.set_function(test->func, wrap_ret2);
        auto f_ret2 = (void(*)())test_ret2.get_ptr();
        test->testeach_args([&] (double a, double b) {
            global_vals[1] = IR::TagVal(b).val;
            global_vals[2] = IR::TagVal(a).val;
            f_ret2();
            return global_vals[0].f64;
        }, "Closure Return 2");
    }

    {
        auto gv_ty = LLVM::get_array_type(test0.T_i8, 10 * sizeof(double));
        auto set_closure_ptr = [&] (auto &wrap, auto &test) {
            wrap.closure_ptr = LLVM::new_global_variable(*test.mod, gv_ty, false,
                                                         llvm::GlobalValue::ExternalLinkage,
                                                         LLVM::get_aggregate_zero(gv_ty),
                                                         "local_vals");
        };
        IR::GenVal *local_vals;

        auto test11 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap11{true};
        set_closure_ptr(wrap11, test11);
        wrap11.add_closure(0, 0);
        test11.set_function(test->func, wrap11);
        auto f11 = (double(*)(double))test11.get_ptr();
        local_vals = (IR::GenVal*)test11.engine.get_symbol("local_vals");
        test->testeach_args([&] (double a, double b) {
            local_vals[0] = IR::TagVal(a).val;
            return f11(b);
        }, "Closure 11");

        auto test12 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap12{true};
        set_closure_ptr(wrap12, test12);
        wrap12.add_closure(1, 0);
        test12.set_function(test->func, wrap12);
        auto f12 = (double(*)(double))test12.get_ptr();
        local_vals = (IR::GenVal*)test12.engine.get_symbol("local_vals");
        test->testeach_args([&] (double a, double b) {
            local_vals[0] = IR::TagVal(b).val;
            return f12(a);
        }, "Closure 12");

        auto test2 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap2{true};
        set_closure_ptr(wrap2, test2);
        wrap2.add_closure(1, 0)
            .add_closure(0, 1);
        test2.set_function(test->func, wrap2);
        auto f2 = (double(*)())test2.get_ptr();
        local_vals = (IR::GenVal*)test2.engine.get_symbol("local_vals");
        test->testeach_args([&] (double a, double b) {
            local_vals[0] = IR::TagVal(b).val;
            local_vals[1] = IR::TagVal(a).val;
            return f2();
        }, "Closure 2");

        auto test_ret = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap_ret{true};
        set_closure_ptr(wrap_ret, test_ret);
        wrap_ret.add_closure(-1, 0)
            .add_closure(0, 1);
        test_ret.set_function(test->func, wrap_ret);
        auto f_ret = (void(*)(double))test_ret.get_ptr();
        local_vals = (IR::GenVal*)test_ret.engine.get_symbol("local_vals");
        test->testeach_args([&] (double a, double b) {
            local_vals[1] = IR::TagVal(a).val;
            f_ret(b);
            return local_vals[0].f64;
        }, "Closure Return");

        auto test_ret2 = test->get_empty_llvm_test();
        LLVM::Codegen::Wrapper wrap_ret2{true};
        set_closure_ptr(wrap_ret2, test_ret2);
        wrap_ret2.add_closure(-1, 0)
            .add_closure(0, 2)
            .add_closure(1, 1);
        test_ret2.set_function(test->func, wrap_ret2);
        auto f_ret2 = (void(*)())test_ret2.get_ptr();
        local_vals = (IR::GenVal*)test_ret2.engine.get_symbol("local_vals");
        test->testeach_args([&] (double a, double b) {
            local_vals[1] = IR::TagVal(b).val;
            local_vals[2] = IR::TagVal(a).val;
            f_ret2();
            return local_vals[0].f64;
        }, "Closure Return 2");
    }
}
