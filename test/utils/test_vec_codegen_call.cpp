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

#include <functional>
#include <string>
#include <type_traits>

using namespace NaCs;

template<typename... Args>
struct CallTester {
    template<typename Arg>
    static IR::Type get_arg()
    {
        STATIC_REQUIRE(std::is_same_v<Arg, bool> || std::is_same_v<Arg, int> ||
                       std::is_same_v<Arg, double>);
        if (std::is_same_v<Arg, bool>)
            return IR::Type::Bool;
        if (std::is_same_v<Arg, int>)
            return IR::Type::Int32;
        return IR::Type::Float64;
    }
    static void test(IR::Builtins f, const std::string &name,
                     std::function<double(Args...)> cb, const std::vector<Args>&... args,
                     bool vecarg=true)
    {
        TestCtx ctx;

        INFO("Testing " + name);

        std::vector<IR::Type> args_tys{get_arg<Args>()...};
        IR::Builder builder(IR::Type::Float64, args_tys);
        int32_t args_idx[3] = {0, 1, 2};
        STATIC_REQUIRE(sizeof...(Args) <= 3);
        builder.createRet(builder.createCall(f, sizeof...(Args), args_idx));
        std::string expected = "Float64 (";
        std::string args_str = "";
        for (size_t i = 0; i < sizeof...(Args); i++) {
            if (i != 0)
                args_str += ", ";
            switch (args_tys[i]) {
            case IR::Type::Bool:
                args_str += "Bool";
                break;
            case IR::Type::Int32:
                args_str += "Int32";
                break;
            case IR::Type::Float64:
                args_str += "Float64";
                break;
            default:
                FAIL("Unknown argument type.");
                break;
            }
            args_str += " %" + std::to_string(i);
        }
        expected += (args_str + ") {\n"
                     "L0:\n"
                     "  Float64 %" + std::to_string(sizeof...(Args)) +
                     " = call " + name + "(" + args_str + ")\n"
                     "  ret Float64 %" + std::to_string(sizeof...(Args)) + "\n"
                     "}");

        test_str_eq(builder.get(), expected);
        auto test = TestVec<double(Args...), true>(cb, args..., ctx, builder.get());
        if (vecarg) {
            test->test_vec_allarg();
        }
    }
};

TEST_CASE("Call") {
    CallTester<double>::test(IR::Builtins::acos, "acos", ::acos, {-1, -0.5, 0, 0.5, 1});
    CallTester<double>::test(IR::Builtins::acosh, "acosh", ::acosh, {2.5, 2, 1.5, 1});
    CallTester<double>::test(IR::Builtins::asin, "asin", ::asin, {-1, -0.5, 0, 0.5, 1});
    CallTester<double>::test(IR::Builtins::asinh, "asinh", ::asinh, {-1, -1.5, 0, 1.5, 1});
    CallTester<double>::test(IR::Builtins::atan, "atan", ::atan, {-1, -1.5, 0, 1.5, 1});
    CallTester<double>::test(IR::Builtins::atanh, "atanh", ::atanh, {-0.9, -0.5, 0, 0.5, 0.9});
    CallTester<double>::test(IR::Builtins::cbrt, "cbrt", ::cbrt, {-10, -5, 0, 5, 10});
    CallTester<double>::test(IR::Builtins::ceil, "ceil", ::ceil,
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::cos, "cos", ::cos,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::cosh, "cosh", ::cosh,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::erf, "erf", ::erf,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::erfc, "erfc", ::erfc,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp, "exp", ::exp,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp10, "exp10", [] (double x) { return pow(10, x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp2, "exp2", [] (double x) { return pow(2, x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::expm1, "expm1", ::expm1,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::abs, "abs", [] (double x) { return abs(x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::floor, "floor", ::floor,
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::gamma, "gamma", ::tgamma,
                             {-10.3, -5.2, -1.1, -0.9, -0.1, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::j0, "j0", ::j0,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10}, false);
    CallTester<double>::test(IR::Builtins::j1, "j1", ::j1,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10}, false);
    CallTester<double>::test(IR::Builtins::lgamma, "lgamma", ::lgamma,
                             {-10.3, -5.2, -1.1, -0.9, -0.1, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log, "log", ::log, {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log10, "log10", ::log10, {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log1p, "log1p", ::log1p, {-0.9, -0.1, 0.1, 4, 9});
    CallTester<double>::test(IR::Builtins::log2, "log2", ::log2, {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::rint, "rint", ::rint,
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::round, "round", ::round,
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::sin, "sin", ::sin,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::sinh, "sinh", ::sinh,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::sqrt, "sqrt", ::sqrt, {0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::tan, "tan", ::tan,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::tanh, "tanh", ::tanh,
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::y0, "y0", ::y0, {0.1, 0.9, 1.1, 5, 10}, false);
    CallTester<double>::test(IR::Builtins::y1, "y1", ::y1, {0.1, 0.9, 1.1, 5, 10}, false);

    CallTester<double,double>::test(IR::Builtins::atan2, "atan2", ::atan2,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::copysign, "copysign", ::copysign,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::fdim, "fdim", ::fdim,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::max, "max", ::fmax,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::min, "min", ::fmin,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::mod, "mod", ::fmod,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::hypot, "hypot", ::hypot,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::pow, "pow", ::pow,
                                    {-10, -5, -1, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::remainder, "remainder", ::remainder,
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 1, 5, 10}
#ifndef NACS_SLEEF_HAS_REMAINDER
                                    , false
#endif
        );

    CallTester<double,double,double>::test(IR::Builtins::fma, "fma", ::fma,
                                           {-10, -5, -1, 0, 1, 5, 10},
                                           {-10, -5, -1, 0, 1, 5, 10},
                                           {-10, -5, -1, 0, 1, 5, 10});

    CallTester<double,int>::test(IR::Builtins::ldexp, "ldexp", ::ldexp,
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {-10, -5, -1, 0, 1, 5, 10});

    CallTester<int,double>::test(IR::Builtins::jn, "jn", ::jn,
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {-10, -5, -1, 0, 1, 5, 10}, false);
    CallTester<int,double>::test(IR::Builtins::yn, "yn", ::yn,
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {1, 5, 10}, false);
}
