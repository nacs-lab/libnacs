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
                     std::function<double(Args...)> cb, const std::vector<Args>&... args)
    {
        auto ctx = IR::ExeContext::get();
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
        TestIR<double(Args...), true>(cb, args..., *ctx, builder.get());
    }
};

TEST_CASE("Call") {
    CallTester<double>::test(IR::Builtins::acos, "acos",
                             static_cast<double(*)(double)>(::acos), {-1, -0.5, 0, 0.5, 1});
    CallTester<double>::test(IR::Builtins::acosh, "acosh",
                             static_cast<double(*)(double)>(::acosh), {2.5, 2, 1.5, 1});
    CallTester<double>::test(IR::Builtins::asin, "asin",
                             static_cast<double(*)(double)>(::asin), {-1, -0.5, 0, 0.5, 1});
    CallTester<double>::test(IR::Builtins::asinh, "asinh",
                             static_cast<double(*)(double)>(::asinh), {-1, -1.5, 0, 1.5, 1});
    CallTester<double>::test(IR::Builtins::atan, "atan",
                             static_cast<double(*)(double)>(::atan), {-1, -1.5, 0, 1.5, 1});
    CallTester<double>::test(IR::Builtins::atanh, "atanh",
                             static_cast<double(*)(double)>(::atanh),
                             {-0.9, -0.5, 0, 0.5, 0.9});
    CallTester<double>::test(IR::Builtins::cbrt, "cbrt",
                             static_cast<double(*)(double)>(::cbrt), {-10, -5, 0, 5, 10});
    CallTester<double>::test(IR::Builtins::ceil, "ceil",
                             static_cast<double(*)(double)>(::ceil),
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::cos, "cos",
                             static_cast<double(*)(double)>(::cos),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::cosh, "cosh",
                             static_cast<double(*)(double)>(::cosh),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::erf, "erf",
                             static_cast<double(*)(double)>(::erf),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::erfc, "erfc",
                             static_cast<double(*)(double)>(::erfc),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp, "exp",
                             static_cast<double(*)(double)>(::exp),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp10, "exp10",
                             [] (double x) { return pow(10, x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::exp2, "exp2",
                             [] (double x) { return pow(2, x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::expm1, "expm1",
                             static_cast<double(*)(double)>(::expm1),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::abs, "abs", [] (double x) { return abs(x); },
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::floor, "floor",
                             static_cast<double(*)(double)>(::floor),
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::gamma, "gamma",
                             static_cast<double(*)(double)>(::tgamma),
                             {-10.3, -5.2, -1.1, -0.9, -0.1, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::j0, "j0",
                             static_cast<double(*)(double)>(::j0),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::j1, "j1",
                             static_cast<double(*)(double)>(::j1),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::lgamma, "lgamma",
                             static_cast<double(*)(double)>(::lgamma),
                             {-10.3, -5.2, -1.1, -0.9, -0.1, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log, "log",
                             static_cast<double(*)(double)>(::log),
                             {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log10, "log10",
                             static_cast<double(*)(double)>(::log10),
                             {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::log1p, "log1p",
                             static_cast<double(*)(double)>(::log1p),
                             {-0.9, -0.1, 0.1, 4, 9});
    CallTester<double>::test(IR::Builtins::log2, "log2",
                             static_cast<double(*)(double)>(::log2),
                             {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::rint, "rint",
                             static_cast<double(*)(double)>(::rint),
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::round, "round",
                             static_cast<double(*)(double)>(::round),
                             {-1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1});
    CallTester<double>::test(IR::Builtins::sin, "sin",
                             static_cast<double(*)(double)>(::sin),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::sinh, "sinh",
                             static_cast<double(*)(double)>(::sinh),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::sqrt, "sqrt",
                             static_cast<double(*)(double)>(::sqrt),
                             {0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::tan, "tan",
                             static_cast<double(*)(double)>(::tan),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::tanh, "tanh",
                             static_cast<double(*)(double)>(::tanh),
                             {-10, -5, -1.1, -0.9, -0.1, 0, 0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::y0, "y0",
                             static_cast<double(*)(double)>(::y0), {0.1, 0.9, 1.1, 5, 10});
    CallTester<double>::test(IR::Builtins::y1, "y1",
                             static_cast<double(*)(double)>(::y1), {0.1, 0.9, 1.1, 5, 10});

    CallTester<double,double>::test(IR::Builtins::atan2, "atan2",
                                    static_cast<double(*)(double,double)>(::atan2),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::copysign, "copysign",
                                    static_cast<double(*)(double,double)>(::copysign),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::fdim, "fdim",
                                    static_cast<double(*)(double,double)>(::fdim),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::max, "max",
                                    static_cast<double(*)(double,double)>(::fmax),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::min, "min",
                                    static_cast<double(*)(double,double)>(::fmin),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::mod, "mod",
                                    static_cast<double(*)(double,double)>(::fmod),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::hypot, "hypot",
                                    static_cast<double(*)(double,double)>(::hypot),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::pow, "pow",
                                    static_cast<double(*)(double,double)>(::pow),
                                    {-10, -5, -1, 1, 5, 10},
                                    {-10, -5, -1, 0, 1, 5, 10});
    CallTester<double,double>::test(IR::Builtins::remainder, "remainder",
                                    static_cast<double(*)(double,double)>(::remainder),
                                    {-10, -5, -1, 0, 1, 5, 10},
                                    {-10, -5, -1, 1, 5, 10});

    CallTester<double,double,double>::test(IR::Builtins::fma, "fma",
                                           static_cast<double(*)(double,double,double)>(::fma),
                                           {-10, -5, -1, 0, 1, 5, 10},
                                           {-10, -5, -1, 0, 1, 5, 10},
                                           {-10, -5, -1, 0, 1, 5, 10});

    CallTester<double,int>::test(IR::Builtins::ldexp, "ldexp",
                                 static_cast<double(*)(double,int)>(::ldexp),
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {-10, -5, -1, 0, 1, 5, 10});

    CallTester<int,double>::test(IR::Builtins::jn, "jn",
                                 static_cast<double(*)(int,double)>(::jn),
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {-10, -5, -1, 0, 1, 5, 10});
    CallTester<int,double>::test(IR::Builtins::yn, "yn",
                                 static_cast<double(*)(int,double)>(::yn),
                                 {-10, -5, -1, 0, 1, 5, 10},
                                 {1, 5, 10});
}
