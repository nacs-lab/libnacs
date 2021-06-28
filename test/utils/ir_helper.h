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

#ifndef __NACS_TEST_IR_HELPER_H__
#define __NACS_TEST_IR_HELPER_H__

#include "winpath_helper.h"

#include "../../lib/nacs-utils/ir.h"
#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/streams.h"

#include <vector>
#include <sstream>
#include <functional>

#include <catch2/catch.hpp>

namespace {

using namespace NaCs;

template<typename T>
static void test_str_eq(T &&v, std::string val)
{
    auto str = sprint(std::forward<T>(v));
    REQUIRE(val == str);
}

namespace detail {

template<typename T1, typename T2>
static bool compare(T1 a, T2 b, bool _approx)
{
    if (_approx)
        return approx(a, b);
    return a == b;
}

template<typename T, size_t... Is>
static void _print_tuple(std::ostream &stm, T &&t, std::index_sequence<Is...>)
{
    auto f = [&] (auto &&v, size_t i) {
                 if (i != 0)
                     stm << ", ";
                 stm << v;
                 return 0;
             };
    (void)f;
    int l[] = { f(std::get<Is>(t), Is)... };
    (void)l;
}

template<typename T>
static void print_tuple(std::ostream &stm, T &&t)
{
    _print_tuple(stm, std::forward<T>(t), std::make_index_sequence<
                 std::tuple_size_v<std::decay_t<T>>>{});
}

template<bool _approx, typename T1, typename T2, typename Tuple>
static void _test_res(T1 res, T2 exp, const char *name, IR::Function &func,
                      Tuple &&args)
{
    if (!compare(res, exp, _approx)) {
        std::stringstream stm;
        stm << name << " test failed on test case: (";
        print_tuple(stm, args);
        stm << ")" << std::endl;
        stm << "Expected: " << (_approx ? "≈ " : "") << exp << std::endl;
        stm << "Got: " << res << std::endl;
        stm << "Code: " << func << std::endl;
        FAIL(stm.str());
    }
    else {
        REQUIRE(true);
    }
}

}

template<bool approx=false, typename T1, typename Res, typename... Args>
static void test_res(Res exp, const char *name, IR::Function &func,
                     T1 res, Args... args)
{
    detail::_test_res<approx>(res, exp, name, func, std::make_tuple(args...));
}

template<typename FT, bool approx>
struct IRTest;

template<typename Res, typename... Args, bool approx>
struct IRTest<Res(Args...), approx> {
    using FT = Res(Args...);
    IR::ExeContext *ir_ctx;
    IR::Function &func;
    std::function<FT> fexp;
    std::tuple<std::vector<Args>...> args;

    IRTest(std::function<FT> fexp, std::vector<Args> ...args,
           IR::ExeContext &ir_ctx, IR::Function &func)
        : ir_ctx(&ir_ctx),
          func(func),
          fexp(fexp),
          args(std::move(args)...)
    {
    }

    void test()
    {
        testeach_args(get_interp_func(), "Interpreted");
    }

    template<typename T1, typename... Args2>
    void test_res(const char *name, T1 res, Args2... args)
    {
        detail::_test_res<approx>(res, fexp(args...), name, func, std::make_tuple(args...));
    }

    auto get_interp_func() -> decltype(ir_ctx->getFunc<FT>(func))
    {
        return ir_ctx->getFunc<FT>(func);
    }

    template<typename F>
    void foreach_args(F &&f)
    {
        _foreach_args<sizeof...(Args)>(std::forward<F>(f));
    }

    template<typename F>
    void testeach_args(F &&f, const char *name)
    {
        foreach_args([&] (Args... args) {
                         test_call(name, f, args...);
                     });
    }

    template<typename F>
    void test_call(const char *name, F &&f, Args... args)
    {
        test_res(name, std::forward<F>(f)(args...), args...);
    }
private:
    template<size_t i, typename F, typename... Args2>
    std::enable_if_t<i == 0> _foreach_args(F &&f, Args2... args2)
    {
        f(args2...);
    }
    template<size_t i, typename F, typename... Args2>
    std::enable_if_t<(i > 0)> _foreach_args(F &&f, Args2... args2)
    {
        auto &_arg = std::get<i - 1>(args);
        for (auto v: _arg) {
            _foreach_args<i - 1>(std::forward<F>(f), v, args2...);
        }
    }
};

template<template<typename, bool> class Test, typename FT, bool approx=false>
struct MkTest;

template<template<typename, bool> class Test,
         typename Res, typename ...Args, bool approx>
struct MkTest<Test, Res(Args...), approx> {
    Test<Res(Args...), approx> test;

    template<typename... TestArgs>
    MkTest(std::function<Res(Args...)> f, std::vector<Args> ...args, TestArgs&&... targs)
        : test(f, std::move(args)..., std::forward<TestArgs>(targs)...)
    {
        test.test();
    }

    constexpr Test<Res(Args...), approx> &operator*()
    {
        return test;
    }

    constexpr const Test<Res(Args...), approx> &operator*() const
    {
        return test;
    }

    constexpr Test<Res(Args...), approx> *operator->()
    {
        return &test;
    }

    constexpr const Test<Res(Args...), approx> *operator->() const
    {
        return &test;
    }
};

template<typename FT, bool approx=false>
using TestIR = MkTest<IRTest, FT, approx>;

}

#endif
