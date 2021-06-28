/*************************************************************************
 *   Copyright (c) 2016 - 2017 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/nacs-utils/ir.h"
#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/timer.h"
#include <iostream>
#include <sstream>
#include <math.h>

#include <catch2/catch.hpp>

using namespace NaCs;

template<typename T> struct type_id;
template<> struct type_id<bool> {
    constexpr static auto value = IR::Type::Bool;
};
template<> struct type_id<int32_t> {
    constexpr static auto value = IR::Type::Int32;
};
template<> struct type_id<double> {
    constexpr static auto value = IR::Type::Float64;
};

template<typename T>
inline constexpr auto type_id_v = type_id<T>::value;

class ArgSetBase {
    virtual bool get_bool(int i) const = 0;
    virtual int32_t get_int32(int i) const = 0;
    virtual double get_float64(int i) const = 0;
public:
    template<typename T> T get(int i) const;
    template<typename Res>
    NACS_NOINLINE void test(Res res, int i, IR::Type t) const
    {
        switch (t) {
        case IR::Type::Bool:
            REQUIRE(res == (Res)get<bool>(i));
            return;
        case IR::Type::Int32:
            REQUIRE(res == (Res)get<int32_t>(i));
            return;
        case IR::Type::Float64:
            REQUIRE(res == (Res)get<double>(i));
            return;
        default:
            FAIL("Invalid type id");
        }
    }
};
template<>
bool ArgSetBase::get<bool>(int i) const
{
    return get_bool(i);
}
template<>
int32_t ArgSetBase::get<int32_t>(int i) const
{
    return get_int32(i);
}
template<>
double ArgSetBase::get<double>(int i) const
{
    return get_float64(i);
}

struct ArgSet1 :ArgSetBase {
    bool get_bool(int i) const override
    {
        return i % 2 == 0;
    }
    int32_t get_int32(int i) const override
    {
        return i * i + i * 2 - 1000000;
    }
    double get_float64(int i) const override
    {
        return (i * 2.5 - 0.4) / (i + 2);
    }
};

struct ArgSet2 :ArgSetBase {
    bool get_bool(int i) const override
    {
        return i % 2 != 0;
    }
    int32_t get_int32(int i) const override
    {
        return i * i + i * 2 + 2000000;
    }
    double get_float64(int i) const override
    {
        return (i * 2.5 - 0.4) / (i + 2) * 1000 + 20;
    }
};

template<typename Ret, typename... Arg, int... I>
static void test_single_arg(IR::ExeContext *exectx, std::integer_sequence<int,I...>,
                            int arg, const std::vector<IR::Type> &ids, const ArgSetBase &argset)
{
    IR::Builder builder(type_id_v<Ret>, ids);
    builder.createRet(arg);
    auto f = exectx->getFunc<Ret(Arg...)>(builder.get());
    argset.test(f(argset.get<Arg>(I)...), arg, ids[arg]);
}

template<typename Ret, typename... Arg, int... I>
static void test_const_ret(IR::ExeContext *exectx, std::integer_sequence<int,I...>,
                           const std::vector<IR::Type> &ids, const ArgSetBase &argset)
{
    IR::Builder builder(type_id_v<Ret>, ids);
    builder.createRet(builder.getConst(IR::TagVal(argset.get<Ret>(-1))));
    auto f = exectx->getFunc<Ret(Arg...)>(builder.get());
    argset.test(f(argset.get<Arg>(I)...), -1, type_id_v<Ret>);
}

template<typename Ret, typename... Arg, int... I>
static void _test_cc_sig(IR::ExeContext *exectx, std::integer_sequence<int,I...> seq,
                         const ArgSetBase &argset)
{
    std::vector<IR::Type> ids = {type_id_v<Arg>...};
    for (size_t i = 0; i < sizeof...(Arg); i++)
        test_single_arg<Ret,Arg...>(exectx, seq, int(i), ids, argset);
    test_const_ret<Ret,Arg...>(exectx, seq, ids, argset);
}

template<typename Ret, typename... Arg>
static void test_cc_sig(IR::ExeContext *exectx)
{
    _test_cc_sig<Ret,Arg...>(exectx, std::make_integer_sequence<int,int(sizeof...(Arg))>(),
                             ArgSet1());
    _test_cc_sig<Ret,Arg...>(exectx, std::make_integer_sequence<int,int(sizeof...(Arg))>(),
                             ArgSet2());
}

template<int nargs, typename... Arg>
struct Tester {
    static inline void test(IR::ExeContext *exectx)
    {
        // Tester<nargs - 1, Arg..., bool>::test(exectx);
        Tester<nargs - 1, Arg..., int32_t>::test(exectx);
        Tester<nargs - 1, Arg..., double>::test(exectx);
    }
};

template<typename... Arg>
struct Tester<-1, Arg...> {
    static inline void test(IR::ExeContext *exectx)
    {
        test_cc_sig<bool, Arg...>(exectx);
        test_cc_sig<int32_t, Arg...>(exectx);
        test_cc_sig<double, Arg...>(exectx);
    }
};

template<typename... ArgPrefix>
static void test_prefix(IR::ExeContext *exectx)
{
    Tester<0,ArgPrefix...>::test(exectx);
    Tester<1,ArgPrefix...>::test(exectx);
    Tester<2,ArgPrefix...>::test(exectx);
    Tester<3,ArgPrefix...>::test(exectx);
    Tester<4,ArgPrefix...>::test(exectx);
}

ANON_TEST_CASE() {
    auto exectx = IR::ExeContext::get();

#if IR_CC_TESTSET == 0
    test_prefix<>(exectx.get());
#elif IR_CC_TESTSET == 1
    // For Win64
    test_prefix<int,double>(exectx.get());
#elif IR_CC_TESTSET == 2
    test_prefix<int,double,double,int>(exectx.get());
#elif IR_CC_TESTSET == 3
    // For Linux x64
    test_prefix<double,int,double,double,int,double,
                double,int,double>(exectx.get());
#elif IR_CC_TESTSET == 4
    test_prefix<double,int,double,int,double,int,
                double,double,double,int,double>(exectx.get());
#elif IR_CC_TESTSET == 5
    // For Linux aarch64
    test_prefix<double,int,double,int,double,int,double,
                int,double,int,double>(exectx.get());
#elif IR_CC_TESTSET == 6
    test_prefix<double,int,double,int,double,int,
                double,int,double,int,double,int,double>(exectx.get());
#elif IR_CC_TESTSET == 7
    // For Linux arm
    test_prefix<double,double,double,int,double,double,double>(exectx.get());
#elif IR_CC_TESTSET == 8
    test_prefix<double,double,int,double,double,int,double,double,double>(exectx.get());
#else
#  error "Unknown test set."
#endif
}
