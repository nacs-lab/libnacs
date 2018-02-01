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

#ifdef NDEBUG
#  undef NDEBUG
#endif

#include <nacs-utils/ir.h>
#include <nacs-utils/number.h>
#include <nacs-utils/timer.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <math.h>

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

template<typename T> static T get_arg1(int i);
template<>
bool get_arg1<bool>(int i)
{
    return i % 2 == 0;
}
template<>
int32_t get_arg1<int32_t>(int i)
{
    return i * i + i * 2 - 1000000;
}
template<>
double get_arg1<double>(int i)
{
    return (i * 2.5 - 0.4) / (i + 2);
}
template<typename Res>
void test_res1(Res res, int i, const std::vector<IR::Type> &ids)
{
    switch (ids[i]) {
    case IR::Type::Bool:
        assert(res == (Res)get_arg1<bool>(i));
        return;
    case IR::Type::Int32:
        assert(res == (Res)get_arg1<int32_t>(i));
        return;
    case IR::Type::Float64:
        assert(res == (Res)get_arg1<double>(i));
        return;
    default:
        assert(0 && "Invalid type id");
    }
}

template<int arg, typename Ret, typename... Arg, int... I>
static inline void test_single_arg(IR::ExeContext *exectx, std::integer_sequence<int,I...>,
                                   const std::vector<IR::Type> &ids)
{
    IR::Builder builder(type_id<Ret>::value, ids);
    builder.createRet(arg);
    auto f = exectx->getFunc<Ret(Arg...)>(builder.get());
    test_res1(f(get_arg1<Arg>(I)...), arg, ids);
}

template<typename Ret, typename... Arg, int... I>
static void _test_cc_sig(IR::ExeContext *exectx, std::integer_sequence<int,I...> seq)
{
    std::vector<IR::Type> ids = {type_id<Arg>::value...};
    // (test_single_arg<I,Ret,Arg...>(exectx, seq, ids), ...);
    int dummy[] = {0, (test_single_arg<I,Ret,Arg...>(exectx, seq, ids), 0)...};
    (void)dummy;
}

template<typename Ret, typename... Arg>
static void test_cc_sig(IR::ExeContext *exectx)
{
    _test_cc_sig<Ret,Arg...>(exectx, std::make_integer_sequence<int,int(sizeof...(Arg))>());
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
        // test_cc_sig<bool, Arg...>(exectx);
        test_cc_sig<int32_t, Arg...>(exectx);
        test_cc_sig<double, Arg...>(exectx);
    }
};

int main()
{
    auto exectx = IR::ExeContext::get();

    Tester<1>::test(exectx.get());
    Tester<2>::test(exectx.get());
    Tester<3>::test(exectx.get());
    Tester<4>::test(exectx.get());
    Tester<5>::test(exectx.get());
    Tester<6>::test(exectx.get());
    Tester<7>::test(exectx.get());
    Tester<8>::test(exectx.get());
    Tester<9>::test(exectx.get());

    return 0;
}
