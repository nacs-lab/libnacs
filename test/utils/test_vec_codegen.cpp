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

#include "codegen_helper.h"

#include "../../lib/utils/number.h"
#include "../../lib/utils/streams.h"
#include "../../lib/utils/timer.h"
#include "../../lib/utils/processor.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <tuple>

#include <math.h>

#if NACS_OS_WINDOWS && (NACS_CPU_X86 || NACS_CPU_X86_64) && \
    !defined(__clang__) && !defined(_MSC_VER)
#  define REF_VEC_ONLY 1
#elif defined(__clang__) || defined(__NO_INLINE__)
// Clang doesn't handle calling convention change with target attribute
// Also the trick below won't wor without inlining
#  define REF_VEC_ONLY 1
#else
#  define REF_VEC_ONLY 0
#endif

static auto &host_info = CPUInfo::get_host();

template<size_t N>
struct vec_runner {
    template<typename F>
    static void run(F &&)
    {}
};

#if NACS_CPU_X86_64 || NACS_CPU_X86
template<>
struct vec_runner<2> {
    template<typename F>
    NACS_NOINLINE __attribute__((target("sse2"), flatten))
    static void run(F &&f)
    {
        f();
    }
};

template<>
struct vec_runner<4> {
    template<typename F>
    NACS_NOINLINE __attribute__((target("avx"), flatten))
    static void run(F &&f)
    {
        f();
    }
};

template<>
struct vec_runner<8> {
    template<typename F>
    NACS_NOINLINE __attribute__((target("avx512f"), flatten))
    static void run(F &&f)
    {
        f();
    }
};
#elif NACS_CPU_AARCH64
template<>
struct vec_runner<2> {
    template<typename F>
    NACS_NOINLINE __attribute__((flatten))
    static void run(F &&f)
    {
        f();
    }
};
#endif

template<typename T>
struct _vec_element {
    using type = T;
};

template<>
struct _vec_element<bool> {
    using type = uint8_t;
};

template<typename T>
using vec_element_t = typename _vec_element<T>::type;

template<typename T, size_t N>
struct _vec_type {
    typedef vec_element_t<T> type __attribute__((vector_size(sizeof(vec_element_t<T>) * N)));
};

template<typename T, size_t N>
using vec_type_t = typename _vec_type<T, N>::type;

template<typename T, size_t N>
struct api_type {
    using type = vec_type_t<T, N>;
};

#if NACS_CPU_X86_64 || NACS_CPU_X86
// Needed on x86 to match the API convention used by sleef.
template<>
struct api_type<int,2> {
    using type = vec_type_t<int, 4>;
};
template<size_t N>
struct api_type<std::enable_if_t<N < 16, bool>,N> {
    using type = vec_type_t<uint8_t, 16>;
};
#endif

template<typename T, size_t N>
using api_type_t = typename api_type<T, N>::type;

template<typename T, size_t N, size_t ...I>
vec_type_t<T, N> _broadcast(std::index_sequence<I...>, T v)
{
    return vec_type_t<T, N>{(I, v)...};
}

template<typename T, size_t N>
vec_type_t<T, N> broadcast(T v)
{
    return _broadcast(std::make_index_sequence<N>{}, v);
}

template<typename FT, bool approx>
struct VecCodegenTest;

template<typename Res, typename... Args, bool approx>
struct VecCodegenTest<Res(Args...), approx> : CodegenTest<Res(Args...), approx> {
    using FT = Res(Args...);
    using SuperT = CodegenTest<FT, approx>;
    using SuperT::CodegenTest;
    using SuperT::foreach_args;
    using SuperT::get_llvm_test;
    using SuperT::test_call;
    using SuperT::test_res;

    void test()
    {
        SuperT::test();
        _test_vec_ret<2>();
        _test_vec_ret<4>();
        _test_vec_ret<8>();
        _test_vec_ret_ref<2>();
        _test_vec_ret_ref<4>();
        _test_vec_ret_ref<8>();
    }
private:
    template<size_t vsize>
    void _test_vec_ret()
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        auto test_vec = get_llvm_test(vec);
        if (!test_vec.f) {
            std::cerr << "Vectorizing return failed." << std::endl;
            abort();
        }
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
#if !REF_VEC_ONLY
        auto fptr = test_vec.get_ptr();
        foreach_args(
            [&] (Args... args) {
                api_type_t<Res, vsize> vres;
                vec_runner<vsize>::run(
                    [&] {
                        auto f = (api_type_t<Res, vsize>(*)(Args...) NACS_VECTORCALL)fptr;
                        vres = f(args...);
                    });
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector return", vres[i], args...);
                }
            });
#endif
    }
    template<size_t vsize>
    void _test_vec_ret_ref()
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        vec.add_ret_ref();
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f) {
            std::cerr << "Vectorizing byref return failed." << std::endl;
            abort();
        }
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
        auto f = (void(*)(Res*, Args...))test_vec.get_ptr();
        foreach_args(
            [&] (Args... args) {
                Res ares[vsize] __attribute__((aligned(64)));
                f(ares, args...);
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector ref return", ares[i], args...);
                }
            });
    }
};

template<typename FT, bool approx=false>
using TestVec = MkTest<VecCodegenTest, FT, approx>;

int main()
{
    TestCtx ctx;

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        TestVec<double(double)>([] (double v) { return v; },
                                {1.2, 4.2}, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Float64, {});
        builder.createRet(builder.getConstFloat(1.1));
        test_str_eq(builder.get(), "Float64 () {\n"
                    "L0:\n"
                    "  ret Float64 1.1\n"
                    "}");
        TestVec<double()>([] { return 1.1; }, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Bool, {});
        builder.createRet(IR::Consts::False);
        test_str_eq(builder.get(), "Bool () {\n"
                    "L0:\n"
                    "  ret Bool false\n"
                    "}");
        TestVec<bool()>([] { return false; }, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Int32, {});
        builder.createRet(builder.getConstInt(42));
        test_str_eq(builder.get(), "Int32 () {\n"
                    "L0:\n"
                    "  ret Int32 42\n"
                    "}");
        TestVec<int()>([] { return 42; }, ctx, builder.get());
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
        TestVec<double(bool, double)>([] (bool b, double v) {
                                          return b ? v : 3.4;
                                      }, {false, true}, {1.3, 4.5}, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        auto val1 = builder.createAdd(builder.getConstFloat(3.4), 0);
        builder.createRet(val1);
        test_str_eq(builder.get(), "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  Float64 %1 = add Float64 3.4, Float64 %0\n"
                    "  ret Float64 %1\n"
                    "}");
        TestVec<double(double)>([] (double v) { return 3.4 + v; },
                                {-1.71, 2.3, 1.3}, ctx, builder.get());
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
        TestVec<double(double, double), true>(
            [] (double v0, double v1) {
                auto v2 = 3.4 + v0;
                return v2 - v2 * v1;
            }, {-1.71, 2.3, 1.3}, {-1.71, 2.3, 1.3}, ctx, builder.get());
    }

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
        TestVec<double(int, int)>(
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
        TestVec<double(int, double)>([] (int i, double v) { return i > v ? v : i; },
                                     {20, -10}, {1.3, 5.6}, ctx, builder.get());
    }

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
        TestVec<double(int), true>(
            [] (int v) { return sin(v) + sin(2 * v); }, {1, 2}, ctx, builder.get());
    }

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        const double data[] = {0.0, 0.1, 0.2, 0.6};
        auto val1 = builder.createInterp(0, 2, 3, sizeof(data) / sizeof(double), data);
        builder.createRet(val1);
        test_str_eq(builder.get(), "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  Float64 %1 = interp [2, (4) +3] (Float64 %0) {0, 0.1, 0.2, 0.6}\n"
                    "  ret Float64 %1\n"
                    "}");
        TestVec<double(double), true>(
            [&] (double x) { return linearInterpolate(x, 2, 3, 4, data); },
            {0.7, 2.3, 3.5, 4.4, 5.5}, ctx, builder.get());
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
        TestVec<double(double, double), true>(
            [] (double t, double len) {
                double v2 = (3.4 + t);
                return v2 - v2 * len;
            }, {2.3, 1.3}, {1.3, 10, 1.0}, ctx, builder.get());
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
        TestVec<double(bool, int, double)>(
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
        TestVec<int(double)>([] (double v) { return (int)v; },
                             {2.3, 2.9, 10}, ctx, builder.get());
    }

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
        TestVec<double(bool, int, double), true>(
            [] (bool b, int i, double v) {
                return b ? cos(i) : sin(v);
            }, {true, false}, {1, 2}, {2.3, 3.8}, ctx, builder.get());
    }

    return 0;
}
