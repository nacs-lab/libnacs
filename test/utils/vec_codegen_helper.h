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

#ifndef __NACS_TEST_VEC_CODEGEN_HELPER_H__
#define __NACS_TEST_VEC_CODEGEN_HELPER_H__

#include "codegen_helper.h"

#include "../../lib/nacs-utils/processor.h"

#include <tuple>

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

template<size_t ...vargi>
static inline constexpr bool num_in_pack(size_t num)
{
    constexpr size_t idxs[] = {vargi...};
    for (size_t i = 0; i < sizeof...(vargi); i++) {
        if (idxs[i] == num) {
            return true;
        }
    }
    return false;
}

template<typename Seq, typename F, size_t vsize, size_t ...vargi>
struct _vec_sig;

template<typename Res, typename... Args, size_t vsize, size_t ...vargi,
         size_t ...I>
struct _vec_sig<std::index_sequence<I...>, Res(Args...), vsize, vargi...> {
    using sarg_tuple_t = std::tuple<Args...>;
    template<size_t i>
    using sarg_t = std::tuple_element_t<i, sarg_tuple_t>;
    template<size_t i>
    using varg_t = api_type_t<sarg_t<i>,vsize>;
    template<size_t i>
    using arg_t = std::conditional_t<num_in_pack<vargi...>(i), varg_t<i>, sarg_t<i>>;
    using func = api_type_t<Res,vsize>(arg_t<I>...) NACS_VECTORCALL;
    using func_ra = api_type_t<Res,vsize>(arg_t<I>&...) NACS_VECTORCALL;
    using func_rr = void(Res&,arg_t<I>...) NACS_VECTORCALL;
    using func_rra = void(Res&,sarg_t<I>&...);
};

template<typename T, size_t vsize>
struct aligned_array {
    T ary[vsize] __attribute__((aligned(64)));
};

template<typename F, size_t vsize, size_t ...vargi>
struct vec_sig;

template<typename Res, typename... Args, size_t vsize, size_t ...vargi>
struct vec_sig<Res(Args...), vsize, vargi...> :
    _vec_sig<std::make_index_sequence<sizeof...(Args)>, Res(Args...), vsize, vargi...> {
};

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
    template<size_t ...options>
    void test_vec_arg()
    {
        _test_vec_arg(std::index_sequence<options...>{},
                      std::index_sequence<>{});
    }
    void test_vec_allarg()
    {
        _test_vec_arg(std::make_index_sequence<sizeof...(Args)>{},
                      std::index_sequence<>{});
    }
private:
    template<size_t vsize>
    void _test_vec_ret()
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        auto test_vec = get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing return failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
#if !REF_VEC_ONLY
        auto fptr = test_vec.get_ptr();
        foreach_args(
            [&] (Args... args) {
                api_type_t<Res, vsize> vres{};
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
        LLVM::Codegen::Wrapper vec_unalign{false};
        vec_unalign.vector_size = vsize;
        vec_unalign.add_ret_ref(sizeof(Res));
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing byref return failed.");
        auto test_vec_unalign = this->get_llvm_test(vec_unalign);
        if (!test_vec_unalign.f)
            FAIL("Vectorizing unaligned byref return failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
        auto f = (void(*)(Res*, Args...))test_vec.get_ptr();
        auto f_unalign = (void(*)(Res*, Args...))test_vec_unalign.get_ptr();
        foreach_args(
            [&] (Args... args) {
                Res ares[vsize * 2] __attribute__((aligned(64)));
                memset(ares, 0, sizeof(ares));
                f(ares, args...);
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector ref return", ares[i], args...);
                }
                for (size_t offset = 0; offset < vsize; offset++) {
                    memset(ares, 0, sizeof(ares));
                    f_unalign(ares + offset, args...);
                    for (size_t i = 0; i < vsize; i++) {
                        test_res("Vector ref return", ares[offset + i], args...);
                    }
                }
            });
    }
    template<size_t i, size_t vsize, size_t ...vargi, typename F, typename... Args2>
    std::enable_if_t<i == 0> _foreach_vecarg(F &&f, Args2&&... args2)
    {
        f(std::forward<Args2>(args2)...);
    }
    template<size_t i, size_t vsize, size_t ...vargi, typename F, typename... Args2>
    std::enable_if_t<(i > 0) && !num_in_pack<vargi...>(i - 1)>
    _foreach_vecarg(F &&f, Args2&&... args2)
    {
        auto &_arg = std::get<i - 1>(this->args);
        for (auto v: _arg) {
            _foreach_vecarg<i - 1, vsize, vargi...>(std::forward<F>(f), v,
                                                    std::forward<Args2>(args2)...);
        }
    }
    template<size_t i, size_t vsize, size_t ...vargi, typename F, typename... Args2>
    std::enable_if_t<(i > 0) && num_in_pack<vargi...>(i - 1)>
    _foreach_vecarg(F &&f, Args2&&... args2)
    {
        using ele_t = typename vec_sig<Res(Args...), vsize, vargi...>::template sarg_t<i - 1>;
        auto &_arg = std::get<i - 1>(this->args);
        auto arg_sz = _arg.size();
        for (size_t j = 0; j < arg_sz; j++) {
            api_type_t<ele_t, vsize> v;
            for (size_t k = 0; k < vsize; k++)
                v[k] = _arg[(k + j) % arg_sz];
            _foreach_vecarg<i - 1, vsize, vargi...>(std::forward<F>(f), v,
                                                    std::forward<Args2>(args2)...);
        }
    }
    template<size_t vsize, size_t ...vargi, typename F>
    void foreach_vecarg(F &&f)
    {
        _foreach_vecarg<sizeof...(Args), vsize, vargi...>(std::forward<F>(f));
    }
    template<size_t I, size_t ...vargi, typename Arg>
    std::enable_if_t<!num_in_pack<vargi...>(I), Arg> get_sarg(Arg arg, size_t)
    {
        return arg;
    }
    template<size_t I, size_t ...vargi, typename Arg>
    auto get_sarg(Arg arg, size_t i)
        -> std::enable_if_t<num_in_pack<vargi...>(I),
                            std::remove_reference_t<decltype(arg[i])>>
    {
        return arg[i];
    }
    template<size_t varg0, size_t ...vargi>
    void _test_vec_arg(std::index_sequence<>, std::index_sequence<varg0, vargi...>)
    {
        _test_vec_arg2<varg0, vargi...>();
    }
    void _test_vec_arg(std::index_sequence<>, std::index_sequence<>)
    {
    }
    template<size_t option0, size_t ...options, size_t ...vargi>
    void _test_vec_arg(std::index_sequence<option0, options...>,
                       std::index_sequence<vargi...>)
    {
        _test_vec_arg(std::index_sequence<options...>{},
                      std::index_sequence<vargi...>{});
        _test_vec_arg(std::index_sequence<options...>{},
                      std::index_sequence<option0, vargi...>{});
    }
    template<size_t ...vargi>
    void _test_vec_arg2()
    {
        _test_vec_arg1<2, vargi...>(std::make_index_sequence<sizeof...(Args)>{});
        _test_vec_arg1<4, vargi...>(std::make_index_sequence<sizeof...(Args)>{});
        _test_vec_arg1<8, vargi...>(std::make_index_sequence<sizeof...(Args)>{});
    }
    template<size_t vsize, size_t ...vargi, size_t ...I>
    void _test_vec_arg1(std::index_sequence<I...> x)
    {
        _test_vec_arg_byval<vsize, vargi...>(x);
        _test_vec_arg_refarg<vsize, vargi...>(x);
        _test_vec_arg_refret<vsize, vargi...>(x);
        _test_vec_arg_refall<vsize, vargi...>(x);
    }
    template<size_t vsize, size_t ...vargi, size_t ...I>
    void _test_vec_arg_byval(std::index_sequence<I...>)
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        int l[] = {(vec.add_vector(vargi), 0)...};
        (void)l;
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing argument failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
#if !REF_VEC_ONLY
        using vs = vec_sig<Res(Args...), vsize, vargi...>;
        auto fptr = test_vec.get_ptr();
        foreach_vecarg<vsize, vargi...>(
            [&] (typename vs::template arg_t<I>... args) {
                api_type_t<Res, vsize> vres{};
                vec_runner<vsize>::run(
                    [&] {
                        auto f = (typename vs::func*)fptr;
                        vres = f(args...);
                    });
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector arguments", vres[i],
                             get_sarg<I, vargi...>(args, i)...);
                }
            });
#endif
    }
    template<size_t vsize, size_t ...vargi, size_t ...I>
    void _test_vec_arg_refarg(std::index_sequence<I...>)
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        int l[] = {(vec.add_vector(vargi), 0)...};
        (void)l;
        for (size_t i = 0; i < sizeof...(Args); i++)
            vec.add_byref(i);
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing argument failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
#if !REF_VEC_ONLY
        using vs = vec_sig<Res(Args...), vsize, vargi...>;
        auto fptr = test_vec.get_ptr();
        foreach_vecarg<vsize, vargi...>(
            [&] (typename vs::template arg_t<I>... args) {
                api_type_t<Res, vsize> vres{};
                vec_runner<vsize>::run(
                    [&] {
                        auto f = (typename vs::func_ra*)fptr;
                        vres = f(args...);
                    });
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector arguments byref", vres[i],
                             get_sarg<I, vargi...>(args, i)...);
                }
            });
#endif
    }
    template<size_t vsize, size_t ...vargi, size_t ...I>
    void _test_vec_arg_refret(std::index_sequence<I...>)
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        LLVM::Codegen::Wrapper vec_unalign{false};
        vec_unalign.vector_size = vsize;
        int l[] = {(vec.add_vector(vargi), 0)..., (vec_unalign.add_vector(vargi), 0)...};
        (void)l;
        vec.add_ret_ref();
        vec_unalign.add_ret_ref(sizeof(Res));
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing argument with byref return failed.");
        auto test_vec_unalign = this->get_llvm_test(vec_unalign);
        if (!test_vec_unalign.f)
            FAIL("Vectorizing argument with unaligned byref return failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
#if !REF_VEC_ONLY
        using vs = vec_sig<Res(Args...), vsize, vargi...>;
        auto fptr = test_vec.get_ptr();
        auto fptr_unalign = test_vec_unalign.get_ptr();
        foreach_vecarg<vsize, vargi...>(
            [&] (typename vs::template arg_t<I>... args) {
                Res ares[vsize * 2] __attribute__((aligned(64)));
                memset(ares, 0, sizeof(ares));
                vec_runner<vsize>::run(
                    [&] {
                        auto f = (typename vs::func_rr*)fptr;
                        f(*ares, args...);
                    });
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector arguments ret ref", ares[i],
                             get_sarg<I, vargi...>(args, i)...);
                }
                for (size_t offset = 0; offset < vsize; offset++) {
                    memset(ares, 0, sizeof(ares));
                    vec_runner<vsize>::run(
                        [&] {
                            auto f_unalign = (typename vs::func_rr*)fptr_unalign;
                            f_unalign(ares[offset], args...);
                        });
                    for (size_t i = 0; i < vsize; i++) {
                        test_res("Vector arguments ret ref", ares[offset + i],
                                 get_sarg<I, vargi...>(args, i)...);
                    }
                }
            });
#endif
    }
    template<size_t vsize, size_t ...vargi, size_t ...I>
    void _test_vec_arg_refall(std::index_sequence<I...>)
    {
        LLVM::Codegen::Wrapper vec{false};
        vec.vector_size = vsize;
        LLVM::Codegen::Wrapper vec_unalign{false};
        vec_unalign.vector_size = vsize;
        int l[] = {(vec.add_vector(vargi), 0)..., (vec_unalign.add_vector(vargi), 0)...,
            (vec.add_byref(I), 0)..., (vec_unalign.add_byref(I, sizeof(Args)), 0)..., };
        (void)l;
        // for (size_t i = 0; i < sizeof...(Args); i++)
        //     vec.add_byref(i);
        vec.add_ret_ref();
        // for (size_t i = 0; i < sizeof...(Args); i++)
        //     vec_unalign.add_byref(i);
        vec_unalign.add_ret_ref(sizeof(Res));
        auto test_vec = this->get_llvm_test(vec);
        if (!test_vec.f)
            FAIL("Vectorizing argument failed.");
        auto test_vec_unalign = this->get_llvm_test(vec_unalign);
        if (!test_vec_unalign.f)
            FAIL("Vectorizing argument failed.");
        if (vsize * 8 > (unsigned)host_info.get_vector_size())
            return;
        using vs = vec_sig<Res(Args...), vsize, vargi...>;
        auto fptr = test_vec.get_ptr();
        auto fptr_unalign = test_vec_unalign.get_ptr();
        foreach_vecarg<vsize, vargi...>(
            [&] (typename vs::template arg_t<I>... args) {
                std::tuple<aligned_array<Args, vsize * 2>...> aargs;
                Res ares[vsize * 2] __attribute__((aligned(64)));
                int init1[] = {(memcpy(std::get<I>(aargs).ary, &args, sizeof(args)), 0)...};
                (void)init1;
                memset(ares, 0, sizeof(ares));
                vec_runner<vsize>::run(
                    [&] {
                        auto f = (typename vs::func_rra*)fptr;
                        f(*ares, *std::get<I>(aargs).ary...);
                    });
                for (size_t i = 0; i < vsize; i++) {
                    test_res("Vector arguments ref all", ares[i],
                             get_sarg<I, vargi...>(args, i)...);
                }
                for (size_t offset = 0; offset < vsize; offset++) {
                    int init2[] = {(memcpy(std::get<I>(aargs).ary + offset,
                                           &args, sizeof(args)), 0)...};
                    (void)init2;
                    memset(ares, 0, sizeof(ares));
                    vec_runner<vsize>::run(
                        [&] {
                            auto f_unalign = (typename vs::func_rra*)fptr_unalign;
                            f_unalign(ares[offset], std::get<I>(aargs).ary[offset]...);
                        });
                    for (size_t i = 0; i < vsize; i++) {
                        test_res("Vector arguments ref all", ares[offset + i],
                                 get_sarg<I, vargi...>(args, i)...);
                    }
                }
            });
    }
};

template<typename FT, bool approx=false>
using TestVec = MkTest<VecCodegenTest, FT, approx>;

#endif
