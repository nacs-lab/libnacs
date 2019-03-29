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

#ifndef __NACS_TEST_CODEGEN_HELPER_H__
#define __NACS_TEST_CODEGEN_HELPER_H__

#include "winpath_helper.h"

// This avoid linking to LLVM
// Without this, LLVM header generates a global variable with
// weak linkage that links to symbol in libLLVM.
#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/compile.h"
#include "../../lib/utils/llvm/execute.h"
#include "../../lib/utils/llvm/utils.h"

#include <vector>
#include <fstream>

namespace {

using namespace NaCs;

struct LLVMTest {
    LLVMTest(llvm::LLVMContext &ll_ctx, LLVM::Exe::Engine &engine, const IR::Function &func)
        : mod(LLVM::new_module("", ll_ctx)),
          engine(engine)
    {
        LLVM::Codegen::Context ctx(mod);
        f = ctx.emit_function(func, "0");
        auto fty = f->getFunctionType();
        assert(!fty->isVarArg());
        for (auto argt: fty->params()) {
            assert(argt == ctx.T_i8 || argt == ctx.T_i32 || argt == ctx.T_f64);
        }
    }
    LLVMTest(llvm::LLVMContext &ll_ctx, LLVM::Exe::Engine &engine, const IR::Function &func,
             const LLVM::Codegen::Wrapper &wrapper)
        : mod(LLVM::new_module("", ll_ctx)),
          engine(engine)
    {
        LLVM::Codegen::Context ctx(mod);
        auto f0 = ctx.emit_function(func, "1", false);
        auto fty0 = f0->getFunctionType();
        assert(!fty0->isVarArg());
        for (auto argt: fty0->params())
            assert(argt == ctx.T_i8 || argt == ctx.T_i32 || argt == ctx.T_f64);
        f = ctx.emit_wrapper(f0, "0", wrapper);
    }
    LLVMTest(LLVMTest&&) = default;
    LLVMTest &operator=(LLVMTest&&) = default;
    LLVMTest &opt()
    {
        LLVM::Compile::optimize(mod);
        return *this;
    }
    LLVMTest &print()
    {
        LLVM::dump(mod);
        return *this;
    }
    void *get_ptr(const char *name=nullptr)
    {
        assert(f);
        llvm::SmallVector<char,0> vec;
        auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(),
                                               mod);
        assert(res);
        if (name && *name) {
            std::fstream stm(name, std::ios_base::out);
            if (stm) {
                stm.write(&vec[0], vec.size());
            }
        }
        auto obj_id = engine.load(&vec[0], vec.size());
        obj_ids.push_back(obj_id);
        assert(obj_id);
        return engine.get_symbol("0");
    }
    ~LLVMTest()
    {
        for (auto id: obj_ids)
            engine.free(id);
        LLVM::delete_module(mod);
    }
    llvm::Module *mod;
    LLVM::Exe::Engine &engine;
    llvm::Function *f = nullptr;
    std::vector<uint64_t> obj_ids;
};

struct TestCtx {
    LLVM::Exe::Engine engine;
    std::unique_ptr<llvm::LLVMContext,void(*)(llvm::LLVMContext*)> llvm;
    std::unique_ptr<IR::ExeContext> ir;
    TestCtx()
        : llvm(LLVM::new_context(), LLVM::delete_context),
          ir(IR::ExeContext::get())
    {
    }
    template<typename... Args>
    LLVMTest get_llvm_test(Args... args)
    {
        return LLVMTest(*llvm, engine, args...);
    }
};

namespace detail {

static bool approx(double a, double b)
{
    double diff = abs(a - b);
    double avg = abs(a + b) / 2;
    return diff < 2e-8 || diff / avg < 2e-8;
}

template<typename T1, typename T2>
static bool compare(T1 a, T2 b, bool _approx)
{
    if (_approx)
        return approx(a, b);
    return a == b;
}

template<typename T, size_t... Is>
static void _print_tuple(T &&t, std::index_sequence<Is...>)
{
    auto f = [] (auto &&v, size_t i) {
                 if (i != 0)
                     std::cerr << ", ";
                 std::cerr << v;
                 return 0;
             };
    (void)f;
    int l[] = { f(std::get<Is>(t), Is)... };
    (void)l;
}

template<typename T>
static void print_tuple(T &&t)
{
    _print_tuple(std::forward<T>(t), std::make_index_sequence<
                 std::tuple_size<std::decay_t<T>>::value>{});
}

template<bool _approx, typename T1, typename T2, typename Tuple>
static void _test_res(T1 res, T2 exp, const char *name, IR::Function &func,
                       Tuple &&args)
{
    if (detail::compare(res, exp, _approx))
        return;
    std::cerr << name << " test failed on test case: (";
    detail::print_tuple(args);
    std::cerr << ")" << std::endl;
    std::cerr << "Expected: " << (_approx ? "≈ " : "") << exp << std::endl;
    std::cerr << "Got: " << res << std::endl;
    std::cerr << "Code: " << func << std::endl;
    abort();
}

template<class Tuple, std::size_t... I>
static inline constexpr decltype(auto)
tie_tuple_impl(Tuple &t, std::index_sequence<I...>)
{
    return std::tie(std::get<I>(t)...);
}

template<class Tuple>
static inline constexpr decltype(auto) tie_tuple(Tuple &t)
{
    return tie_tuple_impl(t, std::make_index_sequence<
                          std::tuple_size<std::decay_t<Tuple>>::value>{});
}

}

template<typename FT>
struct TestSet;

template<typename Res, typename... Args>
struct TestSet<Res(Args...)> {
    Res res;
    std::tuple<Args...> args;
    template<typename Res2, typename... Args2>
    TestSet(Res2 res, Args2 ...args)
        : res(res),
          args(args...)
    {}
    template<typename Res2, typename... Args2>
    TestSet(const TestSet<Res2(Args2...)> &case2)
        : res(case2.res),
          args(case2.args)
    {}
};

template<typename Res, typename... Args>
static TestSet<Res(Args...)> test_set(Res res, Args ...args)
{
    return TestSet<Res(Args...)>(res, args...);
}

template<bool approx=false, typename T1, typename Res, typename... Args>
static void test_res(Res exp, const char *name, IR::Function &func,
                     T1 res, Args... args)
{
    detail::_test_res<approx>(res, exp, name, func, std::make_tuple(args...));
}

template<typename FT, bool approx>
struct CodegenTest;

template<typename Res, typename... Args, bool approx>
struct CodegenTest<Res(Args...), approx> {
    using FT = Res(Args...);
    TestCtx &ctx;
    IR::Function &func;
    std::vector<TestSet<FT>> cases;
    std::unique_ptr<LLVMTest> llvm_test;

    template<typename... Cases>
    CodegenTest(TestCtx &ctx, IR::Function &func, Cases ..._cases)
        : ctx(ctx),
          func(func),
          cases{_cases...}
    {
    }

    void test()
    {
        auto f_interp = get_interp_func();
        for (auto &_case: cases) {
            _test(f_interp, _case, "Interpreted");
        }

        auto f_comp = get_llvm_func();
        for (auto &_case: cases) {
            _test(f_comp, _case, "Compiled");
        }
        test_ref();
    }

    template<typename T1, typename T2, typename... Args2>
    void test_res(T2 exp, const char *name, T1 res, Args2... args)
    {
        detail::_test_res<approx>(res, exp, name, func, std::make_tuple(args...));
    }

    auto get_interp_func() -> decltype(ctx.ir->getFunc<FT>(func))
    {
        return ctx.ir->getFunc<FT>(func);
    }
    FT *get_llvm_func(const char *name=nullptr)
    {
        return (FT*)_get_llvm_test().get_ptr(name);
    }
    template<typename... Args2>
    LLVMTest get_llvm_test(Args2&&... args)
    {
        return ctx.get_llvm_test(func, std::forward<Args2>(args)...);
    }
private:
    void test_ref()
    {
        LLVM::Codegen::Wrapper ret_ref{false};
        ret_ref.add_ret_ref();
        auto test_ret_ref = get_llvm_test(ret_ref);
        auto f_ret_ref = (void(*)(Res&, Args...))test_ret_ref.get_ptr();
        for (auto &_case: cases) {
            Res res;
            auto args = std::tuple_cat(std::tie(res), _case.args);
            applyTuple(f_ret_ref, args);
            detail::_test_res<approx>(res, _case.res, "Ref return", func, _case.args);
        }
        constexpr uint32_t nargs = sizeof...(Args);
        if (!nargs)
            return;
        LLVM::Codegen::Wrapper ref{false};
        for (uint32_t i = 0; i < nargs; i++)
            ref.add_byref(i);
        auto test_ref = get_llvm_test(ref);
        auto f_ref = (Res(*)(const Args&...))test_ref.get_ptr();
        for (auto &_case: cases)
            _test(f_ref, _case, "By ref");
        ref.add_ret_ref();
        auto test_ref2 = get_llvm_test(ref);
        auto f_ref2 = (void(*)(Res&,const Args&...))test_ref2.get_ptr();
        for (auto &_case: cases) {
            Res res;
            auto args = std::tuple_cat(std::tie(res), _case.args);
            applyTuple(f_ref2, args);
            detail::_test_res<approx>(res, _case.res, "Ref return and params", func, _case.args);

            test_ref_alias(f_ref2, _case, std::make_index_sequence<nargs>());
        }
    }
    template<size_t... I, typename F>
    void test_ref_alias(F &&f, const TestSet<FT> &tc, std::index_sequence<I...>)
    {
        int l[] = {(test_ref_alias_i<I>(f, tc), 0)...};
        (void)l;
    }
    template<size_t i, typename F>
    std::enable_if_t<
        std::is_same<std::tuple_element_t<i, std::tuple<Args...>>,
                     Res>::value> test_ref_alias_i(F &&f, const TestSet<FT> &tc)
    {
        auto args = tc.args;
        auto callargs = std::tuple_cat(std::tie(std::get<i>(args)), detail::tie_tuple(args));
        applyTuple(f, callargs);
        detail::_test_res<approx>(std::get<i>(args), tc.res,
                                  "Alias return", func, tc.args);
    }
    template<size_t i, typename F>
    std::enable_if_t<
        !std::is_same<std::tuple_element_t<i, std::tuple<Args...>>,
                      Res>::value> test_ref_alias_i(F&&, const TestSet<FT>&)
    {
    }
    LLVMTest &_get_llvm_test()
    {
        if (!llvm_test)
            llvm_test.reset(new LLVMTest(ctx.get_llvm_test(func)));
        return *llvm_test;
    }
    template<typename F>
    void _test(F &&f, const TestSet<FT> &tc, const char *name)
    {
        auto res = applyTuple(std::forward<F>(f), tc.args);
        detail::_test_res<approx>(res, tc.res, name, func, tc.args);
    }
};

template<typename FT, bool approx=false, typename... Cases>
static CodegenTest<FT, approx> test_codegen(TestCtx &ctx, IR::Function &func, Cases ...cases)
{
    CodegenTest<FT, approx> test{ctx, func, cases...};
    test.test();
    return test;
}

}

#endif
