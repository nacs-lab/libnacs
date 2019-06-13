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

#include "ir_helper.h"

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
        LLVM::Codegen::Context::data_map_t data_map;
        f = ctx.emit_function(func, "0", &data_map);
        populate_data(func, data_map);
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
        LLVM::Codegen::Context::data_map_t data_map;
        auto f0 = ctx.emit_function(func, "1", false, &data_map);
        populate_data(func, data_map);
        auto fty0 = f0->getFunctionType();
        assert(!fty0->isVarArg());
        for (auto argt: fty0->params())
            assert(argt == ctx.T_i8 || argt == ctx.T_i32 || argt == ctx.T_f64);
        f = ctx.emit_wrapper(f0, "0", wrapper);
    }
    LLVMTest(LLVMTest&&) = default;
    LLVMTest &operator=(LLVMTest&&) = delete;
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
        auto cb = [&] (auto &name) {
            auto it = data.find(name);
            if (it == data.end())
                return (uintptr_t)0;
            return (uintptr_t)&it->second[0];
        };
        auto obj_id = engine.load(&vec[0], vec.size(), cb);
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
    std::map<std::string,std::vector<double>> data;
private:
    void populate_data(const IR::Function &func,
                       LLVM::Codegen::Context::data_map_t &data_map)
    {
        for (auto &it: data_map) {
            auto ptr = &func.float_table[it.second.first];
            data[it.first] = std::vector<double>(ptr, ptr + it.second.second);
        }
    }
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

template<typename FT, bool approx>
struct CodegenTest;

template<typename Res, typename... Args, bool approx>
struct CodegenTest<Res(Args...), approx> : IRTest<Res(Args...), approx> {
    using FT = Res(Args...);
    using SuperT = IRTest<FT, approx>;
    using SuperT::foreach_args;
    using SuperT::test_call;
    using SuperT::test_res;
    using SuperT::testeach_args;

    TestCtx &ctx;
    std::unique_ptr<LLVMTest> llvm_test;

    CodegenTest(std::function<FT> fexp, std::vector<Args> ...args,
                TestCtx &ctx, IR::Function &func)
        : SuperT(fexp, std::move(args)..., *ctx.ir, func),
          ctx(ctx)
    {
    }

    void test()
    {
        SuperT::test();
        testeach_args(get_llvm_func(), "Compiled");
        test_ref();
    }

    FT *get_llvm_func(const char *name=nullptr)
    {
        return (FT*)_get_llvm_test().get_ptr(name);
    }
    template<typename... Args2>
    LLVMTest get_llvm_test(Args2&&... args)
    {
        return ctx.get_llvm_test(this->func, std::forward<Args2>(args)...);
    }
private:
    void test_ref()
    {
        LLVM::Codegen::Wrapper ret_ref{false};
        ret_ref.add_ret_ref();
        auto test_ret_ref = get_llvm_test(ret_ref);
        auto f_ret_ref = (void(*)(Res&, Args...))test_ret_ref.get_ptr();
        testeach_args([&] (Args ...args) {
                          Res res;
                          f_ret_ref(res, args...);
                          return res;
                      }, "Ref return");
        constexpr uint32_t nargs = sizeof...(Args);
        if (!nargs)
            return;
        LLVM::Codegen::Wrapper ref{false};
        for (uint32_t i = 0; i < nargs; i++)
            ref.add_byref(i);
        auto test_ref = get_llvm_test(ref);
        testeach_args((Res(*)(const Args&...))test_ref.get_ptr(), "By ref");
        ref.add_ret_ref();
        auto test_ref2 = get_llvm_test(ref);
        auto f_ref2 = (void(*)(Res&,const Args&...))test_ref2.get_ptr();
        foreach_args([&] (Args ...args) {
                         Res res;
                         f_ref2(res, args...);
                         test_res("Ref return and params", res, args...);
                         test_ref_alias(f_ref2, std::make_index_sequence<nargs>(), args...);
                     });
    }
    template<size_t... I, typename F>
    void test_ref_alias(F &&f, std::index_sequence<I...>, Args... args)
    {
        int l[] = {(test_ref_alias_i<I>(f, args...), 0)...};
        (void)l;
    }
    template<size_t i, typename F>
    void _test_ref_alias_i(F &&f, Args... args, Args... args2)
    {
        auto &res = std::get<i>(std::tuple<Args&...>(args...));
        f(res, args...);
        test_res("Alias return", res, args2...);
    }
    template<size_t i, typename F>
    std::enable_if_t<
        std::is_same<std::tuple_element_t<i, std::tuple<Args...>>,
                     Res>::value> test_ref_alias_i(F &&f, Args... args)
    {
        // Copy arguments since we'll override them.
        _test_ref_alias_i<i>(f, args..., args...);
    }
    template<size_t i, typename F>
    std::enable_if_t<
        !std::is_same<std::tuple_element_t<i, std::tuple<Args...>>,
                      Res>::value> test_ref_alias_i(F&&, Args...)
    {
    }
    LLVMTest &_get_llvm_test()
    {
        if (!llvm_test)
            llvm_test.reset(new LLVMTest(ctx.get_llvm_test(this->func)));
        return *llvm_test;
    }
};

template<typename FT, bool approx=false>
using TestCodegen = MkTest<CodegenTest, FT, approx>;

}

#endif
