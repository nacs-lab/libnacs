/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include <nacs-utils/llvm/codegen.h>
#include <nacs-utils/llvm/compile.h>
#include <nacs-utils/llvm/execute.h>
#include <nacs-utils/llvm/utils.h>
#include <nacs-utils/number.h>
#include <nacs-utils/timer.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <math.h>

using namespace NaCs;

template<typename T>
static std::string sprint(T &&v)
{
    std::stringstream stm;
    stm << std::forward<T>(v);
    return stm.str();
}

template<typename T>
static void test_str_eq(T &&v, std::string val)
{
    auto str = sprint(std::forward<T>(v));
    if (val == str)
        return;
    std::cerr << "Test failed:" << std::endl;
    std::cerr << "  Expect: \"" << val << "\"" << std::endl;
    std::cerr << "  Got: \"" << str << "\"" << std::endl;
    abort();
}

struct LLVMTest {
    LLVMTest(llvm::LLVMContext &ll_ctx, LLVM::Exe::Engine &engine, const IR::Function &func)
        : mod(new llvm::Module("", ll_ctx)),
          engine(engine)
    {
        LLVM::Codegen::Context ctx(mod.get());
        f = ctx.emit_function(func, 0);
    }
    LLVMTest(LLVMTest&&) = default;
    LLVMTest &operator=(LLVMTest&&) = default;
    LLVMTest &opt()
    {
        LLVM::Compile::optimize(mod.get());
        return *this;
    }
    LLVMTest &print()
    {
        LLVM::dump(mod.get());
        return *this;
    }
    void *get_ptr()
    {
        llvm::SmallVector<char,0> vec;
        auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(),
                                               mod.get());
        assert(res);
        auto res2 = engine.load(&vec[0], vec.size());
        assert(res2);
        return engine.get_symbol("0");
    }
    std::unique_ptr<llvm::Module> mod;
    LLVM::Exe::Engine &engine;
    llvm::Function *f = nullptr;
};

int main()
{
    LLVM::Exe::Engine engine;
    llvm::LLVMContext llvm_ctx;

    auto gettest = [&] (const IR::Function &func) {
        return LLVMTest(llvm_ctx, engine, func);
    };

    auto exectx = IR::ExeContext::get();

    {
        IR::Builder builder(IR::Type::Float64, {IR::Type::Float64});
        builder.createRet(0);
        test_str_eq(builder.get(),
                    "Float64 (Float64 %0) {\n"
                    "L0:\n"
                    "  ret Float64 %0\n"
                    "}");
        auto f = exectx->getFunc<double(double)>(builder.get());
        assert(f(1.2) == 1.2);
        assert(f(4.2) == 4.2);
        auto f2 = (double(*)(double))gettest(builder.get()).get_ptr();
        assert(f2(1.2) == 1.2);
        assert(f2(4.2) == 4.2);
    }

    {
        IR::Builder builder(IR::Type::Bool, {});
        builder.createRet(IR::Consts::False);
        test_str_eq(builder.get(), "Bool () {\n"
                    "L0:\n"
                    "  ret Bool false\n"
                    "}");
        auto f = exectx->getFunc<bool()>(builder.get());
        assert(f() == false);
        auto f2 = (bool(*)())gettest(builder.get()).get_ptr();
        assert(f2() == false);
    }

    {
        IR::Builder builder(IR::Type::Float64, {});
        builder.createRet(builder.getConstFloat(1.1));
        test_str_eq(builder.get(), "Float64 () {\n"
                    "L0:\n"
                    "  ret Float64 1.1\n"
                    "}");
        auto f = exectx->getFunc<double()>(builder.get());
        assert(f() == 1.1);
        auto f2 = (double(*)())gettest(builder.get()).get_ptr();
        assert(f2() == 1.1);
    }

    {
        IR::Builder builder(IR::Type::Int32, {});
        builder.createRet(builder.getConstInt(42));
        test_str_eq(builder.get(), "Int32 () {\n"
                    "L0:\n"
                    "  ret Int32 42\n"
                    "}");
        auto f = exectx->getFunc<int()>(builder.get());
        assert(f() == 42);
        auto f2 = (int(*)())gettest(builder.get()).get_ptr();
        assert(f2() == 42);
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
        auto f = exectx->getFunc<double(bool, double)>(builder.get());
        assert(f(true, 1.3) == 1.3);
        assert(f(false, 1.3) == 3.4);
        auto f2 = (double(*)(bool, double))gettest(builder.get()).get_ptr();
        assert(f2(true, 1.3) == 1.3);
        assert(f2(false, 1.3) == 3.4);
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
        auto f = exectx->getFunc<double(double, double)>(builder.get());
        assert(f(2.3, 1.3) == -1.71);
        auto f2 = (double(*)(double, double))gettest(builder.get()).get_ptr();
        assert(f2(2.3, 1.3) == -1.71);
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
        auto f = exectx->getFunc<double(int, int)>(builder.get());
        assert(f(3, 2) == 1.5);
        auto f2 = (double(*)(int, int))gettest(builder.get()).get_ptr();
        assert(f2(3, 2) == 1.5);
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
        auto f = exectx->getFunc<double(int, double)>(builder.get());
        assert(f(20, 1.3) == 1.3);
        assert(f(-10, 1.3) == -10);
        auto f2 = (double(*)(int, double))gettest(builder.get()).get_ptr();
        assert(f2(20, 1.3) == 1.3);
        assert(f2(-10, 1.3) == -10);
    }

    {
        IR::Builder builder(IR::Type::Int32,
                            {IR::Type::Int32, IR::Type::Int32});
        auto loop_bb = builder.newBB();
        auto ret_bb = builder.newBB();
        builder.createBr(loop_bb);
        builder.curBB() = loop_bb;
        auto i = builder.createPhi(IR::Type::Int32, 2);
        auto s = builder.createPhi(IR::Type::Int32, 2);
        builder.addPhiInput(i.second, 0, 0);
        builder.addPhiInput(s.second, 0, builder.getConstInt(0));
        auto i2 = builder.createAdd(i.first, builder.getConstInt(1));
        auto s2 = builder.createAdd(s.first, i.first);
        auto cond = builder.createCmp(IR::CmpType::gt, i2, 1);
        builder.createBr(cond, ret_bb, loop_bb);
        builder.addPhiInput(i.second, loop_bb, i2);
        builder.addPhiInput(s.second, loop_bb, s2);
        builder.curBB() = ret_bb;
        builder.createRet(s2);
        test_str_eq(builder.get(), "Int32 (Int32 %0, Int32 %1) {\n"
                    "L0:\n"
                    "  br L1\n"
                    "L1:\n"
                    "  Int32 %2 = phi [ L0: Int32 %0 ], [ L1: Int32 %4 ]\n"
                    "  Int32 %3 = phi [ L0: Int32 0 ], [ L1: Int32 %5 ]\n"
                    "  Int32 %4 = add Int32 %2, Int32 1\n"
                    "  Int32 %5 = add Int32 %3, Int32 %2\n"
                    "  Bool %6 = cmp gt Int32 %4, Int32 %1\n"
                    "  br Bool %6, L2, L1\n"
                    "L2:\n"
                    "  ret Int32 %5\n"
                    "}");
        auto f = exectx->getFunc<int(int, int)>(builder.get());
        assert(f(1, 3) == 6);
        Timer timer;
        assert(f(2, 1000) == 500499);
        timer.print();
        auto f2 = (int(*)(int, int))gettest(builder.get()).get_ptr();
        assert(f2(1, 3) == 6);
        timer.restart();
        assert(f2(2, 1000) == 500499);
        timer.print();
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
        auto f1 = exectx->getFunc<double(int)>(builder.get());
        assert(f1(1) == sin(1) + sin(2));
        assert(f1(2) == sin(4) + sin(2));

        Timer timer;
        for (int i = 0;i < 1000000;i++)
            f1(1);
        timer.print();
        auto f2 = (double(*)(int))gettest(builder.get()).get_ptr();
        timer.restart();
        for (int i = 0;i < 1000000;i++)
            f2(1);
        timer.print();
    }

    {
        const int32_t data[] = {
            3, 2, 7, 50529027, 197379, 2, 3, 0, 1072693248, 3, 0, 1073741824,
            1, 22, 4, 5, -3, 0, 5, 4, -3, 5, 5, 6, -4, 0, 3, 3, 4, 6, 6, 2, 3,
            -3, 1, 2
        };
        IR::Function newfunc((const uint32_t*)data, sizeof(data) / 4);
        test_str_eq(newfunc, "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %5 = sub Float64 1, Float64 %0\n"
                    "  Float64 %4 = mul Float64 1, Float64 %5\n"
                    "  Float64 %6 = mul Float64 2, Float64 %0\n"
                    "  Float64 %3 = add Float64 %4, Float64 %6\n"
                    "  Float64 %2 = fdiv Float64 %3, Float64 1\n"
                    "  ret Float64 %2\n"
                    "}");
        gettest(newfunc).get_ptr();
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
        auto f = exectx->getFunc<double(double)>(builder.get());
        assert(fabs(f(2.3) - linearInterpolate(2.3, 2, 3, 4, data)) < 1e-10);
        assert(fabs(f(3.5) - linearInterpolate(3.5, 2, 3, 4, data)) < 1e-10);
        assert(fabs(f(4.4) - linearInterpolate(4.4, 2, 3, 4, data)) < 1e-10);
        auto f2 = (double(*)(double))gettest(builder.get()).get_ptr();
        assert(fabs(f2(2.3) - linearInterpolate(2.3, 2, 3, 4, data)) < 1e-10);
        assert(fabs(f2(3.5) - linearInterpolate(3.5, 2, 3, 4, data)) < 1e-10);
        assert(fabs(f2(4.4) - linearInterpolate(4.4, 2, 3, 4, data)) < 1e-10);
    }

    {
        const int32_t data[] = {
            3, 2, 4, 50529027, 1, 3, 0, 1072693248, 1, 13, 10, 3, 0, -3, -3, 0, 4, 3, 2, 3, 1,
            1, 2, 4, 0, 1072693248, 0, 1073741824, 0, 1074003968, 0, 1072693248
        };
        IR::Function newfunc((const uint32_t*)data, sizeof(data) / 4);
        test_str_eq(newfunc, "Float64 (Float64 %0, Float64 %1) {\n"
                    "L0:\n"
                    "  Float64 %3 = interp [1, (4) +1] (Float64 %0) {1, 2, 2.5, 1}\n"
                    "  Float64 %2 = add Float64 %3, Float64 %1\n"
                    "  ret Float64 %2\n"
                    "}");
        auto f = (double(*)(double, double))gettest(newfunc)
            // .print()
            // .opt()
            // .print()
            .get_ptr();
        assert(f(1.5, 2) == 4.25);
    }

    return 0;
}