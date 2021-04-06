/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_UTILS_LLVM_CODEGEN_H__
#define __NACS_UTILS_LLVM_CODEGEN_H__

#include "../utils.h"
#include "../ir.h"

#include <llvm/ADT/StringMap.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/MDBuilder.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Type.h>
#include <llvm/IR/Value.h>

#include <functional>

namespace NaCs {
namespace LLVM {
namespace Codegen {

using namespace llvm;

struct Wrapper {
    enum ArgType : uint8_t {
        Closure = 0,
        ByRef = 1 << 0,
        Vector = 1 << 1,
    };
    struct Arg {
        ArgType type{ArgType(0)};
        uint32_t idx{0};
    };
    bool closure = false;
    unsigned vector_size = 1;
    std::map<uint32_t,Arg> arg_map{};
    Wrapper &add_closure(uint32_t arg, uint32_t idx)
    {
        arg_map.emplace(arg, Arg{Closure, idx});
        return *this;
    }
    Wrapper &add_byref(uint32_t arg)
    {
        auto &spec = arg_map[arg];
        spec.type = ArgType(spec.type | ByRef);
        return *this;
    }
    Wrapper &add_ret_ref()
    {
        return add_byref(uint32_t(-1));
    }
    Wrapper &add_vector(uint32_t arg)
    {
        auto &spec = arg_map[arg];
        spec.type = ArgType(spec.type | Vector);
        return *this;
    }
};

class NACS_EXPORT_ Context {
public:
    Context(Module *mod);
    Context(LLVMContext &ctx);
    virtual ~Context();
    virtual Function *emit_wrapper(Function *func, StringRef name, const Wrapper &spec);
    virtual Function *emit_function(const IR::Function &func, StringRef name,
                                    bool _export=true);
    virtual Function *emit_function(const uint8_t *data, size_t sz, StringRef name,
                                    bool _export=true);

    Type *llvm_ty(IR::Type ty) const;
    Type *llvm_argty(IR::Type ty) const;
    Value *emit_const(IR::TagVal c) const;
    Value *emit_convert(IRBuilder<> &builder, IR::Type ty, Value *val) const;

    LLVMContext &get_context() const
    {
        return m_ctx;
    }
    Module *get_module() const
    {
        return m_mod;
    }
    void set_module(Module *mod)
    {
        m_mod = mod;
    }
    // If `use_extern_data` returns `true`, the constant data needed by the compiled code
    // (e.g. for interpolation) will be passed to `add_data`.
    virtual bool use_extern_data() const
    {
        return false;
    }
    virtual void add_extern_data(StringRef, const void*, size_t)
    {
    }
    virtual uintptr_t get_extern_data(StringRef)
    {
        return 0;
    }
    std::function<uintptr_t(const std::string&)> get_extern_resolver()
    {
        if (!use_extern_data())
            return std::function<uintptr_t(const std::string&)>();
        return std::bind(&Context::get_extern_data, this, std::placeholders::_1);
    }

private:
    Value *emit_add(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_sub(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_mul(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_fdiv(IRBuilder<> &builder, Value *val1, Value *val2) const;
    Value *emit_cmp(IRBuilder<> &builder, IR::CmpType cmptyp, Value *val1, Value *val2) const;
#if LLVM_VERSION_MAJOR >= 11
    FunctionCallee ensurePureFunc(StringRef name, FunctionType *ft, bool canread=false) const;
#else
    Constant *ensurePureFunc(StringRef name, FunctionType *ft, bool canread=false) const;
#endif

    Module *m_mod;
    LLVMContext &m_ctx;
public:
    IntegerType *const T_bool;
    IntegerType *const T_i8;
    IntegerType *const T_i32;
    Type *const T_f64;
    FunctionType *const F_f64_f64;
    FunctionType *const F_f64_f64f64;
    FunctionType *const F_f64_f64f64f64;
    FunctionType *const F_f64_f64i32;
    FunctionType *const F_f64_i32f64;
    FunctionType *const F_f64_f64i32pf64;
    UndefValue *const V_undefbool;
private:
    MDBuilder m_mdbuilder;
public:
    MDNode *tbaa_root;
    MDNode *tbaa_const;
};

class NACS_EXPORT_ CachedContext: public Context {
public:
    CachedContext(Module *mod);
    CachedContext(LLVMContext &ctx);
    ~CachedContext();
    Function *emit_function(const uint8_t *data, size_t sz, StringRef name,
                            bool _export=true) override;
    void clear_cache();
private:
    llvm::StringMap<llvm::Function*> m_cache;
};

}
}
}

#endif
