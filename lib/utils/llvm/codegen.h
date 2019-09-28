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

#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/MDBuilder.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Type.h>
#include <llvm/IR/Value.h>

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

class Context {
public:
    using data_map_t = std::map<std::string,std::pair<uint32_t,uint32_t>>;
    NACS_EXPORT(utils) Context(Module *mod);
    NACS_EXPORT(utils) Function *emit_wrapper(Function *func,
                                              StringRef name, const Wrapper &spec);
    // If `data` is not `NULL`, the constant data needed by the compiled code
    // (e.g. for interpolation) will be managed by the caller.
    // The offset (in size of `uint32_t`) and size of the data in the `func`'s constant table
    // will be stored under the symbol name in the map. The caller is responsible
    // for passing that info to the runtime loader so that the data is accessible by the code.
    NACS_EXPORT(utils) Function *emit_function(const IR::Function &func,
                                               StringRef name, bool _export=true,
                                               data_map_t *data=nullptr);
    inline Function *emit_function(const IR::Function &func, StringRef name, data_map_t *data)
    {
        return emit_function(func, name, true, data);
    }

    Type *llvm_ty(IR::Type ty) const;
    Type *llvm_argty(IR::Type ty) const;
    Value *emit_const(IR::TagVal c) const;
    Value *emit_convert(IRBuilder<> &builder, IR::Type ty, Value *val) const;

private:
    Value *emit_add(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_sub(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_mul(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_fdiv(IRBuilder<> &builder, Value *val1, Value *val2) const;
    Value *emit_cmp(IRBuilder<> &builder, IR::CmpType cmptyp, Value *val1, Value *val2) const;
    Constant *ensurePureFunc(StringRef name, FunctionType *ft, bool canread=false) const;

    Module *m_mod;
    LLVMContext &m_ctx;
public:
    IntegerType *T_bool;
    IntegerType *T_i8;
    IntegerType *T_i32;
    Type *T_f64;
    FunctionType *F_f64_f64;
    FunctionType *F_f64_f64f64;
    FunctionType *F_f64_f64f64f64;
    FunctionType *F_f64_f64i32;
    FunctionType *F_f64_i32f64;
    FunctionType *F_f64_f64i32pf64;
private:
    UndefValue *V_undefbool;
    MDBuilder m_mdbuilder;
    MDNode *tbaa_root;
    MDNode *tbaa_const;
    uint64_t m_counter = 0;
};

}
}
}

#endif
