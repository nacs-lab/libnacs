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

class Context {
public:
    NACS_EXPORT(utils) Context(Module *mod);
    Function *emit_function(const IR::Function &func, uint64_t func_id,
                            const std::map<uint32_t,uint32_t> &closure_args)
    {
        return _emit_function(func, func_id, closure_args, true);
    }
    Function *emit_function(const IR::Function &func, uint64_t func_id)
    {
        return _emit_function(func, func_id, {}, true);
    }
private:
    NACS_EXPORT(utils)
    Function *_emit_function(const IR::Function &func, uint64_t func_id,
                             const std::map<uint32_t,uint32_t> &closure_args,
                             bool has_closure);
    Type *llvm_ty(IR::Type ty) const;
    Type *llvm_argty(IR::Type ty) const;
    Value *emit_const(IR::TagVal c) const;
    Value *emit_convert(IRBuilder<> &builder, IR::Type ty, Value *val) const;
    Value *emit_add(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_sub(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_mul(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const;
    Value *emit_fdiv(IRBuilder<> &builder, Value *val1, Value *val2) const;
    Value *emit_cmp(IRBuilder<> &builder, IR::CmpType cmptyp, Value *val1, Value *val2) const;
    Constant *ensurePureFunc(StringRef name, FunctionType *ft, bool canread=false) const;

    Module *m_mod;
    LLVMContext &m_ctx;
    IntegerType *T_bool;
    IntegerType *T_i8;
    IntegerType *T_i32;
    Type *T_f64;
    FunctionType *F_f64_f64;
    FunctionType *F_f64_f64f64;
    FunctionType *F_f64_f64f64f64;
    FunctionType *F_f64_f64i32;
    FunctionType *F_f64_i32f64;
    FunctionType *F_f64_f64f64f64i32pf64;
    MDBuilder m_mdbuilder;
    MDNode *tbaa_root;
    MDNode *tbaa_const;
    uint64_t m_counter = 0;
};

}
}
}

#endif
