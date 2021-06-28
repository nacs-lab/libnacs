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

#ifndef __NACS_UTILS_LLVM_UTILS_H__
#define __NACS_UTILS_LLVM_UTILS_H__

#include "../utils.h"
#include "../ir.h"

#include <llvm/IR/DebugLoc.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Metadata.h>
#include <llvm/IR/Type.h>
#include <llvm/IR/Value.h>
#include <llvm/Support/Debug.h>
#include <llvm/Transforms/Utils/ValueMapper.h>

#define NACS_LLVM_VER (LLVM_VERSION_MAJOR * 10000 + LLVM_VERSION_MINOR * 100 \
                       + LLVM_VERSION_PATCH)

namespace NaCs::LLVM {

NACS_EXPORT(utils) void dump(const llvm::Value *v);
NACS_EXPORT(utils) void dump(const llvm::Type *v);
NACS_EXPORT(utils) void dump(const llvm::Function *f);
NACS_EXPORT(utils) void dump(const llvm::Module *m);
NACS_EXPORT(utils) void dump(const llvm::Metadata *m);
NACS_EXPORT(utils) void dump(const llvm::DebugLoc *dbg);

// The following functions are provided so that the user does not need to link to LLVM.
// This is useful/needed when we statically linking LLVM and cannot
// link to the same version of LLVM as the user.
// See also comments in `lib/nacs-utils/CMakeLists.txt`
NACS_EXPORT(utils) void delete_module(llvm::Module*);
struct module_deleter {
    void operator()(llvm::Module *mod) const
    {
        delete_module(mod);
    }
};
using module_ref = std::unique_ptr<llvm::Module,module_deleter>;
NACS_EXPORT(utils) module_ref new_module(llvm::StringRef, llvm::LLVMContext&);
NACS_EXPORT(utils) void delete_context(llvm::LLVMContext*);
struct context_deleter {
    void operator()(llvm::LLVMContext *ctx) const
    {
        delete_context(ctx);
    }
};
using context_ref = std::unique_ptr<llvm::LLVMContext,context_deleter>;
NACS_EXPORT(utils) context_ref new_context();
NACS_EXPORT(utils) llvm::GlobalVariable *new_global_variable(
    llvm::Module &M, llvm::Type *Ty, bool isConstant,
    llvm::GlobalValue::LinkageTypes Linkage,
    llvm::Constant *Initializer,
    const llvm::Twine &Name = "",
    llvm::GlobalVariable *InsertBefore = nullptr);
NACS_EXPORT(utils) llvm::GlobalVariable *new_global_variable(
    llvm::Type *Ty, bool isConstant,
    llvm::GlobalValue::LinkageTypes Linkage,
    llvm::Constant *Initializer,
    const llvm::Twine &Name = "");
NACS_EXPORT(utils) llvm::ArrayType *get_array_type(llvm::Type *ElementType,
                                                   uint64_t NumElements);
NACS_EXPORT(utils) llvm::ArrayType *get_array_type(llvm::Type *ElementType,
                                                   uint64_t NumElements);
NACS_EXPORT(utils) llvm::ConstantAggregateZero *get_aggregate_zero(llvm::Type *ty);

NACS_EXPORT(utils) IR::Type get_ir_type(llvm::Type*, bool apitype=true);
NACS_EXPORT(utils) llvm::Value *convert_scalar(llvm::IRBuilder<> &builder,
                                               llvm::Type *typ, llvm::Value *val);

/**
 * Clone a function from one module to another alone with all its dependencies.
 * Simplified from the one in julia.
 */
class FunctionMover final : public llvm::ValueMaterializer {
public:
    FunctionMover(llvm::Module *dest);
    ~FunctionMover();
    llvm::Function *clone_function(llvm::Function *F);

private:
    llvm::ValueToValueMapTy m_vmap;
    llvm::Module *m_dest;
    llvm::SmallVector<llvm::Function*, 16> m_lazy_funcs;

    llvm::Function *queue_proto(llvm::Function *F);
    void clone_body(llvm::Function *F);
    void resolve_lazy();
    llvm::Value *clone_proto(llvm::Function *F);

    llvm::Value *materialize(llvm::Value *V) override;
};

}

#endif
