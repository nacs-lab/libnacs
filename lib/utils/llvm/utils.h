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

#include <llvm/IR/DebugLoc.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Metadata.h>
#include <llvm/IR/Type.h>
#include <llvm/IR/Value.h>
#include <llvm/Support/Debug.h>

#define NACS_LLVM_VER (LLVM_VERSION_MAJOR * 10000 + LLVM_VERSION_MINOR * 100 \
                       + LLVM_VERSION_PATCH)

namespace NaCs {
namespace LLVM {

void dump(llvm::Value *v);
void dump(llvm::Type *v);
void dump(llvm::Function *f);
void dump(llvm::Module *m);
void dump(llvm::Metadata *m);
void dump(llvm::DebugLoc *dbg);

// The following functions are provided so that the user does not need to link to LLVM.
// This is useful/needed when we statically linking LLVM and cannot
// link to the same version of LLVM as the user.
// See also comments in `lib/utils/CMakeLists.txt`
llvm::Module *new_module(llvm::StringRef, llvm::LLVMContext&);
void delete_module(llvm::Module*);
llvm::LLVMContext *new_context();
void delete_context(llvm::LLVMContext*);

}
}

#endif
