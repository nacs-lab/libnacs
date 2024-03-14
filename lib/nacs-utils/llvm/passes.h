/*************************************************************************
 *   Copyright (c) 2018 - 2024 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_UTILS_LLVM_PASSES_H__
#define __NACS_UTILS_LLVM_PASSES_H__

#include <llvm/Config/llvm-config.h>

#if LLVM_VERSION_MAJOR <= 16
#  define NACS_ENABLE_LEGACY_PASS 1
#else
#  define NACS_ENABLE_LEGACY_PASS 0
#endif
#if LLVM_VERSION_MAJOR >= 16
#  define NACS_ENABLE_NEW_PASS 1
#else
#  define NACS_ENABLE_NEW_PASS 0
#endif

#if NACS_ENABLE_LEGACY_PASS
#  include <llvm/Pass.h>
#endif
#if NACS_ENABLE_NEW_PASS
#  include <llvm/IR/PassManager.h>
#endif

namespace NaCs::LLVM {

using namespace llvm;

bool fixVectorABI(Function &F);

#if NACS_ENABLE_LEGACY_PASS
Pass *createElimMachOPrefixPass();
Pass *createGlobalRenamePass();
#endif

Pass *createLowerVectorPass();
using resolver_cb_t = std::function<uintptr_t(const std::string&)>;
FunctionPass *createNaCsInstSimplifyPass(const resolver_cb_t &cb=resolver_cb_t());
Pass *createVectorABIPass();

#if NACS_ENABLE_NEW_PASS
struct ElimMachOPrefixPass : PassInfoMixin<ElimMachOPrefixPass> {
    PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};
struct GlobalRenamePass : PassInfoMixin<GlobalRenamePass> {
    PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};
#endif

}

#endif
