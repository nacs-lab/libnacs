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

#include "../macros.h"

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
#  include <llvm/Passes/PassBuilder.h>
#endif

namespace NaCs::LLVM {

using namespace llvm;

using resolver_cb_t = std::function<uintptr_t(const std::string&)>;
bool fixVectorABI(Function &F);
void runGlobalRenamePasses(Module &M);
void runAlwaysInlinerPasses(Module &M);

#if NACS_ENABLE_LEGACY_PASS
Pass *createElimMachOPrefixPass();
Pass *createGlobalRenamePass();
Pass *createLowerVectorPass();
FunctionPass *createNaCsInstSimplifyPass(const resolver_cb_t &cb=resolver_cb_t());
Pass *createVectorABIPass();
#endif

#if NACS_ENABLE_NEW_PASS
struct ElimMachOPrefixPass : PassInfoMixin<ElimMachOPrefixPass> {
    NACS_EXPORT(utils) PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};
struct GlobalRenamePass : PassInfoMixin<GlobalRenamePass> {
    NACS_EXPORT(utils) PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};
struct LowerVectorPass : PassInfoMixin<LowerVectorPass> {
    NACS_EXPORT(utils) PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};
struct NaCsInstSimplifyPass : PassInfoMixin<NaCsInstSimplifyPass> {
    NaCsInstSimplifyPass(const resolver_cb_t &cb=resolver_cb_t())
        : m_cb(cb)
    {
    }
    NACS_EXPORT(utils) PreservedAnalyses run(Function &F, FunctionAnalysisManager &AM);
    static bool isRequired() { return false; }
private:
    resolver_cb_t m_cb;
};
struct VectorABIPass : PassInfoMixin<VectorABIPass> {
    NACS_EXPORT(utils) PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM);
    static bool isRequired() { return true; }
};

struct AnalysisManagers {
    LoopAnalysisManager LAM;
    FunctionAnalysisManager FAM;
    CGSCCAnalysisManager CGAM;
    ModuleAnalysisManager MAM;

    NACS_EXPORT(utils) AnalysisManagers(PassBuilder &PB);
    NACS_EXPORT(utils) AnalysisManagers(TargetMachine &TM, PassBuilder &PB);
    NACS_EXPORT(utils) ~AnalysisManagers();
};

#endif

}

#endif
