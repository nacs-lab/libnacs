/*************************************************************************
 *   Copyright (c) 2024 - 2024 Yichao Yu <yyc1992@gmail.com>             *
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

#include "passes.h"

#include "../utils.h"

#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO/AlwaysInliner.h>

#if NACS_ENABLE_NEW_PASS
#  include <llvm/Analysis/AliasAnalysis.h>
#  include <llvm/Analysis/TargetTransformInfo.h>
#  include <llvm/Target/TargetMachine.h>
#else
#  include <llvm/IR/LegacyPassManager.h>
#endif

namespace NaCs::LLVM {

#if NACS_ENABLE_NEW_PASS
static FunctionAnalysisManager createFAM(TargetMachine &TM)
{
    FunctionAnalysisManager FAM;
    FAM.registerPass([&] {
        AAManager AA;
        TM.registerDefaultAliasAnalyses(AA);
        return AA;
    });
    FAM.registerPass([&] {
        return TargetIRAnalysis(TM.getTargetIRAnalysis());
    });
    FAM.registerPass([&] {
        return TargetLibraryAnalysis(TargetLibraryInfoImpl(TM.getTargetTriple()));
    });
    return FAM;
}

NACS_EXPORT() AnalysisManagers::AnalysisManagers(TargetMachine &TM, PassBuilder &PB)
: LAM(), FAM(createFAM(TM)), CGAM(), MAM()
{
    PB.registerLoopAnalyses(LAM);
    PB.registerFunctionAnalyses(FAM);
    PB.registerCGSCCAnalyses(CGAM);
    PB.registerModuleAnalyses(MAM);
    PB.crossRegisterProxies(LAM, FAM, CGAM, MAM);
}

NACS_EXPORT() AnalysisManagers::AnalysisManagers(PassBuilder &PB)
: LAM(), FAM(), CGAM(), MAM()
{
    PB.registerLoopAnalyses(LAM);
    PB.registerFunctionAnalyses(FAM);
    PB.registerCGSCCAnalyses(CGAM);
    PB.registerModuleAnalyses(MAM);
    PB.crossRegisterProxies(LAM, FAM, CGAM, MAM);
}

NACS_EXPORT() AnalysisManagers::~AnalysisManagers() = default;

NACS_EXPORT() void runAlwaysInlinerPasses(Module &M)
{
    llvm::PassBuilder PB;
    LLVM::AnalysisManagers AM(PB);

    llvm::ModulePassManager MPM;
#  ifndef NDEBUG
    MPM.addPass(llvm::VerifierPass());
#  endif
    MPM.addPass(llvm::AlwaysInlinerPass());
#  ifndef NDEBUG
    MPM.addPass(llvm::VerifierPass());
#  endif
    MPM.run(M, AM.MAM);
}
#else
NACS_EXPORT() void runAlwaysInlinerPasses(Module &M)
{
    llvm::legacy::PassManager PM;
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.add(llvm::createAlwaysInlinerLegacyPass());
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.run(M);
}
#endif

}
