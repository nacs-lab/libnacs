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

#include "compile.h"
#include "compile_p.h"

#include "../utils.h"

#include <llvm/Analysis/TargetTransformInfo.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/MC/MCContext.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Scalar/GVN.h>
#include <llvm/Transforms/Vectorize.h>

namespace NaCs {
namespace LLVM {
namespace Compile {

void addOptimization(legacy::PassManagerBase &pm)
{
    pm.add(createCFGSimplificationPass());
    pm.add(createDeadInstEliminationPass());
    pm.add(createPromoteMemoryToRegisterPass());
    pm.add(createEarlyCSEPass());
    pm.add(createMergePhiPass());
    pm.add(createDeadInstEliminationPass());
    pm.add(createInstructionCombiningPass());

    pm.add(createInstructionCombiningPass()); // Cleanup for scalarrepl.
    pm.add(createSROAPass());                 // Break up aggregate allocas
    pm.add(createInstructionCombiningPass()); // Cleanup for scalarrepl.
    pm.add(createJumpThreadingPass());        // Thread jumps.
    pm.add(createInstructionCombiningPass()); // Combine silly seq's
    pm.add(createReassociatePass());          // Reassociate expressions
    pm.add(createEarlyCSEPass()); //// ****

    pm.add(createLoopIdiomPass()); //// ****
    pm.add(createLoopRotatePass());           // Rotate loops.
    pm.add(createLICMPass());                 // Hoist loop invariants
    pm.add(createLoopUnswitchPass());         // Unswitch loops.
    pm.add(createInstructionCombiningPass());
    pm.add(createIndVarSimplifyPass());       // Canonicalize indvars
    pm.add(createLoopDeletionPass());         // Delete dead loops
    pm.add(createSimpleLoopUnrollPass());     // Unroll small loops

    pm.add(createSROAPass());                 // Break up aggregate allocas
    pm.add(createInstructionCombiningPass()); // Clean up after the unroller
    pm.add(createGVNPass());                  // Remove redundancies
    pm.add(createSCCPPass());                 // Constant prop with SCCP

    pm.add(createSinkingPass()); ////////////// ****
    pm.add(createInstructionSimplifierPass());///////// ****
    pm.add(createInstructionCombiningPass());
    pm.add(createJumpThreadingPass());         // Thread jumps
    pm.add(createDeadStoreEliminationPass());  // Delete dead stores

    // see if all of the constant folding has exposed more loops
    // to simplification and deletion
    // this helps significantly with cleaning up iteration
    pm.add(createCFGSimplificationPass());     // Merge & remove BBs
    pm.add(createLoopIdiomPass());
    pm.add(createLoopDeletionPass());          // Delete dead loops
    pm.add(createJumpThreadingPass());         // Thread jumps
    pm.add(createInstructionCombiningPass());   // Clean up after SLP loop vectorizer

    pm.add(createAggressiveDCEPass());         // Delete dead instructions
    pm.add(createGlobalDCEPass());
    pm.add(createConstantMergePass());
    pm.add(createMergeFunctionsPass());
}

NACS_EXPORT() bool emit_objfile(raw_pwrite_stream &stm, TargetMachine *tgt, Module *M, bool opt)
{
    legacy::PassManager pm;
    pm.add(new TargetLibraryInfoWrapperPass(Triple(tgt->getTargetTriple())));
    pm.add(createTargetTransformInfoWrapperPass(tgt->getTargetIRAnalysis()));
    if (opt)
        addOptimization(pm);
    MCContext *ctx;
    if (tgt->addPassesToEmitMC(pm, ctx, stm))
        return false;
    pm.run(*M);
    return true;
}

// For testing only
NACS_EXPORT() Module *optimize(Module *mod)
{
    legacy::PassManager pm;
    addOptimization(pm);
    pm.run(*mod);
    return mod;
}

static TargetMachine *create_native_target()
{
    InitializeNativeTarget();
    InitializeNativeTargetAsmPrinter();
    InitializeNativeTargetAsmParser();

    TargetOptions options;
    EngineBuilder eb;
    eb.setEngineKind(EngineKind::JIT)
        .setTargetOptions(options)
        .setRelocationModel(Reloc::Static)
        .setOptLevel(CodeGenOpt::Aggressive);
    if (sizeof(void*) > 4)
        eb.setCodeModel(CodeModel::Large);
    Triple TheTriple(sys::getProcessTriple());
    StringMap<bool> HostFeatures;
    sys::getHostCPUFeatures(HostFeatures);
    SmallVector<std::string,10> attr;
    for (auto it = HostFeatures.begin(); it != HostFeatures.end(); it++) {
        if (it->getValue()) {
            attr.append(1, it->getKey().str());
        }
    }
    // Explicitly disabled features need to be added at the end so that
    // they are not reenabled by other features that implies them by default.
    for (auto it = HostFeatures.begin(); it != HostFeatures.end(); it++) {
        if (!it->getValue()) {
            attr.append(1, std::string("-") + it->getKey().str());
        }
    }
    std::string cpu = sys::getHostCPUName();
    return eb.selectTarget(TheTriple, "", cpu, attr);
}

NACS_EXPORT() TargetMachine *get_native_target()
{
    static std::unique_ptr<TargetMachine> tgt(create_native_target());
    return tgt.get();
}

}
}
}
