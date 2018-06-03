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

#include <llvm/IR/LegacyPassManager.h>
#include <llvm/MC/MCContext.h>
#include <llvm/Target/TargetMachine.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Scalar/GVN.h>
#include <llvm/Transforms/Vectorize.h>

namespace NaCs {
namespace LLVM {
namespace Compile {

NACS_EXPORT() bool emit_objfile(raw_pwrite_stream &stm, TargetMachine &tgt, Module *M)
{
    legacy::PassManager pm;
    MCContext *ctx;
    if (tgt.addPassesToEmitMC(pm, ctx, stm))
        return false;
    pm.run(*M);
    return true;
}

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
}

/**
 * Tmp
 */

NACS_EXPORT() Function *optimize(Function *f)
{
    legacy::FunctionPassManager pm(f->getParent());
    addOptimization(pm);
    pm.doInitialization();
    pm.run(*f);
    return f;
}

NACS_EXPORT() Module *optimize(Module *mod)
{
    legacy::PassManager pm;
    addOptimization(pm);
    pm.run(*mod);
    return mod;
}

}
}
}
