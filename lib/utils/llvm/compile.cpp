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
#include "vector_abi.h"
#include "inst_simplify.h"
#include "lower_vector.h"
#include "utils.h"

#include "../utils.h"

#include <llvm/Analysis/TargetTransformInfo.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/IR/Intrinsics.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/MC/MCContext.h>
#include <llvm/Support/TargetSelect.h>
#if LLVM_VERSION_MAJOR >= 7
#  include <llvm/Transforms/InstCombine/InstCombine.h>
#endif
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/IPO/AlwaysInliner.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Scalar/GVN.h>
#include <llvm/Transforms/Utils/ModuleUtils.h>
#if LLVM_VERSION_MAJOR >= 7
#  include <llvm/Transforms/Utils.h>
#endif
#include <llvm/Transforms/Vectorize.h>

#include <mutex>

namespace NaCs {
namespace LLVM {
namespace Compile {

void addOptimization(legacy::PassManagerBase &pm)
{
    pm.add(createCFGSimplificationPass());
    pm.add(createDeadInstEliminationPass());
    pm.add(createEarlyCSEPass());
    pm.add(createDeadInstEliminationPass());
    pm.add(createInstructionCombiningPass());

    pm.add(createVectorABIPass());            // Fix vector ABI
    pm.add(createAlwaysInlinerLegacyPass());  // Respect always_inline
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
    pm.add(createNaCsInstSimplifyPass()); // TODO support external symbol
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

    pm.add(createLowerVectorPass());
    pm.add(createAlwaysInlinerLegacyPass());  // Inlining for lower vector pass

    pm.add(createAggressiveDCEPass());         // Delete dead instructions
    pm.add(createGlobalDCEPass());
    pm.add(createConstantMergePass());
    pm.add(createMergeFunctionsPass());
}

NACS_EXPORT() bool emit_objfile(raw_pwrite_stream &stm, TargetMachine *tgt, Module *M, bool opt)
{
    M->setTargetTriple(tgt->getTargetTriple().getTriple());
    M->setDataLayout(tgt->createDataLayout());
    // By default, if we didn't need executable stack, LLVM emits
    // `.note.GNU-stack` section to disable executable stack.
    // However, we don't use the system linker and don't care about this section
    // and we'd like to disable this section.
    // So far, the only way I can find to do this is to trick LLVM to think we are using
    // the trampoline intrinsics which is done by the code below.
    // This may not be very reliable but should be safe to do so even if it break
    // in the future this should only be a performance (memory consumption) issue
    // rather than a correctness issue.
    auto init_tramp = M->getFunction("llvm.init.trampoline");
    if (!init_tramp || init_tramp->use_empty())
        appendToCompilerUsed(*M, {Intrinsic::getDeclaration(M, Intrinsic::init_trampoline)});
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

NACS_EXPORT() bool emit_objfile(SmallVectorImpl<char> &vec, TargetMachine *tgt,
                                Module *M, bool opt)
{
    llvm::raw_svector_ostream stm(vec);
    return emit_objfile(stm, tgt, M, opt);
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
    static std::once_flag flag;
    std::call_once(flag, [] {
            InitializeNativeTarget();
            InitializeNativeTargetAsmPrinter();
            InitializeNativeTargetAsmParser();
        });

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

NACS_EXPORT() std::unique_ptr<TargetMachine> create_target(StringRef triple,
                                                           StringRef cpu, StringRef features)
{
    static std::once_flag flag;
    std::call_once(flag, [] {
            InitializeAllTargets();
            InitializeAllTargetInfos();
            InitializeAllTargetMCs();
            InitializeAllAsmPrinters();
            InitializeAllAsmParsers();
        });
    Triple TheTriple(triple);
    TargetOptions options;
    EngineBuilder eb;
    eb.setEngineKind(EngineKind::JIT)
        .setTargetOptions(options)
        .setRelocationModel(Reloc::Static)
        .setOptLevel(CodeGenOpt::Aggressive);
    if (TheTriple.isArch64Bit())
        eb.setCodeModel(CodeModel::Large);
    SmallVector<StringRef,16> attr_sr;
    features.split(attr_sr, ",", -1, false);
    SmallVector<std::string,16> attr;
    for (auto a : attr_sr)
        attr.push_back(a.str());
    std::string ec;
    eb.setErrorStr(&ec);
    auto tgt = eb.selectTarget(TheTriple, "", cpu, attr);
    if (!tgt)
        throw std::runtime_error(ec);
    return std::unique_ptr<TargetMachine>(tgt);
}

}
}
}
