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

#include "passes.h"
#include "utils.h"

#include "../utils.h"

#include <llvm/Analysis/TargetTransformInfo.h>
#include <llvm/Analysis/TargetLibraryInfo.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/IR/Intrinsics.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/MC/MCContext.h>
#include <llvm/Support/Host.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Transforms/InstCombine/InstCombine.h>
#include <llvm/Transforms/IPO/AlwaysInliner.h>
#include <llvm/Transforms/Scalar/GVN.h>
#include <llvm/Transforms/Utils/ModuleUtils.h>
#include <llvm/Transforms/Utils.h>
#include <llvm/Transforms/Vectorize.h>

#if NACS_ENABLE_NEW_PASS
#  include <llvm/Transforms/IPO/ConstantMerge.h>
#  include <llvm/Transforms/IPO/GlobalDCE.h>
#  include <llvm/Transforms/IPO/MergeFunctions.h>
#  include <llvm/Transforms/Scalar/ADCE.h>
#  include <llvm/Transforms/Scalar/DCE.h>
#  include <llvm/Transforms/Scalar/DeadStoreElimination.h>
#  include <llvm/Transforms/Scalar/EarlyCSE.h>
#  include <llvm/Transforms/Scalar/JumpThreading.h>
#  include <llvm/Transforms/Scalar/IndVarSimplify.h>
#  include <llvm/Transforms/Scalar/LICM.h>
#  include <llvm/Transforms/Scalar/LoopDeletion.h>
#  include <llvm/Transforms/Scalar/LoopIdiomRecognize.h>
#  include <llvm/Transforms/Scalar/LoopRotation.h>
#  include <llvm/Transforms/Scalar/LoopUnrollPass.h>
#  include <llvm/Transforms/Scalar/Reassociate.h>
#  include <llvm/Transforms/Scalar/SCCP.h>
#  include <llvm/Transforms/Scalar/Sink.h>
#  include <llvm/Transforms/Scalar/SROA.h>
#  include <llvm/Transforms/Scalar/SimplifyCFG.h>
#  include <llvm/Transforms/Scalar/SimpleLoopUnswitch.h>
#else
#  include <llvm/Transforms/IPO.h>
#  include <llvm/Transforms/Scalar.h>
#endif

#include <mutex>

namespace NaCs::LLVM::Compile {

#if NACS_ENABLE_NEW_PASS
static void addOptimization(ModulePassManager &MPM)
{
#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
    {
        FunctionPassManager FPM;
        FPM.addPass(SimplifyCFGPass());
        FPM.addPass(DCEPass());
        FPM.addPass(EarlyCSEPass());
        FPM.addPass(DCEPass());
        FPM.addPass(InstCombinePass());
        MPM.addPass(createModuleToFunctionPassAdaptor(std::move(FPM)));
    }

#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
    MPM.addPass(VectorABIPass());            // Fix vector ABI
#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
    MPM.addPass(AlwaysInlinerPass());        // Respect always_inline
    {
        FunctionPassManager FPM;
        FPM.addPass(InstCombinePass());      // Cleanup for scalarrepl.
        FPM.addPass(SROAPass(SROAOptions::PreserveCFG)); // Break up aggregate allocas
        FPM.addPass(InstCombinePass());      // Cleanup for scalarrepl.
        FPM.addPass(JumpThreadingPass());    // Thread jumps.
        FPM.addPass(InstCombinePass());      // Combine silly seq's
        FPM.addPass(ReassociatePass());      // Reassociate expressions
        FPM.addPass(EarlyCSEPass());         //// ****
        {
            LoopPassManager LPM;
            LPM.addPass(LoopIdiomRecognizePass()); //// ****
            LPM.addPass(LoopRotatePass());         // Rotate loops.
            LPM.addPass(LICMPass(LICMOptions()));  // Hoist loop invariants
            LPM.addPass(SimpleLoopUnswitchPass()); // Unswitch loops.

            FPM.addPass(createFunctionToLoopPassAdaptor(std::move(LPM), true));
        }
        FPM.addPass(InstCombinePass());      // Combine silly seq's

        {
            LoopPassManager LPM;
            LPM.addPass(IndVarSimplifyPass());     // Canonicalize indvars
            LPM.addPass(LoopDeletionPass());       // Delete dead loops

            FPM.addPass(createFunctionToLoopPassAdaptor(std::move(LPM)));
        }

        FPM.addPass(LoopUnrollPass());       // Unroll small loops
        FPM.addPass(SROAPass(SROAOptions::PreserveCFG)); // Break up aggregate allocas

        FPM.addPass(InstCombinePass());     // Clean up after the unroller
        FPM.addPass(GVNPass());             // Remove redundancies
        FPM.addPass(SCCPPass());            // Constant prop with SCCP
        FPM.addPass(SinkingPass()); ////////////// ****

#ifndef NDEBUG
        FPM.addPass(VerifierPass());
#endif
        FPM.addPass(NaCsInstSimplifyPass()); // TODO support external symbol
#ifndef NDEBUG
        FPM.addPass(VerifierPass());
#endif

        FPM.addPass(InstCombinePass());
        FPM.addPass(JumpThreadingPass());    // Thread jumps.
        FPM.addPass(DSEPass());              // Delete dead stores

        // see if all of the constant folding has exposed more loops
        // to simplification and deletion
        // this helps significantly with cleaning up iteration
        FPM.addPass(SimplifyCFGPass());      // Merge & remove BBs

        {
            LoopPassManager LPM;
            LPM.addPass(LoopIdiomRecognizePass());
            LPM.addPass(LoopDeletionPass());       // Delete dead loops

            FPM.addPass(createFunctionToLoopPassAdaptor(std::move(LPM)));
        }

        FPM.addPass(JumpThreadingPass());    // Thread jumps.
        FPM.addPass(InstCombinePass());      // Clean up after SLP loop vectorizer

        MPM.addPass(createModuleToFunctionPassAdaptor(std::move(FPM)));
    }

#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
    MPM.addPass(LowerVectorPass());
#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
    MPM.addPass(AlwaysInlinerPass());  // Inlining for lower vector pass

    {
        FunctionPassManager FPM;
        FPM.addPass(ADCEPass());         // Delete dead instructions

        MPM.addPass(createModuleToFunctionPassAdaptor(std::move(FPM)));
    }
    MPM.addPass(GlobalDCEPass());
    MPM.addPass(ConstantMergePass());
    MPM.addPass(MergeFunctionsPass());
#ifndef NDEBUG
    MPM.addPass(VerifierPass());
#endif
}
#else
static void addOptimization(legacy::PassManagerBase &pm)
{
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
    pm.add(createCFGSimplificationPass());
    pm.add(createDeadCodeEliminationPass());
    pm.add(createEarlyCSEPass());
    pm.add(createDeadCodeEliminationPass());
    pm.add(createInstructionCombiningPass());

#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
    pm.add(createVectorABIPass());            // Fix vector ABI
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
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
#if LLVM_VERSION_MAJOR < 15
    // Removed from legacy pass manager
    pm.add(createLoopUnswitchPass());         // Unswitch loops.
#endif
    pm.add(createInstructionCombiningPass());
    pm.add(createIndVarSimplifyPass());       // Canonicalize indvars
    pm.add(createLoopDeletionPass());         // Delete dead loops
    pm.add(createSimpleLoopUnrollPass());     // Unroll small loops

    pm.add(createSROAPass());                 // Break up aggregate allocas
    pm.add(createInstructionCombiningPass()); // Clean up after the unroller
    pm.add(createGVNPass());                  // Remove redundancies
    pm.add(createSCCPPass());                 // Constant prop with SCCP

    pm.add(createSinkingPass()); ////////////// ****
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
    pm.add(createNaCsInstSimplifyPass()); // TODO support external symbol
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
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

#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
    pm.add(createLowerVectorPass());
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
    pm.add(createAlwaysInlinerLegacyPass());  // Inlining for lower vector pass

    pm.add(createAggressiveDCEPass());         // Delete dead instructions
    pm.add(createGlobalDCEPass());
    pm.add(createConstantMergePass());
    pm.add(createMergeFunctionsPass());
#ifndef NDEBUG
    pm.add(createVerifierPass());
#endif
}
#endif

NACS_EXPORT() bool emit_objfile(raw_pwrite_stream &stm, TargetMachine *tgt, Module *M, bool opt)
{
    auto &triple = tgt->getTargetTriple();
    M->setTargetTriple(triple.getTriple());
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
#if 0
    // Disable for now since the verifier doesn't like this.
    auto init_tramp = M->getFunction("llvm.init.trampoline");
    if (!init_tramp || init_tramp->use_empty())
        appendToCompilerUsed(*M, {Intrinsic::getDeclaration(M, Intrinsic::init_trampoline)});
#endif

#if NACS_ENABLE_NEW_PASS
    PassBuilder PB;
    AnalysisManagers AM(*tgt, PB);
    ModulePassManager MPM;
    if (triple.getObjectFormat() == Triple::ObjectFormatType::MachO)
        MPM.addPass(ElimMachOPrefixPass());
    if (opt)
        addOptimization(MPM);
    MPM.run(*M, AM.MAM);

    legacy::PassManager pm;
    MCContext *ctx;
    if (tgt->addPassesToEmitMC(pm, ctx, stm))
        return false;
    pm.run(*M);
#else
    legacy::PassManager pm;
    pm.add(new TargetLibraryInfoWrapperPass(Triple(tgt->getTargetTriple())));
    pm.add(createTargetTransformInfoWrapperPass(tgt->getTargetIRAnalysis()));
    if (triple.getObjectFormat() == Triple::ObjectFormatType::MachO)
        pm.add(createElimMachOPrefixPass());
    if (opt)
        addOptimization(pm);
    MCContext *ctx;
    if (tgt->addPassesToEmitMC(pm, ctx, stm))
        return false;
    pm.run(*M);
#endif
    return true;
}

NACS_EXPORT() bool emit_objfile(SmallVectorImpl<char> &vec, TargetMachine *tgt,
                                Module *M, bool opt)
{
    llvm::raw_svector_ostream stm(vec);
    return emit_objfile(stm, tgt, M, opt);
}

// For testing only (to avoid directly linking to libLLVM in test).
// Does not need to be efficient...
NACS_EXPORT() bool emit_objfile(std::vector<char> &vec, TargetMachine *tgt,
                                Module *M, bool opt)
{
    llvm::SmallVector<char,0> svec;
    auto res = emit_objfile(svec, tgt, M, opt);
    vec.resize(svec.size());
    memcpy(&vec[0], &svec[0], svec.size());
    return res;
}

// For testing only
NACS_EXPORT() Module *optimize(Module *mod)
{
#if NACS_ENABLE_NEW_PASS
    PassBuilder PB;
    AnalysisManagers AM(PB);
    ModulePassManager MPM;
    addOptimization(MPM);
    MPM.run(*mod, AM.MAM);
#else
    legacy::PassManager pm;
    addOptimization(pm);
    pm.run(*mod);
#endif
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
    auto cpu = sys::getHostCPUName();
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
