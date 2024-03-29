/*************************************************************************
 *   Copyright (c) 2019 - 2024 Yichao Yu <yyc1992@gmail.com>             *
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

#define DEBUG_TYPE "lower_vector"

#include "passes.h"
#include "codegen_p.h"
#include "utils.h"

#include <llvm/ADT/SmallVector.h>
#include <llvm/ADT/SmallSet.h>
#include <llvm/IR/InstIterator.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Module.h>

namespace NaCs::LLVM {

using namespace llvm;

// This pass lowers the LLVM vector instructions, including vectorized llvm libm intrinsics
// and other vectorized operators, to function calls that our runtime loader can understand.

namespace {

#if LLVM_VERSION_MAJOR < 11
#  define FixedVectorType VectorType
#endif

struct LowerVector {
    bool run(Module &M);
    void replace_frem(BinaryOperator *binop);
    bool process_function(Function &F);
    static void replace_function(Function &F, StringRef new_name);
    static bool process_intrinsic(Function &F, const std::string &name);

    SmallSet<Function*, 4> m_tofix{};
};

void LowerVector::replace_function(Function &F, StringRef new_name)
{
    auto _new_f = F.getParent()->getOrInsertFunction(new_name, F.getFunctionType(),
                                                     F.getAttributes());
    auto new_f = _new_f.getCallee();
    F.replaceAllUsesWith(new_f);
    fixVectorABI(*cast<Function>(new_f));
}

bool LowerVector::process_intrinsic(Function &F, const std::string &name)
{
    auto fname = F.getName();
    if (fname == "llvm." + name + ".v2f64") {
        replace_function(F, name + ".2");
        return true;
    }
    else if (fname == "llvm." + name + ".v4f64") {
        replace_function(F, name + ".4");
        return true;
    }
    else if (fname == "llvm." + name + ".v8f64") {
        replace_function(F, name + ".8");
        return true;
    }
    return false;
}

void LowerVector::replace_frem(BinaryOperator *binop)
{
    auto ty = cast<FixedVectorType>(binop->getType());
    auto nele = ty->getNumElements();
    std::string name = "fmod.";
    // No need to go through general purpose number -> string conversion.
    if (nele == 2) {
        name += "2";
    }
    else if (nele == 4) {
        name += "4";
    }
    else if (nele == 8) {
        name += "8";
    }
    auto fty = FunctionType::get(ty, {ty, ty}, false);
    auto f = Codegen::ensurePureExtern(binop->getModule(), fty, name);
    // The vector ABI fix looks at the caller of vector functions on windows
    // in order to set the calling convension on the call instruction.
    // We should do it after all the instructions have be fixed.
#if LLVM_VERSION_MAJOR >= 11
    m_tofix.insert(cast<Function>(f.getCallee()));
#else
    m_tofix.insert(cast<Function>(f));
#endif
    IRBuilder<> builder(binop);
    auto call = builder.CreateCall(f, {binop->getOperand(0), binop->getOperand(1)});
    binop->replaceAllUsesWith(call);
}

bool LowerVector::process_function(Function &F)
{
    if (process_intrinsic(F, "sin") || process_intrinsic(F, "cos") ||
        process_intrinsic(F, "pow") || process_intrinsic(F, "exp") ||
        process_intrinsic(F, "exp2") || process_intrinsic(F, "log") ||
        process_intrinsic(F, "log10") || process_intrinsic(F, "log2"))
        return true;
    // TODO check rounding intrinsics
    if (F.empty())
        return false;
    SmallVector<Instruction*, 32> to_remove;
    bool changed = false;
    for (auto &I: instructions(F)) {
        auto binop = dyn_cast<BinaryOperator>(&I);
        if (!binop || binop->getOpcode() != Instruction::FRem)
            continue;
        auto ty = binop->getType();
        if (!ty->isVectorTy())
            continue;
        // No SVE support for now.
        auto vecty = cast<FixedVectorType>(ty);
        if (!vecty->getElementType()->isDoubleTy())
            continue;
        switch (vecty->getNumElements()) {
        default:
            continue;
        case 2:
        case 4:
        case 8:
            replace_frem(binop);
            changed = true;
            to_remove.push_back(binop);
            binop->dropAllReferences();
        }
    }
    for (auto I: to_remove)
        I->eraseFromParent();
    return changed;
}

bool LowerVector::run(Module &M)
{
    bool changed = false;
    for (auto &f: M)
        changed |= process_function(f);
    for (auto f: m_tofix)
        fixVectorABI(*f);
    return changed;
}

#if NACS_ENABLE_LEGACY_PASS
struct LegacyLowerVectorPass : public ModulePass {
    static char ID;
    LegacyLowerVectorPass()
        : ModulePass(ID)
    {}

private:
    bool runOnModule(Module &M) override
    {
        LowerVector lv;
        return lv.run(M);
    }
};

char LegacyLowerVectorPass::ID = 0;
static RegisterPass<LegacyLowerVectorPass> X("LowerVector",
                                       "Lower Vector Instructions and Calls Pass",
                                       false /* Only looks at CFG */,
                                       false /* Analysis Pass */);
#endif

}

#if NACS_ENABLE_LEGACY_PASS
NACS_EXPORT() Pass *createLowerVectorPass()
{
    return new LegacyLowerVectorPass();
}
#endif

#if NACS_ENABLE_NEW_PASS
NACS_EXPORT() PreservedAnalyses
LowerVectorPass::run(Module &M, ModuleAnalysisManager &AM)
{
    LowerVector lv;
    if (lv.run(M))
        return PreservedAnalyses::allInSet<CFGAnalyses>();
    return PreservedAnalyses::all();
}
#endif

}
