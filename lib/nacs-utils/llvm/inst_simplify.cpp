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

#include "inst_simplify.h"

#include "../number.h"
#include "execute.h"

#include <llvm/ADT/DepthFirstIterator.h>
#include <llvm/Analysis/AssumptionCache.h>
#include <llvm/Analysis/ConstantFolding.h>
#include <llvm/Analysis/InstructionSimplify.h>
#include <llvm/Analysis/OptimizationRemarkEmitter.h>
#include <llvm/Analysis/TargetLibraryInfo.h>
#include <llvm/InitializePasses.h>
#include <llvm/IR/BasicBlock.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/Module.h>
#include <llvm/Transforms/Utils/Local.h>


namespace NaCs::LLVM {

template<typename... Args>
static double call_ptr(uintptr_t ptr, Args... args)
{
    return ((double(*)(Args...))ptr)(args...);
}

static Value *simplify_inst(Instruction *I, const SimplifyQuery &SQ,
                            const resolver_cb_t &cb, OptimizationRemarkEmitter *ORE)
{
#if LLVM_VERSION_MAJOR >= 15
    auto v = simplifyInstruction(I, SQ, ORE);
#else
    auto v = SimplifyInstruction(I, SQ, ORE);
#endif
    if (v) {
        if (!isa<Instruction>(v))
            return v;
        I = cast<Instruction>(v);
    }
    // All the functions we are intersted in returns `double`
    if (!I->getType()->isDoubleTy())
        return v;
    if (!isa<CallInst>(I))
        return v;
    auto call = cast<CallInst>(I);
    auto callee = call->getCalledFunction();
    if (!callee)
        return v;
    auto get_i32 = [&] (Value *v, auto *res) {
        auto ci = dyn_cast<ConstantInt>(v);
        if (!ci)
            return false;
        if (std::is_signed_v<decltype(*res)>) {
            *res = ci->getSExtValue();
        }
        else {
            *res = ci->getZExtValue();
        }
        return true;
    };
    auto get_f64 = [&] (Value *v, double *res) {
        auto cf = dyn_cast<ConstantFP>(v);
        if (!cf)
            return false;
        *res = cf->getValueAPF().convertToDouble();
        return true;
    };
    auto clear_v = [&] {
        if (v) {
            I->eraseFromParent();
        }
    };
    auto name = callee->getName();
    auto &llvm_ctx = callee->getContext();
    if (name == "interp" && call->arg_size() == 3) {
        double x;
        uint32_t npoints;
        if (!get_f64(call->getArgOperand(0), &x) ||
            !get_i32(call->getArgOperand(1), &npoints))
            return v;
        auto pv = dyn_cast<Constant>(call->getArgOperand(2));
        if (!pv)
            return v;
        APInt offset;
        GlobalValue *gval;
        if (!IsConstantOffsetFromGlobal(pv, gval, offset, SQ.DL))
            return v;
        auto gv = dyn_cast<GlobalVariable>(gval);
        if (!gv)
            return v;
        const char *data;
        if (gv->hasInitializer()) {
            auto init = dyn_cast<ConstantDataSequential>(gv->getInitializer());
            if (!init)
                return v;
            data = init->getRawDataValues().data();
        }
        else {
            auto gname = gv->getName();
            if (!cb)
                return v;
            data = (const char*)cb(gname.str());
            if (!data) {
                return v;
            }
        }
        const double *points = (const double*)(data + (int32_t)offset.getLimitedValue());
        clear_v();
        return ConstantFP::get(Type::getDoubleTy(llvm_ctx),
                               linearInterpolate(x, npoints, points));
    }
    auto call_sym = [&] (auto... args) -> Value* {
        if (auto addr = Exe::Resolver::resolve_ir_sym(name.str())) {
            clear_v();
            return ConstantFP::get(Type::getDoubleTy(llvm_ctx),
                                   call_ptr(addr, args...));
        }
        return v;
    };
    if (call->arg_size() == 1) {
        double x;
        if (!get_f64(call->getArgOperand(0), &x))
            return v;
        return call_sym(x);
    }
    else if (call->arg_size() == 2) {
        double xf;
        if (!get_f64(call->getArgOperand(0), &xf)) {
            int xi;
            if (!get_i32(call->getArgOperand(0), &xi))
                return v;
            double y;
            if (!get_f64(call->getArgOperand(1), &y))
                return v;
            return call_sym(xi, y);
        }
        double yf;
        if (!get_f64(call->getArgOperand(1), &yf)) {
            int yi;
            if (!get_i32(call->getArgOperand(1), &yi))
                return v;
            return call_sym(xf, yi);
        }
        return call_sym(xf, yf);
    }
    return v;
}

NACS_EXPORT() bool instSimplify(Function &F, const SimplifyQuery &SQ,
                                const resolver_cb_t &cb, OptimizationRemarkEmitter *ORE)
{
    SmallPtrSet<const Instruction*, 8> S1, S2, *to_simplify = &S1, *next = &S2;
    bool changed = false;

    do {
        for (auto BB: depth_first(&F.getEntryBlock())) {
            // The iterator must be incremented before the loop body since we are
            // deleting the instruction.
            for (BasicBlock::iterator BI = BB->begin(), BE = BB->end(); BI != BE;) {
                Instruction *I = &*BI++;
                // The first time through the loop `to_simplify` is empty and we try to
                // simplify all instructions. On later iterations `to_simplify` is not
                // empty and we only bother simplifying instructions that are in it.
                if (!to_simplify->empty() && !to_simplify->count(I))
                    continue;

                // Don't waste time simplifying unused instructions.
                if (!I->use_empty()) {
                    if (Value *V = simplify_inst(I, SQ, cb, ORE)) {
                        // Mark all uses for resimplification next time round the loop.
                        for (User *U : I->users())
                            next->insert(cast<Instruction>(U));
                        I->replaceAllUsesWith(V);
                        changed = true;
                    }
                }
                if (RecursivelyDeleteTriviallyDeadInstructions(I, SQ.TLI)) {
                    // RecursivelyDeleteTriviallyDeadInstruction can remove more than one
                    // instruction, so simply incrementing the iterator does not work.
                    // When instructions get deleted re-iterate instead.
                    BI = BB->begin();
                    BE = BB->end();
                    changed = true;
                }
            }
        }

        // Place the list of instructions to simplify on the next loop iteration
        // into `to_simplify`.
        std::swap(to_simplify, next);
        next->clear();
    } while (!to_simplify->empty());

    return changed;
}

NACS_EXPORT() bool instSimplify(Function &F, const resolver_cb_t &cb,
                                OptimizationRemarkEmitter *ORE)
{
    const DataLayout &DL = F.getParent()->getDataLayout();
    return instSimplify(F, SimplifyQuery(DL), cb, ORE);
}

namespace {
struct NaCsInstSimplify : public FunctionPass {
    static char ID;
    resolver_cb_t m_cb;
    NaCsInstSimplify(const resolver_cb_t &cb=resolver_cb_t()) :
        FunctionPass(ID),
        m_cb(cb)
    {
        initializeDominatorTreeWrapperPassPass(*PassRegistry::getPassRegistry());
        initializeAssumptionCacheTrackerPass(*PassRegistry::getPassRegistry());
        initializeTargetLibraryInfoWrapperPassPass(*PassRegistry::getPassRegistry());
        initializeOptimizationRemarkEmitterWrapperPassPass(*PassRegistry::getPassRegistry());
    }

    void getAnalysisUsage(AnalysisUsage &AU) const override
    {
        AU.setPreservesCFG();
        AU.addRequired<DominatorTreeWrapperPass>();
        AU.addRequired<AssumptionCacheTracker>();
        AU.addRequired<TargetLibraryInfoWrapperPass>();
        AU.addRequired<OptimizationRemarkEmitterWrapperPass>();
    }

    /// runOnFunction - Remove instructions that simplify.
    bool runOnFunction(Function &F) override
    {
        if (skipFunction(F))
            return false;

        const DominatorTree *DT =
            &getAnalysis<DominatorTreeWrapperPass>().getDomTree();
        const TargetLibraryInfo *TLI =
            &getAnalysis<TargetLibraryInfoWrapperPass>().getTLI(F);
        AssumptionCache *AC =
            &getAnalysis<AssumptionCacheTracker>().getAssumptionCache(F);
        OptimizationRemarkEmitter *ORE =
            &getAnalysis<OptimizationRemarkEmitterWrapperPass>().getORE();
        const DataLayout &DL = F.getParent()->getDataLayout();
        const SimplifyQuery SQ(DL, TLI, DT, AC);
        return instSimplify(F, SQ, m_cb, ORE);
    }
};
}

char NaCsInstSimplify::ID = 0;
// Public interface to the simplify instructions pass.
NACS_EXPORT() FunctionPass *createNaCsInstSimplifyPass(const resolver_cb_t &cb)
{
    return new NaCsInstSimplify(cb);
}

}
