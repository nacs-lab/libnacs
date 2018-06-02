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

#define DEBUG_TYPE "merge_phi"

#include "codegen_p.h"
#include "utils.h"

#include <llvm/ADT/SmallVector.h>
#include <llvm/IR/ConstantFolder.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Value.h>

namespace NaCs {
namespace LLVM {

using namespace llvm;

/**
 * Combine `phi` nodes with `cmp`, `select` and possibly other `phi` nodes.
 * This helps optimize our simple and ugly lowering of `phi` node.
 * We lower `phi` node to something like,
 *
 * %prev_bb = load i32, i32 *%bb_num
 * %phi0_input0 = load @T0, @T0 *%slot0
 * %phi0_input1 = load @T0, @T0 *%slot1
 * %phi0_input2 = load @T0, @T0 *%slot2
 * ; ...
 * %phi0_flag0 = icmp eq i32 %prev_bb, <constant 0>
 * %phi0_flag1 = icmp eq i32 %prev_bb, <constant 1>
 * ; ...
 * %phi0_select0 = select i1 %phi0_flag0, %phi0_input1, %phi0_input0
 * %phi0_select1 = select i1 %phi0_flag1, %phi0_input2, %phi0_select0
 * ; ...
 *
 * The `%phi*_input*` may be converted before used as input to `select` in which case
 * this pass will not optimize it very well. After `mem2reg` (and `earlycse`), this becomes
 *
 * %prev_bb = phi i32 ... ; constants only, this is the mem2reg of the prev_bb load
 * %phi0_input0 = phi @T0 ... ; these are the mem2reg of the input load
 * ; ...
 * ; code below are not changed other than cse
 * %phi0_flag0 = icmp eq i32 %prev_bb, <constant 0>
 * %phi0_flag1 = icmp eq i32 %prev_bb, <constant 1>
 * ; ...
 * %phi0_select0 = select i1 %phi0_flag0, @T0 %phi0_input1, @T0 %phi0_input0
 * %phi0_select1 = select i1 %phi0_flag1, @T0 %phi0_input2, @T0 %phi0_select0
 * ; ...
 *
 * This pass does three optimizations to collapsed the code.
 * 1. Merge `cmp` of `phi` node and constant. The `phi` node can only have constant inputs.
 *    The `cmp` is converted into a `phi` node with the inputs being constant foldings
 *    of the original `phi` inputs with the constant `cmp`, i.e. this turns,
 *
 *    %i = phi @T [ @C0, %L1 ], [ @C1, %L2 ]
 *    %f = icmp COND @T %i, @C2
 *
 *    into
 *
 *    %i = phi @T [ @C0, %L1 ], [ @C1, %L2 ]
 *    %f = phi i1 [ icmp COND (@T @C0, @T @C2), %L1 ], [ icmp COND (@T @C1, @T @C2), %L2 ]
 *
 *    This converts all the `%phi*_flag*` to `phi` nodes.
 * 2. Merge `select` of `i1` `phi` node with constant int inputs.
 *    The inputs of the `phi` node must be constant or dominate the first non-phi instruction
 *    in the BB. The `select` is converted into a `phi` node with the inputs being
 *    constant foldings of the original `phi` inputs with the `select`, i.e. this turns
 *
 *    %f = phi i1 [ true, %L1 ], [ false, %L2 ]
 *    %v = select i1 %f, @T %v1, @T %v2
 *
 *    into
 *
 *    %f = phi i1 [ true, %L1 ], [ false, %L2 ]
 *    %v = phi @T [ %v1, %L1 ], [ %v2, %L2 ]
 *
 *    This merges the phi node from the first pass with the `%phi*_select*`s.
 * 3. Merge `phi` nodes with `phi` node in the same BB as input, i.e. this turns
 *
 *    %p1 = phi @T [ %v1, %L1 ], [ %v2, %L2 ]
 *    %p2 = phi @T [ %p1, %L1 ], [ %v3, %L2 ]
 *
 *    into
 *
 *    %p1 = phi @T [ %v1, %L1 ], [ %v2, %L2 ]
 *    %p2 = phi @T [ %v1, %L1 ], [ %v3, %L2 ]
 *
 *    This shrink down the long `phi` chain generated by the passes above.
 */

struct MergePhi : public BasicBlockPass {
    static char ID;
    MergePhi() : BasicBlockPass(ID)
    {}

private:
    bool doInitialization(Function &F) override;
    bool runOnBasicBlock(BasicBlock &bb) override;

    bool checkPhiBool(PHINode *phi) const;
    bool processPhiBool(PHINode *phi, Instruction *first_non_phi) const;
    bool mergePhiSelect(BasicBlock &bb) const;

    bool checkPhiCmp(PHINode *phi) const;
    bool processPhiCmp(PHINode *phi, Instruction *first_non_phi) const;
    bool mergePhiCmp(BasicBlock &bb) const;

    bool mergePhiPhi(BasicBlock &bb) const;

    Type *T_bool;
    ConstantFolder m_folder;
};

bool MergePhi::doInitialization(Function &F)
{
    T_bool = Type::getInt1Ty(F.getContext());
    return false;
}

bool MergePhi::checkPhiBool(PHINode *phi) const
{
    if (phi->getType() != T_bool)
        return false;
    for (auto &op: phi->incoming_values()) {
        if (!isa<ConstantInt>(op.get())) {
            return false;
        }
    }
    return true;
}

bool MergePhi::processPhiBool(PHINode *phi, Instruction *first_non_phi) const
{
    bool changed = false;
    auto check_select_op = [&] (Value *op) {
        if (isa<Constant>(op) || isa<Argument>(op))
            return true;
        auto inst = dyn_cast<Instruction>(op);
        if (!inst)
            return false;
        if (inst->getParent() != phi->getParent())
            return true;
        if (isa<PHINode>(inst))
            return true;
        return false;
    };
    SmallVector<std::pair<SelectInst*,PHINode*>,16> replace;
    for (auto &use: phi->uses()) {
        auto select = dyn_cast<SelectInst>(use.getUser());
        if (!select)
            continue;
        if (use.getOperandNo() != 0)
            continue;
        if (select->getParent() != phi->getParent())
            continue;
        assert(select->getCondition() == phi);
        if (!check_select_op(select->getTrueValue()) ||
            !check_select_op(select->getFalseValue())) {
            continue;
        }
        changed = true;
        auto nincoming = phi->getNumIncomingValues();
        PHINode *newphi = PHINode::Create(select->getType(), nincoming,
                                          "", first_non_phi);
        for (unsigned i = 0; i < nincoming; i++) {
            auto val = cast<ConstantInt>(phi->getIncomingValue(i))->isZero() ?
                select->getFalseValue() : select->getTrueValue();
            newphi->addIncoming(val, phi->getIncomingBlock(i));
        }
        replace.emplace_back(select, newphi);
    }
    for (auto p: replace) {
        p.first->replaceAllUsesWith(p.second);
        p.first->eraseFromParent();
    }
    return changed;
}

bool MergePhi::mergePhiSelect(BasicBlock &bb) const
{
    SmallVector<PHINode*,16> phis;
    Instruction *first_non_phi = nullptr;
    for (auto &I: bb) {
        if (auto phi = dyn_cast<PHINode>(&I)) {
            if (checkPhiBool(phi))
                phis.push_back(phi);
            continue;
        }
        first_non_phi = &I;
        break;
    }
    if (phis.empty())
        return false;
    assert(first_non_phi);
    bool changed = false;
    for (auto phi: phis)
        changed = changed | processPhiBool(phi, first_non_phi);
    return changed;
}

bool MergePhi::checkPhiCmp(PHINode *phi) const
{
    for (auto &op: phi->incoming_values()) {
        if (!isa<Constant>(op.get())) {
            return false;
        }
    }
    return true;
}

bool MergePhi::processPhiCmp(PHINode *phi, Instruction *first_non_phi) const
{
    bool changed = false;
    SmallVector<std::pair<CmpInst*,PHINode*>,16> replace;
    for (auto &use: phi->uses()) {
        auto cmp = dyn_cast<CmpInst>(use.getUser());
        if (!cmp)
            continue;
        if (cmp->getParent() != phi->getParent())
            continue;
        const bool isfp = cmp->isFPPredicate();
        auto create_cmp = [&] (CmpInst::Predicate pred, Constant *op1, Constant *op2) {
            if (isfp) {
                return m_folder.CreateFCmp(pred, op1, op2);
            }
            return m_folder.CreateICmp(pred, op1, op2);
        };
        auto opno = use.getOperandNo();
        assert(opno == 0 || opno == 1);
        auto otherop = dyn_cast<Constant>(cmp->getOperand(1 - opno));
        if (!otherop)
            continue;
        changed = true;
        auto nincoming = phi->getNumIncomingValues();
        PHINode *newphi = PHINode::Create(T_bool, nincoming, "", first_non_phi);
        for (unsigned i = 0; i < nincoming; i++) {
            auto val = opno == 0 ?
                create_cmp(cmp->getPredicate(),
                           cast<Constant>(phi->getIncomingValue(i)), otherop) :
                create_cmp(cmp->getPredicate(), otherop,
                           cast<Constant>(phi->getIncomingValue(i)));
            newphi->addIncoming(val, phi->getIncomingBlock(i));
        }
        replace.emplace_back(cmp, newphi);
    }
    for (auto p: replace) {
        p.first->replaceAllUsesWith(p.second);
        p.first->eraseFromParent();
    }
    return changed;
}

bool MergePhi::mergePhiCmp(BasicBlock &bb) const
{
    SmallVector<PHINode*,16> phis;
    Instruction *first_non_phi = nullptr;
    for (auto &I: bb) {
        if (auto phi = dyn_cast<PHINode>(&I)) {
            if (checkPhiCmp(phi))
                phis.push_back(phi);
            continue;
        }
        first_non_phi = &I;
        break;
    }
    if (phis.empty())
        return false;
    assert(first_non_phi);
    bool changed = false;
    for (auto phi: phis)
        changed = changed | processPhiCmp(phi, first_non_phi);
    return changed;
}

bool MergePhi::mergePhiPhi(BasicBlock &bb) const
{
    bool changed = false;
    for (auto &I: bb) {
        auto phi = dyn_cast<PHINode>(&I);
        if (!phi)
            break;
        auto nincoming = phi->getNumIncomingValues();
        for (unsigned i = 0; i < nincoming; i++) {
            auto iv = dyn_cast<PHINode>(phi->getIncomingValue(i));
            if (!iv || iv->getParent() != phi->getParent())
                continue;
            phi->setIncomingValue(i, iv->getIncomingValueForBlock(phi->getIncomingBlock(i)));
        }
    }
    return changed;
}

bool MergePhi::runOnBasicBlock(BasicBlock &bb)
{
    bool changed = mergePhiCmp(bb);
    changed = mergePhiSelect(bb) | changed;
    changed = mergePhiPhi(bb) | changed;
    return changed;
}

char MergePhi::ID = 0;
static RegisterPass<MergePhi> X("MergePhi", "Merge phi node, select and cmp",
                                false /* Only looks at CFG */,
                                false /* Analysis Pass */);

Pass *createMergePhiPass()
{
    return new MergePhi();
}

}
}
