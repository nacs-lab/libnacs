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
 * Combine `phi` node that has only constant values with `select` that uses them.
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
    Constant *V_true;
    Constant *V_false;
    ConstantFolder m_folder;
};

bool MergePhi::doInitialization(Function &F)
{
    T_bool = Type::getInt1Ty(F.getContext());
    V_true = ConstantInt::get(T_bool, 1);
    V_false = ConstantInt::get(T_bool, 0);
    return false;
}

bool MergePhi::checkPhiBool(PHINode *phi) const
{
    if (phi->getType() != T_bool)
        return false;
    for (auto &op: phi->incoming_values()) {
        if (!isa<Constant>(op.get())) {
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
            auto val = m_folder.CreateSelect(cast<Constant>(phi->getIncomingValue(i)),
                                             cast<Constant>(select->getTrueValue()),
                                             cast<Constant>(select->getFalseValue()));
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
