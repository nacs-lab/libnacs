/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "vectorize.h"

#include <llvm/ADT/SmallPtrSet.h>
#include <llvm/ADT/SmallSet.h>
#include <llvm/Analysis/ValueTracking.h>
#include <llvm/IR/BasicBlock.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/Transforms/Utils/ValueMapper.h>

namespace NaCs {
namespace LLVM {

namespace {

template<typename T>
static Function *createFunction(const Function &F, StringRef name, unsigned vec_size,
                                const T &vec_args, bool _export)
{
    auto orig_ft = F.getFunctionType();
    if (orig_ft->isVarArg())
        return nullptr;
    SmallVector<Type*, 8> fsig(orig_ft->param_begin(), orig_ft->param_end());
    for (auto arg: vec_args) {
        auto argt = fsig[arg];
        if (!VectorType::isValidElementType(argt))
            return nullptr;
        fsig[arg] = VectorType::get(argt, vec_size);
    }
    auto orig_rt = orig_ft->getReturnType();
    Type *rt;
    if (orig_rt->isVoidTy()) {
        rt = orig_rt;
    }
    else if (!VectorType::isValidElementType(orig_rt)) {
        return nullptr;
    }
    else {
        rt = VectorType::get(orig_rt, vec_size);
    }
    auto ftype = FunctionType::get(rt, fsig, false);

    auto f = Function::Create(ftype, _export ? GlobalValue::ExternalLinkage :
                              GlobalValue::PrivateLinkage, name,
                              const_cast<Module*>(F.getParent()));
    if (_export)
        f->setVisibility(GlobalValue::ProtectedVisibility);
    f->setAttributes(AttributeList::get(F.getContext(), F.getAttributes().getFnAttributes(),
                                        AttributeSet(), {}));
    f->addFnAttr(Attribute::AlwaysInline);
    return f;
}

static void trivialVectorize(const Function &F, Function *new_f, unsigned vec_size)
{
    auto orig_ft = F.getFunctionType();
    auto new_args = new_f->arg_begin();
    auto nargs = orig_ft->getNumParams();
    SmallVector<Value*, 8> fargs(nargs);
    for (unsigned i = 0; i < nargs; i++) {
        Value *new_arg = &*(new_args + i);
        auto argty = orig_ft->getParamType(i);
        fargs[i] = new_arg->getType() == argty ? new_arg : UndefValue::get(argty);
    }
    BasicBlock *b0 = BasicBlock::Create(new_f->getContext(), "top", new_f);
    IRBuilder<> builder(b0);
    auto orig_res = builder.CreateCall(const_cast<Function*>(&F), fargs);
    auto res = builder.CreateVectorSplat(vec_size, orig_res);
    builder.CreateRet(res);
}

static bool isTriviallyVectorizable(const Function &F, const SmallVectorImpl<unsigned> &vec_args)
{
    if (vec_args.empty())
        return true;
    auto args = F.arg_begin();
    for (auto i: vec_args) {
        if (!(args + i)->use_empty()) {
            return false;
        }
    }
    return true;
}

struct Vectorizer {
    Vectorizer(const Function &_F, const SmallVectorImpl<unsigned> &_vec_args)
        : F(_F),
          vec_args()
    {
        vec_args.insert(_vec_args.begin(), _vec_args.end());
    }

    Function *vectorize(StringRef name, unsigned vec_size, bool _export)
    {
        auto &entry = F.getEntryBlock();
        // FIXME: Single BB for now.
        if (!isa<ReturnInst>(entry.getTerminator()))
            return nullptr;
        if (!scan_bb(entry))
            return nullptr;

        auto new_f = createFunction(F, name, vec_size, vec_args, _export);
        BasicBlock *b0 = BasicBlock::Create(new_f->getContext(), "top", new_f);
        IRBuilder<> builder(b0);

        ValueToValueMapTy vmap;
        ValueMapper mapper(vmap);
        SmallPtrSet<Value*, 16> new_vec_vals{};

        auto arg_it = new_f->arg_begin();
        for (auto &arg: F.args()) {
            auto new_arg = &*(arg_it++);
            if (vec_args.count(arg.getArgNo()))
                new_vec_vals.insert(new_arg);
            vmap[&arg] = new_arg;
        }

        auto map_val =
            [&] (const Value *val, bool vec) -> Value* {
                auto it = vmap.find(val);
                bool found = it != vmap.end();
                assert(found || !isa<Instruction>(val));
                Value *new_val = found ? (Value*)it->second : const_cast<Value*>(val);
                if (!vec || new_vec_vals.count(new_val))
                    return new_val;
                return builder.CreateVectorSplat(vec_size, new_val);
            };

        for (auto &inst: entry) {
            if (inst.isTerminator()) {
                auto &ret = cast<ReturnInst>(inst);
                if (auto orig_val = ret.getReturnValue()) {
                    builder.CreateRet(map_val(orig_val, true));
                }
                else {
                    builder.CreateRetVoid();
                }
                break;
            }
            if (!vec_insts.count(&inst)) {
                auto new_inst = inst.clone();
                vmap[&inst] = new_inst;
                mapper.remapInstruction(*new_inst);
                builder.Insert(new_inst);
                continue;
            }
            auto add_vec_inst =
                [&] (Value *val) {
                    vmap[&inst] = val;
                    new_vec_vals.insert(val);
                };
            // Logic copied from `InnerLoopVectorizer::widenInstruction`
            // in llvm `LoopVectorize.cpp`
            if (auto *binop = dyn_cast<BinaryOperator>(&inst)) {
                auto *a = map_val(binop->getOperand(0), true);
                auto *b = map_val(binop->getOperand(1), true);
                auto *v = builder.CreateBinOp(binop->getOpcode(), a, b);
                if (auto *new_binop = dyn_cast<BinaryOperator>(v))
                    new_binop->copyIRFlags(binop);
                add_vec_inst(v);
            }
            else if (auto *cmp = dyn_cast<CmpInst>(&inst)) {
                Value *a = map_val(cmp->getOperand(0), true);
                Value *b = map_val(cmp->getOperand(1), true);
                Value *c = nullptr;
                if (inst.getOpcode() == Instruction::FCmp) {
                    // Propagate fast math flags.
                    IRBuilder<>::FastMathFlagGuard FMFG(builder);
                    builder.setFastMathFlags(cmp->getFastMathFlags());
                    c = builder.CreateFCmp(cmp->getPredicate(), a, b);
                }
                else {
                    assert(isa<ICmpInst>(cmp));
                    c = builder.CreateICmp(cmp->getPredicate(), a, b);
                }
                add_vec_inst(c);
            }
            else if (auto *gep = dyn_cast<GetElementPtrInst>(&inst)) {
                auto *ptr = map_val(gep->getPointerOperand(), false);
                SmallVector<Value*, 4> indices;
                for (auto &u: gep->indices())
                    indices.push_back(map_val(u.get(), false));
                auto elty = gep->getSourceElementType();
                add_vec_inst(gep->isInBounds() ?
                             builder.CreateInBoundsGEP(elty, ptr, indices) :
                             builder.CreateGEP(elty, ptr, indices));
            }
            else if (auto *select = dyn_cast<SelectInst>(&inst)) {
                auto *cond = map_val(select->getCondition(), false);
                Value *trueop = map_val(select->getTrueValue(), true);
                Value *falseop = map_val(select->getFalseValue(), true);
                add_vec_inst(builder.CreateSelect(cond, trueop, falseop));
            }
            else if (auto *cast = dyn_cast<CastInst>(&inst)) {
                Type *dest_ty = VectorType::get(cast->getType(), vec_size);
                Value *a = map_val(cast->getOperand(0), true);
                add_vec_inst(builder.CreateCast(cast->getOpcode(), a, dest_ty));
            }
            else {
                llvm_unreachable("Unhandled instruction!");
            }
        }

        return new_f;
    }

private:
    bool scan_bb(const BasicBlock &bb)
    {
        // FIXME: Ignoring phi node for now
        for (auto it = bb.getFirstNonPHI()->getIterator(); ; ++it) {
            auto *inst = &*it;
            // **This is the loop termination condition.**
            if (inst->isTerminator()) {
                // FIXME: only allow `ret` for now
                // TODO: br, switch
                if (isa<ReturnInst>(inst))
                    break;
                return false;
            }
            if (!isSafeToSpeculativelyExecute(inst))
                return false;
            bool has_vec = false;
            for (auto *op: inst->operand_values()) {
                if (auto arg = dyn_cast<Argument>(op)) {
                    if (vec_args.count(arg->getArgNo())) {
                        has_vec = true;
                        break;
                    }
                }
                else if (auto inst_op = dyn_cast<Instruction>(op)) {
                    if (vec_insts.count(inst_op)) {
                        has_vec = true;
                        break;
                    }
                }
            }
            if (!has_vec)
                continue;
            vec_insts.insert(inst);
            if (isa<BinaryOperator>(inst) || isa<ICmpInst>(inst) || isa<FCmpInst>(inst) ||
                isa<GetElementPtrInst>(inst) || isa<SelectInst>(inst) || isa<CastInst>(inst))
                continue;
#if LLVM_VERSION_MAJOR >= 8
            // Disable for now since IRBuilder doesn't have support for this yet...
            // // LLVM 8 have got its first unary operator
            // if (isa<UnaryOperator>(inst))
            //     continue;
#endif
            // FIXME: No call allowed for now.
            return false;
        }
        return true;
    }

    const Function &F;
    SmallSet<unsigned, 8> vec_args{};
    SmallPtrSet<const Instruction*, 16> vec_insts{};
};

}

Function *vectorizeFunction(const Function &F, StringRef name, unsigned vec_size,
                            const SmallVectorImpl<unsigned> &vec_args, bool _export)
{
    if (isTriviallyVectorizable(F, vec_args)) {
        auto new_f = createFunction(F, name, vec_size, vec_args, _export);
        if (!new_f)
            return nullptr;
        trivialVectorize(F, new_f, vec_size);
        return new_f;
    }
    Vectorizer vectorizer(F, vec_args);
    return vectorizer.vectorize(name, vec_size, _export);
}

}
}
