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

#include <llvm/IR/BasicBlock.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/IRBuilder.h>

namespace NaCs {
namespace LLVM {

namespace {

static Function *createFunction(const Function &F, StringRef name, unsigned vec_size,
                                const SmallVectorImpl<unsigned> &vec_args, bool _export)
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
    if (!VectorType::isValidElementType(orig_rt))
        return nullptr;
    auto rt = VectorType::get(orig_rt, vec_size);
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
    return nullptr;
}

}
}
