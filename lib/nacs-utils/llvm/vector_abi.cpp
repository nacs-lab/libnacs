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

#define DEBUG_TYPE "vector_abi"

#include "passes.h"
#include "utils.h"

#if LLVM_VERSION_MAJOR >= 16
#  include <llvm/TargetParser/Triple.h>
#else
#  include <llvm/ADT/Triple.h>
#endif
#include <llvm/ADT/SmallVector.h>
#include <llvm/IR/CallingConv.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Module.h>

namespace NaCs::LLVM {

using namespace llvm;

// This pass fixes the ABI for functions that takes vector as argument/return value.
// It assumes that the function signature used in the IR are direct mapping of the C
// signature.
// We do this transformation here instead of in the frontend when the code are emitted
// to make sure that the code we emit there can be used for any target.

// 1. x86_64/i386, both return value and arguments
//    For windows use `x86_vectorcallcc` and pass as is.
//    For size less than 16 extend to 16 and pass as vector.
//    For longer vector, pass as is.
//    The transformation for smaller than 16 is to be consistent with sleef API/ABI.
// 2. Linux/apple aarch64
//    * Return value
//      Pass as is.
//    * Arguments
//      For size <= 4, ext and bitcast to i32
//      For longer vector (size <= 16), round size up to 8/16.

namespace {

#if LLVM_VERSION_MAJOR < 11
#  define FixedVectorType VectorType
#endif

struct VectorABI {
    static bool run_on_function(Function &F);
    static bool run_on_module(Module &M);

private:
    static bool handle_x86_func(const DataLayout &DL, FixedVectorType *T_v128,
                                Function &f, bool is_win);
    static bool handle_aarch64_func(const DataLayout &DL, IntegerType *T_i32,
                                    Function &f);

    // Rename the current function and set to private linkage
    // Create a new function with signature `new_fty` that inherits
    // the linkage of the old function and returns it.
    static Function *clone_to_api(Function &F, FunctionType *new_fty);
    static Value *set_nele_vector(IRBuilder<> &builder, Value *v, unsigned nele);
    static Value *cast_vector(IRBuilder<> &builder, FixedVectorType *ty, Value *v,
                              const DataLayout &DL);
    static bool is_vector_func(Function &F, bool export_only);
    static bool x86_abi(Module &M, bool is_win);
    static bool aarch64_abi(Module &M);
};

Function *VectorABI::clone_to_api(Function &F, FunctionType *new_fty)
{
    auto old_name = F.getName().str();
    auto lt = F.getLinkage();
    auto vis = F.getVisibility();
    auto dll = F.getDLLStorageClass();
    auto attr = AttributeList::get(F.getContext(), getFnAttrs(F),
                                   AttributeSet(), {});
    F.setVisibility(GlobalValue::DefaultVisibility);
    F.setDLLStorageClass(GlobalValue::DefaultStorageClass);
    F.setLinkage(GlobalValue::PrivateLinkage);
    F.setName(old_name + ".vec");
    F.removeFnAttr(Attribute::NoInline);
    F.addFnAttr(Attribute::AlwaysInline);
    auto newf = Function::Create(new_fty, lt, old_name, F.getParent());
    newf->setVisibility(vis);
    newf->setDLLStorageClass(dll);
    newf->setAttributes(attr);
    newf->setCallingConv(F.getCallingConv());
    return newf;
}

Value *VectorABI::set_nele_vector(IRBuilder<> &builder, Value *v, unsigned nele)
{
    auto ty = cast<FixedVectorType>(v->getType());
    auto old_nele = ty->getNumElements();
#if LLVM_VERSION_MAJOR >= 11
    SmallVector<int, 16> mask(nele);
    for (unsigned i = 0; i < nele; i++)
        mask[i] = int(i < old_nele ? i : old_nele);
#else
    SmallVector<uint32_t, 16> mask(nele);
    for (unsigned i = 0; i < nele; i++)
        mask[i] = i < old_nele ? i : old_nele;
#endif
    return builder.CreateShuffleVector(v, UndefValue::get(ty), mask);
}

Value *VectorABI::cast_vector(IRBuilder<> &builder, FixedVectorType *ty, Value *v,
                              const DataLayout &DL)
{
    auto old_ty = cast<FixedVectorType>(v->getType());
    auto old_elsz = DL.getTypeSizeInBits(old_ty->getElementType());
    auto old_sz = old_ty->getNumElements() * old_elsz;
    auto elsz = DL.getTypeSizeInBits(ty->getElementType());
    auto sz = ty->getNumElements() * elsz;
    if (sz % old_elsz == 0) {
        v = set_nele_vector(builder, v, sz / old_elsz);
        if (v->getType() != ty)
            v = builder.CreateBitCast(v, ty);
        return v;
    }
    else if (old_sz % elsz == 0) {
        v = builder.CreateBitCast(v, FixedVectorType::get(ty->getElementType(), old_sz / elsz));
        return set_nele_vector(builder, v, ty->getNumElements());
    }
    v = builder.CreateBitCast(v, FixedVectorType::get(Type::getInt8Ty(ty->getContext()),
                                                      old_sz / 8));
    v = set_nele_vector(builder, v, sz);
    return builder.CreateBitCast(v, ty);
}

// Ignore global alias.... We don't emit them and I don't want to think about them....
bool VectorABI::is_vector_func(Function &F, bool export_only)
{
    // LLVM knows how to deal with it so let it be...
    if (F.isIntrinsic())
        return false;
    // Already handled.
    if (F.hasFnAttribute("nacs.vector_abi"))
        return false;
    if (export_only) {
        auto link = F.getLinkage();
        // We use a list of known non-exported linkage here since
        // this is purely an optimization to reduce the work and it is always safe
        // to treat every function as exported.
        if (link == GlobalValue::PrivateLinkage || link == GlobalValue::InternalLinkage) {
            return false;
        }
    }
    auto ft = F.getFunctionType();
    // We don't support vararg function.
    // I don't expect to get them but if we do, just ignore them instead of generating some
    // invalid IR...
    if (ft->isVarArg())
        return false;
    if (ft->getReturnType()->isVectorTy())
        return true;
    for (auto t: ft->params()) {
        if (t->isVectorTy()) {
            return true;
        }
    }
    return false;
}

static bool fix_call_cc(Function &f)
{
    bool changed = false;
    for (const auto &use: f.uses()) {
        // We can replace this with CallBase on later version of LLVM.
        // We don't use any other kind of call instructions though so it doesn't
        // really matter...
        if (auto inst = dyn_cast<CallInst>(use.getUser())) {
            if (inst->getCalledFunction() == &f &&
                inst->getCallingConv() != CallingConv::X86_VectorCall) {
                inst->setCallingConv(CallingConv::X86_VectorCall);
                changed = true;
            }
        }
    }
    return changed;
}

bool VectorABI::handle_x86_func(const DataLayout &DL, FixedVectorType *T_v128,
                                Function &f, bool is_win32)
{
    // On win32 we need to fix the caller, which may be created
    // after the function itself has been fixed.
    // Therefore we need to apply this fix again even if we see a processed function.
    if (is_win32 && f.getCallingConv() == CallingConv::X86_VectorCall)
        return fix_call_cc(f);
    if (!is_vector_func(f, !is_win32))
        return false;
    if (is_win32) {
        // Mark function processed.
        f.addFnAttr("nacs.vector_abi", "");
        f.setCallingConv(CallingConv::X86_VectorCall);
        fix_call_cc(f);
        auto link = f.getLinkage();
        if (link == GlobalValue::PrivateLinkage ||
            link == GlobalValue::InternalLinkage) {
            return true;
        }
    }
    // For size less than 16 extend to 16 and pass as vector.
    // For longer vector, pass as is.
    auto ft = f.getFunctionType();
    Type *ret = ft->getReturnType();
    bool fix_ret = false;
    SmallVector<Type*, 8> argts;
    SmallVector<unsigned, 8> tofix;
    if (ret->isVectorTy()) {
        if (DL.getTypeSizeInBits(ret) < 128) {
            fix_ret = true;
            ret = T_v128;
        }
    }
    for (auto t: ft->params()) {
        if (!t->isVectorTy() || DL.getTypeSizeInBits(t) >= 128) {
            argts.push_back(t);
            continue;
        }
        tofix.push_back(argts.size());
        argts.push_back(T_v128);
    }
    if (tofix.empty() && !fix_ret)
        return is_win32;
    // Mark function processed. Already marked on windows.
    if (!is_win32)
        f.addFnAttr("nacs.vector_abi", "");
    auto new_fty = FunctionType::get(ret, argts, false);
    auto newf = clone_to_api(f, new_fty);
    if (f.empty()) {
        // Declaration of external vector ABI functions
        BasicBlock *b0 = BasicBlock::Create(f.getContext(), "top", &f);
        IRBuilder<> builder(b0);
        SmallVector<Value*, 8> args(argts.size());
        uint32_t fix_i = 0;
        for (auto &arg: f.args()) {
            auto argno = arg.getArgNo();
            if (fix_i >= tofix.size() || tofix[fix_i] != argno) {
                args[argno] = &arg;
                continue;
            }
            args[argno] = cast_vector(builder, T_v128, &arg, DL);
            fix_i++;
        }
        Value *res = builder.CreateCall(newf, args);
        if (fix_ret)
            res = cast_vector(builder, cast<FixedVectorType>(ft->getReturnType()), res, DL);
        if (res->getType()->isVoidTy()) {
            builder.CreateRetVoid();
        }
        else {
            builder.CreateRet(res);
        }
        return true;
    }
    BasicBlock *b0 = BasicBlock::Create(newf->getContext(), "top", newf);
    IRBuilder<> builder(b0);
    SmallVector<Value*, 8> args(argts.size());
    uint32_t fix_i = 0;
    for (auto &arg: newf->args()) {
        auto argno = arg.getArgNo();
        if (fix_i >= tofix.size() || tofix[fix_i] != argno) {
            args[argno] = &arg;
            continue;
        }
        args[argno] = cast_vector(builder, cast<FixedVectorType>(ft->getParamType(argno)),
                                  &arg, DL);
        fix_i++;
    }
    Value *res = builder.CreateCall(&f, args);
    if (fix_ret)
        res = cast_vector(builder, T_v128, res, DL);
    if (res->getType()->isVoidTy()) {
        builder.CreateRetVoid();
    }
    else {
        builder.CreateRet(res);
    }
    return true;
}

bool VectorABI::x86_abi(Module &M, bool is_win32)
{
    auto &DL = M.getDataLayout();
    auto T_v128 = FixedVectorType::get(Type::getInt32Ty(M.getContext()), 4);
    bool changed = false;
    for (auto &f: M)
        changed |= handle_x86_func(DL, T_v128, f, is_win32);
    return changed;
}

bool VectorABI::handle_aarch64_func(const DataLayout &DL, IntegerType *T_i32,
                                    Function &f)
{
    if (!is_vector_func(f, true))
        return false;
    // For size <= 4, ext and bitcast to i32
    // For longer vector (size <= 16), round size up to 8/16.
    // This is different from the calling convention used by gcc
    // since it uses the GP registers to pass small vector types
    // (at least for <32bit).
    // However, clang passes these in the vector registers
    // so using them directly in C declaration causes a mismatch between
    // gcc and clang code.
    // AAPCS64 does not seem to cover this explicitly
    // (the vector types mentioned in it has a minimum size of 8 bytes)
    // so I'm not sure which one is more correct.
    // Since GCC can be convinced to pass things through the vector registers
    // by extending it to a larger type, but it's harder to convince clang
    // to pass vector in GP registers, let's just use the vector registers
    // to pass everything.
    auto ft = f.getFunctionType();
    SmallVector<Type*, 8> argts;
    SmallVector<unsigned, 8> tofix;
    auto to_vector_abi = [&] (Type *t) -> Type* {
        if (!t->isVectorTy())
            return t;
        auto argsz = DL.getTypeSizeInBits(t);
        // Check if big enough & power of two
        if (argsz < 64 || (argsz & (argsz - 1)) != 0) {
            // Now find the next power of 2
            auto sz = 1 << (32 - __builtin_clz(argsz));
            if (sz < 64)
                sz = 64;
            auto ele = cast<FixedVectorType>(t)->getElementType();
            auto elesz = DL.getTypeSizeInBits(ele);
            if (sz % elesz == 0) {
                return FixedVectorType::get(ele, sz / elesz);
            }
            else {
                return FixedVectorType::get(Type::getInt8Ty(ele->getContext()),
                                            sz / 8);
            }
        }
        return t;
    };
    Type *orig_ret = ft->getReturnType();
    Type *ret = to_vector_abi(orig_ret);
    for (auto t: ft->params()) {
        auto newt = to_vector_abi(t);
        if (newt != t)
            tofix.push_back(argts.size());
        argts.push_back(newt);
    }
    if (tofix.empty() && ret == ft->getReturnType())
        return false;
    // Mark function processed.
    f.addFnAttr("nacs.vector_abi", "");
    auto new_fty = FunctionType::get(ret, argts, false);
    auto newf = clone_to_api(f, new_fty);
    if (f.empty()) {
        // Declaration of external vector ABI functions
        BasicBlock *b0 = BasicBlock::Create(f.getContext(), "top", &f);
        IRBuilder<> builder(b0);
        SmallVector<Value*, 8> args(argts.size());
        uint32_t fix_i = 0;
        for (auto &arg: f.args()) {
            auto argno = arg.getArgNo();
            if (fix_i >= tofix.size() || tofix[fix_i] != argno) {
                args[argno] = &arg;
                continue;
            }
            args[argno] = cast_vector(builder, cast<FixedVectorType>(argts[argno]), &arg, DL);
            fix_i++;
        }
        Value *res = builder.CreateCall(newf, args);
        auto res_ty = res->getType();
        if (res_ty->isVoidTy()) {
            builder.CreateRetVoid();
        }
        else if (res_ty == orig_ret) {
            builder.CreateRet(res);
        }
        else {
            builder.CreateRet(cast_vector(builder, cast<FixedVectorType>(orig_ret),
                                          res, DL));
        }
        return true;
    }
    BasicBlock *b0 = BasicBlock::Create(newf->getContext(), "top", newf);
    IRBuilder<> builder(b0);
    SmallVector<Value*, 8> args(argts.size());
    uint32_t fix_i = 0;
    for (auto &arg: newf->args()) {
        auto argno = arg.getArgNo();
        if (fix_i >= tofix.size() || tofix[fix_i] != argno) {
            args[argno] = &arg;
            continue;
        }
        args[argno] = cast_vector(builder, cast<FixedVectorType>(ft->getParamType(argno)),
                                  &arg, DL);
        fix_i++;
    }
    Value *res = builder.CreateCall(&f, args);
    auto res_ty = res->getType();
    if (res->getType()->isVoidTy()) {
        builder.CreateRetVoid();
    }
    else if (res_ty == ret) {
        builder.CreateRet(res);
    }
    else {
        builder.CreateRet(cast_vector(builder, cast<FixedVectorType>(ret), res, DL));
    }
    return true;
}

bool VectorABI::aarch64_abi(Module &M)
{
    auto &DL = M.getDataLayout();
    auto T_i32 = Type::getInt32Ty(M.getContext());
    bool changed = false;
    for (auto &f: M)
        changed |= handle_aarch64_func(DL, T_i32, f);
    return changed;
}

bool VectorABI::run_on_function(Function &F)
{
    auto *M = F.getParent();
    assert(M && "Function must have a parent");
    Triple triple(M->getTargetTriple());
    auto &DL = M->getDataLayout();
    switch (triple.getArch()) {
    case Triple::x86:
    case Triple::x86_64:
        return handle_x86_func(DL, FixedVectorType::get(Type::getInt32Ty(M->getContext()), 4),
                               F, triple.getOS() == Triple::Win32);
    case Triple::aarch64:
    case Triple::aarch64_be:
        return handle_aarch64_func(DL, Type::getInt32Ty(M->getContext()), F);
    default:
        return false;
    }
}

bool VectorABI::run_on_module(Module &M)
{
    Triple triple(M.getTargetTriple());
    switch (triple.getArch()) {
    case Triple::x86:
    case Triple::x86_64:
        return x86_abi(M, triple.getOS() == Triple::Win32);
    case Triple::aarch64:
    case Triple::aarch64_be:
        return aarch64_abi(M);
    default:
        return false;
    }
}

#if NACS_ENABLE_LEGACY_PASS
struct LegacyVectorABIPass : public ModulePass {
    static char ID;
    LegacyVectorABIPass()
        : ModulePass(ID)
    {}

private:
    bool runOnModule(Module &M) override
    {
        return VectorABI::run_on_module(M);
    }
};

char LegacyVectorABIPass::ID = 0;
static RegisterPass<LegacyVectorABIPass> X("VectorABI", "Fixing Vector ABI Pass",
                                           false /* Only looks at CFG */,
                                           false /* Analysis Pass */);
#endif

}

NACS_EXPORT() bool fixVectorABI(Function &F)
{
    return VectorABI::run_on_function(F);
}

#if NACS_ENABLE_LEGACY_PASS
NACS_EXPORT() Pass *createVectorABIPass()
{
    return new LegacyVectorABIPass();
}
#endif

#if NACS_ENABLE_NEW_PASS
NACS_EXPORT() PreservedAnalyses
VectorABIPass::run(Module &M, ModuleAnalysisManager &AM)
{
    if (VectorABI::run_on_module(M))
        return PreservedAnalyses::allInSet<CFGAnalyses>();
    return PreservedAnalyses::all();
}
#endif

}
