/*************************************************************************
 *   Copyright (c) 2019 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "env.h"

#include "../nacs-utils/llvm/analysis.h"
#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/llvm/compile.h"
#include "../nacs-utils/llvm/execute.h"
#include "../nacs-utils/llvm/inst_simplify.h"
#include "../nacs-utils/llvm/utils.h"

#include <stdexcept>

#include <llvm/IR/Constants.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Utils/Cloning.h>

namespace NaCs::Seq {

NACS_EXPORT() IR::Type Var::type() const
{
    if (is_const())
        return get_const().typ;
    if (is_extern())
        return get_extern().first;
    assert(is_call());
    auto f = get_callee();
    if (!f.is_llvm)
        return f.var->type();
    return LLVM::get_ir_type(f.llvm->getReturnType());
}

NACS_EXPORT() llvm::Type *Var::llvm_type() const
{
    if (is_call()) {
        auto f = get_callee();
        if (!f.is_llvm)
            return f.var->llvm_type();
        return f.llvm->getReturnType();
    }
    assert(is_const() || is_extern());
    auto ir_type = is_const() ? get_const().typ : get_extern().first;
    return m_env.cg_context()->llvm_ty(ir_type);
}

NACS_INTERNAL void Var::fill_args(llvm::ArrayRef<Arg> args, int nfreeargs)
{
    size_t n = args.size();
    m_args.resize(n);
    for (size_t i = 0; i < n; i++)
        m_args[i] = args[i];
    m_n_freeargs = nfreeargs;
}

void Var::assign_call(Var *func, llvm::ArrayRef<Arg> args, int nfreeargs)
{
    fill_args(args, nfreeargs);
    m_is_const = false;
    m_is_extern = false;
    m_value.func.is_llvm = false;
    m_value.func.var = func;
}

void Var::assign_call(llvm::Function *func, llvm::ArrayRef<Arg> args, int nfreeargs)
{
    fill_args(args, nfreeargs);
    m_is_const = false;
    m_is_extern = false;
    m_value.func.is_llvm = true;
    m_value.func.llvm = func;
}

NACS_EXPORT() std::ostream &operator<<(std::ostream &stm, const Arg &arg)
{
    if (arg.is_const()) {
        stm << arg.get_const();
    }
    else if (arg.is_var()) {
        stm << '%' << arg.get_var()->varid();
    }
    else if (arg.is_arg()) {
        stm << "arg[" << arg.get_arg() << ']';
    }
    else {
        stm << "<undef>";
    }
    return stm;
}

NACS_EXPORT() void Var::print(std::ostream &stm, bool newline) const
{
    auto id = varid();
    if (m_extern_ref > 0)
        stm << '*';
    if (id < 0) {
        stm << "%<invalid_ref>";
        if (newline)
            stm << std::endl;
        return;
    }
    stm << '%' << id;
    if (nfreeargs()) {
        stm << '(';
        for (int i = 0; i < nfreeargs(); i++) {
            if (i != 0)
                stm << ", ";
            stm << "arg[" << i << ']';
        }
        stm << ')';
    }
    stm << " = ";
    if (is_const()) {
        stm << get_const();
        if (newline)
            stm << std::endl;
        return;
    }
    else if (is_extern()) {
        stm << get_extern().first << " extern(0x"
            << std::hex << get_extern().second << ')' << std::dec;
        if (newline)
            stm << std::endl;
        return;
    }
    auto f = get_callee();
    stm << type();
    if (f.is_llvm) {
        stm << " llvm:";
        llvm::raw_os_ostream lstm(stm);
        f.llvm->printAsOperand(lstm, false);
    }
    else {
        stm << " %" << f.var->varid();
    }
    stm << '(';
    auto nargs = args().size();
    for (size_t i = 0; i < nargs; i++) {
        if (i != 0)
            stm << ", ";
        stm << args()[i];
    }
    stm << ')';
    if (newline) {
        stm << std::endl;
    }
}

NACS_EXPORT() bool Var::argument_unused(int idx) const
{
    for (auto &arg: args()) {
        if (arg.is_arg() && arg.get_arg() == idx) {
            return false;
        }
    }
    return true;
}

bool Var::inline_callee()
{
    assert(is_call());
    auto f = get_callee();
    bool changed = false;
    if (f.is_llvm) {
        if (auto c = LLVM::Analysis::returns_const(*f.llvm)) {
            if (auto ci = llvm::dyn_cast<llvm::ConstantInt>(c)) {
                if (ci->getBitWidth() == 1) {
                    // bool
                    assign_const(ci->getZExtValue() != 0);
                }
                else if (ci->getBitWidth() <= 32) {
                    assign_const((int32_t)ci->getSExtValue());
                }
                else {
                    assign_const(ci->getValue().signedRoundToDouble());
                }
                return true;
            }
            else if (auto cf = llvm::dyn_cast<llvm::ConstantFP>(c)) {
                assign_const(cf->getValueAPF().convertToDouble());
                return true;
            }
        }
        else if (auto argl = LLVM::Analysis::returns_argument(*f.llvm)) {
            auto arg = args()[argl->getArgNo()];
            if (arg.is_const()) {
                // There might be an argument conversion.
                assign_const(arg.get_const().convert(type()));
                return true;
            }
            else if (arg.is_var()) {
                assign_var(arg.get_var());
                f = get_callee();
                changed = true;
                goto non_llvm;
            }
            // After a single round of optimization there should be no caller
            // of this function anymore so a forward of argument isn't very useful
            // for optimization (since such info is only useful for the caller).
            // The original caller of this (if there were any)
            // will get the same information after inlining the LLVM function.
            // The only caller will be external user. If necessary,
            // a check can be added by the external user to optimize for this case.
        }
        return false;
    }
non_llvm:
    if (auto copy = f.var->get_assigned_var()) {
        f.var = copy;
        changed = true;
    }
    if (f.var->is_const()) {
        assign_const(f.var->get_const());
        return true;
    }
    else if (f.var->is_extern()) {
        // Maintain the pointer identity of the `Var*` with external parameter.
        if (!args().empty() || nfreeargs() > 0) {
            assign_var(f.var);
            return true;
        }
        return changed;
    }
    assert(f.var->is_call());
    // Now we should actually inline the called function into the caller so that
    // we can directly deal with calling LLVM functions.
    // We don't want this to result in repeated calculation though, so if the
    // callee is a standalone computable variable (i.e. has no free arguments)
    // we'll just make this a copy of that
    // (the user might eliminate the need for this variable later).
    if (f.var->nfreeargs() == 0) {
        changed |= !args().empty();
        assign_var(f.var);
        return changed;
    }
    // The callee has free arguments, we'll just make a copy of the llvm function
    // and inline the arguments.
    m_value.func = f.var->get_callee();
    auto new_args = f.var->m_args;
    for (auto &arg: new_args) {
        if (!arg.is_arg())
            continue;
        auto argi = arg.get_arg();
        if (argi >= (ssize_t)args().size())
            throw std::out_of_range("Missing argument.");
        arg = args()[argi];
    }
    m_args = new_args;
    return true;
}

static llvm::Type *merge_type(llvm::Type *t1, llvm::Type *t2)
{
    if (t1 == t2)
        return t1;
    if (t1->isIntegerTy()) {
        if (t2->isIntegerTy()) {
            return (t1->getPrimitiveSizeInBits() > t2->getPrimitiveSizeInBits() ?
                    t1 : t2);
        }
        else if (t2->isFloatingPointTy()) {
            return t2;
        }
    }
    else if (t1->isFloatingPointTy()) {
        if (t2->isIntegerTy()) {
            return t1;
        }
        else if (t2->isFloatingPointTy()) {
            return (t1->getPrimitiveSizeInBits() > t2->getPrimitiveSizeInBits() ?
                    t1 : t2);
        }
    }
    throw std::invalid_argument("Argument type mismatch.");
}

NACS_INTERNAL void Var::optimize_llvmf(llvm::Function *f)
{
    llvm::legacy::FunctionPassManager FPM(f->getParent());
#ifndef NDEBUG
    FPM.add(llvm::createVerifierPass());
#endif
    FPM.add(llvm::createSCCPPass());
    FPM.add(LLVM::createNaCsInstSimplifyPass(m_env.cg_context()->get_extern_resolver()));
    FPM.add(llvm::createCFGSimplificationPass());
    FPM.add(llvm::createEarlyCSEPass());
    FPM.add(llvm::createCFGSimplificationPass());
#ifndef NDEBUG
    FPM.add(llvm::createVerifierPass());
#endif

    FPM.doInitialization();
    FPM.run(*f);
    FPM.doFinalization();
}

NACS_INTERNAL bool Var::reduce_args()
{
    assert(m_in_env && !m_env.m_varid_dirty);
    assert(get_callee().is_llvm);
    // Collect information about arguments
    struct ArgsInfo {
        bool has_const{false};
        std::map<int,size_t> arg_idxs;
        std::map<Var*,size_t> var_idxs;
        llvm::SmallVector<std::pair<size_t,ssize_t>,8> skip;
    } info;
    size_t newi = 0;
    size_t nargs = args().size();
    for (size_t i = 0; i < nargs; i++, newi++) {
        auto &arg = args()[i];
        if (arg.is_var()) {
            auto argv = arg.get_var();
            auto it = info.var_idxs.find(argv);
            if (it != info.var_idxs.end()) {
                info.skip.emplace_back(i, it->second);
                newi--;
                continue;
            }
            info.var_idxs.emplace(argv, newi);
        }
        else if (arg.is_arg()) {
            auto argi = arg.get_arg();
            auto it = info.arg_idxs.find(argi);
            if (it != info.arg_idxs.end()) {
                info.skip.emplace_back(i, it->second);
                newi--;
                continue;
            }
            info.arg_idxs.emplace(argi, newi);
        }
        else if (arg.is_const()) {
            info.has_const = true;
            info.skip.emplace_back(i, -1);
            newi--;
        }
    }
    if (info.skip.empty() && !info.has_const)
        return false;
    auto oldf = get_callee().llvm;
    auto oldfty = oldf->getFunctionType();
    assert(!oldfty->isVarArg());
    assert(oldfty->getNumParams() == nargs);
    auto newnargs = newi;
    auto rt = oldfty->getReturnType();
    llvm::SmallVector<llvm::Type*, 8> fsig(newnargs);
    {
        size_t i = 0;
        size_t newi = 0;
        for (auto skip: info.skip) {
            for (; i < skip.first; i++, newi++)
                fsig[newi] = oldfty->getParamType(i);
            assert(skip.second < (ssize_t)i);
            if (skip.second >= 0)
                fsig[skip.second] = merge_type(fsig[skip.second], oldfty->getParamType(i));
            i++;
        }
        for (; i < nargs; i++, newi++) {
            fsig[newi] = oldfty->getParamType(i);
        }
    }
    auto fty = llvm::FunctionType::get(rt, fsig, false);
    auto llvm_mod = m_env.llvm_module();
    auto &llvm_ctx = llvm_mod->getContext();
    auto &cgctx = *m_env.cg_context();
    auto f = llvm::Function::Create(fty, llvm::GlobalValue::ExternalLinkage,
                                    oldf->getName() + ".l", llvm_mod);
    f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
    f->setAttributes(llvm::AttributeList::get(llvm_ctx,
                                              oldf->getAttributes().getFnAttributes(),
                                              {}, {}));
    auto *b0 = llvm::BasicBlock::Create(llvm_ctx, "top", f);
    llvm::IRBuilder<> builder(b0);
    llvm::SmallVector<llvm::Value*, 8> call_args(nargs);
    decltype(m_args) newargs(newnargs);
    try {
        auto argit = f->arg_begin();
        size_t i = 0;
        size_t newi = 0;
        for (auto skip: info.skip) {
            for (; i < skip.first; i++, newi++) {
                newargs[newi] = args()[i];
                call_args[i] = &*(argit + newi);
            }
            assert(skip.second < (ssize_t)i);
            if (skip.second < 0) {
                call_args[i] = cgctx.emit_const(args()[i].get_const());
            }
            else {
                call_args[i] = &*(argit + skip.second);
            }
            i++;
        }
        for (; i < nargs; i++, newi++) {
            newargs[newi] = args()[i];
            call_args[i] = &*(argit + newi);
        }
        for (size_t i = 0; i < nargs; i++) {
            // Hard code this for now (assume T_i8 are booleans)
            auto ty = oldfty->getParamType(i);
            auto v = call_args[i];
            if (ty == cgctx.T_i8)
                v = LLVM::convert_scalar(builder, cgctx.T_bool, v);
            call_args[i] = LLVM::convert_scalar(builder, ty, v);
        }
    }
    catch (...) {
        builder.ClearInsertionPoint();
        f->eraseFromParent();
        throw;
    }
    auto res = builder.CreateCall(oldf, call_args);
    builder.CreateRet(res);
    llvm::InlineFunctionInfo IFI;
#if LLVM_VERSION_MAJOR >= 11
    llvm::InlineFunction(*res, IFI);
#else
    llvm::InlineFunction(res, IFI);
#endif

    optimize_llvmf(f);

    assign_call(f, newargs, nfreeargs());
    return true;
}

bool Var::optimize_call()
{
    assert(is_call());
    // The caller should inline the callee before this
    if (!get_callee().is_llvm)
        return false;
    bool changed = false;
    if (!args().empty()) {
        changed |= reduce_args();
    }
    else {
        optimize_llvmf(get_callee().llvm);
    }
    changed |= inline_callee();
    // Optimized to constant or copy of another variable.
    if (!is_call() || !get_callee().is_llvm)
        return changed;
    if (!args().empty())
        return changed;
    // Zero arguments: try compiling and running the function.
    auto f = get_callee().llvm;
    // We can run the function without side effects.
    if (!f->isSpeculatable() && !(f->doesNotThrow() && f->onlyReadsMemory()))
        return changed;
    // Make sure that we understand the type
    IR::Type typ = LLVM::get_ir_type(f->getReturnType());
    if (!uint8_t(typ))
        return changed;
    auto llvm_mod = m_env.llvm_module();
    auto &llvm_ctx = llvm_mod->getContext();
    llvm::Module temp_mod("", llvm_ctx);
    // Copy the function and its dependencies to a temporary module for compilation.
    LLVM::FunctionMover mover(&temp_mod);
    auto newf = mover.clone_function(f);
    // The passes below (in `emit_objfile`) may recreate the functions
    // so the function handle may not be valid anymore.
    auto new_name = newf->getName().str();
    llvm::SmallVector<char,0> vec;
    // Emit object file, no need to optimize since we are only calling it once.
    auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(),
                                           &temp_mod, false);
    // Allow compilation failure here
    if (!res)
        return changed;
    LLVM::Exe::Engine engine;
    auto obj_id = engine.load(&vec[0], vec.size(), m_env.cg_context()->get_extern_resolver());
    if (!obj_id)
        return changed;
    auto sym = engine.get_symbol(new_name);
    if (!sym) {
        engine.free(obj_id);
        return changed;
    }
    IR::TagVal val;
    switch (typ) {
    case IR::Type::Bool:
        val = ((bool(*)())sym)();
        break;
    case IR::Type::Int32:
        val = ((int(*)())sym)();
        break;
    case IR::Type::Float64:
        val = ((double(*)())sym)();
        break;
    default:
        return changed;
    }
    engine.free(obj_id);
    assign_const(val);
    return true;
}

}
