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
#include "../nacs-utils/llvm/global_rename.h"
#include "../nacs-utils/llvm/utils.h"
#include "../nacs-utils/streams.h"

#include <llvm/ADT/StringMap.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/Utils/Cloning.h>

#include <vector>

namespace NaCs::Seq {

namespace {
struct DefaultCGContext: LLVM::Codegen::CachedContext {
    DefaultCGContext(llvm::LLVMContext &llvm_ctx)
        : LLVM::Codegen::CachedContext(llvm_ctx)
    {
    }
private:
    bool use_extern_data() const override
    {
        return true;
    }
    void add_extern_data(llvm::StringRef name, const void *_data, size_t size) override
    {
        auto data = (const double*)_data;
        m_data[name] = std::vector<double>{data, data + size / sizeof(double)};
    }
    uintptr_t get_extern_data(llvm::StringRef name) override
    {
        auto it = m_data.find(name);
        if (it != m_data.end())
            return (uintptr_t)&it->second[0];
        return 0;
    }
    llvm::StringMap<std::vector<double>> m_data;
};
}

NACS_EXPORT_ Env::Env(std::unique_ptr<LLVM::Codegen::Context> cgctx)
    : m_llvm_mod(new llvm::Module("", cgctx->get_context())),
      m_cgctx(std::move(cgctx))
{
    m_cgctx->set_module(m_llvm_mod.get());
}

NACS_EXPORT_ Env::Env(llvm::LLVMContext &llvm_ctx)
    : Env(std::make_unique<DefaultCGContext>(llvm_ctx))
{
}

NACS_EXPORT() Env::~Env()
{
    while (m_vars) {
        m_vars->release();
    }
}

NACS_INTERNAL Var *Env::new_var()
{
    auto var = new Var(*this);
    if (m_vars)
        m_vars->m_prev = &var->m_next;
    var->m_prev = &m_vars;
    var->m_next = m_vars;
    m_vars = var;
    m_varid_dirty = true;
    return var;
}

NACS_EXPORT() Var *Env::new_const(IR::TagVal c)
{
    auto var = new_var();
    var->assign_const(c);
    return var;
}

static void check_args(llvm::ArrayRef<Arg> args, int nfreeargs)
{
    for (auto &arg: args) {
        if (arg.is_arg() && arg.get_arg() >= nfreeargs) {
            throw std::out_of_range("Argument index out of bound.");
        }
        else if (arg.is_var() && arg.get_var()->nfreeargs() > 0) {
            throw std::invalid_argument("Argument is not fixed.");
        }
    }
}

NACS_EXPORT() Var *Env::new_call(Var *func, llvm::ArrayRef<Arg> args, int nfreeargs)
{
    // Allow the callee to have 0 argument
    // to be consistent with what we allow in optimiztion
    if (func->nfreeargs() > 0 && func->nfreeargs() != (ssize_t)args.size())
        throw std::invalid_argument("Argument number mismatch");
    check_args(args, nfreeargs);
    // There are significantly fewer things to check since we assume the callee
    // is already verified.
    auto var = new_var();
    var->assign_call(func, args, nfreeargs);
    return var;
}

static bool check_type(llvm::Type *ty)
{
    return uint8_t(LLVM::get_ir_type(ty)) > 0;
}

NACS_INTERNAL Var *Env::_new_call(llvm::Function *func, llvm::ArrayRef<Arg> args,
                                  int nfreeargs)
{
    if (func->arg_size() != args.size())
        throw std::invalid_argument("Argument number mismatch");
    check_args(args, nfreeargs);
    auto var = new_var();
    var->assign_call(func, args, nfreeargs);
    return var;
}

NACS_EXPORT() Var *Env::new_call(llvm::Function *func, llvm::ArrayRef<Arg> args,
                                 int nfreeargs)
{
    auto fty = func->getFunctionType();
    if (!check_type(fty->getReturnType()))
        throw std::invalid_argument("Invalid return type.");
    if (func->isVarArg())
        throw std::invalid_argument("Vararg function not supported.");
    for (auto t: fty->params()) {
        if (!check_type(t)) {
            throw std::invalid_argument("Invalid argumentt type.");
        }
    }
    return _new_call(func, args, nfreeargs);
}

NACS_EXPORT() Var *Env::new_call(const IR::Function &func, llvm::ArrayRef<Arg> args,
                                 int nfreeargs)
{
    return _new_call(cg_context()->emit_function(func, "f"), args, nfreeargs);
}

NACS_EXPORT() Var *Env::new_extern(std::pair<IR::Type,uint64_t> ext)
{
    auto var = new_var();
    var->assign_extern(ext);
    return var;
}

NACS_EXPORT() void Env::compute_varid()
{
    int nvars = 0;
    for (auto var NACS_UNUSED: *this)
        nvars++;
    for (auto var: *this)
        var->m_varid = --nvars;
    m_varid_dirty = false;
    assert(nvars == 0);
}

NACS_EXPORT() void Env::compute_varuse()
{
    // Collect all roots.
    llvm::SmallVector<Var*, 16> roots;
    for (auto var: *this) {
        var->m_used = false;
        if (var->m_extern_ref > 0) {
            roots.push_back(var);
        }
    }
    auto ref_empty = [] (Var *var) {
        if (!var->args().empty())
            return false;
        if (!var->is_call())
            return true;
        return var->get_callee().is_llvm;
    };
    auto get_ref = [] (Var *var, ssize_t idx) -> Var* {
        if (idx < 0) {
            if (var->is_call()) {
                auto f = var->get_callee();
                if (!f.is_llvm) {
                    return assume(f.var);
                }
            }
            return nullptr;
        }
        auto arg = var->args()[idx];
        return arg.is_var() ? arg.get_var() : nullptr;
    };

    // Mark starting from roots
    // Use a manual stack for better optimization and less actual stack usage so that
    // we don't need to worry about stack overflow.
    // All items on the stack is inbound (see `push` below) `pop` does not need bounds check
    llvm::SmallVector<std::pair<Var*,size_t>, 16> stack;
    for (auto var: roots) {
        if (ref_empty(var))
            continue;
        ssize_t idx = -1;
        auto next = [&] {
            // Next argument
            if (++idx < (ssize_t)var->args().size())
                return true;
            // Done
            if (stack.empty())
                return false;
            std::tie(var, idx) = stack.pop_back_val();
            return true;
        };
        auto iterate = [&] {
            auto new_var = get_ref(var, idx);
            if (!new_var || new_var->m_used)
                return next();
            new_var->m_used = true;
            // If this is a externally and internally referenced variable,
            // it'll be scanned from the root so we can stop here.
            if (new_var->m_extern_ref > 0 || ref_empty(new_var))
                return next();
            // Save the next one to be processed to the stack
            // If we are already at the last one, we don't need to do anything.
            if (idx + 1 < (ssize_t)var->args().size())
                stack.emplace_back(var, idx + 1);
            var = new_var;
            idx = -1;
            return true;
        };
        while (iterate()) {
        }
    }

    m_varuse_dirty = false;
}

NACS_EXPORT() void Env::gc()
{
    compute_varuse();
    for (auto _var = m_vars; _var;) {
        auto var = _var;
        _var = _var->m_next;
        if (var->m_extern_ref > 0 || var->m_used)
            continue;
        var->release();
        m_varid_dirty = true;
    }
    // Remove unused variables should not affect varuse.
    assert(!m_varuse_dirty);
}

bool Env::optimize_local()
{
    assert(!m_varid_dirty);
    if (!m_vars)
        return false;
    bool changed = false;
    // 0: unhandled
    // 1: working on it
    // 2: done
    llvm::SmallVector<uint8_t, 64> var_states(m_vars->m_varid + 1, 0); // nvars
    llvm::SmallVector<std::pair<Var*,ssize_t>, 16> stack;
    ssize_t idx;
    Var *var;
    auto next = [&] {
        // Next argument
        if (++idx <= (ssize_t)var->args().size())
            return true;
        var_states[var->m_varid] = 2;
        // Done
        if (stack.empty())
            return false;
        std::tie(var, idx) = stack.pop_back_val();
        return true;
    };
    auto visit_field = [&] (Var *new_var) {
        auto &vs = var_states[new_var->m_varid];
        if (vs == 2)
            return next();
        if (!new_var->is_call()) {
            vs = 2;
            return next();
        }
        if (vs == 1)
            throw std::runtime_error("Dependency loop detected.");
        vs = 1;
        if (idx + 1 <= (ssize_t)var->args().size()) {
            stack.emplace_back(var, idx + 1);
        }
        else {
            var_states[var->m_varid] = 2;
        }
        var = new_var;
        idx = -2;
        return true;
    };
    auto iterate = [&] {
        if (idx == -2) {
            // Optimize callee
            assert(var->is_call());
            auto f = var->get_callee();
            if (!f.is_llvm)
                return visit_field(f.var);
            return next();
        }
        else if (idx == -1) {
            // Check callee
            changed |= var->inline_callee();
            return next();
        }
        else if (idx < (ssize_t)var->args().size()) {
            // Optimize arguments
            // This must be a call since there are at least 1 argument.
            assert(var->is_call());
            auto f = var->get_callee();
            if (f.is_llvm) {
                if (LLVM::Analysis::argument_unused(*f.llvm, idx)) {
                    // Unused, ignore.
                    changed |= !var->args()[idx].is_const();
                    var->set_arg(idx, Arg::create_const(false));
                    return next();
                }
            }
            else {
                if (f.var->argument_unused(idx)) {
                    // Unused, ignore.
                    changed |= !var->args()[idx].is_const();
                    var->set_arg(idx, Arg::create_const(false));
                    return next();
                }
            }
            auto &arg = var->args()[idx];
            if (!arg.is_var())
                return next();
            auto arg_var = arg.get_var();
            if (arg_var->is_const()) {
                var->set_arg(idx, Arg::create_const(arg_var->get_const()));
                changed = true;
                return next();
            }
            return visit_field(arg_var);
        }
        else {
            // Check arguments
            if (!var->is_call())
                return next();
            bool has_arg = false;
            // There could be constant variable argument again
            // since the argument var optimization might have optimized
            // the variable to a constant
            // Normalize that first.
            for (size_t i = 0; i < var->args().size(); i++) {
                auto &arg = var->args()[i];
                if (arg.is_arg()) {
                    has_arg = true;
                    continue;
                }
                if (!arg.is_var())
                    continue;
                auto arg_var = arg.get_var();
                if (arg_var->is_const()) {
                    var->set_arg(i, Arg::create_const(arg_var->get_const()));
                    changed = true;
                }
                else if (auto copy = arg_var->get_assigned_var()) {
                    var->set_arg(i, Arg::create_var(copy));
                    changed = true;
                }
            }
            // We don't want to surprise the user with mismatch argument number.
            // However, we allowe it in the special case where none of the arguments are used.
            if (!has_arg && var->nfreeargs() != 0) {
                var->m_n_freeargs = 0;
                changed = true;
            }
            // Now all of the non-const arguments are used by the function
            // and all of the constant arguments are stored as constant `Arg`
            // (instead of `Arg` of constant `Var`)
            // We can execute the function iff all the arguments are constant.
            if (var->optimize_call()) {
                changed = true;
                if (var->is_call()) {
                    // If `Var::optimize_call` returns `true`,
                    // we have changed the LLVM function
                    // and should check the arguments property again
                    // if it is still a call.
                    // See the `test_return_arg_indirect` test case.
                    idx = -1;
                    return next();
                }
            }
            return next();
        }
    };
    for (auto root: *this) {
        if (root->m_extern_ref <= 0)
            continue;
        if (!root->is_call()) {
            var_states[root->m_varid] = 2;
            continue;
        }
        auto &vs = var_states[root->m_varid];
        // We are just starting a cycle so there shouldn't be nothing in progress.
        assert(vs != 1);
        // Already processed
        if (vs == 2)
            continue;
        var = root;
        idx = -2;
        vs = 1;
        while (iterate()) {
        }
    }
    return changed;
}

NACS_EXPORT() void Env::optimize()
{
    if (!m_vars)
        return;

    // Initialize variable ID eagerly so the optimization functions don't need to check
    if (m_varid_dirty)
        compute_varid();
    optimize_local();
    gc();
    // TODO: global dependency optimization
    finalize_vars();
}

void Env::finalize_vars()
{
    // We don't use any LLVM global info for optimization so we only need to do the LLVM
    // global DCE once.
    for (auto &go: llvm_module()->global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    auto type_matches = [&] (Var *var) {
        auto f = var->get_callee();
        assert(f.is_llvm);
        auto args = var->args();
        uint32_t nargs = args.size();
        auto fty = f.llvm->getFunctionType();
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            if (arg.is_arg())
                continue;
            // Constants should have been optimized out.
            assert(arg.is_var());
            if (fty->getParamType(i) != arg.get_var()->llvm_type()) {
                return false;
            }
        }
        return true;
    };
    // Make sure the argument type matches the LLVM type
    for (auto var: *this) {
        if (!var->is_call())
            continue;
        auto f = var->get_callee();
        // There can be assignments
        if (!f.is_llvm) {
            assert(var->args().empty());
            continue;
        }
        else if (type_matches(var)) {
            f.llvm->setLinkage(llvm::GlobalValue::ExternalLinkage);
            f.llvm->setVisibility(llvm::GlobalValue::ProtectedVisibility);
            continue;
        }
        auto args = var->args();
        uint32_t nargs = args.size();
        llvm::SmallVector<llvm::Type*,8> arg_types(nargs);
        auto oldfty = f.llvm->getFunctionType();
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            if (arg.is_arg()) {
                arg_types[i] = oldfty->getParamType(i);
            }
            else {
                // Constants should have been optimized out.
                assert(arg.is_var());
                arg_types[i] = arg.get_var()->llvm_type();
            }
        }
        auto rt = oldfty->getReturnType();
        auto llvm_mod = llvm_module();
        auto &llvm_ctx = llvm_mod->getContext();
        auto fty = llvm::FunctionType::get(rt, arg_types, false);
        auto newf = llvm::Function::Create(fty, llvm::GlobalValue::ExternalLinkage,
                                           f.llvm->getName() + ".t", llvm_mod);
        newf->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        newf->setAttributes(llvm::AttributeList::get(llvm_ctx,
                                                     f.llvm->getAttributes().getFnAttributes(),
                                                     {}, {}));
        auto *b0 = llvm::BasicBlock::Create(llvm_ctx, "top", newf);
        llvm::IRBuilder<> builder(b0);
        llvm::SmallVector<llvm::Value*,8> call_args(nargs);
        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = args[i];
            llvm::Value *larg = newf->getArg(i);
            if (!arg.is_arg()) {
                // Constants should have been optimized out.
                assert(arg.is_var());
                if (oldfty->getParamType(i) != larg->getType()) {
                    if (arg.get_var()->type() == IR::Type::Bool ||
                        oldfty->getParamType(i) == m_cgctx->T_i8)
                        larg = LLVM::convert_scalar(builder, m_cgctx->T_bool, larg);
                    larg = LLVM::convert_scalar(builder, oldfty->getParamType(i), larg);
                }
            }
            call_args[i] = larg;
        }
        auto call = builder.CreateCall(f.llvm, call_args);
        builder.CreateRet(call);
        llvm::InlineFunctionInfo IFI;
#if LLVM_VERSION_MAJOR >= 11
        llvm::InlineFunction(*call, IFI);
#else
        llvm::InlineFunction(call, IFI);
#endif
        var->assign_call(newf, args, var->nfreeargs());
    }
    llvm::legacy::PassManager PM;
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.add(llvm::createGlobalDCEPass());
    // Shorten all global names since we don't care what they are
    // and this should slightly reduce the compiled binary size.
    PM.add(LLVM::createGlobalRenamePass());
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif

    PM.run(*llvm_module());
}

NACS_EXPORT() void Env::print(std::ostream &stm) const
{
    llvm::SmallVector<Var*, 32> vars(begin(), end());
    if (vars.empty()) {
        stm << "Variables: <empty>" << std::endl;
    }
    else {
        stm << "Variables: <" << vars.size() << ">" << std::endl;
    }
    for (auto it = vars.rbegin(), end = vars.rend(); it != end; ++it) {
        stm << "  ";
        (*it)->print(stm);
        stm << std::endl;
    }
    stm << std::endl;
    if (!llvm_module()) {
        stm << "LLVM: <null>" << std::endl;
    }
    else {
        stm << "LLVM:";
        string_ostream sstm;
        llvm::raw_os_ostream lstm(sstm);
        llvm_module()->print(lstm, nullptr);
        auto str = sstm.get_buf();
        auto p = &str[0];
        auto end = p + str.size();
        for (; p < end && *p == '\n'; p++) {
        }
        if (p < end) {
            stm << std::endl;
            while (p < end) {
                auto lend = (char*)memchr(p, '\n', end - p);
                stm.write("  ", 2);
                if (!lend) {
                    stm.write(p, end - p);
                    stm.put('\n');
                    break;
                }
                stm.write(p, lend + 1 - p);
                p = lend + 1;
            }
        }
        else {
            stm << " <empty>" << std::endl;
        }
    }
}

}
