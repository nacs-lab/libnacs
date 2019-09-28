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

#include "../utils/llvm/codegen.h"
#include "../utils/llvm/utils.h"
#include "../utils/streams.h"

#include <llvm/ADT/StringMap.h>
#include <llvm/Support/raw_os_ostream.h>

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
