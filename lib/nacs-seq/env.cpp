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

#include "../nacs-utils/llvm/utils.h"

namespace NaCs::Seq {

NACS_EXPORT_ Env::Env(llvm::LLVMContext &llvm_ctx)
    : m_llvm_mod(new llvm::Module("", llvm_ctx))
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
    if (func->arg_size() != args.size())
        throw std::invalid_argument("Argument number mismatch");
    check_args(args, nfreeargs);
    auto var = new_var();
    var->assign_call(func, args, nfreeargs);
    return var;
}

NACS_EXPORT() Var *Env::new_extern(std::pair<IR::Type,uint64_t> ext)
{
    auto var = new_var();
    var->assign_extern(ext);
    return var;
}

}
