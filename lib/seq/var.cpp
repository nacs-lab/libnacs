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

#include "../utils/llvm/utils.h"

#include <stdexcept>

namespace NaCs::Seq {

NACS_EXPORT() IR::Type Var::type() const
{
    if (is_const())
        return m_value.constant.typ;
    if (is_extern())
        return m_value.ext.first;
    assert(is_call());
    auto f = get_callee();
    if (!f.is_llvm)
        return f.var->type();
    return LLVM::get_ir_type(f.llvm->getReturnType());
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

}
