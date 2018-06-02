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

#include "codegen.h"
#include "utils.h"

namespace NaCs {
namespace LLVM {
namespace Codegen {

Context::Context(Module *mod)
    : m_mod(mod),
      m_ctx(mod->getContext()),
      T_bool(Type::getInt1Ty(m_ctx)),
      T_i32(Type::getInt32Ty(m_ctx)),
      T_isz(sizeof(void*) == 8 ? Type::getInt64Ty(m_ctx) : T_i32),
      T_f64(Type::getDoubleTy(m_ctx)),
      F_f64_f64(FunctionType::get(T_f64, {T_f64}, false)),
      F_f64_f64f64(FunctionType::get(T_f64, {T_f64, T_f64}, false)),
      F_f64_f64f64f64(FunctionType::get(T_f64, {T_f64, T_f64, T_f64}, false)),
      F_f64_f64i32(FunctionType::get(T_f64, {T_f64, T_i32}, false)),
      F_f64_i32f64(FunctionType::get(T_f64, {T_i32, T_f64}, false))
{
}

}
}
}
