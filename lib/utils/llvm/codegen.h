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

#include "../utils.h"

#include <llvm/IR/Module.h>
#include <llvm/IR/Type.h>

namespace NaCs {
namespace LLVM {
namespace Codegen {

using namespace llvm;

class NACS_EXPORT(utils) Context {
public:
    Context(Module *mod);
private:
    Module *m_mod;
    LLVMContext &m_ctx;
    IntegerType *T_bool;
    IntegerType *T_i32;
    IntegerType *T_isz;
    Type *T_f64;
    FunctionType *F_f64_f64;
    FunctionType *F_f64_f64f64;
    FunctionType *F_f64_f64f64f64;
    FunctionType *F_f64_f64i32;
    FunctionType *F_f64_i32f64;
};

}
}
}
