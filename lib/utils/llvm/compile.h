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
#include "../ir.h"

#include <llvm/IR/Module.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Target/TargetMachine.h>

namespace NaCs {
namespace LLVM {
namespace Compile {

using namespace llvm;

Module *optimize(Module *mod);
bool emit_objfile(raw_pwrite_stream &stm, TargetMachine *tgt, Module *M, bool opt=true);
TargetMachine *get_native_target();

}
}
}
