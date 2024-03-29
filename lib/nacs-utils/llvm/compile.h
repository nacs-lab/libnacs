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

#ifndef __NACS_UTILS_LLVM_COMPILE_H__
#define __NACS_UTILS_LLVM_COMPILE_H__

#include "../utils.h"
#include "../ir.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/ADT/SmallVector.h>
#include <llvm/IR/Module.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Target/TargetMachine.h>

namespace NaCs::LLVM::Compile {

using namespace llvm;

NACS_EXPORT(utils) Module *optimize(Module *mod);
NACS_EXPORT(utils) bool emit_objfile(raw_pwrite_stream &stm, TargetMachine *tgt, Module *M, bool opt=true);
NACS_EXPORT(utils) bool emit_objfile(SmallVectorImpl<char> &vec, TargetMachine *tgt, Module *M, bool opt=true);
// For testing only (to avoid directly linking to libLLVM in test).
NACS_EXPORT(utils) bool emit_objfile(std::vector<char> &vec, TargetMachine *tgt, Module *M, bool opt=true);
NACS_EXPORT(utils) TargetMachine *get_native_target();
NACS_EXPORT(utils) std::unique_ptr<TargetMachine> create_target(StringRef triple, StringRef cpu,
                                                                StringRef features);

}

#endif
