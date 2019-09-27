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

#ifndef __NACS_UTILS_LLVM_INST_SIMPLIFY_H__
#define __NACS_UTILS_LLVM_INST_SIMPLIFY_H__

#include <llvm/IR/InstrTypes.h> // Required to workaround bug in `InstructionSimplify.h`
#include <llvm/Analysis/InstructionSimplify.h>
#include <llvm/IR/Function.h>
#include <llvm/Pass.h>

namespace NaCs {
namespace LLVM {

using namespace llvm;

using resolver_cb_t = std::function<uintptr_t(const std::string&)>;

bool instSimplify(Function &F, const SimplifyQuery &SQ,
                  const resolver_cb_t &cb=resolver_cb_t(),
                  OptimizationRemarkEmitter *ORE=nullptr);
bool instSimplify(Function &F, const resolver_cb_t &cb=resolver_cb_t(),
                  OptimizationRemarkEmitter *ORE=nullptr);
FunctionPass *createNaCsInstSimplifyPass(const resolver_cb_t &cb=resolver_cb_t());

}
}

#endif
