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

#include "analysis.h"
#include "utils.h"

#include <llvm/IR/Instructions.h>

namespace NaCs::LLVM::Analysis {

// We could do something fancier like detecting no-op uses of the argument
// However, that's a better job for DCE and our input shouldn't have a lot of (any?)
// dead code anyway.
NACS_EXPORT(utils) bool argument_unused(const Function &f, unsigned argno)
{
    auto args = f.arg_begin();
    return (args + argno)->use_empty();
}

NACS_EXPORT(utils) Constant *returns_const(const Function &f)
{
    if (!f.isSpeculatable() && !(f.doesNotThrow() && f.onlyReadsMemory()))
        return nullptr;
    Constant *ret_val = nullptr;
    for (const BasicBlock &bb: f) {
        auto ret_inst = dyn_cast<ReturnInst>(bb.getTerminator());
        if (!ret_inst)
            continue;
        auto new_ret_val = dyn_cast<Constant>(ret_inst->getReturnValue());
        // Not returning a Constant.
        if (!new_ret_val)
            return nullptr;
        // Returning a different Constant.
        if (ret_val && ret_val != new_ret_val)
            return nullptr;
        ret_val = new_ret_val;
    }
    return ret_val;
}

NACS_EXPORT(utils) Argument *returns_argument(const Function &f)
{
    if (!f.isSpeculatable() && !(f.doesNotThrow() && f.onlyReadsMemory()))
        return nullptr;
    Argument *ret_val = nullptr;
    for (const BasicBlock &bb: f) {
        auto ret_inst = dyn_cast<ReturnInst>(bb.getTerminator());
        if (!ret_inst)
            continue;
        auto new_ret_val = dyn_cast<Argument>(ret_inst->getReturnValue());
        // Not returning an argument.
        if (!new_ret_val)
            return nullptr;
        // Returning a different argument
        if (ret_val && ret_val != new_ret_val)
            return nullptr;
        ret_val = new_ret_val;
    }
    return ret_val;
}

}
