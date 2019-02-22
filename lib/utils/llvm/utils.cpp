/*************************************************************************
 *   Copyright (c) 2018 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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
#include "utils.h"

namespace NaCs {
namespace LLVM {

using namespace llvm;

NACS_EXPORT() void dump(Value *v)
{
    v->print(dbgs(), true);
    dbgs() << "\n";
}

NACS_EXPORT() void dump(Type *v)
{
    v->print(dbgs(), true);
    dbgs() << "\n";
}

NACS_EXPORT() void dump(Function *f)
{
    f->print(dbgs(), nullptr, false, true);
}

NACS_EXPORT() void dump(Module *m)
{
    m->print(dbgs(), nullptr);
}

NACS_EXPORT() void dump(Metadata *m)
{
    m->print(dbgs());
    dbgs() << "\n";
}

NACS_EXPORT() void dump(DebugLoc *dbg)
{
    dbg->print(dbgs());
    dbgs() << "\n";
}

NACS_EXPORT() Module *new_module(StringRef name, LLVMContext &ctx)
{
    return new Module(name, ctx);
}

NACS_EXPORT() void delete_module(Module *mod)
{
    delete mod;
}

NACS_EXPORT() LLVMContext *new_context()
{
    return new LLVMContext;
}

NACS_EXPORT() void delete_context(LLVMContext *ctx)
{
    delete ctx;
}

}
}
