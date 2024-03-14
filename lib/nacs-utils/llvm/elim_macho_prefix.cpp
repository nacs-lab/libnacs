/*************************************************************************
 *   Copyright (c) 2021 - 2024 Yichao Yu <yyc1992@gmail.com>             *
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

#define DEBUG_TYPE "elim_macho_prefix"

#include "passes.h"
#include "utils.h"

namespace NaCs::LLVM {

using namespace llvm;

// This pass adds a `\01` prefix in front of all the exported globals.
// This prevents the backends from adding a `_` prefix to the symbol names
// and makes sure the symbol name in the LLVM module is the same name one can use
// in get_symbol from RTLD.

namespace {

static bool elimMachOPrefix(Module &M)
{
    bool changed = false;
    for (auto &g: M.global_values()) {
        // We can't rename declarations since they are referring to external
        // objects. We don't need to rename local objects since
        // their name aren't significan't and won't appear in the
        // final binary anyway.
        if (g.isDeclaration() || g.hasLocalLinkage())
            continue;
        auto old_name = g.getName();
        // Try not to apply the fix twice.
        // This allow the module to be compiled again
        // (even though the function names has changed.)
        // Real users shouldn't really need this but this is used in the test.
        if (!old_name.empty() && old_name[0] == '\x01')
            continue;
        std::string new_name = "\x01" + old_name.str();
        assert(!M.getNamedValue(new_name));
        changed = true;
        g.setName(new_name);
    }
    return changed;
}

#if NACS_ENABLE_LEGACY_PASS
struct LegacyElimMachOPrefixPass : public ModulePass {
    static char ID;
    LegacyElimMachOPrefixPass()
        : ModulePass(ID)
    {}

private:
    bool runOnModule(Module &M) override
    {
        return elimMachOPrefix(M);
    }
};

char LegacyElimMachOPrefixPass::ID = 0;
static RegisterPass<LegacyElimMachOPrefixPass> X("ElimMachOPrefix",
                                           "Prevent underscore prefix",
                                           false /* Only looks at CFG */,
                                           false /* Analysis Pass */);
#endif

}

#if NACS_ENABLE_LEGACY_PASS
NACS_EXPORT() Pass *createElimMachOPrefixPass()
{
    return new LegacyElimMachOPrefixPass();
}
#endif

#if NACS_ENABLE_NEW_PASS
NACS_EXPORT() PreservedAnalyses
ElimMachOPrefixPass::run(Module &M, ModuleAnalysisManager &AM)
{
    if (elimMachOPrefix(M))
        return PreservedAnalyses::allInSet<CFGAnalyses>();
    return PreservedAnalyses::all();
}
#endif

}
