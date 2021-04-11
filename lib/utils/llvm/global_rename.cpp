/*************************************************************************
 *   Copyright (c) 2020 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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

#define DEBUG_TYPE "global_rename"

#include "global_rename.h"
#include "utils.h"

namespace NaCs::LLVM {

using namespace llvm;

// This pass renames all the exported globals assuming that the exact name of the names
// to be unimportant. User is expected to use reference to the `GlobalValue`s to access
// the globals later (including obtaining the new name)

namespace {

struct GlobalRenamePass : public ModulePass {
    static char ID;
    GlobalRenamePass()
        : ModulePass(ID)
    {}

private:
    bool runOnModule(Module &M) override;
};

class NameCounter {
public:
    const char *next_name()
    {
        static constexpr char digits[] =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        static_assert(sizeof(digits) == 65, "");
        for (unsigned i = 0; i < sizeof(counters); i++) {
            if (counters[i] < 64) {
                counters[i] += 1;
                name[i] = digits[counters[i] - 1];
                break;
            }
            counters[i] = 1;
            name[i] = digits[0];
        }
        return name;
    }
private:
    uint8_t counters[8] = {};
    char name[8] = {};
};

bool GlobalRenamePass::runOnModule(Module &M)
{
    bool changed = false;
    NameCounter counter;
    for (auto &g: M.global_values()) {
        // We can't rename declarations since they are referring to external
        // objects. We don't need to rename local objects since
        // their name aren't significan't and won't appear in the
        // final binary anyway.
        if (g.isDeclaration() || g.hasLocalLinkage())
            continue;
        auto new_name = counter.next_name();
        if (g.getName() == new_name)
            continue;
        changed = true;
        if (auto old_g = M.getNamedValue(new_name))
            old_g->setName(".");
        g.setName(new_name);
    }
    return changed;
}

char GlobalRenamePass::ID = 0;
static RegisterPass<GlobalRenamePass> X("GlobalRename", "Rename all exported globals",
                                        false /* Only looks at CFG */,
                                        false /* Analysis Pass */);

}

NACS_EXPORT() Pass *createGlobalRenamePass()
{
    return new GlobalRenamePass();
}

}
