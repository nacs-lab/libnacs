/*************************************************************************
 *   Copyright (c) 2020 - 2024 Yichao Yu <yyc1992@gmail.com>             *
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

#include "passes.h"
#include "utils.h"

#include <llvm/IR/Verifier.h>

#if NACS_ENABLE_NEW_PASS
#  include <llvm/Transforms/IPO/GlobalDCE.h>
#else
#  include <llvm/IR/LegacyPassManager.h>
#  include <llvm/Transforms/IPO.h>
#endif

namespace NaCs::LLVM {

using namespace llvm;

// This pass renames all the exported globals assuming that the exact name of the names
// to be unimportant. User is expected to use reference to the `GlobalValue`s to access
// the globals later (including obtaining the new name)

namespace {

class NameCounter {
public:
    const char *next_name()
    {
        // Do not include `.` in the list below since that can poentially cause names
        // starting with `llvm.` to be generated.
        // We may also want to avoid lib functions in the future.
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

bool renameGlobal(Module &M)
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
    retry:
        auto new_name = counter.next_name();
        if (g.getName() == new_name)
            continue;
        if (auto old_g = M.getNamedValue(new_name)) {
            // Don't rename declaration
            if (old_g->isDeclaration())
                goto retry;
            // Either a local function or a global that we want to rename
            // Just set the name to something random to avoid conflict.
            old_g->setName(".");
        }
        changed = true;
        g.setName(new_name);
    }
    return changed;
}

#if NACS_ENABLE_LEGACY_PASS
struct LegacyGlobalRenamePass : public ModulePass {
    static char ID;
    LegacyGlobalRenamePass()
        : ModulePass(ID)
    {}

private:
    bool runOnModule(Module &M) override
    {
        return renameGlobal(M);
    }
};

char LegacyGlobalRenamePass::ID = 0;
static RegisterPass<LegacyGlobalRenamePass> X("GlobalRename", "Rename all exported globals",
                                        false /* Only looks at CFG */,
                                        false /* Analysis Pass */);
#endif

}

#if NACS_ENABLE_LEGACY_PASS
NACS_EXPORT() Pass *createGlobalRenamePass()
{
    return new LegacyGlobalRenamePass();
}
#endif

#if NACS_ENABLE_NEW_PASS
NACS_EXPORT() PreservedAnalyses
GlobalRenamePass::run(Module &M, ModuleAnalysisManager &AM)
{
    if (renameGlobal(M))
        return PreservedAnalyses::allInSet<CFGAnalyses>();
    return PreservedAnalyses::all();
}

NACS_EXPORT() void runGlobalRenamePasses(Module &M)
{
    llvm::PassBuilder PB;
    LLVM::AnalysisManagers AM(PB);

    llvm::ModulePassManager MPM;
#  ifndef NDEBUG
    MPM.addPass(llvm::VerifierPass());
#  endif
    MPM.addPass(llvm::GlobalDCEPass());
    // Shorten all global names since we don't care what they are
    // and this should slightly reduce the compiled binary size.
    MPM.addPass(LLVM::GlobalRenamePass());
#  ifndef NDEBUG
    MPM.addPass(llvm::VerifierPass());
#  endif
    MPM.run(M, AM.MAM);
}
#else
NACS_EXPORT() void runGlobalRenamePasses(Module &M)
{
    llvm::legacy::PassManager PM;
#  ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#  endif
    PM.add(llvm::createGlobalDCEPass());
    // Shorten all global names since we don't care what they are
    // and this should slightly reduce the compiled binary size.
    PM.add(LLVM::createGlobalRenamePass());
#  ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#  endif
    PM.run(M);
}
#endif

}
