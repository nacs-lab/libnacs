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

#include "utils.h"
#include "cpu_p.h"

#include <llvm/ADT/Triple.h>
#include <llvm/Support/Host.h>
#include <llvm/Transforms/Utils/Cloning.h>

namespace NaCs {
namespace LLVM {

using namespace llvm;

NACS_EXPORT() void dump(const Value *v)
{
    v->print(dbgs(), true);
    dbgs() << "\n";
}

NACS_EXPORT() void dump(const Type *v)
{
    v->print(dbgs(), true);
    dbgs() << "\n";
}

NACS_EXPORT() void dump(const Function *f)
{
    f->print(dbgs(), nullptr, false, true);
}

NACS_EXPORT() void dump(const Module *m)
{
    m->print(dbgs(), nullptr);
}

NACS_EXPORT() void dump(const Metadata *m)
{
    m->print(dbgs());
    dbgs() << "\n";
}

NACS_EXPORT() void dump(const DebugLoc *dbg)
{
    dbg->print(dbgs());
    dbgs() << "\n";
}

NACS_EXPORT() module_ref new_module(StringRef name, LLVMContext &ctx)
{
    return module_ref(new Module(name, ctx));
}

NACS_EXPORT() void delete_module(Module *mod)
{
    delete mod;
}

NACS_EXPORT() context_ref new_context()
{
    return context_ref(new LLVMContext);
}

NACS_EXPORT() void delete_context(LLVMContext *ctx)
{
    delete ctx;
}

NACS_EXPORT() IR::Type get_ir_type(llvm::Type *typ, bool apitype)
{
    if (typ->isDoubleTy())
        return IR::Type::Float64;
    if (typ->isIntegerTy()) {
        auto bits = typ->getPrimitiveSizeInBits();
        if (bits == 32) {
            return IR::Type::Int32;
        }
        else if ((apitype && bits == 8) || (!apitype && bits == 1)) {
            return IR::Type::Bool;
        }
    }
    return IR::Type::_Min;
}

NACS_EXPORT_ FunctionMover::FunctionMover(Module *dest)
    : ValueMaterializer(), m_vmap(), m_dest(dest), m_lazy_funcs(0)
{
}

NACS_EXPORT() FunctionMover::~FunctionMover()
{
}

NACS_INTERNAL Function *FunctionMover::queue_proto(Function *F)
{
    assert(!F->isDeclaration());
    auto newf = Function::Create(F->getFunctionType(), Function::ExternalLinkage,
                                 F->getName(), m_dest);
    m_lazy_funcs.push_back(F);
    m_vmap[F] = newf;
    return newf;
}

NACS_INTERNAL void FunctionMover::clone_body(Function *F)
{
    auto newf = (Function*)(Value*)m_vmap[F];
    assert(newf);

    auto destit = newf->arg_begin();
    for (auto &arg: F->args()) {
        destit->setName(arg.getName()); // Copy the name over...
        m_vmap[&arg] = &*(destit++); // Add mapping to m_vmap
    }

    SmallVector<ReturnInst*, 8> Returns;
    CloneFunctionInto(newf, F, m_vmap, true, Returns, "", nullptr, nullptr, this);
}

NACS_EXPORT() Function *FunctionMover::clone_function(Function *F)
{
    auto newf = (Function*)MapValue(F, m_vmap, RF_None, nullptr, this);
    resolve_lazy();
    return newf;
}

NACS_INTERNAL void FunctionMover::resolve_lazy()
{
    while (!m_lazy_funcs.empty()) {
        clone_body(m_lazy_funcs.pop_back_val());
    }
}

NACS_INTERNAL Value *FunctionMover::clone_proto(Function *F)
{
    auto newf = m_dest->getFunction(F->getName());
    if (!newf) {
        newf = Function::Create(F->getFunctionType(), Function::ExternalLinkage,
                                F->getName(), m_dest);
        // FunctionType does not include any attributes. Copy them over manually
        // as codegen may make decisions based on the presence of certain attributes
        newf->copyAttributesFrom(F);
    }
    return newf;
}

NACS_INTERNAL Value *FunctionMover::materialize(Value *V)
{
    if (auto F = dyn_cast<Function>(V)) {
        if (F->getParent() == m_dest)
            return F;
        if (F->isIntrinsic() || F->isDeclaration())
            return clone_proto(F);
        return queue_proto(F);
    }
    else if (auto GV = dyn_cast<GlobalVariable>(V)) {
        if (auto oldGV = m_dest->getGlobalVariable(GV->getName()))
            return oldGV;
        auto newGV = new GlobalVariable(*m_dest,
                                        GV->getType()->getElementType(),
                                        GV->isConstant(),
                                        GV->getLinkage(),
                                        NULL,
                                        GV->getName(),
                                        NULL,
                                        GV->getThreadLocalMode(),
                                        GV->getType()->getPointerAddressSpace());
        newGV->copyAttributesFrom(GV);
        if (GV->isDeclaration())
            return newGV;
        if (GV->hasInitializer()) {
            Value *C = MapValue(GV->getInitializer(), m_vmap, RF_None, NULL, this);
            newGV->setInitializer(cast<Constant>(C));
        }
        return newGV;
    }
    return nullptr;
};

const std::string &get_cpu_arch()
{
    static const std::string arch = Triple(sys::getProcessTriple()).getArchName();
    return arch;
}

const std::string &get_cpu_name()
{
    static const std::string name = sys::getHostCPUName().str();
    return name;
}

const std::string &get_cpu_features()
{
    static const std::string features =
        [] {
            StringMap<bool> HostFeatures;
            sys::getHostCPUFeatures(HostFeatures);
            std::string attr;
            for (auto &ele: HostFeatures) {
                if (ele.getValue()) {
                    if (!attr.empty()) {
                        attr.append(",+");
                    }
                    else {
                        attr.append("+");
                    }
                    attr.append(ele.getKey().str());
                }
            }
            // Explicitly disabled features need to be added at the end so that
            // they are not reenabled by other features that implies them by default.
            for (auto &ele: HostFeatures) {
                if (!ele.getValue()) {
                    if (!attr.empty()) {
                        attr.append(",-");
                    }
                    else {
                        attr.append("-");
                    }
                    attr.append(ele.getKey().str());
                }
            }
            return attr;
        }();
    return features;
}

}
}
