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

#include <llvm/ExecutionEngine/RuntimeDyld.h>
#include <llvm/ExecutionEngine/SectionMemoryManager.h>

#include <memory>

namespace NaCs {
namespace LLVM {
namespace Exe {

using namespace llvm;

class Resolver : public JITSymbolResolver {
public:
    Resolver() = default;
private:
    JITSymbol findSymbolInLogicalDylib(const std::string&) override;
    JITSymbol findSymbol(const std::string &name) override;
    uintptr_t find_extern(const std::string &name);
};

class MemMgr : public SectionMemoryManager {
};

class Engine {
public:
    Engine();
    ~Engine();
    bool load(const object::ObjectFile &O);
    bool load(const char *p, size_t len);
    void *get_symbol(StringRef name);
private:
    void reset_dyld();
    MemMgr m_memmgr;
    Resolver m_resolver;
    std::unique_ptr<RuntimeDyld> m_dyld;
};

}
}
}
