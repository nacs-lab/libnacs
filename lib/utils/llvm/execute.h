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

#ifndef __NACS_UTILS_LLVM_EXECUTE_H__
#define __NACS_UTILS_LLVM_EXECUTE_H__

#include <llvm/ADT/SmallVector.h>
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

class AllocTracker {
public:
    AllocTracker(size_t blocksz)
        : m_blocksz(blocksz)
    {
    }
    template<typename BlockAlloc>
    void *alloc(size_t sz, BlockAlloc &&block_alloc);
    template<typename BlockFree>
    void free(void *ptr, size_t sz, BlockFree &&block_free);

private:
    template<typename BlockFree>
    void free_real(void *ptr, size_t sz, BlockFree &&block_free,
                   SmallVector<std::map<void*,size_t>::iterator,2> replaces={});

    const size_t m_blocksz;
    std::map<void*,size_t> m_freelist;
    std::map<void*,size_t> m_blocks;
};

class MemMgr : public RuntimeDyld::MemoryManager {
public:
    MemMgr();
    ~MemMgr();
    uint64_t new_group();
    void free_group(uint64_t id);

private:
    uint8_t *allocateCodeSection(uintptr_t sz, unsigned align,
                                 unsigned sid, StringRef sname) override;
    uint8_t *allocateDataSection(uintptr_t sz, unsigned align,
                                 unsigned sid, StringRef sname, bool ro) override;
    void registerEHFrames(uint8_t *addr, uint64_t load_addr, size_t sz) override;
    void deregisterEHFrames() override;
    bool finalizeMemory(std::string *ErrMsg) override;
    SectionMemoryManager tmp;

    uint64_t m_grp_cnt = 0;
};

class Engine {
public:
    Engine();
    ~Engine();
    uint64_t load(const object::ObjectFile &O);
    uint64_t load(const char *p, size_t len);
    void *get_symbol(StringRef name);
    void free(uint64_t id);
private:
    void reset_dyld();
    MemMgr m_memmgr;
    Resolver m_resolver;
    std::unique_ptr<RuntimeDyld> m_dyld;
};

}
}
}

#endif
