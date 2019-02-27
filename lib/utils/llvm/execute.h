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

#include "../mem.h"

#include <llvm/ADT/SmallVector.h>
#include <llvm/ExecutionEngine/RuntimeDyld.h>
#include <llvm/ExecutionEngine/SectionMemoryManager.h>

#include <memory>
#include <utility>

namespace NaCs {
namespace LLVM {
namespace Exe {

using namespace llvm;

#if LLVM_VERSION_MAJOR >= 7
using ResolverBase = LegacyJITSymbolResolver;
#else
using ResolverBase = JITSymbolResolver;
#endif

class Resolver : public ResolverBase {
public:
    Resolver() = default;
private:
    JITSymbol findSymbolInLogicalDylib(const std::string&) override;
    JITSymbol findSymbol(const std::string &name) override;
    uintptr_t find_extern(const std::string &name);
};

/**
 * Track small allocations out of bigger blocks.
 * The minimum block size is configurable at construction time and the size of
 * each block will be a integer multiple of the min size.
 */
class AllocTracker {
public:
    AllocTracker(size_t blocksz)
        : m_blocksz(blocksz)
    {
    }

    // Allocate a piece of memory of size `sz` and alignment `align`.
    // If one cannot be found in the freelist `block_alloc` will be called
    // with the desired block size to allocate a block of memory.
    // It is guaranteed that an allocation won't cross multiple blocks even if `block_alloc`
    // returns two blocks that are adjacent in memory.
    // The return value of `block_alloc` is assumed to be properly aligned.
    template<typename BlockAlloc>
    void *alloc(size_t sz, size_t align, BlockAlloc &&block_alloc);
    template<typename BlockAlloc>
    void *alloc(size_t sz, BlockAlloc &&block_alloc)
    {
        return alloc(sz, 1, std::forward<BlockAlloc>(block_alloc));
    }
    // Freeing a piece of memory at `ptr` and size `sz`.
    // If this frees up a whole block, `block_free` will be called to free the empty block.
    template<typename BlockFree>
    void free(void *ptr, size_t sz, BlockFree &&block_free);
    size_t block_size() const
    {
        return m_blocksz;
    }
    template<typename BlockFree>
    void free_all(BlockFree &&block_free, bool free_cache=true);

private:
    // Similar to the public API `free` but `ptr` and `sz` are the address and the size
    // of the now free'd memory after merging with other free pieces **in the same block**
    // `replaces` may contain `m_freelist` iterators that should be replaced/removed
    // since they are merged into a single piece of memory.
    template<typename BlockFree>
    void free_real(void *ptr, size_t sz, BlockFree &&block_free,
                   SmallVector<std::map<void*,size_t>::iterator,2> replaces={});

    // Minimum block size
    const size_t m_blocksz;
    // Map of free memory pieces.
    // The address (key of the map) is the **end** address of the free space.
    // This way we can allocate from the start of the range without changing the key.
    std::map<void*,size_t> m_freelist;
    // Allocated blocks. The address in this one is the start address.
    std::map<void*,size_t> m_blocks;
    // Possibly cache one block of the minimum size.
    void *m_lastblock = nullptr;
};

class RWAllocator {
public:
    RWAllocator(size_t block_size);
    void *alloc(size_t size, size_t align);
    void free(void *ptr, size_t size);
    ~RWAllocator();

private:
    AllocTracker m_tracker;
};

struct ROAlloc {
    void *wraddr;
    void *rtaddr;
};
template<bool exec>
class ROAllocator {
public:
    virtual ROAlloc alloc(size_t size, size_t align) = 0;
    virtual void finalize() = 0;
    virtual void free(ROAlloc alloc, size_t size) = 0;
    virtual ~ROAllocator()
    {}
};

class MemMgr : public RuntimeDyld::MemoryManager {
public:
    MemMgr();
    ~MemMgr();
    uint64_t new_group();
    void free_group(uint64_t id);

private:
    struct Alloc {
        enum Type : uint8_t {
            RW,
            RO,
            RX
        };
        void *rtaddr;
        size_t size;
        Type type;
        void *wraddr = nullptr;
    };
    uint8_t *allocateCodeSection(uintptr_t sz, unsigned align,
                                 unsigned sid, StringRef sname) override;
    uint8_t *allocateDataSection(uintptr_t sz, unsigned align,
                                 unsigned sid, StringRef sname, bool ro) override;
    void notifyObjectLoaded(RuntimeDyld &dyld, const object::ObjectFile &obj) override;
    void registerEHFrames(uint8_t *addr, uint64_t load_addr, size_t sz) override;
    void deregisterEHFrames() override;
    bool finalizeMemory(std::string *ErrMsg) override;

    size_t m_blocksz;
    DualMap m_dualmap{false};
    MemWriter m_memwriter{false};
    RWAllocator m_rwalloc;
    std::unique_ptr<ROAllocator<false>> m_roalloc;
    std::unique_ptr<ROAllocator<true>> m_rxalloc;

    SmallVector<Alloc,4> *m_cur_allocs = nullptr;
    std::map<uint64_t,SmallVector<Alloc,4>> m_allocs;

    uint64_t m_grp_cnt = 0;
};

class Engine {
public:
    NACS_EXPORT(utils) Engine();
    NACS_EXPORT(utils) ~Engine();
    NACS_EXPORT(utils) uint64_t load(const object::ObjectFile &O);
    NACS_EXPORT(utils) uint64_t load(const char *p, size_t len);
    // Must be called before `free`ing any other files after the `load`.
    NACS_EXPORT(utils) void *get_symbol(StringRef name);
    // This does not invalidates pointers from any other files
    // but it does clears all symbols.
    NACS_EXPORT(utils) void free(uint64_t id);
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
