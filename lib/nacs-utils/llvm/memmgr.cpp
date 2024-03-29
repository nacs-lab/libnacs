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

#include "../utils.h"
#include "../number.h"
#include "../fd_utils.h"

#include "execute.h"

namespace NaCs::LLVM::Exe {

template<typename BlockAlloc>
void *AllocTracker::alloc(size_t sz, size_t align, BlockAlloc &&block_alloc)
{
    // Simple first fit allocator
    auto fallback = m_freelist.end();
    for (auto I = m_freelist.begin(), E = m_freelist.end(); I != E; I++) {
        auto freesz = I->second;
        if (freesz < sz)
            continue;
        auto start = (void*)((char*)I->first - freesz);
        auto padding = getPadding((uintptr_t)start, align);
        if (!padding) {
            if (freesz > sz) {
                I->second = freesz - sz;
                return start;
            }
            m_freelist.erase(I);
            return start;
        }
        if (padding + sz < freesz)
            continue;
        fallback = I;
    }
    if (fallback != m_freelist.end()) {
        auto freesz = fallback->second;
        auto start = (void*)((char*)fallback->first - freesz);
        auto padding = getPadding((uintptr_t)start, align);
        auto end_pad = freesz - padding - sz;
        if (end_pad != 0) {
            fallback->second = end_pad;
        }
        else {
            m_freelist.erase(fallback);
        }
        auto res = (void*)((char*)start + padding);
        m_freelist[res] = padding;
        return res;
    }
    auto blksz = alignTo(sz, m_blocksz);
    void *ptr;
    if (likely(blksz == m_blocksz) && m_lastblock) {
        ptr = m_lastblock;
        m_lastblock = nullptr;
    }
    else {
        ptr = block_alloc(blksz);
    }
    m_blocks[ptr] = blksz;
    auto szleft = blksz - sz;
    if (likely(szleft > 0)) {
        auto endptr = (void*)((char*)ptr + blksz);
        m_freelist[endptr] = szleft;
    }
    return ptr;
}

template<typename BlockFree>
void AllocTracker::free(void *ptr, size_t sz, BlockFree &&block_free)
{
    auto I = m_freelist.lower_bound(ptr);
    if (I == m_freelist.end()) {
        free_real(ptr, sz, std::forward<BlockFree>(block_free));
        return;
    }
    SmallVector<std::map<void*,size_t>::iterator,2> replaces;
    if (I->first == ptr) {
        if (m_blocks.find(ptr) == m_blocks.end()) {
            // Do not merge allocations from different blocks.
            replaces.push_back(I);
            ptr = (void*)((char*)I->first - I->second);
            sz = I->second + sz;
        }
        ++I;
        if (I == m_freelist.end()) {
            free_real(ptr, sz, std::forward<BlockFree>(block_free),
                      std::move(replaces));
            return;
        }
    }
    if ((char*)I->first - I->second == (char*)ptr + sz) {
        if (m_blocks.find((void*)((char*)I->first - I->second)) == m_blocks.end()) {
            // Do not merge allocations from different blocks.
            replaces.push_back(I);
            sz = I->second + sz;
        }
    }
    free_real(ptr, sz, std::forward<BlockFree>(block_free),
              std::move(replaces));
}

template<typename BlockFree>
void AllocTracker::free_real(void *ptr, size_t sz, BlockFree &&block_free,
                             SmallVector<std::map<void*,size_t>::iterator,2> replaces)
{
    // First check if we can free a block
    if (sz >= m_blocksz) {
        // This can only happen if the memory is large enough.
        auto I = m_blocks.find(ptr);
        if (I != m_blocks.end() && I->second == sz) {
            if (likely(sz == m_blocksz) && !m_lastblock) {
                m_lastblock = ptr;
            }
            else {
                block_free(ptr, sz);
            }
            m_blocks.erase(I);
            for (auto fi: replaces)
                m_freelist.erase(fi);
            return;
        }
    }
    auto endptr = (void*)((char*)ptr + sz);
    bool found = false;
    for (auto fi: replaces) {
        if (!found && fi->first == endptr) {
            // Reuse if possible
            fi->second = sz;
            continue;
        }
        m_freelist.erase(fi);
    }
    if (!found) {
        m_freelist[endptr] = sz;
    }
}

template<typename BlockFree>
void AllocTracker::free_all(BlockFree &&block_free, bool free_cache)
{
    if (free_cache && m_lastblock) {
        block_free(m_lastblock, m_blocksz);
        m_lastblock = nullptr;
    }
    for (auto block: m_blocks)
        block_free(block.first, block.second);
    m_blocks.clear();
    m_freelist.clear();
}

inline RWAllocator::RWAllocator(size_t block_size)
    : m_tracker(block_size)
{
}

void *RWAllocator::alloc(size_t size, size_t align)
{
    return m_tracker.alloc(size, align, [&] (size_t bsize) {
            return mapAnonPage(bsize, Prot::RW);
        });
}

void RWAllocator::free(void *ptr, size_t size)
{
    m_tracker.free(ptr, size, [&] (void *ptr, size_t bsize) {
            return unmapPage(ptr, bsize);
        });
}

RWAllocator::~RWAllocator()
{
    m_tracker.free_all(unmapPage);
}

template<bool exec>
class DualMapAllocator : public ROAllocator<exec> {
public:
    DualMapAllocator(DualMap &dual_map, size_t block_size)
        : m_dualmap(dual_map),
          m_tracker(block_size)
    {
        if (!m_dualmap.init()) {
            throw std::runtime_error("Failed to initialize DualMap.");
        }
    }
    ROAlloc alloc(size_t size, size_t align) override
    {
        bool new_block = false;
        auto rtaddr = m_tracker.alloc(size, align, [&] (size_t bsize) {
                new_block = true;
                auto new_alloc = m_dualmap.alloc(bsize, exec);
                auto res = m_block_infos.insert({bsize + (char*)new_alloc.first,
                            BlockInfo{bsize, new_alloc.second, nullptr}});
                assert(res.second);
                m_cur_blocks.push_back(res.first);
                return new_alloc.first;
            });
        if (new_block)
            return {rtaddr, rtaddr};
        auto info_it = m_block_infos.upper_bound(rtaddr);
        assert(info_it != m_block_infos.end());
        auto &info = info_it->second;
        auto offset = (char*)rtaddr - ((char*)info_it->first - info_it->second.size);
        assert(offset >= 0);
        if (!info.wraddr)
            return {rtaddr, rtaddr};
        if (info.wraddr == (void*)-1)
            info.wraddr = m_dualmap.remap_wraddr(info.id, info.size);
        return {(char*)info.wraddr + offset, rtaddr};
    }
    void finalize() override
    {
        for (auto info_it: m_cur_blocks) {
            auto &info = info_it->second;
            assert(!info.wraddr || info.wraddr == (void*)-1);
            // already handled
            if (info.wraddr)
                continue;
            if (!protectPage((char*)info_it->first - info.size, info.size,
                             exec ? Prot::RX : Prot::RO))
                checkErrno(-1, "Cannot set page protection");
            info.wraddr = (void*)-1;
        }
        m_cur_blocks.clear();
    }
    void free(ROAlloc alloc, size_t size) override
    {
        m_tracker.free(alloc.rtaddr, size, [&] (void *ptr, size_t bsize) {
                auto it = free_block(ptr, bsize);
                m_block_infos.erase(it);
            });
    }
    ~DualMapAllocator() override
    {
        m_tracker.free_all([&] (void *ptr, size_t bsize) {
                free_block(ptr, bsize);
            });
    }

private:
    struct BlockInfo {
        // Block size
        size_t size;
        // ID for dual map allocation
        uintptr_t id;
        // Writeable address
        // 0 means that this is a new allocation and the runtime address is still writable
        // -1 means that the runtime address is not writable anymore but the writable
        // address haven't be allocated yet.
        void *wraddr;
    };
    typename std::map<void*,BlockInfo>::iterator free_block(void *ptr, size_t bsize)
    {
        void *key = (char*)ptr + bsize;
        auto it = m_block_infos.find(key);
        assert(it != m_block_infos.end());
        auto wraddr = it->second.wraddr;
        if (wraddr == (void*)-1)
            wraddr = nullptr;
        m_dualmap.free(ptr, it->second.id, bsize, wraddr);
        return it;
    }
    DualMap &m_dualmap;
    AllocTracker m_tracker;
    std::map<void*,BlockInfo> m_block_infos; // keys are the end addresses
    SmallVector<typename std::map<void*,BlockInfo>::iterator,3> m_cur_blocks;
};

template<bool exec>
class MemWriterAllocator : public ROAllocator<exec> {
public:
    MemWriterAllocator(MemWriter &mem_writer, size_t block_size)
        : m_writer(mem_writer),
          m_tracker(block_size),
          m_wrtracker(block_size)
    {
        if (!m_writer.init()) {
            throw std::runtime_error("Failed to initialize MemWriter.");
        }
    }
    ROAlloc alloc(size_t size, size_t align) override
    {
        bool new_block = false;
        auto rtaddr = m_tracker.alloc(size, align, [&] (size_t bsize) {
                new_block = true;
                auto ptr = mapAnonPage(bsize, Prot::RW);
                m_new_blocks.emplace_back(ptr, bsize);
                return ptr;
            });
        if (new_block)
            return {rtaddr, rtaddr};
        for (auto nb: m_new_blocks) {
            if (nb.first <= rtaddr && (char*)nb.first + nb.second > (char*)rtaddr) {
                return {rtaddr, rtaddr};
            }
        }
        auto wraddr = alloc_wraddr(size, align);
        m_allocs.push_back({wraddr, rtaddr, size});
        return {wraddr, rtaddr};
    }
    void finalize() override
    {
        for (auto nb: m_new_blocks) {
            if (!protectPage(nb.first, nb.second, exec ? Prot::RX : Prot::RO)) {
                checkErrno(-1, "Cannot set page protection");
            }
        }
        m_new_blocks.clear();
        for (auto alloc: m_allocs)
            m_writer.write(alloc.rtaddr, alloc.wraddr, alloc.size);
        m_allocs.clear();
    }
    void free(ROAlloc alloc, size_t size) override
    {
        m_tracker.free(alloc.rtaddr, size, unmapPage);
    }
    ~MemWriterAllocator() override
    {
        m_tracker.free_all(unmapPage);
        m_wrtracker.free_all(unmapPage);
    }

private:
    void *alloc_wraddr(size_t size, size_t align)
    {
        return m_wrtracker.alloc(size, align, [&] (size_t size) {
                return mapAnonPage(size, Prot::RW);
            });
    }
    struct Alloc {
        void *wraddr;
        void *rtaddr;
        size_t size;
    };
    MemWriter &m_writer;
    AllocTracker m_tracker;
    AllocTracker m_wrtracker;
    SmallVector<std::pair<void*,size_t>,2> m_new_blocks;
    SmallVector<Alloc,3> m_allocs;
};

MemMgr::MemMgr()
    : m_blocksz(alignTo(8 * 1024 * 1024, page_size)),
      m_rwalloc(m_blocksz)
{
    if (m_memwriter.init()) {
        m_roalloc.reset(new MemWriterAllocator<false>(m_memwriter, m_blocksz));
        m_rxalloc.reset(new MemWriterAllocator<true>(m_memwriter, m_blocksz));
    }
    else {
        m_roalloc.reset(new DualMapAllocator<false>(m_dualmap, m_blocksz));
        m_rxalloc.reset(new DualMapAllocator<true>(m_dualmap, m_blocksz));
    }
}

MemMgr::~MemMgr()
{
}

uint64_t MemMgr::new_group()
{
    auto gid = ++m_grp_cnt;
    m_cur_allocs = &m_allocs[gid];
    return gid;
}

void MemMgr::free_group(uint64_t id)
{
    auto it = m_allocs.find(id);
    if (it == m_allocs.end())
        return;
    for (auto alloc: it->second) {
        switch (alloc.type) {
        case Alloc::RW:
            m_rwalloc.free(alloc.rtaddr, alloc.size);
            break;
        case Alloc::RO:
            m_roalloc->free({alloc.wraddr, alloc.rtaddr}, alloc.size);
            break;
        case Alloc::RX:
            m_rxalloc->free({alloc.wraddr, alloc.rtaddr}, alloc.size);
            break;
        default:
            abort();
        }
    }
    m_allocs.erase(it);
}

uint8_t *MemMgr::allocateCodeSection(uintptr_t sz, unsigned align, unsigned, StringRef)
{
    assert(m_cur_allocs);
    auto res = m_rxalloc->alloc(sz, align);
    m_cur_allocs->push_back(Alloc{res.rtaddr, sz, Alloc::RX, res.wraddr});
    return (uint8_t*)res.wraddr;
}

uint8_t *MemMgr::allocateDataSection(uintptr_t sz, unsigned align,
                                     unsigned, StringRef, bool ro)
{
    assert(m_cur_allocs);
    if (!ro) {
        auto ptr = m_rwalloc.alloc(sz, align);
        m_cur_allocs->push_back(Alloc{ptr, sz, Alloc::RW});
        return (uint8_t*)ptr;
    }
    auto res = m_roalloc->alloc(sz, align);
    m_cur_allocs->push_back(Alloc{res.rtaddr, sz, Alloc::RO, res.wraddr});
    return (uint8_t*)res.wraddr;
}

void MemMgr::registerEHFrames(uint8_t*, uint64_t, size_t)
{
}

void MemMgr::deregisterEHFrames()
{
}

bool MemMgr::finalizeMemory(std::string*)
{
    assert(m_cur_allocs);
    m_roalloc->finalize();
    m_rxalloc->finalize();

    for (auto &alloc: *m_cur_allocs) {
        // ensure the mapped pages are consistent
        sys::Memory::InvalidateInstructionCache(alloc.rtaddr, alloc.size);
        if (alloc.rtaddr != alloc.wraddr && alloc.wraddr) {
            sys::Memory::InvalidateInstructionCache(alloc.wraddr, alloc.size);
        }
    }
    return false;
}

void MemMgr::notifyObjectLoaded(RuntimeDyld &dyld, const object::ObjectFile&)
{
    for (auto &alloc: *m_cur_allocs) {
        if (alloc.rtaddr == alloc.wraddr || !alloc.wraddr)
            continue;
        dyld.mapSectionAddress(alloc.wraddr, (uintptr_t)alloc.rtaddr);
    }
}

}
