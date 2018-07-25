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

#include "execute.h"

namespace NaCs {
namespace LLVM {
namespace Exe {

/**
 */

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
    auto ptr = block_alloc(blksz);
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
            block_free(ptr, sz);
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

MemMgr::MemMgr()
{
}

MemMgr::~MemMgr()
{
}

uint64_t MemMgr::new_group()
{
    return ++m_grp_cnt;
}

void MemMgr::free_group(uint64_t id)
{
    (void)id;
}

uint8_t *MemMgr::allocateCodeSection(uintptr_t sz, unsigned align,
                                     unsigned sid, StringRef sname)
{
    if (!m_grp_cnt)
        abort();
    return tmp.allocateCodeSection(sz, align, sid, sname);
}

uint8_t *MemMgr::allocateDataSection(uintptr_t sz, unsigned align,
                                     unsigned sid, StringRef sname, bool ro)
{
    if (!m_grp_cnt)
        abort();
    return tmp.allocateDataSection(sz, align, sid, sname, ro);
}

void MemMgr::registerEHFrames(uint8_t*, uint64_t, size_t)
{
}

void MemMgr::deregisterEHFrames()
{
}

bool MemMgr::finalizeMemory(std::string *ErrMsg)
{
    if (!m_grp_cnt)
        abort();
    return tmp.finalizeMemory(ErrMsg);
}

}
}
}
