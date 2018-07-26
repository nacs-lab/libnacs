/*************************************************************************
 *   Copyright (c) 2014 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include <cstring>

#if NACS_OS_WINDOWS
#  include <windows.h>
#else
#  include <sys/mman.h>
#endif

#ifndef __NACS_UTILS_MEM_H__
#define __NACS_UTILS_MEM_H__

namespace NaCs {

namespace Mem {

template<typename T> static inline T
load_unalign(const void *p, size_t idx=0)
{
    T v;
    std::memcpy(&v, (const char*)p + sizeof(T) * idx, sizeof(T));
    return v;
}

template<typename T> static inline void
store_unalign(void *p, T v, size_t idx=0)
{
    std::memcpy((const char*)p + sizeof(T) * idx, &v, sizeof(T));
}

template<typename T> static inline T
read(volatile void *addr, size_t idx=0)
{
    return ((volatile T*)addr)[idx];
}

template<typename T, typename T2> static inline void
write(volatile void *addr, T2 val, size_t idx=0)
{
    ((volatile T*)addr)[idx] = val;
}
}

NACS_EXPORT(utils) extern const size_t page_size;
void *mapPhyAddr(void*, size_t);
void *getPhyAddr(void*);
#if NACS_OS_WINDOWS
enum class Prot : int {
    RW = PAGE_READWRITE,
    RX = PAGE_EXECUTE,
    RO = PAGE_READONLY
};
#else
enum class Prot : int {
    RW = PROT_READ | PROT_WRITE,
    RX = PROT_READ | PROT_EXEC,
    RO = PROT_READ
};
#endif
// Map anonymous pages of `size`. Guaranteed to be aligned to page size.
void *mapAnonPage(size_t size, Prot prot);
// Unmap pages.
// Windows does not support partially unmap a previously allocated region.
// If partially freeing the memory is needed, use `decommitPage` instead.
void unmapPage(void *ptr, size_t size);
// Change the protection on the page
bool protectPage(void *ptr, size_t size, Prot prot);

// Mark the page as unused. This does not free the address range but allows the OS
// to reused the physical memory. On Unix, the page is still accessible after the call.
void decommitPage(void *ptr, size_t size);
// Recommit a previously decommited page and set the page protection to `prot`.
// This is required on windows after a page is decommited.
bool recommitPage(void *ptr, size_t size, Prot prot);

class DualMap {
public:
    DualMap(bool do_init)
    {
        if (do_init) {
            init();
        }
    }
    void init();
private:
#if !NACS_OS_WINDOWS
    int m_fd = -1;
    size_t m_offset = 0;
    size_t m_size = 0;
#endif
};

// Pre-allocate memory for a number of objects so that they can be allocated quickly.
template<typename T, size_t n>
class SmallAllocator {
    struct alignas(T) Ceil {
        char buff[sizeof(T)];
    };
    Ceil m_buff[n];
    static constexpr auto nbits32 = (n + 31) / 32;
    uint32_t m_bits[nbits32];
    size_t m_cur = 0;
    Ceil *find_ceil()
    {
        if (n == 0)
            return nullptr;
        const auto old_cur = m_cur;
        auto cur = old_cur;
        do {
            auto bits = m_bits[cur];
            auto idx = __builtin_ffs(bits);
            if (idx) {
                idx -= 1;
                m_bits[cur] = bits & ~(1 << idx);
                m_cur = cur;
                return &m_buff[idx + cur * 32];
            }
            cur++;
            if (cur == nbits32) {
                cur = 0;
            }
        } while (cur != old_cur);
        return nullptr;
    }
    NACS_INLINE Ceil *alloc_mem()
    {
        if (auto cached = find_ceil())
            return cached;
        return new Ceil;
    }
    NACS_INLINE void free_mem(Ceil *mem)
    {
        if (n && m_buff <= mem && mem < m_buff + n) {
            auto idx = mem - m_buff;
            m_bits[idx / 32] |= 1 << (idx % 32);
        }
        else {
            delete mem;
        }
    }
public:
    SmallAllocator()
    {
        memset(&m_bits, 0xff, sizeof(m_bits));
        if (auto nleft = n % 32) {
            m_bits[nbits32 - 1] = (1 << nleft) - 1;
        }
    }
    template<typename... Args>
    NACS_INLINE T *alloc(Args&&... args)
    {
        auto mem = alloc_mem();
        return new (&mem->buff) T(std::forward<Args>(args)...);
    }
    NACS_INLINE void free(T *p)
    {
        p->~T();
        free_mem(reinterpret_cast<Ceil*>(p));
    }
};

}

#endif
