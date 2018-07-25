/*************************************************************************
 *   Copyright (c) 2015 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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

#include "mem.h"
#include "fd_utils.h"

#if NACS_OS_WINDOWS
#  include <windows.h>
#else
#  include <unistd.h>
#endif

namespace NaCs {

#if NACS_OS_WINDOWS
static size_t win_get_pagesz()
{
    SYSTEM_INFO systemInfo;
    GetSystemInfo(&systemInfo);
    return systemInfo.dwPageSize;
}
NACS_EXPORT() extern const size_t page_size = win_get_pagesz();
#else
NACS_EXPORT() extern const size_t page_size = sysconf(_SC_PAGESIZE);
#endif

NACS_EXPORT() void*
mapPhyAddr(void *phy_addr, size_t len)
{
#if NACS_OS_LINUX
    static int fd = open("/dev/mem", O_RDWR | O_SYNC);
    return mapFile(fd, off_t(phy_addr), len);
#else
    (void)phy_addr;
    (void)len;
    return nullptr;
#endif
}

NACS_EXPORT() void*
getPhyAddr(void *virt_addr)
{
#if NACS_OS_LINUX
    static int page_map = open("/proc/self/pagemap", O_RDONLY);
    static uint32_t page_size = getpagesize();
    uintptr_t virt_offset = uintptr_t(virt_addr) % page_size;
    uintptr_t virt_pfn = uintptr_t(virt_addr) / page_size;

    uint64_t page_info;
    auto res = pread(page_map, &page_info, sizeof(page_info),
                     virt_pfn * sizeof(uint64_t));
    (void)res;
    auto phy_pfn = uintptr_t(page_info & ((1ll << 55) - 1));
    auto phy_page = phy_pfn * page_size;
    return (void*)(phy_page + virt_offset);
#else
    (void)virt_addr;
    return nullptr;
#endif
}

NACS_EXPORT() void *mapAnonPage(size_t size, Prot prot)
{
#if NACS_OS_WINDOWS
    void *mem = VirtualAlloc(NULL, size, MEM_COMMIT, (int)prot);
#else
    void *mem = mmap(nullptr, size, (int)prot,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (mem == MAP_FAILED)
        mem = nullptr;
#endif
    return mem;
}

NACS_EXPORT() void unmapPage(void *ptr, size_t size)
{
#if NACS_OS_WINDOWS
    VirtualFree(ptr, size, MEM_DECOMMIT);
#else
    munmap(ptr, size);
#endif
}

NACS_EXPORT() bool protectPage(void *ptr, size_t size, Prot prot)
{
#if NACS_OS_WINDOWS
    DWORD old_prot;
    return VirtualProtect(ptr, size, (DWORD)prot, &old_prot);
#else
    return mprotect(ptr, size, (int)prot) == 0;
#endif
}

}
