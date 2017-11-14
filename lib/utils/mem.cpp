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

#include <sys/mman.h>
#include <unistd.h>

namespace NaCs {

NACS_EXPORT void*
mapPhyAddr(void *phy_addr, size_t len)
{
    static int fd = open("/dev/mem", O_RDWR | O_SYNC);
    return mapFile(fd, off_t(phy_addr), len);
}

NACS_EXPORT void*
getPhyAddr(void *virt_addr)
{
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
}

}
