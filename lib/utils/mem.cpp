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

#include "fd_utils.h"
#include "mem.h"
#include "number.h"

#if !NACS_OS_WINDOWS
#  include <unistd.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  if NACS_OS_DARWIN && !defined(MAP_ANONYMOUS)
#    define MAP_ANONYMOUS MAP_ANON
#  endif
#endif
#if NACS_OS_LINUX
#  include <sys/syscall.h>
#elif NACS_OS_FREEBSD
#  include <sys/types.h>
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
    (void)size;
    VirtualFree(ptr, 0, MEM_RELEASE);
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

void decommitPage(void *ptr, size_t size)
{
#if NACS_OS_WINDOWS
    VirtualFree(ptr, size, MEM_DECOMMIT);
#else
    madvise(ptr, size, MADV_DONTNEED);
#endif
}

bool recommitPage(void *ptr, size_t size, Prot prot)
{
#if NACS_OS_WINDOWS
    return VirtualAlloc(ptr, size, MEM_COMMIT, (int)prot);
#else
    return mprotect(ptr, size, (int)prot) == 0;
#endif
}

#if NACS_OS_WINDOWS
NACS_EXPORT() bool DualMap::init()
{
    return true;
}

NACS_EXPORT() std::pair<void*,uintptr_t> DualMap::alloc(size_t size, bool exec)
{
    // As far as I can tell `CreateFileMapping` cannot be resized on windows.
    // Also, creating big file mapping and then map pieces of it seems to
    // consume too much global resources. Therefore, we use each file mapping
    // as a block on windows
    DWORD file_mode = exec ? PAGE_EXECUTE_READWRITE : PAGE_READWRITE;
    HANDLE hdl = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
                                   file_mode, 0, size, NULL);
    // We set the maximum permissions for this to the maximum for this file, and then
    // VirtualProtect, such that the debugger can still access these
    // pages and set breakpoints if it wants to.
    DWORD map_mode = FILE_MAP_ALL_ACCESS | (exec ? FILE_MAP_EXECUTE : 0);
    void *addr = MapViewOfFile(hdl, map_mode, 0, 0, size);
    if (!addr)
        throw std::runtime_error("Cannot map initial view.");
    // Return the initially read-write page.
    VirtualProtect(addr, size, PAGE_READWRITE, &file_mode);
    return std::make_pair(addr, (uintptr_t)hdl);
}

NACS_EXPORT() void *DualMap::remap_wraddr(uintptr_t id, size_t size)
{
    void *addr = MapViewOfFile((HANDLE)id, FILE_MAP_ALL_ACCESS,
                               0, 0, size);
    if (!addr)
        throw std::runtime_error("Cannot map RW view.");
    return addr;
}

NACS_EXPORT() void DualMap::free(void *ptr, uintptr_t id, size_t, void *wraddr)
{
    if (wraddr)
        UnmapViewOfFile(wraddr);
    UnmapViewOfFile(ptr);
    CloseHandle((HANDLE)id);
}
#else
bool DualMap::checkFdOrClose(int fd)
{
    if (fd == -1)
        return false;
    fcntl(fd, F_SETFD, FD_CLOEXEC);
    fchmod(fd, S_IRWXU);
    if (ftruncate(fd, m_region_sz) != 0) {
        close(fd);
        return false;
    }
    // This can fail due to `noexec` mount option ....
    void *ptr = mmap(nullptr, page_size, PROT_READ | PROT_EXEC, MAP_SHARED, fd, 0);
    if (ptr == MAP_FAILED) {
        close(fd);
        return false;
    }
    munmap(ptr, page_size);
    return true;
}

NACS_EXPORT() bool DualMap::init()
{
    if (m_fd != -1)
        return true;
    // Linux and FreeBSD can create an anonymous fd without touching the
    // file system.
#  ifdef __NR_memfd_create
    m_fd = (int)syscall(__NR_memfd_create, "nacs-utils", 0);
    if (checkFdOrClose(m_fd))
        return true;
#  endif
#  if NACS_OS_FREEBSD
    m_fd = shm_open(SHM_ANON, O_RDWR, S_IRWXU);
    if (checkFdOrClose(m_fd))
        return true;
#  endif
    char shm_name[] = "nacs-utils-0123456789-0123456789/tmp///";
    pid_t pid = getpid();
#  if !NACS_OS_DARWIN
    // `shm_open` can't be mapped exec on mac
    do {
        snprintf(shm_name, sizeof(shm_name),
                 "nacs-utils-%d-%d", (int)pid, rand());
        m_fd = shm_open(shm_name, O_RDWR | O_CREAT | O_EXCL, S_IRWXU);
        if (checkFdOrClose(m_fd)) {
            shm_unlink(shm_name);
            return true;
        }
    } while (errno == EEXIST);
#  endif
    FILE *tmpf = tmpfile();
    if (tmpf) {
        m_fd = dup(fileno(tmpf));
        fclose(tmpf);
        if (checkFdOrClose(m_fd)) {
            return true;
        }
    }
    snprintf(shm_name, sizeof(shm_name), "/tmp/nacs-utils-%d-XXXXXX", (int)pid);
    m_fd = mkstemp(shm_name);
    if (checkFdOrClose(m_fd)) {
        unlink(shm_name);
        return true;
    }
    m_fd = -1;
    return false;
}

NACS_EXPORT() std::pair<void*,uintptr_t> DualMap::alloc(size_t size, bool)
{
    // Find an free region first
    for (auto it = m_freeregions.begin(), end = m_freeregions.end(); it != end; ++it) {
        auto rsz = it->second;
        if (unlikely(rsz < size))
            continue;
        auto offset = it->first - rsz;
        if (rsz == size) {
            m_freeregions.erase(it);
        }
        else {
            it->second = rsz - size;
        }
        return std::make_pair(remap_wraddr(offset, size), offset);
    }
    // Allocate from the end and make sure the file is large enough.
    auto offset = m_maxoffset;
    m_maxoffset += size;
    if (unlikely(m_maxoffset > m_filesize)) {
        m_filesize = alignTo(m_maxoffset, m_region_sz);
        checkErrno(ftruncate(m_fd, m_filesize));
    }
    return std::make_pair(remap_wraddr(offset, size), offset);
}

NACS_EXPORT() void *DualMap::remap_wraddr(uintptr_t id, size_t size)
{
    void *addr = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_SHARED, m_fd, id);
    if (addr == MAP_FAILED)
        throw std::runtime_error("Failed to map RW view.");
    return addr;
}

NACS_EXPORT() void DualMap::free(void *ptr, uintptr_t id, size_t size, void *wraddr)
{
    // We don't need anything from the pointers so just free them directly.
    if (wraddr)
        munmap(wraddr, size);
    munmap(ptr, size);
    // Now record that the offset region is free.
    auto it = m_freeregions.lower_bound(id);
    auto end = id + size;
    if (it == m_freeregions.end()) {
        // Cannot merge with any free regions in the list.
        // Check if it can be merged with the free region at the end.
        if (end == m_maxoffset)
            goto free_end;
        m_freeregions.insert(it, std::make_pair(end, size));
        goto free_mem;
    }
    if (it->first == id) {
        // We found a previous one to merge with.
        id = id - it->second;
        size = it->second + size;
        it = m_freeregions.erase(it);
        if (it == m_freeregions.end()) {
            // Check if it can be merged with the free region at the end.
            if (end == m_maxoffset)
                goto free_end;
            m_freeregions.insert(it, std::make_pair(end, size));
            goto free_mem;
        }
    }
    // Now `it->first` is guaranteed to be greater than id,
    // meaning that we definitely won't merge with the empty region at the end.
    if (it->first - it->second == end) {
        it->second += size;
    }
    else {
        m_freeregions.insert(it, std::make_pair(end, size));
    }
    goto free_mem;
free_end:
    m_maxoffset = id;
    if (m_filesize - m_maxoffset > m_region_sz) {
        m_filesize = alignTo(m_maxoffset, m_region_sz);
        ftruncate(m_fd, m_filesize);
        return;
    }
free_mem:
#  if NACS_OS_LINUX
    // Free up the memory
    if (fallocate(m_fd, FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE, id, size) == -1) {
        // ignore error
        fallocate(m_fd, FALLOC_FL_ZERO_RANGE | FALLOC_FL_KEEP_SIZE, id, size);
    }
#  else
    return;
#  endif
}
#endif

}
