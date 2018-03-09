/*************************************************************************
 *   Copyright (c) 2013 - 2014 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef _NACS_UTILS_FD_UTILS_H_
#define _NACS_UTILS_FD_UTILS_H_

#include "utils.h"
#include "macros.h"

#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <assert.h>
#if NACS_OS_LINUX
#  include <sys/eventfd.h>
#endif

#include <stdexcept>
#include <system_error>

namespace NaCs {

void *mapFile(const char *name, off_t offset, size_t len);
void *mapFile(int fd, off_t offset, size_t len);
bool sendFD(int sock, int fd);
int recvFD(int sock);
bool fdSetCloexec(int fd, bool cloexec);
bool fdSetNonBlock(int fd, bool nonblock);

#if !NACS_OS_WINDOWS
class FLock {
    int m_fd;
public:
    NACS_INLINE
    FLock(int fd) : m_fd(fd)
    {
        if (fd < 0) {
            throw std::runtime_error("Invalid FD.");
        }
    }
    NACS_INLINE
    FLock(const char *fname) :
        FLock(open(fname, O_RDWR | O_CREAT, 0644))
    {}
    NACS_INLINE void
    lock()
    {
        if (flock(m_fd, LOCK_EX) == -1) {
            throw std::runtime_error("Failed to acquire lock.");
        }
    }
    NACS_INLINE void
    unlock() noexcept
    {
        flock(m_fd, LOCK_UN);
    }
};
#endif

template<typename T, typename... Arg>
static inline T checkErrno(T res, Arg&&... arg)
{
    if (res == -1)
        throw std::system_error(errno, std::system_category(), std::forward<Arg>(arg)...);
    return res;
}

#if NACS_OS_LINUX
// Helper functions for event fd.

// Open a eventfd and check for error
static inline int openEvent(int init, int flags)
{
    return checkErrno(eventfd(init, flags), "Cannot open eventfd.");
}

// Write certain number of events (default to 1) to the event fd.
// Retry if aborted by signal and return if the write succeeded.
static inline bool writeEvent(int fd, uint64_t numevt=1)
{
    while (true) {
        auto res = write(fd, &numevt, sizeof(uint64_t));
        if (res == -1) {
            if (errno == EINTR)
                continue;
            return false;
        }
        assert(res == 8);
        return true;
    }
}

// Read events from the event fd.
// Retry if aborted by signal and return the number of events read.
static inline uint64_t readEvent(int fd)
{
    while (true) {
        uint64_t numevt;
        auto res = read(fd, &numevt, sizeof(uint64_t));
        if (res == -1) {
            if (errno == EINTR)
                continue;
            return 0;
        }
        assert(res == 8);
        return numevt;
    }
}
#endif

}

#endif
