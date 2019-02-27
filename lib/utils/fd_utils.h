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
#include "errors.h"

#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <assert.h>
#if NACS_OS_LINUX
#  include <sys/eventfd.h>
#endif

namespace NaCs {

/**
 * Unix only functions that manipulate files
 */

// offset does **NOT** have to be integer multiple of page size.
NACS_EXPORT(utils) void *mapFile(const char *name, off_t offset, size_t len);
NACS_EXPORT(utils) void *mapFile(int fd, off_t offset, size_t len);
NACS_EXPORT(utils) bool sendFD(int sock, int fd);
NACS_EXPORT(utils) int recvFD(int sock);
NACS_EXPORT(utils) bool fdSetCloexec(int fd, bool cloexec);
NACS_EXPORT(utils) bool fdSetNonBlock(int fd, bool nonblock);

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
