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

#include "fd_utils.h"
#include "log.h"
#include "number.h"
#include "macros.h"

#ifndef NACS_OS_WINDOWS_
#  include <sys/socket.h>
#  include <unistd.h>
#  include <sys/mman.h>
#  include <fcntl.h>
#endif

namespace NaCs {

NACS_EXPORT void*
mapFile(int fd, off_t offset, size_t len)
{
#ifndef NACS_OS_WINDOWS_
    static const auto page_size = sysconf(_SC_PAGESIZE);
    NACS_RET_IF_FAIL(fd >= 0, nullptr);
    off_t start = offset - offset % page_size;
    off_t end = alignTo(offset + len, page_size);

    void *base = mmap(nullptr, end - start, PROT_READ | PROT_WRITE,
                      MAP_SHARED, fd, start);
    NACS_RET_IF_FAIL(base != (void*)-1, nullptr);

    return (char*)base + offset - start;
#else
    (void)fd;
    (void)offset;
    (void)len;
    return nullptr;
#endif
}

NACS_EXPORT void*
mapFile(const char *name, off_t offset, size_t len)
{
    NACS_RET_IF_FAIL(name, nullptr);
    int fd = open(name, O_RDWR | O_SYNC);
    auto res = mapFile(fd, offset, len);
    close(fd);
    return res;
}

NACS_EXPORT bool
sendFD(int sock, int fd)
{
#ifndef NACS_OS_WINDOWS_
    NACS_RET_IF_FAIL(fd >= 0 && sock >= 0, false);
    char buf = 0;
    struct iovec iov = {
        .iov_base = &buf,
        .iov_len = 1
    };
    union {
        struct cmsghdr cmsghdr;
        char control[CMSG_SPACE(sizeof(int))];
    } cmsgu;
    memset(&cmsgu, 0, sizeof(cmsgu));
    struct msghdr msg = {
        .msg_name = nullptr,
        .msg_namelen = 0,
        .msg_iov = &iov,
        .msg_iovlen = 1,
        .msg_control = cmsgu.control,
        .msg_controllen = sizeof(cmsgu.control),
        .msg_flags = 0,
    };
    struct cmsghdr *cmsg = CMSG_FIRSTHDR(&msg);
    cmsg->cmsg_len = CMSG_LEN(sizeof(int));
    cmsg->cmsg_level = SOL_SOCKET;
    cmsg->cmsg_type = SCM_RIGHTS;
    memcpy(CMSG_DATA(cmsg), &fd, sizeof(int));
    return sendmsg(sock, &msg, 0) >= 0;
#else
    (void)sock;
    (void)fd;
    return false;
#endif
}

NACS_EXPORT int
recvFD(int sock)
{
#ifndef NACS_OS_WINDOWS_
    NACS_RET_IF_FAIL(sock >= 0, -1);
    char buf = 0;
    struct iovec iov = {
        .iov_base = &buf,
        .iov_len = 1
    };
    union {
        struct cmsghdr cmsghdr;
        char control[CMSG_SPACE(sizeof(int))];
    } cmsgu;
    memset(&cmsgu, 0, sizeof(cmsgu));
    struct msghdr msg = {
        .msg_name = nullptr,
        .msg_namelen = 0,
        .msg_iov = &iov,
        .msg_iovlen = 1,
        .msg_control = cmsgu.control,
        .msg_controllen = sizeof(cmsgu.control),
        .msg_flags = 0,
    };
    NACS_RET_IF_FAIL(recvmsg(sock, &msg, 0) >= 0, -1);
    struct cmsghdr *cmsg = CMSG_FIRSTHDR(&msg);
    NACS_RET_IF_FAIL(cmsg && cmsg->cmsg_len == CMSG_LEN(sizeof(int)) &&
                     cmsg->cmsg_level == SOL_SOCKET &&
                     cmsg->cmsg_type == SCM_RIGHTS, -1);
    int fd;
    memcpy(&fd, CMSG_DATA(cmsg), sizeof(int));
    return fd;
#else
    (void)sock;
    return false;
#endif
}

NACS_EXPORT bool
fdSetCloexec(int fd, bool cloexec)
{
#ifndef NACS_OS_WINDOWS_
    long flags;
    flags = fcntl(fd, F_GETFD, 0);
    if (flags == -1) {
        return false;
    }
    if (cloexec) {
        flags |= FD_CLOEXEC;
    } else {
        flags &= ~FD_CLOEXEC;
    }
    return fcntl(fd, F_SETFD, flags) != -1;
#else
    (void)fd;
    (void)cloexec;
    return false;
#endif
}

}

NACS_EXPORT bool
fdSetNonBlock(int fd, bool nonblock)
{
#ifndef NACS_OS_WINDOWS_
    long flags;
    flags = fcntl(fd, F_GETFL, 0);
    if (flags == -1) {
        return false;
    }
    if (nonblock) {
        nonblock |= O_NONBLOCK;
    } else {
        nonblock &= ~O_NONBLOCK;
    }
    return fcntl(fd, F_SETFL, flags) != -1;
#else
    (void)fd;
    (void)nonblock;
    return false;
#endif
}
