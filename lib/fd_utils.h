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

NACS_BEGIN_DECLS

void *nacsMapFile(const char *name, off_t offset, size_t len);
bool nacsSendFD(int sock, int fd);
int nacsRecvFD(int sock);
bool nacsFDSetCloexec(int fd, bool cloexec);
bool nacsFDSetNonBlock(int fd, bool nonblock);

NACS_END_DECLS

#endif
