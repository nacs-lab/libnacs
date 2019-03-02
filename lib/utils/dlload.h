/*************************************************************************
 *   Copyright (c) 2013 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef _NACS_UTILS_DLLOAD_H_
#define _NACS_UTILS_DLLOAD_H_

#include "utils.h"

namespace NaCs {

namespace DL {

enum Flags {
    LOCAL = 1,
    GLOBAL = 2,
    LAZY = 4,
    NOW = 8,
    /* Linux/glibc and MacOS X: */
    NODELETE = 16,
    NOLOAD = 32,
    /* Linux/glibc: */
    DEEPBIND = 64,
    /* MacOS X 10.5+: */
    FIRST = 128
};

NACS_EXPORT(utils) void *open(const char *filename, int flags);
NACS_EXPORT(utils) bool close(void *handle);
NACS_EXPORT(utils) void *sym(void *handle, const char *symbol);

}

}

#endif
