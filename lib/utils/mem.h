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
read(volatile void *addr)
{
    return *(volatile T*)addr;
}

template<typename T, typename T2> static inline void
write(volatile void *addr, T2 val)
{
    *(volatile T*)addr = val;
}
}

void *mapPhyAddr(void*, size_t);
void *getPhyAddr(void*);

}

#endif
