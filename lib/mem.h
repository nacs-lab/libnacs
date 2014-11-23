/*************************************************************************
 *   Copyright (c) 2014 - 2014 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_UTILS_MEM_H__
#define __NACS_UTILS_MEM_H__

#define mem_barrier() asm volatile("" ::: "memory")

#define MEMDEF_RW(type, vtype, suffix)                  \
    typedef volatile type vtype;                        \
    static NACS_INLINE type                             \
    mem_read##suffix(volatile void *addr)               \
    {                                                   \
        return *(vtype*)addr;                           \
    }                                                   \
    static NACS_INLINE void                             \
    mem_write##suffix(volatile void *addr, vtype val)   \
    {                                                   \
        *(vtype*)addr = val;                            \
    }

MEMDEF_RW(uint8_t, vuint8, 8)
MEMDEF_RW(uint16_t, vuint16, 16)
MEMDEF_RW(uint32_t, vuint32, 32)
MEMDEF_RW(uint64_t, vuint64, 64)

#undef MEMDEF_RW

#endif
