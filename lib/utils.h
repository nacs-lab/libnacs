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

#include "macros.h"

#ifndef __NACS_UTILS_UTILS_H__
#define __NACS_UTILS_UTILS_H__

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
// GNU extension
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef __cplusplus
template<typename T>
static NACS_INLINE T
nacsSetBit(T orig, uint8_t bit, bool val)
{
    if (val) {
        return orig | (static_cast<T>(1) << bit);
    } else {
        return orig & ~(static_cast<T>(1) << bit);
    }
}
#else
#define nacsSetBit(orig, bit, val)              \
    ((val) ?                                    \
     ((orig) | (((typeof(orig))1) << (bit))) :  \
     ((orig) & ~(((typeof(orig))1) << (bit))))
#endif

#endif
