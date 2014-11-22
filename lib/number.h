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

#ifndef _NACS_UTILS_NUMBER_H_
#define _NACS_UTILS_NUMBER_H_

#include "utils.h"

#ifdef __cplusplus

template <typename T1, typename T2>
NACS_INLINE static auto
nacsMax(const T1 &a, const T2 &b) -> decltype((a > b) ? a : b)
{
    return (a > b) ? a : b;
}

template <typename T1, typename T2>
NACS_INLINE static auto
nacsMin(const T1 &a, const T2 &b) -> decltype((a < b) ? a : b)
{
    return (a < b) ? a : b;
}

template <typename T>
NACS_INLINE static T
nacsAbs(const T &a)
{
    return (a > 0) ? a : -a;
}
#else
#define nacsMax(a, b)                           \
    ({                                          \
        typeof(a) _a = (a);                     \
        typeof(b) _b = (b);                     \
        (_a > _b) ? _a : _b;                    \
    })
#define nacsMin(a, b)                           \
    ({                                          \
        typeof(a) _a = (a);                     \
        typeof(b) _b = (b);                     \
        (_a < _b) ? _a : _b;                    \
    })
#define nacsAbs(a)                              \
    ({                                          \
        typeof(a) _a = (a);                     \
        (_a > 0) ? _a : -_a;                    \
    })
#endif
#define nacsBound(a, b, c) nacsMax(a, nacsMin(b, c))
#define nacsLimit(v, l) nacsBound(0, v, l)

NACS_INLINE static uintptr_t
nacsGetPadding(uintptr_t len, uintptr_t align)
{
    uintptr_t left;
    if ((left = len % align)) {
        return align - left;
    }
    return 0;
}

NACS_INLINE static uintptr_t
nacsAlignTo(uintptr_t len, uintptr_t align)
{
    return len + nacsGetPadding(len, align);
}

#endif
