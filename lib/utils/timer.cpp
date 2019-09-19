/*************************************************************************
 *   Copyright (c) 2013 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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

#include "timer.h"

#include <time.h>

#include <vector>

namespace NaCs {

NACS_EXPORT() uint64_t getTime()
{
    timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return uint64_t(t.tv_sec) * 1000 * 1000 * 1000 + t.tv_nsec;
}

NACS_EXPORT() uint64_t getRes()
{
    timespec t;
    clock_getres(CLOCK_MONOTONIC, &t);
    return uint64_t(t.tv_sec) * 1000 * 1000 * 1000 + t.tv_nsec;
}

#ifndef CLOCK_MONOTONIC_COARSE
#  define CLOCK_MONOTONIC_COARSE CLOCK_MONOTONIC
#endif

NACS_EXPORT() uint64_t getCoarseTime()
{
    timespec t;
    clock_gettime(CLOCK_MONOTONIC_COARSE, &t);
    return uint64_t(t.tv_sec) * 1000 * 1000 * 1000 + t.tv_nsec;
}

NACS_EXPORT() uint64_t getCoarseRes()
{
    timespec t;
    clock_getres(CLOCK_MONOTONIC_COARSE, &t);
    return uint64_t(t.tv_sec) * 1000 * 1000 * 1000 + t.tv_nsec;
}

NACS_EXPORT() uint64_t getElapse(uint64_t prev)
{
    return getTime() - prev;
}

}
