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

#include "timer.h"
#include "number.h"
#include <time.h>
#include <pthread.h>

#include <vector>

#ifdef CLOCK_THREAD_CPUTIME_ID
#  define CLOCK_ID CLOCK_THREAD_CPUTIME_ID
#else
#  define CLOCK_ID CLOCK_MONOTONIC
#endif

NACS_EXPORT uint64_t
nacsGetTime()
{
    struct timespec time_spec;
    clock_gettime(CLOCK_ID, &time_spec);
    return (((uint64_t)time_spec.tv_sec) * 1000000000ull +
            (uint64_t)time_spec.tv_nsec);
}

NACS_EXPORT uint64_t
nacsGetElapse(uint64_t prev)
{
    return nacsGetTime() - prev;
}

static thread_local std::vector<uint64_t> nacs_tics;

NACS_EXPORT void
nacsTic()
{
    nacs_tics.push_back(0);
    auto &back = nacs_tics.back();
    back = nacsGetTime();
}

NACS_EXPORT uint64_t
nacsToc()
{
    uint64_t cur_time = nacsGetTime();
    if (!nacs_tics.size()) {
        return 0;
    }
    uint64_t old_time = nacs_tics.back();
    nacs_tics.pop_back();
    return cur_time - old_time;
}
