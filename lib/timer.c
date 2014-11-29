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

typedef struct {
    size_t alloc;
    size_t num;
    uint64_t tics[0];
} NaCsTics;

static NaCsTics*
nacsTicsResize(NaCsTics *tics, size_t new_size)
{
    new_size = nacsAlignTo(new_size, 16);
    if (!tics || new_size > tics->alloc || new_size * 2 < tics->alloc) {
        tics = (NaCsTics*)realloc(tics, (sizeof(NaCsTics) +
                                         sizeof(uint64_t) * new_size));
        tics->alloc = new_size;
    }
    return tics;
}

static pthread_key_t nacs_tics_key;
static pthread_once_t nacs_tics_key_once = PTHREAD_ONCE_INIT;

static void
nacsMakeTicsKey()
{
    pthread_key_create(&nacs_tics_key, free);
}

static NaCsTics*
nacsGetTics()
{
    pthread_once(&nacs_tics_key_once, nacsMakeTicsKey);
    return (NaCsTics*)pthread_getspecific(nacs_tics_key);
}

static void
nacsSetTics(NaCsTics *tics)
{
    pthread_once(&nacs_tics_key_once, nacsMakeTicsKey);
    pthread_setspecific(nacs_tics_key, tics);
}

NACS_EXPORT void
nacsTic()
{
    NaCsTics *old_tics = nacsGetTics();
    size_t num = old_tics ? old_tics->num : 0;
    NaCsTics *tics = nacsTicsResize(old_tics, num + 1);
    tics->num = num + 1;
    if (tics != old_tics) {
        nacsSetTics(tics);
    }
    tics->tics[num] = nacsGetTime();
}

NACS_EXPORT uint64_t
nacsToc()
{
    uint64_t cur_time = nacsGetTime();
    NaCsTics *old_tics = nacsGetTics();
    if (!old_tics || !old_tics->num) {
        return 0;
    }
    old_tics->num--;
    uint64_t old_time = old_tics->tics[old_tics->num];
    NaCsTics *tics = nacsTicsResize(old_tics, old_tics->num);
    if (tics != old_tics) {
        nacsSetTics(tics);
    }
    return cur_time - old_time;
}
