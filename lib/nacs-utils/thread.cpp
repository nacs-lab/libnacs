/*************************************************************************
 *   Copyright (c) 2015 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "thread.h"

#include <pthread.h>
#include <sched.h>

namespace NaCs::Thread {

NACS_EXPORT() bool pin(int cpu)
{
#if NACS_OS_LINUX
    auto self = pthread_self();

    int ret;
    if (likely(cpu < CPU_SETSIZE)) {
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(cpu, &cpuset);
        ret = pthread_setaffinity_np(self, sizeof(cpu_set_t), &cpuset);
    }
    else {
        auto cpuset = CPU_ALLOC(cpu + 1);
        auto setsize = CPU_ALLOC_SIZE(cpu + 1);
        CPU_ZERO_S(setsize, cpuset);
        CPU_SET_S(cpu, setsize, cpuset);
        ret = pthread_setaffinity_np(self, setsize, cpuset);
        CPU_FREE(cpuset);
    }
    return ret == 0;
#else
    return false;
#endif
}
}
