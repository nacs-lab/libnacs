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

#if NACS_OS_LINUX
#  include <sys/ioctl.h>
#  include <linux/perf_event.h>
#  include <linux/hw_breakpoint.h>
#  include <asm/unistd.h>
#endif

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

#if NACS_OS_LINUX
namespace {

static inline int perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                                   int cpu, int group_fd, unsigned long flags)
{
    return (int)syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
}

}

PerfCounter::PerfCounter(Type type)
{
    struct perf_event_attr pe = {};
    pe.size = sizeof(pe);
    pe.disabled = 1;
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;
    switch (type) {
    case CPUCycles:
        pe.type = PERF_TYPE_HARDWARE;
        pe.config = PERF_COUNT_HW_CPU_CYCLES;
        break;
    case CPUInsts:
        pe.type = PERF_TYPE_HARDWARE;
        pe.config = PERF_COUNT_HW_INSTRUCTIONS;
        break;
    default:
        return;
    }
    m_fd = perf_event_open(&pe, 0, -1, -1, 0);
}

void PerfCounter::start(bool reset)
{
    if (m_fd == -1)
        return;
    if (reset)
        ioctl(m_fd, PERF_EVENT_IOC_RESET, 0);
    ioctl(m_fd, PERF_EVENT_IOC_ENABLE, 0);
}

int64_t PerfCounter::finish(bool stop)
{
    if (m_fd == -1)
        return 0;
    if (stop)
        ioctl(m_fd, PERF_EVENT_IOC_DISABLE, 0);
    long long res;
    read(m_fd, &res, sizeof(res));
    return res;
}
#else
PerfCounter::PerfCounter(Type)
{
}

void PerfCounter::start(bool reset)
{
}

int64_t PerfCounter::finish(bool stop)
{
    return 0;
}
#endif

}
