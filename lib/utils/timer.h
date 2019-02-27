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

#ifndef _NACS_UTILS_TIMER_H_
#define _NACS_UTILS_TIMER_H_

#include "log.h"
#include "macros.h"

#include <inttypes.h>

#ifndef PRIu64
#  if NACS_OS_WINDOWS
#    define PRIu64 "I64u"
#  else
#    define PRIu64 "llu"
#  endif
#endif
#define PRTime PRIu64

#define TICKS_PER_SECOND (1000000000ll)
#define TICKS_PER_US (1000ll)
#define TICKS_PER_MS (1000000ll)

namespace NaCs {

NACS_EXPORT(utils) uint64_t getTime();
NACS_EXPORT(utils) uint64_t getRes();
NACS_EXPORT(utils) uint64_t getCoarseTime();
NACS_EXPORT(utils) uint64_t getCoarseRes();
NACS_EXPORT(utils) uint64_t getElapse(uint64_t prev);

static inline void
printTime(uint64_t time)
{
    Log::log("Time: %" PRIu64 "\n", time);
}

struct Timer {
    uint64_t elapsed(bool restart=false)
    {
        auto new_time = getTime();
        auto res = new_time - m_start_time;
        if (restart)
            m_start_time = new_time;
        return res;
    }
    void print(bool restart=false)
    {
        printTime(elapsed(restart));
    }
    void restart()
    {
        elapsed(true);
    }
private:
    uint64_t m_start_time = getTime();
};

}

#endif
