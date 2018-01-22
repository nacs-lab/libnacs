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

#include <chrono>
#include <vector>

namespace NaCs {

NACS_EXPORT(utils) uint64_t
getTime()
{
    using namespace std::chrono;
    return time_point_cast<nanoseconds>(high_resolution_clock::now())
        .time_since_epoch().count();
}

NACS_EXPORT(utils) uint64_t
getElapse(uint64_t prev)
{
    return getTime() - prev;
}

static thread_local std::vector<uint64_t> tics_list;

NACS_EXPORT(utils) void
tic()
{
    tics_list.push_back(0);
    auto &back = tics_list.back();
    back = getTime();
}

NACS_EXPORT(utils) uint64_t
toc()
{
    uint64_t cur_time = getTime();
    if (!tics_list.size()) {
        return 0;
    }
    uint64_t old_time = tics_list.back();
    tics_list.pop_back();
    return cur_time - old_time;
}

}
