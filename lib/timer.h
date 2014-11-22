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
#include <inttypes.h>

NACS_BEGIN_DECLS

#define PRTime PRIu64
uint64_t nacsGetTime();
uint64_t nacsGetElapse(uint64_t prev);
void nacsTic();
uint64_t nacsToc();

#define nacsPrintTime(time) nacsForceLog("Time: %" PRTime "\n", time)
#define nacsPrintElapse(prev) nacsPrintTime(nacsGetElapse(prev))
#define nacsPrintToc() nacsPrintTime(nacsToc())

#define TICKS_PER_SECOND (1000000000ll)
#define TICKS_PER_US (1000ll)
#define TICKS_PER_MS (1000000ll)

NACS_END_DECLS

#endif
