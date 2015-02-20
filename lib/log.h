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

#ifndef _NACS_UTILS_LOG_H_
#define _NACS_UTILS_LOG_H_

#include "utils.h"
#include <stdarg.h>

typedef enum {
    NACS_LOG_DEBUG,
    NACS_LOG_INFO,
    NACS_LOG_WARN,
    NACS_LOG_ERROR,
    NACS_LOG_FORCE
} NaCsLogLevel;

extern NaCsLogLevel nacs_log_level;

static NACS_INLINE bool
nacsCheckLogLevel(unsigned level)
{
    return nacsUnlikely(level <= NACS_LOG_FORCE && level >= nacs_log_level);
}

void nacsSetLog(FILE *log_f);
FILE *nacsGetLog();

__attribute__((format(printf, 3, 4)))
void _nacsLog(NaCsLogLevel level, const char *func, const char *fmt, ...);

__attribute__((format(printf, 3, 0)))
void _nacsLogV(NaCsLogLevel level, const char *func,
               const char *fmt, va_list ap);

#define __nacsLog(__level, fmt, args...)                                \
    do {                                                                \
        unsigned level = (__level);                                     \
        if (!nacsCheckLogLevel(level)) {                                \
            break;                                                      \
        }                                                               \
        _nacsLog((NaCsLogLevel)level, __FUNCTION__, fmt, ##args);       \
    } while (0)

#define __nacsLogV(__level, fmt, ap)                            \
    do {                                                        \
        unsigned level = (__level);                             \
        if (!nacsCheckLogLevel(level)) {                        \
            break;                                              \
        }                                                       \
        _nacsLogV((NaCsLogLevel)level, __FUNCTION__, fmt, ap);  \
    } while (0)

#define nacsDebug(fmt, args...)                 \
    __nacsLog(NACS_LOG_DEBUG, fmt, ##args)
#define nacsInfo(fmt, args...)                  \
    __nacsLog(NACS_LOG_INFO, fmt, ##args)
#define nacsWarn(fmt, args...)                  \
    __nacsLog(NACS_LOG_WARN, fmt, ##args)
#define nacsError(fmt, args...)                 \
    __nacsLog(NACS_LOG_ERROR, fmt, ##args)
#define nacsLog(fmt, args...)                   \
    __nacsLog(NACS_LOG_FORCE, fmt, ##args)

#define nacsDebugV(fmt, args...)                \
    __nacsLogV(NACS_LOG_DEBUG, fmt, ##args)
#define nacsInfoV(fmt, args...)                 \
    __nacsLogV(NACS_LOG_INFO, fmt, ##args)
#define nacsWarnV(fmt, args...)                 \
    __nacsLogV(NACS_LOG_WARN, fmt, ##args)
#define nacsErrorV(fmt, args...)                \
    __nacsLogV(NACS_LOG_ERROR, fmt, ##args)
#define nacsLogV(fmt, args...)                  \
    __nacsLogV(NACS_LOG_FORCE, fmt, ##args)

void nacsBacktrace();

#endif
