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

#ifndef _NACS_UTILS_LOG_H_
#define _NACS_UTILS_LOG_H_

#include "utils.h"
#include <stdarg.h>

namespace NaCs {
namespace Log {

typedef enum {
    DEBUG,
    INFO,
    WARN,
    ERROR,
    FORCE
} Level;

extern Level level;

static NACS_INLINE bool checkLevel(unsigned _level)
{
    return NaCs::unlikely(_level <= FORCE && _level >= level);
}

void setLog(FILE *log_f);
FILE *getLog();

__attribute__((format(printf, 3, 4)))
void _log(Level level, const char *func, const char *fmt, ...);

__attribute__((format(printf, 3, 0)))
void _logV(Level level, const char *func, const char *fmt, va_list ap);

} // Log

void backtrace();

} // NaCs

#define __nacsLog(__level, fmt, args...)                        \
    do {                                                        \
        auto level = (NaCs::Log::Level)(__level);               \
        if (!NaCs::Log::checkLevel(level))                      \
            break;                                              \
        NaCs::Log::_log(level, __FUNCTION__, fmt, ##args);      \
    } while (0)

#define __nacsLogV(__level, fmt, ap)                    \
    do {                                                \
        auto level = (NaCs::Log::Level)(__level);       \
        if (!NaCs::Log::checkLevel(level))              \
            break;                                      \
        NaCs::Log::_logV(level, __FUNCTION__, fmt, ap); \
    } while (0)

#define nacsDebug(fmt, args...)                 \
    __nacsLog(NaCs::Log::DEBUG, fmt, ##args)
#define nacsInfo(fmt, args...)                  \
    __nacsLog(NaCs::Log::INFO, fmt, ##args)
#define nacsWarn(fmt, args...)                  \
    __nacsLog(NaCs::Log::WARN, fmt, ##args)
#define nacsError(fmt, args...)                 \
    __nacsLog(NaCs::Log::ERROR, fmt, ##args)
#define nacsLog(fmt, args...)                   \
    __nacsLog(NaCs::Log::FORCE, fmt, ##args)

#define nacsDebugV(fmt, args...)                \
    __nacsLogV(NaCs::Log::DEBUG, fmt, ##args)
#define nacsInfoV(fmt, args...)                 \
    __nacsLogV(NaCs::Log::INFO, fmt, ##args)
#define nacsWarnV(fmt, args...)                 \
    __nacsLogV(NaCs::Log::WARN, fmt, ##args)
#define nacsErrorV(fmt, args...)                \
    __nacsLogV(NaCs::Log::ERROR, fmt, ##args)
#define nacsLogV(fmt, args...)                  \
    __nacsLogV(NaCs::Log::FORCE, fmt, ##args)

#endif
