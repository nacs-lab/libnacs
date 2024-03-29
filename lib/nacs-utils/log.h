/*************************************************************************
 *   Copyright (c) 2013 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include <functional>

#include <stdarg.h>

namespace NaCs {
namespace Log {

// Do not use this on windows since the compiler doesn't seem to like the
// version of `inttypes.h` it finds...
#if !NACS_OS_WINDOWS
#  define NACS_PRINTF_FORMA_ATTR(a, b) __attribute__((format(printf, a, b)))
#else
#  define NACS_PRINTF_FORMA_ATTR(a, b)
#endif

typedef enum {
    Debug,
    Info,
    Warn,
    Error,
    Force
} Level;

extern NACS_EXPORT(utils) Level level;

using cb_t = std::function<void(Level,const char*,const char*)>;

NACS_EXPORT(utils) void pushLogger(cb_t cb);
NACS_EXPORT(utils) void popLogger();
NACS_EXPORT(utils) bool printPID();
NACS_EXPORT(utils) void printPID(bool b);

NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(3, 4)
void _log(Level level, const char *func, const char *fmt, ...);

NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(3, 0)
void _logV(Level level, const char *func, const char *fmt, va_list ap);

NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 2) void info(const char *fmt, ...);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 2) void warn(const char *fmt, ...);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 2) void error(const char *fmt, ...);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 2) void log(const char *fmt, ...);

NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 0) void infoV(const char *fmt, va_list ap);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 0) void warnV(const char *fmt, va_list ap);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 0) void errorV(const char *fmt, va_list ap);
NACS_EXPORT(utils) NACS_PRINTF_FORMA_ATTR(1, 0) void logV(const char *fmt, va_list ap);

} // Log

NACS_EXPORT(utils) void backtrace();

} // NaCs

#define nacsDbg(fmt, args...)                                           \
    NaCs::Log::_log(NaCs::Log::Debug, __FUNCTION__, fmt, ##args)

#define nacsDbgV(fmt, ap)                                       \
    NaCs::Log::_logV(NaCs::Log::Debug, __FUNCTION__, fmt, ap)

#endif
