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

#include "log.h"
#include "number.h"

#include <mutex>

#include <stdarg.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#if !NACS_OS_WINDOWS
#  include <execinfo.h>
#endif

namespace NaCs {
namespace Log {

NACS_EXPORT() Level level = [] {
    auto env = getenv("NACS_LOG");
    if (!env)
        return Info;
    if (strcasecmp(env, "debug") == 0)
        return Debug;
    else if (strcasecmp(env, "info") == 0)
        return Info;
    else if (strcasecmp(env, "warn") == 0 || strcasecmp(env, "warning") == 0)
        return Warn;
    else if (strcasecmp(env, "error") == 0)
        return Error;
    else if (strcasecmp(env, "none") == 0)
        return Force;
    return Info;
}();

static FILE *log_f = stderr;

NACS_EXPORT() FILE *getLog()
{
    return log_f;
}

NACS_EXPORT() void setLog(FILE *f)
{
    log_f = f ? f : stderr;
}

static bool print_pid = true;
NACS_EXPORT() bool printPID()
{
    return print_pid;
}

NACS_EXPORT() void printPID(bool b)
{
    print_pid = b;
}

static NACS_INLINE bool checkLevel(unsigned _level)
{
    return _level <= Force && _level >= level;
}

NACS_EXPORT() void _logV(Level level, const char *func, const char *fmt, va_list ap)
{
    NACS_RET_IF_FAIL(checkLevel(level));
    static const char *log_prefixes[] = {
        [Debug] = "Debug",
        [Info] = "Info",
        [Warn] = "Warn",
        [Error] = "Error",
        // [Force] = "Log",
    };

    static std::mutex log_lock;
    {
        std::lock_guard<std::mutex> lk(log_lock);
        if (print_pid) {
            int pid = getpid();
            if (level == Force) {
                fprintf(log_f, "%d: ", pid);
            }
            else if (func) {
                fprintf(log_f, "%s-%d %s ", log_prefixes[(int)level], pid, func);
            }
            else {
                fprintf(log_f, "%s-%d ", log_prefixes[(int)level], pid);
            }
        }
        else if (level == Force) {
        }
        else if (func) {
            fprintf(log_f, "%s: %s ", log_prefixes[(int)level], func);
        }
        else {
            fprintf(log_f, "%s: ", log_prefixes[(int)level]);
        }
        vfprintf(log_f, fmt, ap);
    }
    fflush(log_f);
}

NACS_EXPORT() void _log(Level level, const char *func, const char *fmt, ...)
{
    NACS_RET_IF_FAIL(checkLevel(level));
    va_list ap;
    va_start(ap, fmt);
    _logV(level, func, fmt, ap);
    va_end(ap);
}

NACS_EXPORT() void infoV(const char *fmt, va_list ap)
{
    _logV(Info, nullptr, fmt, ap);
}
NACS_EXPORT() void warnV(const char *fmt, va_list ap)
{
    _logV(Warn, nullptr, fmt, ap);
}
NACS_EXPORT() void errorV(const char *fmt, va_list ap)
{
    _logV(Error, nullptr, fmt, ap);
}
NACS_EXPORT() void logV(const char *fmt, va_list ap)
{
    _logV(Force, nullptr, fmt, ap);
}

NACS_EXPORT() void info(const char *fmt, ...)
{
    NACS_RET_IF_FAIL(checkLevel(Info));
    va_list ap;
    va_start(ap, fmt);
    infoV(fmt, ap);
    va_end(ap);
}

NACS_EXPORT() void warn(const char *fmt, ...)
{
    NACS_RET_IF_FAIL(checkLevel(Warn));
    va_list ap;
    va_start(ap, fmt);
    warnV(fmt, ap);
    va_end(ap);
}

NACS_EXPORT() void error(const char *fmt, ...)
{
    NACS_RET_IF_FAIL(checkLevel(Error));
    va_list ap;
    va_start(ap, fmt);
    errorV(fmt, ap);
    va_end(ap);
}

NACS_EXPORT() void log(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    logV(fmt, ap);
    va_end(ap);
}

} // Log

NACS_EXPORT() void backtrace()
{
#if !NACS_OS_WINDOWS
    void *buff[1024];
    int size = ::backtrace(buff, 1024);
    int fd = fileno(Log::log_f);
    if (fd == -1)
        fd = STDERR_FILENO;
    backtrace_symbols_fd(buff, size, fd);
#endif
}

} // Nacs
