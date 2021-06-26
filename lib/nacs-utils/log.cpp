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
#include <string>
#include <vector>

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#if !NACS_OS_WINDOWS
#  include <execinfo.h>
#endif

namespace NaCs {
namespace Log {

static thread_local std::vector<cb_t> loggers;

NACS_EXPORT() void pushLogger(cb_t cb)
{
    loggers.push_back(cb);
}

NACS_EXPORT() void popLogger()
{
    if (loggers.empty())
        return;
    loggers.pop_back();
}

static cb_t get_logger()
{
    if (loggers.empty())
        return cb_t();
    return loggers.back();
}

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
    auto logger = get_logger();
    if (logger) {
#if NACS_OS_LINUX
        char *str = nullptr;
        if (vasprintf(&str, fmt, ap) == -1) {
            logger(Error, __func__, "Unable to allocate memory for message\n");
            return;
        }
#else
        va_list aq;
        va_copy(aq, ap);
        auto size = vsnprintf(nullptr, 0, fmt, aq);
        va_end(aq);
        // size doesn't include the NUL byte at the end.
        char *str = (char*)malloc(size + 1);
        auto size2 = vsnprintf(str, size + 1, fmt, ap);
        assert(size == size2);
        (void)size2;
#endif
        logger(level, func, str);
        free(str);
        return;
    }

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
                fprintf(stderr, "%d: ", pid);
            }
            else if (func) {
                fprintf(stderr, "%s-%d %s ", log_prefixes[(int)level], pid, func);
            }
            else {
                fprintf(stderr, "%s-%d ", log_prefixes[(int)level], pid);
            }
        }
        else if (level == Force) {
        }
        else if (func) {
            fprintf(stderr, "%s: %s ", log_prefixes[(int)level], func);
        }
        else {
            fprintf(stderr, "%s: ", log_prefixes[(int)level]);
        }
        vfprintf(stderr, fmt, ap);
    }
    fflush(stderr);
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
    auto logger = Log::get_logger();
    if (!logger) {
        backtrace_symbols_fd(buff, size, STDERR_FILENO);
        return;
    }
    auto strs = backtrace_symbols(buff, size);
    std::string linebuff; // Buffer to add EOL
    for (int i = 0; i < size; i++) {
        linebuff = strs[i];
        linebuff += '\n';
        logger(Log::Force, nullptr, linebuff.data());
    }
    free(strs);
#endif
}

} // Nacs
