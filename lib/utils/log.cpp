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
        return INFO;
    if (strcasecmp(env, "debug") == 0)
        return DEBUG;
    else if (strcasecmp(env, "info") == 0)
        return INFO;
    else if (strcasecmp(env, "warn") == 0 || strcasecmp(env, "warning") == 0)
        return WARN;
    else if (strcasecmp(env, "error") == 0)
        return ERROR;
    else if (strcasecmp(env, "none") == 0)
        return FORCE;
    return INFO;
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

NACS_EXPORT() void _logV(Level level, const char *func, const char *fmt, va_list ap)
{
    NACS_RET_IF_FAIL(level >= level && ((int)level) >= 0 && level <= FORCE);
    static const char *log_prefixes[] = {
        [DEBUG] = "Debug-",
        [INFO] = "Info-",
        [WARN] = "Warn-",
        [ERROR] = "Error-",
        // [FORCE] = "Log-",
    };

    int pid = getpid();

    static std::mutex log_lock;
    {
        std::lock_guard<std::mutex> lk(log_lock);
        if (level == FORCE) {
            fprintf(log_f, "%d: ", pid);
        } else {
            fprintf(log_f, "%s%d %s ", log_prefixes[(int)level], pid, func);
        }
        vfprintf(log_f, fmt, ap);
    }
    fflush(log_f);
}

NACS_EXPORT() void _log(Level level, const char *func, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    _logV(level, func, fmt, ap);
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
