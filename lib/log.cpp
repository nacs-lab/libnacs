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
#include <unistd.h>
#include <stdarg.h>
#include <mutex>

// TODO lock
NACS_EXPORT NaCsLogLevel nacs_log_level = NACS_LOG_ERROR;
static FILE *logf = stderr;

NACS_EXPORT void
setLog(FILE *f)
{
    logf = f ? f : stderr;
}

NACS_EXPORT void
_nacsLogV(NaCsLogLevel level, const char *fname, int line, const char *func,
          const char *fmt, va_list ap)
{
    NACS_RET_IF_FAIL(level >= nacs_log_level && ((int)level) >= 0 &&
                     level <= NACS_LOG_FORCE);
    static const char *log_prefixes[] = {
        [NACS_LOG_DEBUG] = "nacsDebug-",
        [NACS_LOG_INFO] = "nacsInfo-",
        [NACS_LOG_WARN] = "nacsWarn-",
        [NACS_LOG_ERROR] = "nacsError-",
        [NACS_LOG_FORCE] = "nacsLog-",
    };

    static std::mutex log_lock;

    std::lock_guard<std::mutex> lk(log_lock);

    fprintf(logf, "%s%d (%s:%d) %s ", log_prefixes[(int)level], getpid(),
            fname, line, func);
    vfprintf(logf, fmt, ap);
}

NACS_EXPORT void
_nacsLog(NaCsLogLevel level, const char *fname, int line, const char *func,
         const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    _nacsLogV(level, fname, line, func, fmt, ap);
    va_end(ap);
}

/* NACS_EXPORT void */
/* nacsBacktrace() */
/* { */
/* #ifdef NACS_ENABLE_BACKTRACE */
/*     void *buff[1024]; */
/*     size_t size = backtrace(buff, 1024); */
/*     backtrace_symbols_fd(buff, size, STDERR_FILENO); */
/* #endif */
/* } */
