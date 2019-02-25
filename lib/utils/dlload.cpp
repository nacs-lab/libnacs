/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

// The implementation is mostely borrowed from julia.

#include "dlload.h"
#if NACS_OS_WINDOWS
#  include <windows.h>
#  include <direct.h>
#else
#  include <unistd.h>
#  include <dlfcn.h>
#endif

namespace NaCs {

namespace DL {

NACS_EXPORT() void *open(const char *filename, Flags flags)
{
#if NACS_OS_WINDOWS
    size_t len = MultiByteToWideChar(CP_UTF8, 0, filename, -1, nullptr, 0);
    if (!len)
        return nullptr;
    WCHAR *wfilename = (WCHAR*)alloca(len * sizeof(WCHAR));
    if (!MultiByteToWideChar(CP_UTF8, 0, filename, -1, wfilename, len))
        return nullptr;
    return LoadLibraryExW(wfilename, nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
#else
#  define CONVERT_FLAG(flags, FLAG) (flags & FLAG ? RTLD_ ## FLAG : 0)
    dlerror(); /* Reset error status. */
    return dlopen(filename, CONVERT_FLAG(flags, LOCAL) | CONVERT_FLAG(flags, GLOBAL) |
#  ifdef RTLD_NODELETE
                  CONVERT_FLAG(flags, NODELETE) |
#  endif
#  ifdef RTLD_NOLOAD
                  CONVERT_FLAG(flags, NOLOAD) |
#  endif
#  ifdef RTLD_DEEPBIND
                  CONVERT_FLAG(flags, DEEPBIND) |
#  endif
#  ifdef RTLD_FIRST
                  CONVERT_FLAG(flags, FIRST) |
#  endif
                  (flags & NOW ? RTLD_NOW : RTLD_LAZY));
#endif
}

NACS_EXPORT() bool close(void *handle)
{
#if NACS_OS_WINDOWS
    if (!handle)
        return false;
    return FreeLibrary((HMODULE)handle);
#else
    dlerror(); /* Reset error status. */
    if (!handle)
        return false;
    return dlclose(handle) == 0;
#endif
}

NACS_EXPORT() void *sym(void *handle, const char *symbol)
{
#if NACS_OS_WINDOWS
    return GetProcAddress((HMODULE)handle, symbol);
#else
    dlerror(); /* Reset error status. */
    return dlsym(handle, symbol);
#endif
}

}

}
