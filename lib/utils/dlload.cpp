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

#include <string>

namespace NaCs {

namespace DL {

#if NACS_OS_WINDOWS
static HMODULE get_libnacs_handle()
{
    static HMODULE hdl =
        [] {
            HMODULE hdl;
            GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                               GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                               (LPCWSTR)(uintptr_t)(&open),
                               (HMODULE*)&hdl);
            return hdl;
        }();
    return hdl;
}

static inline bool is_dir_separator(wchar_t c)
{
    return c == L'/' || c == L'\\';
}

static const std::wstring &get_libnacs_dir()
{
    static const std::wstring dir =
        [] {
            auto hdl = get_libnacs_handle();
            DWORD len = GetModuleFileNameW(hdl, nullptr, 0);
            if (!len)
                return std::wstring();
            std::wstring res(len, 0);
            GetModuleFileNameW(hdl, &res[0], len);
            while (len-- > 0) {
                if (is_dir_separator(res[len])) {
                    res.resize(len);
                    return res;
                }
            }
            // Shouldn't happen...
            res = L"C:\\";
            return res;
        }();
    return dir;
}

// This implementation is copied from dotnet/coreclr (MIT License)
static bool is_full_path(wchar_t *path, DWORD len)
{
    // There is no way to specify a fixed path with one character (or less).
    if (len < 2)
        return false;
    // There is no valid way to specify a relative path with two initial slashes or
    // "\?" as "?" isn't valid for drive relative paths and "\??\" is equivalent to "\\?\"
    if (is_dir_separator(path[0]))
        return path[1] == L'?' || is_dir_separator(path[1]);
    return len >= 3 && (path[1] == L':') && is_dir_separator(path[2]);
}

NACS_EXPORT() void *open(const char *filename, int)
{
    DWORD len = MultiByteToWideChar(CP_UTF8, 0, filename, -1, nullptr, 0);
    if (!len)
        return nullptr;
    WCHAR *wfilename = (WCHAR*)alloca(len * sizeof(WCHAR));
    if (!MultiByteToWideChar(CP_UTF8, 0, filename, -1, wfilename, len))
        return nullptr;
    if (is_full_path(wfilename, len - 1))
        return LoadLibraryExW(wfilename, nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
    if (auto hdl = LoadLibraryExW(wfilename, nullptr, 0))
        return hdl;
    auto path = get_libnacs_dir() + wfilename;
    return LoadLibraryExW(path.data(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
}

NACS_EXPORT() bool close(void *handle)
{
    if (!handle)
        return false;
    return FreeLibrary((HMODULE)handle);
}

static void *_sym(HMODULE handle, const char *symbol)
{
    return (void*)GetProcAddress((HMODULE)handle, symbol);
}

NACS_EXPORT() void *sym(void *handle, const char *symbol)
{
    if (handle)
        return _sym((HMODULE)handle, symbol);
    if (auto res = _sym(get_libnacs_handle(), symbol))
        return res;
    static auto exe_hdl = GetModuleHandle(nullptr);
    if (auto res = _sym(exe_hdl, symbol))
        return res;
    static auto ntdll_hdl = GetModuleHandle("ntdll.dll");
    if (auto res = _sym(ntdll_hdl, symbol))
        return res;
    return nullptr;
}
#else
static void* get_libnacs_handle()
{
    static void* hdl =
        [] {
            Dl_info info;
            if (dladdr((void*)(uintptr_t)&open, &info) && info.dli_fname) {
                if (auto hdl = dlopen(info.dli_fname, RTLD_NOW)) {
                    return hdl;
                }
            }
            return dlopen(nullptr, RTLD_NOW);
        }();
    return hdl;
}

#if NACS_OS_DARWIN
static char const *const extensions[] = {".dylib"};
#else
static char const *const extensions[] = {".so"};
#endif
static constexpr int n_extensions = sizeof(extensions) / sizeof(char*);

static bool endswith_extension(const char *path)
{
    if (!path)
        return false;
    size_t len = strlen(path);
    for (int i = 0; i < n_extensions; i++) {
        const char *ext = extensions[i];
        size_t extlen = strlen(ext);
        // Skip version extensions if present
        size_t j = len - 1;
        while (j > 0) {
            if (path[j] != '.' && (path[j] < '0' || path[j] > '9'))
                break;
            j--;
        }
        if (j < extlen)
            continue;
        if ((j == len - 1 || path[j + 1] == '.') &&
            memcmp(ext, path + j - extlen + 1, extlen) == 0) {
            return true;
        }
    }
    return false;
}

static void *_open(const char *filename, int flags)
{
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
}

NACS_EXPORT() void *open(const char *filename, int flags)
{
    if (auto hdl = _open(filename, flags))
        return hdl;
    if (!endswith_extension(filename)) {
        std::string name_str(filename);
        for (int i = 0; i < n_extensions; i++) {
            if (auto hdl = _open((name_str + extensions[i]).c_str(), flags)) {
                return hdl;
            }
        }
    }
    return nullptr;
}

NACS_EXPORT() bool close(void *handle)
{
    dlerror(); /* Reset error status. */
    if (!handle)
        return false;
    return dlclose(handle) == 0;
}

NACS_EXPORT() void *sym(void *handle, const char *symbol)
{
    if (!handle)
        handle = get_libnacs_handle();
    dlerror(); /* Reset error status. */
    return dlsym(handle, symbol);
}
#endif

}

}
