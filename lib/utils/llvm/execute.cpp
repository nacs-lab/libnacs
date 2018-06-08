/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "execute.h"
#include "../number.h"

#if NACS_OS_WINDOWS
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

namespace NaCs {
namespace LLVM {
namespace Exe {

#if NACS_OS_WINDOWS
static HMODULE get_self_hdl()
{
    HMODULE handle;
    if (!GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS,
                            (LPCWSTR)uintptr_t(get_self_hdl), &handle))
        return 0;
    return handle;
}
#endif

JITSymbol Resolver::findSymbolInLogicalDylib(const std::string&)
{
    return nullptr;
}

JITSymbol Resolver::findSymbol(const std::string &name)
{
    if (auto addr = find_extern(name))
        return {addr, JITSymbolFlags::Exported};
    return nullptr;
}

uintptr_t Resolver::find_extern(const std::string &name)
{
    if (name == "interp")
        return (uintptr_t)linearInterpolate;
#if NACS_OS_WINDOWS
    static auto self_hdl = get_self_hdl();
    if (uintptr_t ptr = (uintptr_t)GetProcAddress(self_hdl, name.c_str()))
        return ptr;
    static auto exe_hdl = GetModuleHandle(nullptr);
    if (uintptr_t ptr = (uintptr_t)GetProcAddress(exe_hdl, name.c_str()))
        return ptr;
    static auto ntdll_hdl = GetModuleHandle("ntdll.dll");
    if (uintptr_t ptr = (uintptr_t)GetProcAddress(ntdll_hdl, name.c_str()))
        return ptr;
    return 0;
#else
    static auto self_hdl = dlopen(nullptr, RTLD_LAZY);
    return (uintptr_t)dlsym(self_hdl, name.c_str());
#endif
}

NACS_EXPORT() Engine::Engine()
{
}

}
}
}
