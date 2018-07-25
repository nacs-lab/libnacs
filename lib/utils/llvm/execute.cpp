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
#include "../ir_p.h"

#if NACS_OS_WINDOWS
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

#include <llvm/Support/MemoryBuffer.h>

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
#ifndef NACS_HAS_EXP10
    if (name == "exp10")
        return (uintptr_t)nacs_exp10;
#endif
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
    reset_dyld();
}

NACS_EXPORT() Engine::~Engine()
{
}

void Engine::reset_dyld()
{
    m_dyld.reset(new RuntimeDyld(m_memmgr, m_resolver));
}

NACS_EXPORT() bool Engine::load(const object::ObjectFile &obj)
{
    if (!m_dyld->loadObject(obj))
        return false;
    m_dyld->finalizeWithMemoryManagerLocking();
    return true;
}

NACS_EXPORT() bool Engine::load(const char *p, size_t len)
{
    MemoryBufferRef buff(StringRef(p, len), "");
    auto obj = object::ObjectFile::createObjectFile(buff);
    if (!obj)
        return false;
    return load(*(*obj));
}

NACS_EXPORT() void *Engine::get_symbol(StringRef name)
{
    return (void*)(uintptr_t)m_dyld->getSymbol(name).getAddress();
}

}
}
}