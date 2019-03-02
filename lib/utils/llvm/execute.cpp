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
#include "../dlload.h"

#if NACS_OS_WINDOWS
#  include <windows.h>
#endif

#include <llvm/Support/MemoryBuffer.h>

namespace NaCs {
namespace LLVM {
namespace Exe {

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
        return (uintptr_t)static_cast<double(*)(double, uint32_t, const double*)>(
            linearInterpolate);
    if (auto addr = DL::sym(IR::get_openlibm_handle(), name.c_str()))
        return (uintptr_t)addr;
#ifndef NACS_HAS_EXP10
    if (name == "exp10")
        return (uintptr_t)nacs_exp10;
#endif
    if (auto addr = DL::sym(nullptr, name.c_str()))
        return (uintptr_t)addr;
#if NACS_OS_WINDOWS
    static auto exe_hdl = GetModuleHandle(nullptr);
    if (uintptr_t ptr = (uintptr_t)GetProcAddress(exe_hdl, name.c_str()))
        return ptr;
    static auto ntdll_hdl = GetModuleHandle("ntdll.dll");
    if (uintptr_t ptr = (uintptr_t)GetProcAddress(ntdll_hdl, name.c_str()))
        return ptr;
#endif
    return 0;
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

NACS_EXPORT() uint64_t Engine::load(const object::ObjectFile &obj)
{
    if (!m_dyld)
        reset_dyld();
    auto id = m_memmgr.new_group();
    if (!m_dyld->loadObject(obj)) {
        m_memmgr.free_group(id);
        return 0;
    }
    m_dyld->finalizeWithMemoryManagerLocking();
    return id;
}

NACS_EXPORT() uint64_t Engine::load(const char *p, size_t len)
{
    MemoryBufferRef buff(StringRef(p, len), "");
    auto obj = object::ObjectFile::createObjectFile(buff);
    if (!obj)
        return 0;
    return load(*(*obj));
}

NACS_EXPORT() void *Engine::get_symbol(StringRef name)
{
    return (void*)(uintptr_t)m_dyld->getSymbol(name).getAddress();
}

NACS_EXPORT() void Engine::free(uint64_t id)
{
    m_memmgr.free_group(id);
    // Reusing the memory could really confuse the relocator so we need to reset the dyld
    // when we free any memory.
    m_dyld.reset(nullptr);
}

}
}
}
