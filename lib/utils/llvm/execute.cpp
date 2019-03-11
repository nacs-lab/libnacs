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
#include "../processor.h"

#include <llvm/Support/MemoryBuffer.h>

namespace NaCs {
namespace LLVM {
namespace Exe {

// Clang is not happy if we redeclare the function even if it's in a function local scope
// since it moves all the declarations to the global scope...
// Give them different names to work around that.
// Thanks to http://zwizwa.be/-/c/20100825-142132
#define __asm_sym(var_suffix, name)                     \
    ([] {                                               \
         extern void asm_sym ## var_suffix() asm(name); \
         return asm_sym ## var_suffix;                  \
     }())
// Indirection to make sure `__COUNTER__` is expanded correctly.
#define _asm_sym(var_suffix, name) __asm_sym(var_suffix, name)
#define asm_sym(name) _asm_sym(__COUNTER__, name)

#ifdef ENABLE_SIMD
#  if NACS_CPU_X86 || NACS_CPU_X86_64
#    if NACS_OS_WINDOWS && NACS_CPU_X86_64
#      define VEC_SUFFIX(x) "@@" #x
#    elif NACS_OS_WINDOWS
#      error "SIMD not supported on Win32"
#    else
#      define VEC_SUFFIX(x)
#    endif
#  elif NACS_CPU_AARCH64
#  else
#    error "SIMD not supported"
#  endif
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
        return (uintptr_t)static_cast<double(*)(double, uint32_t, const double*)>(
            linearInterpolate);
    if (auto openlibm_hdl = IR::get_openlibm_handle())
        if (auto addr = DL::sym(openlibm_hdl, name.c_str()))
            return (uintptr_t)addr;
#ifndef NACS_HAS_EXP10
    if (name == "exp10")
        return (uintptr_t)nacs_exp10;
#endif
#ifdef ENABLE_SIMD
#  if NACS_CPU_X86 || NACS_CPU_X86_64
    // We need this to convience gcc to find the correct symbol on windows.
    static const auto &host_info = CPUInfo::get_host();
    if (name == "interp.2") {
        auto addr = asm_sym("_ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd" VEC_SUFFIX(32));
        if (host_info.test_feature(X86::Feature::avx2))
            addr = asm_sym("_ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd" VEC_SUFFIX(32));
        return (uintptr_t)addr;
    }
    else if (name == "interp.4") {
        auto addr = asm_sym("_ZN4NaCs22linearInterpolate4_avxEDv4_djPKd" VEC_SUFFIX(48));
        if (host_info.test_feature(X86::Feature::avx2))
            addr = asm_sym("_ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd" VEC_SUFFIX(48));
        return (uintptr_t)addr;
    }
#    if !NACS_OS_WINDOWS
    else if (name == "interp.8") {
        __m512d(*addr)(__m512d, uint32_t, const double*);
        addr = linearInterpolate8_avx512f;
        return (uintptr_t)addr;
    }
#    endif
#  elif NACS_CPU_AARCH64
    if (name == "interp.2") {
        return (uintptr_t)static_cast<float64x2_t(*)(float64x2_t, uint32_t, const double*)>(
            linearInterpolate2_asimd);
    }
    else if (name == "interp.2x2") {
        return (uintptr_t)static_cast<float64x2x2_t(*)(float64x2x2_t, uint32_t, const double*)>(
            linearInterpolate4_asimd);
    }
    else if (name == "interp.2x4") {
        return (uintptr_t)static_cast<float64x2x4_t(*)(float64x2x4_t, uint32_t, const double*)>(
            linearInterpolate8_asimd);
    }
#  endif
#endif
    return (uintptr_t)DL::sym(nullptr, name.c_str());
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
