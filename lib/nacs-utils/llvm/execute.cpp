/*************************************************************************
 *   Copyright (c) 2018 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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

#include "config.h"

#include <llvm/Support/MemoryBuffer.h>

#if NACS_OS_WINDOWS
#  include "nacs_lib_names.h"
#endif

namespace NaCs::LLVM::Exe {

#define __asm_sym_real(var_suffix, name, suffix, type)                  \
    ([] {                                                               \
        extern auto asm_sym ## var_suffix type asm(name) suffix;        \
        return (void(*)())asm_sym ## var_suffix;                        \
    }())
#define asm_sym_real(var_suffix, name, suffix, type...)                 \
    NACS_SWITCH(type, __asm_sym_real(var_suffix, name, suffix, type),   \
                __asm_sym_real(var_suffix, name, suffix, () -> void))

// Clang is not happy if we redeclare the function even if it's in a function local scope
// since it moves all the declarations to the global scope...
// Give them different names to work around that.
// Thanks to http://zwizwa.be/-/c/20100825-142132
#define __asm_sym(var_suffix, name, type...)    \
    asm_sym_real(var_suffix, name, , ##type)
// Indirection to make sure `__COUNTER__` is expanded correctly.
#define _asm_sym(var_suffix, name, type...) __asm_sym(var_suffix, name, ##type)
#define asm_sym(name, type...) _asm_sym(__COUNTER__, name, ##type)
#define __asm_sym_w(var_suffix, name, type...)                          \
    asm_sym_real(var_suffix, name, __attribute__((weak)), ##type)
// Indirection to make sure `__COUNTER__` is expanded correctly.
#define _asm_sym_w(var_suffix, name, type...) __asm_sym_w(var_suffix, name, ##type)
#define asm_sym_w(name, type...) _asm_sym_w(__COUNTER__, name, ##type)

#ifdef ENABLE_SIMD
#  if NACS_CPU_X86 || NACS_CPU_X86_64
#    if NACS_OS_WINDOWS && NACS_CPU_X86_64
#      define VEC_SUFFIX(x) "@@" #x
#    elif NACS_OS_WINDOWS
#      error "SIMD not supported on Win32"
#    else
#      define VEC_SUFFIX(x)
#    endif
#    define _check_sleef(host_info, sleef_sym_f, var, sym, s1, s2, s3) do { \
        if (var == #sym ".2") {                                         \
            if (host_info.test_feature(X86::Feature::avx2))             \
                return sleef_sym_f(#sym "d2", "avx2128" VEC_SUFFIX(s1)); \
            if (host_info.test_feature(X86::Feature::sse41))            \
                return sleef_sym_f(#sym "d2", "sse4" VEC_SUFFIX(s1));   \
            return sleef_sym_f(#sym "d2", "sse2" VEC_SUFFIX(s1));       \
        }                                                               \
        if (var == #sym ".4") {                                         \
            if (host_info.test_feature(X86::Feature::avx2) &&           \
                host_info.test_feature(X86::Feature::fma))              \
                return sleef_sym_f(#sym "d4", "avx2" VEC_SUFFIX(s2));   \
            if (host_info.test_feature(X86::Feature::fma4))             \
                return sleef_sym_f(#sym "d4", "fma4" VEC_SUFFIX(s2));   \
            return sleef_sym_f(#sym "d4", "avx" VEC_SUFFIX(s2));        \
        }                                                               \
        if (var == #sym ".8") {                                         \
            return sleef_sym_f(#sym "d8", "avx512f" VEC_SUFFIX(s3));    \
        }                                                               \
    } while (0)
#    define _check_sleef_d(host_info, sleef_sym_f, var, sym)    \
    _check_sleef(host_info, sleef_sym_f, var, sym, 16, 32, 64)
#    define _check_sleef_dd(host_info, sleef_sym_f, var, sym)   \
    _check_sleef(host_info, sleef_sym_f, var, sym, 32, 64, 128)
#    define _check_sleef_ddd(host_info, sleef_sym_f, var, sym)  \
    _check_sleef(host_info, sleef_sym_f, var, sym, 48, 96, 192)
#    define _check_sleef_di(host_info, sleef_sym_f, var, sym)   \
    _check_sleef(host_info, sleef_sym_f, var, sym, 32, 48, 96)
#  elif NACS_CPU_AARCH64
#    define _check_sleef_d(host_info, sleef_sym_f, var, sym) do {       \
        if (var == #sym ".2") {                                         \
            return sleef_sym_f(#sym "d2", "advsimd");                   \
        }                                                               \
    } while (0)
#    define _check_sleef_dd _check_sleef_d
#    define _check_sleef_ddd _check_sleef_d
#    define _check_sleef_di _check_sleef_d
#  else
#    error "SIMD not supported"
#  endif
#  if NACS_OS_WINDOWS
// Weak symbol doesn't seem to work on windows.
// For release build, non-existing symbol returns a non-null address of some wrapper.
static void *libsleef_handle(void)
{
    static void *hdl = DL::open(NACS_SLEEF_NAME, DL::GLOBAL | DL::LAZY);
    return hdl;
}
#    define sleef_asm_sym_w(name)                               \
    ([] {                                                       \
        static auto sym = DL::sym(libsleef_handle(), name);     \
        return sym;                                             \
    } ())
#  else
#    define sleef_asm_sym_w(name) asm_sym_w(name)
#  endif
#  define sleef_sym_u(prefix, suffix)                                   \
    ([] {                                                               \
        if (auto addr = sleef_asm_sym_w("Sleef_" prefix "_u35" suffix)) \
            return (uintptr_t)addr;                                     \
        return (uintptr_t)asm_sym("Sleef_" prefix "_u10" suffix);       \
    } ())
#  define sleef_sym_u15(prefix, suffix)                                 \
    ([] {                                                               \
        if (auto addr = sleef_asm_sym_w("Sleef_" prefix "_u35" suffix)) \
            return (uintptr_t)addr;                                     \
        return (uintptr_t)asm_sym("Sleef_" prefix "_u15" suffix);       \
    } ())
#  define sleef_sym_u05(prefix, suffix)                                 \
    ([] {                                                               \
        if (auto addr = sleef_asm_sym_w("Sleef_" prefix "_u35" suffix)) \
            return (uintptr_t)addr;                                     \
        return (uintptr_t)asm_sym("Sleef_" prefix "_u05" suffix);       \
    } ())
#  define sleef_sym(prefix, suffix)                                     \
    ([] {                                                               \
        if (auto addr = sleef_asm_sym_w("Sleef_" prefix "_u35" suffix)) \
            return (uintptr_t)addr;                                     \
        if (auto addr = sleef_asm_sym_w("Sleef_" prefix "_u05" suffix)) \
            return (uintptr_t)addr;                                     \
        return (uintptr_t)asm_sym("Sleef_" prefix "_" suffix);          \
    } ())
#  define check_sleef_u_d(host_info, var, sym)          \
    _check_sleef_d(host_info, sleef_sym_u, var, sym)
#  define check_sleef_u15_d(host_info, var, sym)        \
    _check_sleef_d(host_info, sleef_sym_u15, var, sym)
#  define check_sleef_d(host_info, var, sym)            \
    _check_sleef_d(host_info, sleef_sym, var, sym)
#  define check_sleef_u_dd(host_info, var, sym)         \
    _check_sleef_dd(host_info, sleef_sym_u, var, sym)
#  define check_sleef_u05_dd(host_info, var, sym)       \
    _check_sleef_dd(host_info, sleef_sym_u05, var, sym)
#  define check_sleef_dd(host_info, var, sym)           \
    _check_sleef_dd(host_info, sleef_sym, var, sym)
#  define check_sleef_ddd(host_info, var, sym)          \
    _check_sleef_ddd(host_info, sleef_sym, var, sym)
#  define check_sleef_di(host_info, var, sym)           \
    _check_sleef_di(host_info, sleef_sym, var, sym)
#endif

class Resolver::SetCB {
public:
    SetCB(Resolver &resolver, const cb_t *cb)
        : m_resolver(resolver)
    {
        m_resolver.m_cb = cb;
    }

    ~SetCB()
    {
        m_resolver.m_cb = nullptr;
    }
private:
    Resolver &m_resolver;
};

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

#if NACS_OS_WINDOWS && NACS_CPU_X86_64
static size_t find_vector_suffix(const std::string &name)
{
    size_t len = name.size();
    for (size_t i = len; i >= 2; i--) {
        char c = name[i - 1];
        if (c >= '0' && c <= '9')
            continue;
        if (c != '@' || name[i - 2] != '@')
            return 0;
        return i - 2;
    }
    return 0;
}
#endif

NACS_EXPORT() uintptr_t Resolver::resolve_ir_sym(const std::string &orig_name)
{
#if NACS_OS_WINDOWS && NACS_CPU_X86_64
    std::string devec_name;
    const std::string *pname = &orig_name;
    if (auto prefix_len = find_vector_suffix(orig_name)) {
        devec_name = orig_name.substr(0, prefix_len);
        pname = &devec_name;
    }
    const std::string &name = *pname;
#else
    const std::string &name = orig_name;
#endif
    if (name == "interp")
        return (uintptr_t)static_cast<double(*)(double, uint32_t, const double*)>(
            linearInterpolate);
    if (auto openlibm_hdl = IR::get_openlibm_handle())
        if (auto addr = DL::sym(openlibm_hdl, orig_name.c_str()))
            return (uintptr_t)addr;
#ifndef NACS_HAS_EXP10
    if (name == "exp10")
        return (uintptr_t)nacs_exp10;
#endif
#ifdef ENABLE_SIMD
    static const auto &host_info = CPUInfo::get_host();
#  if NACS_CPU_X86 || NACS_CPU_X86_64
    // We need this to convience gcc to find the correct symbol on windows.
    if (name == "interp.2") {
        auto addr = asm_sym("_ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd" VEC_SUFFIX(32),
                            (__m128d, uint32_t, const double*) -> __m128d);
        if (host_info.test_feature(X86::Feature::avx2))
            addr = asm_sym("_ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd" VEC_SUFFIX(32),
                           (__m128d, uint32_t, const double*) -> __m128d);
        return (uintptr_t)addr;
    }
    else if (name == "interp.4") {
        auto addr = asm_sym("_ZN4NaCs22linearInterpolate4_avxEDv4_djPKd" VEC_SUFFIX(48),
                            (__m256d, uint32_t, const double*) -> __m256d);
        if (host_info.test_feature(X86::Feature::avx2))
            addr = asm_sym("_ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd" VEC_SUFFIX(48),
                           (__m256d, uint32_t, const double*) -> __m256d);
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
    (void)host_info;
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
    check_sleef_u_d(host_info, name, acos);
    check_sleef_u_d(host_info, name, acosh);
    check_sleef_u_d(host_info, name, asin);
    check_sleef_u_d(host_info, name, asinh);
    check_sleef_u_d(host_info, name, atan);
    check_sleef_u_d(host_info, name, atanh);
    check_sleef_u_d(host_info, name, cbrt);
    check_sleef_d(host_info, name, ceil);
    check_sleef_u_d(host_info, name, cos);
    check_sleef_u_d(host_info, name, cosh);
    check_sleef_u_d(host_info, name, erf);
    check_sleef_u15_d(host_info, name, erfc);
    check_sleef_u_d(host_info, name, exp);
    check_sleef_u_d(host_info, name, exp10);
    check_sleef_u_d(host_info, name, exp2);
    check_sleef_u_d(host_info, name, expm1);
    check_sleef_d(host_info, name, fabs);
    check_sleef_d(host_info, name, floor);
    check_sleef_u_d(host_info, name, tgamma);
    check_sleef_u_d(host_info, name, lgamma);
    check_sleef_u_d(host_info, name, log);
    check_sleef_u_d(host_info, name, log10);
    check_sleef_u_d(host_info, name, log1p);
    check_sleef_u_d(host_info, name, log2);
    check_sleef_d(host_info, name, rint);
    check_sleef_d(host_info, name, round);
    check_sleef_u_d(host_info, name, sin);
    check_sleef_u_d(host_info, name, sinh);
    check_sleef_d(host_info, name, sqrt);
    check_sleef_u_d(host_info, name, tan);
    check_sleef_u_d(host_info, name, tanh);

    check_sleef_u_dd(host_info, name, atan2);
    check_sleef_dd(host_info, name, copysign);
    check_sleef_dd(host_info, name, fdim);
    check_sleef_dd(host_info, name, fmax);
    check_sleef_dd(host_info, name, fmin);
    check_sleef_dd(host_info, name, fmod);
    check_sleef_u05_dd(host_info, name, hypot);
    check_sleef_u_dd(host_info, name, pow);

    check_sleef_ddd(host_info, name, fma);
    check_sleef_di(host_info, name, ldexp);
#endif
    return (uintptr_t)DL::sym(nullptr, orig_name.c_str());
}

uintptr_t Resolver::find_extern(const std::string &name)
{
    if (m_cb) {
        if (auto ptr = (*m_cb)(name)) {
            return ptr;
        }
    }
    return resolve_ir_sym(name);
}

NACS_EXPORT() Engine::Engine()
{
    reset_dyld();
}

NACS_EXPORT() Engine::~Engine()
{
}

NACS_EXPORT() void Engine::reset_dyld()
{
    m_dyld.reset(new RuntimeDyld(m_memmgr, m_resolver));
}

NACS_EXPORT() uint64_t Engine::load(const object::ObjectFile &obj, const Resolver::cb_t &cb)
{
    if (!m_dyld)
        reset_dyld();
    m_errstr.clear();
    auto id = m_memmgr.new_group();
    Resolver::SetCB setter(m_resolver, cb ? &cb : nullptr);
    // If the loading errors, we should reset the loader so that it does not remember
    // any wrong state. Calling the `free(id)` function accomplishes this.
    if (!m_dyld->loadObject(obj) || m_dyld->hasError()) {
        m_errstr = m_dyld->getErrorString().str();
        free(id);
        return 0;
    }
    m_dyld->finalizeWithMemoryManagerLocking();
    if (m_dyld->hasError()) {
        m_errstr = m_dyld->getErrorString().str();
        free(id);
        return 0;
    }
    return id;
}

NACS_EXPORT() uint64_t Engine::load(const char *p, size_t len, const Resolver::cb_t &cb)
{
    MemoryBufferRef buff(StringRef(p, len), "");
    auto obj = object::ObjectFile::createObjectFile(buff);
    if (!obj)
        return 0;
    return load(*(*obj), cb);
}

NACS_EXPORT() void *Engine::get_symbol(StringRef name) const
{
    return (void*)(uintptr_t)m_dyld->getSymbol(name).getAddress();
}

NACS_EXPORT() void Engine::free(uint64_t id)
{
    m_memmgr.free_group(id);
    // Reusing the memory could really confuse the relocator so we need to reset the dyld
    // when we free any memory.
    // This is also used to clear any sideeffect when the loading process errored.
    m_dyld.reset(nullptr);
}

}
