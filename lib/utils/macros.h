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

#ifndef _NACS_UTILS_MACROS_H_
#define _NACS_UTILS_MACROS_H_

/**
 * \file macros.h
 * \author Yichao Yu <yyc1992@gmail.com>
 * \brief Definitions of some useful macros.
 */

#define NACS_OS_FREEBSD 0
#define NACS_OS_LINUX 0
#define NACS_OS_WINDOWS 0
#define NACS_OS_DARWIN 0

#if defined(__FreeBSD__)
#  undef NACS_OS_FREEBSD
#  define NACS_OS_FREEBSD 1
#elif defined(__linux__)
#  undef NACS_OS_LINUX
#  define NACS_OS_LINUX 1
#elif defined(_WIN32) || defined(_WIN64)
#  undef NACS_OS_WINDOWS
#  define NACS_OS_WINDOWS 1
#elif defined(__APPLE__) && defined(__MACH__)
#  undef NACS_OS_DARWIN
#  define NACS_OS_DARWIN 1
#endif

#define NACS_CPU_X86_64 0
#define NACS_CPU_X86 0
#define NACS_CPU_AARCH64 0
#define NACS_CPU_AARCH32 0
#define NACS_CPU_PPC64 0
#define NACS_CPU_PPC32 0

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || \
    defined(__x86_64) || defined(_M_X64) || defined(_M_AMD64)
#  undef NACS_CPU_X86_64
#  define NACS_CPU_X86_64 1
#elif defined(i386) || defined(__i386) || defined(__i386__) || defined(_M_IX86) || defined(_X86_)
#  undef NACS_CPU_X86
#  define NACS_CPU_X86 1
#elif defined(__aarch64__)
#  undef NACS_CPU_AARCH64
#  define NACS_CPU_AARCH64 1
#elif defined(__arm__) || defined(_M_ARM)
#  undef NACS_CPU_AARCH32
#  define NACS_CPU_AARCH32 1
#elif defined(__PPC64__)
#  undef NACS_CPU_PPC64
#  define NACS_CPU_PPC64 1
#elif defined(_ARCH_PPC)
#  undef NACS_CPU_PPC32
#  define NACS_CPU_PPC32 1
#endif

#ifdef __has_builtin
#  define nacs_has_builtin(x) __has_builtin(x)
#else
#  define nacs_has_builtin(x) 0
#endif

/** \defgroup nacs_switch Macros for detecting empty arguments
 * \brief Used to implement function overloading and default arguments in c.
 *
 * The idea of this implementation is borrowed from the following
 * [post](https://gustedt.wordpress.com/2010/06/08/detect-empty-macro-arguments)
 * and is modified in order to fit with our usage.
 * @{
 */

#define __NACS_USE_3(_1, _2, _3, ...) _3
// Test if args has one comma
#define __NACS_HAS_COMMA1(ret_true, ret_false, args...) \
    __NACS_USE_3(args, ret_true, ret_false)
// Convert parentheses to comma, used to check if the next character is "("
#define __NACS_CONVERT_PAREN(...) ,
// Check if arg starts with "("
#define __NACS_IS_PAREN_(ret_true, ret_false, arg)                      \
    __NACS_HAS_COMMA1(ret_true, ret_false, __NACS_CONVERT_PAREN arg)
// Extra layer just to make sure more evaluation (if any) is done than the
// separator path.
#define __NACS_IS_PAREN(ret_true, ret_false, arg)       \
    __NACS_IS_PAREN_(ret_true, ret_false, arg)
// Check if arg is not empty and does not start with "("
// Will not work if arg has comma or is the name of a function like macro
#define __NACS_IS_SEP(ret_true, ret_false, arg)                         \
    __NACS_HAS_COMMA1(ret_false, ret_true, __NACS_CONVERT_PAREN arg())
#define __NACS_IS_EMPTY_PAREN_TRUE(ret_true, ret_false, arg) ret_false
#define __NACS_IS_EMPTY_PAREN_FALSE(ret_true, ret_false, arg)   \
    __NACS_IS_SEP(ret_false, ret_true, arg)

/**
 * \brief Test if \param arg is empty.
 * Evaluate to \param ret_true if it is not empty,
 * evaluate to \param ret_false otherwise.
 *
 * NOTE: this may not work if \param arg is a macro that can be evaluated to a
 * comma separate list without parentheses or is the name of
 * a function like macro.
 */
#define NACS_SWITCH(arg, ret_true, ret_false)           \
    __NACS_IS_PAREN(__NACS_IS_EMPTY_PAREN_TRUE,         \
                    __NACS_IS_EMPTY_PAREN_FALSE, arg)   \
    (ret_false, ret_true, arg)

/**
 * Evaluate to \param def if \param v is empty and to \param v otherwise.
 * \sa NACS_SWITCH for restrictions.
 **/
#define NACS_DEFAULT(v, def) NACS_SWITCH(v, v, def)

/** @} */

#define nacsMakeVersion(a, b, c...)                     \
    ((a) << 16 | (b) << 8 | (NACS_DEFAULT(c, 0)))
#ifdef __GNUC__
#  define NACS_GCC_VERSION                                              \
    nacsMakeVersion(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
#else
#  define NACS_GCC_VERSION 0
#endif
#define NACS_CHECK_GCC_VERSION(args...)         \
    (NACS_GCC_VERSION >= nacsMakeVersion(args))

/**
 * \brief always inline the function.
 *
 * Should only be used for small functions
 */
#define NACS_INLINE __attribute__((always_inline)) inline
#define NACS_NOINLINE __attribute__((noinline))

/**
 * \brief Export symbol.
 */
#if NACS_OS_WINDOWS
#  define NACS_EXPORT_SWITCH_REAL(s)                                    \
    NACS_SWITCH(s, __declspec(dllimport), __declspec(dllexport))
#  define NACS_EXPORT_SWITCH(lib)                       \
    NACS_EXPORT_SWITCH_REAL(NACS_EXPORT_LIB_##lib())
#  define NACS_EXPORT(lib...)                                           \
    NACS_SWITCH(lib, NACS_EXPORT_SWITCH(lib), __declspec(dllexport))
#  define NACS_PROTECTED(lib...) NACS_EXPORT(lib)
#  define NACS_INTERNAL
#else
#  define NACS_EXPORT(...) __attribute__((visibility("default")))
#  define NACS_PROTECTED(...) __attribute__((visibility("protected")))
#  define NACS_INTERNAL __attribute__((visibility("internal")))
#endif

/**
 * \def NACS_BEGIN_DECLS
 * \brief start declaring c-linked functions.
 */
/**
 * \def NACS_END_DECLS
 * \brief end declaring c-linked functions.
 */
// For public c headers
#ifdef __cplusplus
#  define NACS_BEGIN_DECLS extern "C" {
#  define NACS_END_DECLS }
#else
#  define NACS_BEGIN_DECLS
#  define NACS_END_DECLS
#endif

/**
 * Suppress unused parameter warning on variable \param x.
 */
#define NACS_UNUSED __attribute__((unused))

#ifdef __GNUC__
#  define NACS_NORETURN __attribute__((noreturn))
#elif defined(_MSC_VER)
#  define NACS_NORETURN __declspec(noreturn)
#else
#  define NACS_NORETURN
#endif

/**
 * cast the \param member pointer \param ptr of a structure \param type to
 * the containing structure.
 */
#define nacsContainerOf(ptr, type, member)              \
    ((type*)(((char*)(ptr)) - offsetof(type, member)))

#define NACS_RET_IF_FAIL(exp, val...) do {              \
        if (!NaCs::likely(exp)) {                       \
            return (NACS_DEFAULT(val, (void)0));        \
        }                                               \
    } while (0)

#if (defined(__x86_64__) || defined(__x86_64) || defined(__i386) || defined(__i386__)) && \
    !defined(__clang__)
#  define NACS_PACKED __attribute__((__packed__, gcc_struct))
#else
#  define NACS_PACKED __attribute__((__packed__))
#endif

#endif
