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

#ifndef _NACS_UTILS_MACROS_H_
#define _NACS_UTILS_MACROS_H_

/**
 * \file macros.h
 * \author Yichao Yu <yyc1992@gmail.com>
 * \brief Definitions of some useful macros.
 */

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

/**
 * Evaluate to _\param f if \param v is empty and to \param f otherwise.
 * \sa NACS_SWITCH for restrictions.
 **/
#define NACS_SWITCH_(v, f) NACS_SWITCH(v, f, _##f)

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
 * \brief Export symbol.
 */
#define NACS_EXPORT __attribute__((visibility("default")))

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
 * \brief always inline the function.
 *
 * Should only be used for small functions
 */
#define NACS_INLINE __attribute__((always_inline)) inline

/**
 * Suppress unused parameter warning on variable \param x.
 */
#define NACS_UNUSED(x) ((void)(x))

/**
 * cast the \param member pointer \param ptr of a structure \param type to
 * the containing structure.
 */
#define nacsContainerOf(ptr, type, member)              \
    ((type*)(((void*)(ptr)) - offsetof(type, member)))

/**
 * Tell the compiler that \param exp is likely to be \param var.
 */
#if NACS_CHECK_GCC_VERSION(3, 0)
#  define nacsExpect(exp, var) __builtin_expect(exp, var)
#else
#  define nacsExpect(exp, var) (exp)
#endif

/**
 * Tell the compiler that \param x is likely to be true.
 */
#define nacsLikely(x) nacsExpect(!!(x), 1)

/**
 * Tell the compiler that \param x is likely to be false.
 */
#define nacsUnlikely(x) nacsExpect(!!(x), 0)

#define NACS_RET_IF_FAIL(exp, val...) do {              \
        if (!nacsLikely(exp)) {                         \
            return (NACS_DEFAULT(val, (void)0));        \
        }                                               \
    } while (0)

#endif
