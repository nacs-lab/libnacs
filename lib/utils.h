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

#include "macros.h"

#ifndef __NACS_UTILS_UTILS_H__
#define __NACS_UTILS_UTILS_H__

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
// GNU extension
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef __cplusplus

#include <utility>

template<typename T>
static NACS_INLINE T
nacsSetBit(T orig, uint8_t bit, bool val)
{
    if (val) {
        return orig | (static_cast<T>(1) << bit);
    } else {
        return orig & ~(static_cast<T>(1) << bit);
    }
}

namespace NaCs {

template<typename T>
class ScopeSwap {
    T m_orig_val;
    T *m_var;
    ScopeSwap(const ScopeSwap&) = delete;
    void operator=(const ScopeSwap&) = delete;
public:
    ScopeSwap(T &var, T&&new_val);
    ScopeSwap(ScopeSwap &&other);
    ~ScopeSwap();
};

template<typename T, typename T2>
static inline auto
make_scope_swap(T &var, T2&&new_val)
{
    return ScopeSwap<T>(var, std::forward<T2>(new_val));
}

template<typename T>
ScopeSwap<T>::ScopeSwap(ScopeSwap &&other)
    : m_orig_val(std::move(other.m_orig_val)),
      m_var(other.m_var)
{
    other.m_var = nullptr;
}

template<typename T>
ScopeSwap<T>::ScopeSwap(T &var, T&&new_val)
    : m_orig_val(new_val),
      m_var(&var)
{
    std::swap(m_orig_val, *m_var);
}

template<typename T>
ScopeSwap<T>::~ScopeSwap()
{
    if (m_var) {
        std::swap(m_orig_val, *m_var);
    }
}

}

#else
#define nacsSetBit(orig, bit, val)              \
    ((val) ?                                    \
     ((orig) | (((typeof(orig))1) << uint8_t(bit))) :  \
     ((orig) & ~(((typeof(orig))1) << uint8_t(bit))))
#endif

#endif
