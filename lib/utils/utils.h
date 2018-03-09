/*************************************************************************
 *   Copyright (c) 2013 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <utility>
#include <type_traits>
#include <tuple>
#include <string>
#include <map>

namespace NaCs {

/**
 * Tell the compiler that \param val is likely to be \param exp.
 */
template<typename T1, typename T2>
static NACS_INLINE T1 expect(T1 val, T2 exp)
{
#if NACS_CHECK_GCC_VERSION(3, 0)
    return __builtin_expect(val, exp);
#else
    (void)exp;
    return val;
#endif
}

/**
 * Tell the compiler that \param x is likely to be true.
 */
template<typename T>
static NACS_INLINE bool likely(T x)
{
    return expect(bool(x), true);
}

/**
 * Tell the compiler that \param x is likely to be false.
 */
template<typename T>
static NACS_INLINE bool unlikely(T x)
{
    return expect(bool(x), false);
}

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

struct CDeleter {
    template<typename T>
    void operator()(T *p) {
        free((void*)p);
    }
};

template<typename T1, typename T2>
static constexpr bool isBaseOf =
    std::is_base_of<T1, std::decay_t<T2>>::value;

// std::experimental::apply
template<class F, class Tuple, std::size_t... I>
static inline constexpr decltype(auto)
_applyTupleImpl(F &&f, Tuple &&t, std::index_sequence<I...>)
{
    return std::forward<F>(f)(std::get<I>(std::forward<Tuple>(t))...);
}

template<class F, class Tuple>
static inline constexpr decltype(auto) applyTuple(F&& f, Tuple&& t)
{
    return _applyTupleImpl(std::forward<F>(f), std::forward<Tuple>(t),
                           std::make_index_sequence<
                           std::tuple_size<std::decay_t<Tuple>>::value>{});
}

template<typename T>
struct StrCache {
    struct CacheEntry {
        size_t age;
        size_t sz;
        T v;
    };
    StrCache(size_t lim)
        : m_lim(lim),
          m_cache{},
          m_id(0),
          m_totalsz(0)
    {
    }
    const T *get(const std::string &key) const
    {
        auto it = m_cache.find(key);
        if (it == m_cache.end())
            return nullptr;
        return &it->second.v;
    }
    void ejectOldest()
    {
        auto it = m_cache.begin();
        if (it == m_cache.end()) {
            m_totalsz = 0;
            return;
        }
        for (auto it2 = it;it2 != m_cache.end();++it2) {
            if (it2->second.age < it->second.age) {
                it = it2;
            }
        }
        m_totalsz -= it->second.sz;
        m_cache.erase(it);
    }
    const T *set(const std::string &key, T &&v)
    {
        if (const T *r = get(key))
            return r;
        size_t sz = v.cacheSize();
        if (sz >= m_lim)
            return nullptr;
        while (sz + m_totalsz > m_lim)
            ejectOldest();
        auto &entry = m_cache[key];
        entry.age = m_id++;
        entry.sz = sz;
        entry.v = std::move(v);
        m_totalsz += sz;
        return &entry.v;
    }

private:
    size_t m_lim;
    std::map<std::string,CacheEntry> m_cache;
    size_t m_id;
    size_t m_totalsz;
};

}

#endif
