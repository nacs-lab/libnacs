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

#if NACS_CPU_X86_64 || NACS_CPU_X86
#  include <immintrin.h>
#endif

namespace NaCs {

/**
 * Debug helper
 */
NACS_EXPORT(utils) void breakpoint(const void *p=nullptr);

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

#if nacs_has_builtin(__builtin_assume)
template<typename T>
static NACS_INLINE T assume(T v)
{
    __builtin_assume(bool(v));
    return v;
}
#elif defined(__GNUC__)
template<typename T>
static NACS_INLINE T assume(T v)
{
    if (!bool(v))
        __builtin_unreachable();
    return v;
}
#else
template<typename T>
static NACS_INLINE T assume(T v)
{
    return v;
}
#endif

namespace CPU {
#ifdef __MIC__
static NACS_INLINE void pause()
{
    _mm_delay_64(100);
}
static NACS_INLINE void wake()
{
}
#elif NACS_CPU_X86_64 || NACS_CPU_X86  /* !__MIC__ */
static NACS_INLINE void pause()
{
    _mm_pause();
}
static NACS_INLINE void wake()
{
}
#elif NACS_CPU_AARCH64 || (NACS_CPU_AARCH32 && __ARM_ARCH >= 7)
static NACS_INLINE void pause()
{
    __asm__ volatile ("wfe" ::: "memory");
}
static NACS_INLINE void wake()
{
    __asm__ volatile ("sev" ::: "memory");
}
#else
static NACS_INLINE void pause()
{
}
static NACS_INLINE void wake()
{
}
#endif
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
static constexpr bool isBaseOf = std::is_base_of_v<T1, std::decay_t<T2>>;

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

// Library Fundamentals TS v2. Most likely replaced by concept in C++20.
namespace detail {
template <class Default, class AlwaysVoid,
          template<class...> class Op, class... Args>
struct detector {
    using value_t = std::false_type;
    using type = Default;
};

template <class Default, template<class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
    using value_t = std::true_type;
    using type = Op<Args...>;
};

struct nonesuch {
    ~nonesuch() = delete;
    nonesuch(nonesuch const&) = delete;
    void operator=(nonesuch const&) = delete;
};

} // namespace detail

template <template<class...> class Op, class... Args>
using is_detected = typename detail::detector<detail::nonesuch, void, Op, Args...>::value_t;

template <template<class...> class Op, class... Args>
using detected_t = typename detail::detector<detail::nonesuch, void, Op, Args...>::type;

template <class Default, template<class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template< template<class...> class Op, class... Args >
constexpr inline bool is_detected_v = is_detected<Op, Args...>::value;
template< class Default, template<class...> class Op, class... Args >
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;
template <class Expected, template<class...> class Op, class... Args>
constexpr inline bool is_detected_exact_v = is_detected_exact<Expected, Op, Args...>::value;

template <class To, template<class...> class Op, class... Args>
using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;
template <class To, template<class...> class Op, class... Args>
constexpr inline bool is_detected_convertible_v =
    is_detected_convertible<To, Op, Args...>::value;

namespace detail {
// Detect if method exist.
// See https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function
template<typename T,typename V>
using has_postvisit_t =
    decltype(std::declval<std::remove_reference_t<T>>().postvisit(std::declval<V>()));
template<typename T,typename V>
inline constexpr auto has_postvisit_v = is_detected_v<has_postvisit_t,T,V>;

template<typename T,bool>
struct postvisit_caller {
    template<typename V>
    static void call(T&&, V&&)
    {
    }
};

template<typename T>
struct postvisit_caller<T,true> {
    template<typename V>
    static void call(T &&visitor, V &&v)
    {
        visitor.postvisit(std::forward<V>(v));
    }
};

template<typename T,typename V>
void postvisit(T &&visitor, V &&v)
{
    postvisit_caller<T,has_postvisit_v<T,V>>::call(std::forward<T>(visitor),
                                                   std::forward<V>(v));
}

}

template<typename T>
struct FieldIterator {
    // FieldIterator(T);
    // T get();
    // T parent();
    // bool is_end();
    // FieldIterator &operator++();
};

template<typename T, typename Visitor>
void visit_dfs(T v, Visitor &&visitor);

template<typename T,typename Iterator,template<typename...> class Vector>
struct DFSVisitor {
    // Process an object before scanning. Return `true` if it should be scanned.
    // `postvisit` will be called on this visit IFF this returns `true`.
    // bool previsit(T); // Required
    // void postvisit(T); // Optional

private:
    using iterator = Iterator;
    Vector<Iterator> m_stack;
    template<typename T2, typename Visitor> friend void visit_dfs(T2 v, Visitor &&visitor);
};

template<typename T, typename Visitor>
void visit_dfs(T v, Visitor &&visitor)
{
    using Iterator = typename std::decay_t<Visitor>::iterator;
    // Use a manual stack for better optimization and less actual stack usage so that
    // we don't need to worry about stack overflow.
    // If postvisit isn't used,
    // all items on the stack are inbound (see `push` below)
    // and `pop` does not need bounds check
    auto &stack = visitor.m_stack;
    Iterator it(v);
    constexpr bool has_postvisit = detail::has_postvisit_v<Visitor,T>;
    auto next = [&] {
        ++it;
        // Next argument
        if (!it.is_end())
            return true;
        detail::postvisit(visitor, it.parent());
        // Done
        if (stack.empty())
            return false;
        it = stack.back();
        stack.pop_back();
        if (has_postvisit) {
            while (it.is_end()) {
                detail::postvisit(visitor, it.parent());
                if (stack.empty())
                    return false;
                it = stack.back();
                stack.pop_back();
            }
        }
        return true;
    };
    auto iterate = [&] {
        auto new_v = it.get();
        if (!new_v || !visitor.previsit(new_v))
            return next();
        // Save the next one to be processed to the stack
        // If we are already at the last one and don't have postvisit,
        // we don't need to do anything.
        ++it;
        if (has_postvisit || !it.is_end())
            stack.push_back(std::move(it));
        it = Iterator(new_v);
        return true;
    };
    while (iterate()) {
    }
}

}

#endif
