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

#ifndef _NACS_UTILS_NUMBER_H_
#define _NACS_UTILS_NUMBER_H_

#include <cmath>
#include <type_traits>
#include <utility>

namespace NaCs {

template<typename T1, typename T2>
static inline constexpr auto
max(T1 &&a, T2 &&b)
{
    return (a > b) ? a : b;
}

template<typename First, typename... Rest>
static inline constexpr auto
max(First &&first, Rest&&... rest)
{
    return max(std::forward<First>(first),
               max(std::forward<Rest>(rest)...));
}

template<typename T1, typename T2>
static inline constexpr auto
min(T1 &&a, T2 &&b)
{
    return (a < b) ? a : b;
}

template<typename First, typename... Rest>
static inline constexpr auto
min(First &&first, Rest&&... rest)
{
    return min(std::forward<First>(first),
               min(std::forward<Rest>(rest)...));
}

template<typename T1, typename T2, typename T3>
static inline constexpr auto
bound(T1 &&v1, T2 &&v2, T3 &&v3)
{
    return max(std::forward<T1>(v1), min(std::forward<T2>(v2),
                                         std::forward<T3>(v3)));
}

template<typename T>
static inline auto
square(const T &a)
{
    return a * a;
}

template<typename T1, typename T2>
static inline constexpr auto
limit(T1 &&v1, T2 &&v2)
{
    return bound(0, std::forward<T1>(v1), std::forward<T2>(v2));
}

template<typename T1, typename T2>
static inline auto
getPadding(T1 len, T2 align) -> decltype(len + align)
{
    if (auto left = len % align) {
        return align - left;
    }
    return 0;
}

template<typename T1, typename T2>
static inline auto
alignTo(T1 len, T2 align)
{
    return len + getPadding(len, align);
}

template<typename T, int upper, int lower>
static constexpr inline std::enable_if_t<(upper == lower), int>
getBits(T)
{
    static_assert(std::is_unsigned<T>(), "");
    return lower;
}

template<typename T, int upper=sizeof(T) * 8, int lower=0>
static constexpr inline std::enable_if_t<(upper > lower), int>
getBits(T v)
{
    static_assert(std::is_unsigned<T>(), "");
    constexpr int mid = (upper + lower) / 2;
    constexpr T thresh = 1 << mid;
    if (v >= thresh) {
        return getBits<T, upper, mid + 1>(v);
    } else {
        return getBits<T, mid, lower>(v);
    }
}

}

#endif
