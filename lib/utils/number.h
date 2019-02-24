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

#ifndef _NACS_UTILS_NUMBER_H_
#define _NACS_UTILS_NUMBER_H_

#include "utils.h"

#include <stdint.h>

#include <cmath>

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
static inline constexpr auto
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

template<typename T, uint8_t upper, uint8_t lower>
static constexpr inline std::enable_if_t<(upper == lower), uint8_t>
getBits(T)
{
    static_assert(std::is_unsigned<T>(), "");
    return lower;
}

template<typename T, uint8_t upper=sizeof(T) * 8, uint8_t lower=0>
static constexpr inline std::enable_if_t<(upper > lower), uint8_t>
getBits(T v)
{
    static_assert(std::is_unsigned<T>(), "");
    constexpr uint8_t mid = (upper + lower) / 2;
    constexpr T thresh = static_cast<T>(1) << mid;
    if (v >= thresh) {
        return getBits<T, upper, mid + 1>(v);
    } else {
        return getBits<T, mid, lower>(v);
    }
}

template<typename T>
static constexpr inline T
setBit(T orig, uint8_t bit, bool val)
{
    if (val) {
        return orig | (static_cast<T>(1) << bit);
    } else {
        return orig & ~(static_cast<T>(1) << bit);
    }
}

template<typename T>
struct _SumSingle {
    template<typename T2>
    constexpr inline auto
    operator()(T2 &&v)
    {
        return v;
    }
};

static constexpr struct {
    template<typename First>
    inline constexpr auto
    operator()(First &&first) const
    {
        return _SumSingle<std::decay_t<First> >()(std::forward<First>(first));
    }
    template<typename First, typename... Rest>
    inline constexpr auto
    operator()(First &&first, Rest&&... rest) const
    {
        return (operator()(std::forward<First>(first)) +
                operator()(std::forward<Rest>(rest)...));
    }
} sumAll{};

template<typename... Args>
struct _SumSingle<std::tuple<Args...> > {
    template<typename Tuple>
    constexpr inline auto
    operator()(Tuple &&tuple)
    {
        return applyTuple(sumAll, std::forward<Tuple>(tuple));
    }
};

NACS_EXPORT(utils) double linearInterpolate(double x, uint32_t npoints, const double *points);
NACS_EXPORT(utils) double linearInterpolate(double x, double x0, double dx,
                                            uint32_t npoints, const double *points);

#if defined(ENABLE_SIMD) && (NACS_CPU_X86 || NACS_CPU_X86_64)
NACS_EXPORT(utils) __m128d linearInterpolate2_sse2(__m128d x, uint32_t npoints,
                                                   const double *points);
NACS_EXPORT(utils) __m128d linearInterpolate2_sse2(__m128d x, __m128d x0, __m128d dx,
                                                   uint32_t npoints, const double *points);
NACS_EXPORT(utils) __m128d linearInterpolate2_avx2(__m128d x, uint32_t npoints,
                                                   const double *points);
NACS_EXPORT(utils) __m128d linearInterpolate2_avx2(__m128d x, __m128d x0, __m128d dx,
                                                   uint32_t npoints, const double *points);
NACS_EXPORT(utils) __m256d linearInterpolate4_avx(__m256d x, uint32_t npoints,
                                                  const double *points);
NACS_EXPORT(utils) __m256d linearInterpolate4_avx(__m256d x, __m256d x0, __m256d dx,
                                                  uint32_t npoints, const double *points);
NACS_EXPORT(utils) __m256d linearInterpolate4_avx2(__m256d x, uint32_t npoints,
                                                   const double *points);
NACS_EXPORT(utils) __m256d linearInterpolate4_avx2(__m256d x, __m256d x0, __m256d dx,
                                                   uint32_t npoints, const double *points);
#endif

}

#endif
