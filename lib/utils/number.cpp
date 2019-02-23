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

#include "number.h"

#if defined(ENABLE_SIMD) && (NACS_CPU_X86 || NACS_CPU_X86_64)
#  include <immintrin.h>
extern "C" {
// including sleef.h does not include all the function prototypes unless
// the feature is enabled for the whole file.
// We'll just pull the definitions we need instead so that we could do dispatch ourselves.
typedef struct {
  __m128d x, y;
} nacs_m128d_2;

typedef struct {
  __m256d x, y;
} nacs_m256d_2;

typedef struct {
  __m512d x, y;
} nacs_m512d_2;

// The name sleef here is just a placeholder to make sure this expand to dllimport
// Any name other than utils should work here.
NACS_EXPORT(sleef) nacs_m128d_2 Sleef_modfd2_sse2(__m128d);
NACS_EXPORT(sleef) nacs_m128d_2 Sleef_modfd2_sse4(__m128d);
NACS_EXPORT(sleef) nacs_m128d_2 Sleef_modfd2_avx2128(__m128d);

NACS_EXPORT(sleef) nacs_m256d_2 Sleef_modfd4_avx(__m256d);
NACS_EXPORT(sleef) nacs_m256d_2 Sleef_modfd4_avx2(__m256d);

NACS_EXPORT(sleef) nacs_m512d_2 Sleef_modfd8_avx512f(__m512d);

}
#endif

namespace NaCs {

NACS_EXPORT() double linearInterpolate(double x, uint32_t npoints, const double *points)
{
    if (unlikely(x <= 0)) {
        return points[0];
    }
    else if (unlikely(x >= 1)) {
        return points[npoints - 1];
    }
    x = x * (npoints - 1);
    double lof = 0.0;
    x = modf(x, &lof);
    uint32_t lo = (uint32_t)lof;
    double vlo = points[lo];
    if (x == 0)
        return vlo;
    double vhi = points[lo + 1];
    return x * vhi + (1 - x) * vlo;
}

NACS_EXPORT() double linearInterpolate(double x, double x0, double dx,
                                       uint32_t npoints, const double *points)
{
    return linearInterpolate((x - x0) / dx, npoints, points);
}

#if defined(ENABLE_SIMD) && (NACS_CPU_X86 || NACS_CPU_X86_64)

template<int id> static inline nacs_m128d_2 modfd2(__m128d);

template<>
nacs_m128d_2 modfd2<0>(__m128d v)
{
    return Sleef_modfd2_sse2(v);
}

template<>
nacs_m128d_2 modfd2<1>(__m128d v)
{
    return Sleef_modfd2_sse4(v);
}

template<>
nacs_m128d_2 modfd2<2>(__m128d v)
{
    return Sleef_modfd2_avx2128(v);
}

template<int id>
static inline __m128d linearInterpolate2(__m128d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m128)(x > 0);
    auto ovr_ok = (__m128)(x < 1);
    auto ok = _mm_and_ps(ovr_ok, und_ok);
    x = x * (npoints - 1);
    x = (__m128d)_mm_and_ps((__m128)x, ok);
    auto modres = modfd2<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto idx0 = (int)lof[0];
    auto idx1 = (int)lof[1];
    auto vlo = (__m128d){points[idx0], points[idx1]};
    auto vhi = (__m128d){points[idx0 + 1], points[idx1 + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = (__m128d)_mm_and_ps((__m128)res, ok);
    auto und_res = _mm_andnot_ps((__m128)(__m128d){points[0], points[0]}, und_ok);
    auto ovr_res = _mm_andnot_ps((__m128)(__m128d){points[npoints - 1],
                                                   points[npoints - 1]}, ovr_ok);
    res = (__m128d)_mm_or_ps((__m128)res, und_res);
    res = (__m128d)_mm_or_ps((__m128)res, ovr_res);
    return res;
}

template<int id>
static inline __m128d linearInterpolate2(__m128d x, __m128d x0, __m128d dx,
                                         uint32_t npoints, const double *points)
{
    return linearInterpolate2<id>((x - x0) / dx, npoints, points);
}

__attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2<0>(x, npoints, points);
}

__attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2<0>(x, x0, dx, npoints, points);
}

__attribute__((target("sse4.1")))
NACS_EXPORT() __m128d linearInterpolate2_sse4(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2<1>(x, npoints, points);
}

__attribute__((target("sse4.1")))
NACS_EXPORT() __m128d linearInterpolate2_sse4(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2<1>(x, x0, dx, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2<2>(x, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2<2>(x, x0, dx, npoints, points);
}

template<int id> static inline nacs_m256d_2 modfd4(__m256d);

template<>
nacs_m256d_2 modfd4<0>(__m256d v)
{
    return Sleef_modfd4_avx(v);
}

template<>
nacs_m256d_2 modfd4<1>(__m256d v)
{
    return Sleef_modfd4_avx2(v);
}

template<int id> __attribute__((target("avx")))
static inline __m256d linearInterpolate4(__m256d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m256)(x > 0);
    auto ovr_ok = (__m256)(x < 1);
    auto ok = _mm256_and_ps(ovr_ok, und_ok);
    x = x * (npoints - 1);
    x = (__m256d)_mm256_and_ps((__m256)x, ok);
    auto modres = modfd4<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto idx0 = (int)lof[0];
    auto idx1 = (int)lof[1];
    auto idx2 = (int)lof[2];
    auto idx3 = (int)lof[3];
    auto vlo = (__m256d){points[idx0], points[idx1], points[idx2], points[idx3]};
    auto vhi = (__m256d){points[idx0 + 1], points[idx1 + 1],
                         points[idx2 + 1], points[idx3 + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = (__m256d)_mm256_and_ps((__m256)res, ok);
    auto und_res = _mm256_andnot_ps((__m256)_mm256_broadcast_sd(points), und_ok);
    auto ovr_res = _mm256_andnot_ps((__m256)_mm256_broadcast_sd(&points[npoints - 1]), ovr_ok);
    res = (__m256d)_mm256_or_ps((__m256)res, und_res);
    res = (__m256d)_mm256_or_ps((__m256)res, ovr_res);
    return res;
}

template<int id>
static inline __m256d linearInterpolate4(__m256d x, __m256d x0, __m256d dx,
                                         uint32_t npoints, const double *points)
{
    return linearInterpolate4<id>((x - x0) / dx, npoints, points);
}

__attribute__((target("avx")))
NACS_EXPORT() __m256d linearInterpolate4_avx(__m256d x, uint32_t npoints, const double *points)
{
    return linearInterpolate4<0>(x, npoints, points);
}

__attribute__((target("avx")))
NACS_EXPORT() __m256d linearInterpolate4_avx(__m256d x, __m256d x0, __m256d dx,
                                             uint32_t npoints, const double *points)
{
    return linearInterpolate4<0>(x, x0, dx, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m256d linearInterpolate4_avx2(__m256d x, uint32_t npoints, const double *points)
{
    return linearInterpolate4<1>(x, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m256d linearInterpolate4_avx2(__m256d x, __m256d x0, __m256d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate4<1>(x, x0, dx, npoints, points);
}

#endif

}
