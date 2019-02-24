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
  __m256d x, y;
} nacs_m256d_2;

typedef struct {
  __m512d x, y;
} nacs_m512d_2;

// The name sleef here is just a placeholder to make sure this expand to dllimport
// Any name other than utils should work here.
NACS_EXPORT(sleef) nacs_m256d_2 Sleef_modfd4_avx(__m256d);
NACS_EXPORT(sleef) nacs_m256d_2 Sleef_modfd4_avx2(__m256d);

NACS_EXPORT(sleef) nacs_m512d_2 Sleef_modfd8_avx512f(__m512d);

}
#endif

namespace NaCs {

__attribute__((always_inline))
static inline double _linearInterpolate(double x, uint32_t npoints, const double *points)
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

NACS_EXPORT() double linearInterpolate(double x, uint32_t npoints, const double *points)
{
    return _linearInterpolate(x, npoints, points);
}

NACS_EXPORT() double linearInterpolate(double x, double x0, double dx,
                                       uint32_t npoints, const double *points)
{
    return _linearInterpolate((x - x0) / dx, npoints, points);
}

#if defined(ENABLE_SIMD) && (NACS_CPU_X86 || NACS_CPU_X86_64)

__attribute__((always_inline))
static inline __m128d linearInterpolate2(__m128d x, uint32_t npoints, const double *points)
{
    // From benchmark, on a Skylake CPU, this is actually faster than a vectorized version.
    // (This has higher instruction count but a even higher IPC)
    return __m128d{_linearInterpolate(x[0], npoints, points),
            _linearInterpolate(x[1], npoints, points)};
}

__attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2(x, npoints, points);
}

__attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2((x - x0) / dx, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2(x, npoints, points);
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2((x - x0) / dx, npoints, points);
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

typedef int v4si __attribute__((__vector_size__(16)));

template<int id> __attribute__((target("avx"), always_inline, flatten))
static inline __m256d linearInterpolate4(__m256d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m256)(x > 0);
    auto ovr_ok = (__m256)(x < 1);
    x = x * (npoints - 1);
    x = (__m256d)_mm256_and_ps((__m256)x, _mm256_and_ps(ovr_ok, und_ok));
    auto modres = modfd4<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = (v4si)_mm256_cvtpd_epi32(lof);
    auto vlo = (__m256d){points[lo[0]], points[lo[1]], points[lo[2]], points[lo[3]]};
    auto vhi = (__m256d){points[lo[0] + 1], points[lo[1] + 1],
                         points[lo[2] + 1], points[lo[3] + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd(_mm256_broadcast_sd(points), res, (__m256d)und_ok);
    res = _mm256_blendv_pd(_mm256_broadcast_sd(&points[npoints - 1]), res, (__m256d)ovr_ok);
    return res;
}

template<int id> __attribute__((target("avx"), always_inline, flatten))
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
    auto und_ok = (__m256)(x > 0);
    auto ovr_ok = (__m256)(x < 1);
    x = x * (npoints - 1);
    auto modres = modfd4<1>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = _mm256_cvtpd_epi32(lof);
    auto ok = (__m256d)_mm256_and_ps(ovr_ok, und_ok);
    auto vlo = _mm256_mask_i32gather_pd(_mm256_undefined_pd(), points, lo, ok, 1);
    auto vhi = _mm256_mask_i32gather_pd(_mm256_undefined_pd(), (points + 1), lo, ok, 1);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd(_mm256_broadcast_sd(points), res, (__m256d)und_ok);
    res = _mm256_blendv_pd(_mm256_broadcast_sd(&points[npoints - 1]), res, (__m256d)ovr_ok);
    return res;
}

__attribute__((target("avx2,fma")))
NACS_EXPORT() __m256d linearInterpolate4_avx2(__m256d x, __m256d x0, __m256d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate4_avx2((x - x0) / dx, npoints, points);
}

template<int id> static inline nacs_m512d_2 modfd8(__m512d);

template<>
nacs_m512d_2 modfd8<0>(__m512d v)
{
    return Sleef_modfd8_avx512f(v);
}

typedef int v8si __attribute__((__vector_size__(32)));

template<int id> __attribute__((target("avx512f,avx512dq"), always_inline, flatten))
static inline __m512d linearInterpolate8(__m512d x, uint32_t npoints, const double *points)
{
    auto und_ok = _mm512_cmp_pd_mask(x, _mm512_set1_pd(0), _CMP_GT_OS);
    auto ovr_ok = _mm512_cmp_pd_mask(x, _mm512_set1_pd(1), _CMP_LT_OS);
    x = x * (npoints - 1);
    auto modres = modfd8<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = (v8si)_mm512_cvtpd_epi32(lof);
    __mmask8 ok = und_ok & ovr_ok;
    auto vlo = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, (__m256i)lo, points, 1);
    auto vhi = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, (__m256i)lo,
                                        (points + 1), 1);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm512_mask_blend_pd(und_ok, _mm512_set1_pd(points[0]), res);
    res = _mm512_mask_blend_pd(ovr_ok, _mm512_set1_pd(points[npoints - 1]), res);
    return res;
}

template<int id> __attribute__((target("avx512f,avx512dq"), always_inline, flatten))
static inline __m512d linearInterpolate8(__m512d x, __m512d x0, __m512d dx,
                                         uint32_t npoints, const double *points)
{
    return linearInterpolate8<id>((x - x0) / dx, npoints, points);
}

__attribute__((target("avx512f,avx512dq")))
NACS_EXPORT() __m512d linearInterpolate8_avx512f(__m512d x, uint32_t npoints,
                                                 const double *points)
{
    return linearInterpolate8<0>(x, npoints, points);
}

__attribute__((target("avx512f,avx512dq")))
NACS_EXPORT() __m512d linearInterpolate8_avx512f(__m512d x, __m512d x0, __m512d dx,
                                                 uint32_t npoints, const double *points)
{
    return linearInterpolate8<0>(x, x0, dx, npoints, points);
}

#endif

}
