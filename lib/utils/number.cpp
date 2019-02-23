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
inline nacs_m128d_2 modfd2<0>(__m128d v)
{
    return Sleef_modfd2_sse2(v);
}

template<>
inline nacs_m128d_2 modfd2<1>(__m128d v)
{
    return Sleef_modfd2_sse4(v);
}

template<>
inline nacs_m128d_2 modfd2<2>(__m128d v)
{
    return Sleef_modfd2_avx2128(v);
}

template<int id> static inline __m128d selectd2(__m128, __m128d, __m128d);

template<>
__m128d selectd2<0>(__m128 mask, __m128d v0, __m128d v1)
{
    return (__m128d)_mm_or_ps(_mm_and_ps(mask, (__m128)v0), _mm_andnot_ps(mask, (__m128)v1));
}

template<> __attribute__((target("sse4.1")))
inline __m128d selectd2<1>(__m128 mask, __m128d v0, __m128d v1)
{
    return _mm_blendv_pd((__m128d)mask, v0, v1);
}

template<> __attribute__((target("avx")))
inline __m128d selectd2<2>(__m128 mask, __m128d v0, __m128d v1)
{
    return _mm_blendv_pd((__m128d)mask, v0, v1);
}

typedef int v4si __attribute__((__vector_size__(16)));

template<int id> __attribute__((always_inline, flatten))
static inline __m128d linearInterpolate2(__m128d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m128)(x > 0);
    auto ovr_ok = (__m128)(x < 1);
    x = x * (npoints - 1);
    x = (__m128d)_mm_and_ps((__m128)x, _mm_and_ps(ovr_ok, und_ok));
    auto modres = modfd2<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = (v4si)_mm_cvtpd_epi32(lof);
    auto vlo = (__m128d){points[lo[0]], points[lo[1]]};
    auto vhi = (__m128d){points[lo[0] + 1], points[lo[1] + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = selectd2<id>(und_ok, res, (__m128d){points[0], points[0]});
    res = selectd2<id>(ovr_ok, res, (__m128d){points[npoints - 1], points[npoints - 1]});
    return res;
}

template<int id> __attribute__((always_inline, flatten))
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
    res = _mm256_blendv_pd((__m256d)und_ok, res, _mm256_broadcast_sd(points));
    res = _mm256_blendv_pd((__m256d)ovr_ok, res, _mm256_broadcast_sd(&points[npoints - 1]));
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
    x = (__m256d)_mm256_and_ps((__m256)x, _mm256_and_ps(ovr_ok, und_ok));
    auto modres = modfd4<1>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = _mm256_cvtpd_epi32(lof);
    auto vlo = _mm256_i32gather_pd(points, lo, 1);
    auto vhi = _mm256_i32gather_pd(points + 1, lo, 1);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd((__m256d)und_ok, res, _mm256_broadcast_sd(points));
    res = _mm256_blendv_pd((__m256d)ovr_ok, res, _mm256_broadcast_sd(&points[npoints - 1]));
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
    __mmask8 ok = und_ok & ovr_ok;
    x = x * (npoints - 1);
    auto modres = modfd8<id>(x);
    x = modres.x;
    auto lof = modres.y;
    auto lo = (v8si)_mm512_cvtpd_epi32(lof);
    auto vlo = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, (__m256i)lo, points, 1);
    auto vhi = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, (__m256i)lo,
                                        (points + 1), 1);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm512_mask_blend_pd(und_ok, res, _mm512_set1_pd(points[0]));
    res = _mm512_mask_blend_pd(ovr_ok, res, _mm512_set1_pd(points[npoints - 1]));
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
