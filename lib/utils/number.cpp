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

namespace NaCs {

__attribute__((always_inline))
static inline double _linearInterpolate(double x, uint32_t npoints, const double *points)
{
    if (unlikely(x <= 0))
        return points[0];
    if (unlikely(x >= 1))
        return points[npoints - 1];
    x = x * (npoints - 1);
    int lo = (int)x;
    x = x - lo;
    double vlo = points[lo];
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

#if NACS_CPU_X86 || NACS_CPU_X86_64

__attribute__((always_inline))
static inline __m128d linearInterpolate2(__m128d x, uint32_t npoints, const double *points)
{
    // From benchmark, on a Skylake CPU, this is actually faster than a vectorized version.
    // (This has higher instruction count but a even higher IPC)
    return __m128d{_linearInterpolate(x[0], npoints, points),
            _linearInterpolate(x[1], npoints, points)};
}

__attribute__((target("avx2,fma"), always_inline))
static inline __m128d _linearInterpolate2_avx2(__m128d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m128)(x > 0);
    auto ovr_ok = (__m128)(x < 1);
    x = x * (npoints - 1);
    auto lo = _mm_cvttpd_epi32(x);
    x = x - _mm_cvtepi32_pd(lo);
    auto ok = (__m128d)_mm_and_ps(ovr_ok, und_ok);
    auto vlo = _mm_mask_i32gather_pd(_mm_undefined_pd(), points, lo, ok, 8);
    auto vhi = _mm_mask_i32gather_pd(_mm_undefined_pd(), (points + 1), lo, ok, 8);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm_blendv_pd(_mm_set1_pd(points[0]), res, (__m128d)und_ok);
    res = _mm_blendv_pd(_mm_set1_pd(points[npoints - 1]), res, (__m128d)ovr_ok);
    return res;
}

typedef int v4si __attribute__((__vector_size__(16)));

__attribute__((target("avx"), always_inline, flatten))
static inline __m256d linearInterpolate4(__m256d x, uint32_t npoints, const double *points)
{
    auto und_ok = (__m256)(x > 0);
    auto ovr_ok = (__m256)(x < 1);
    x = x * (npoints - 1);
    x = (__m256d)_mm256_and_ps((__m256)x, _mm256_and_ps(ovr_ok, und_ok));
    auto lo = (v4si)_mm256_cvttpd_epi32(x);
    x = x - _mm256_cvtepi32_pd((__m128i)lo);
    auto vlo = (__m256d){points[lo[0]], points[lo[1]], points[lo[2]], points[lo[3]]};
    auto vhi = (__m256d){points[lo[0] + 1], points[lo[1] + 1],
                         points[lo[2] + 1], points[lo[3] + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd(_mm256_broadcast_sd(points), res, (__m256d)und_ok);
    res = _mm256_blendv_pd(_mm256_broadcast_sd(&points[npoints - 1]), res, (__m256d)ovr_ok);
    return res;
}

__attribute__((target("avx"), always_inline, flatten))
static inline __m256d linearInterpolate4(__m256d x, __m256d x0, __m256d dx,
                                         uint32_t npoints, const double *points)
{
    return linearInterpolate4((x - x0) / dx, npoints, points);
}

__attribute__((target("avx2,fma"), always_inline))
static inline __m256d _linearInterpolate4_avx2(__m256d x, uint32_t npoints,
                                               const double *points)
{
    auto und_ok = (__m256)(x > 0);
    auto ovr_ok = (__m256)(x < 1);
    x = x * (npoints - 1);
    auto lo = _mm256_cvttpd_epi32(x);
    x = x - _mm256_cvtepi32_pd(lo);
    auto ok = (__m256d)_mm256_and_ps(ovr_ok, und_ok);
    auto vlo = _mm256_mask_i32gather_pd(_mm256_undefined_pd(), points, lo, ok, 8);
    auto vhi = _mm256_mask_i32gather_pd(_mm256_undefined_pd(), (points + 1), lo, ok, 8);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd(_mm256_broadcast_sd(points), res, (__m256d)und_ok);
    res = _mm256_blendv_pd(_mm256_broadcast_sd(&points[npoints - 1]), res, (__m256d)ovr_ok);
    return res;
}

__attribute__((target("avx512f,avx512dq"), always_inline, flatten))
static inline __m512d linearInterpolate8(__m512d x, uint32_t npoints, const double *points)
{
    auto und_ok = _mm512_cmp_pd_mask(x, _mm512_set1_pd(0), _CMP_GT_OS);
    auto ovr_ok = _mm512_cmp_pd_mask(x, _mm512_set1_pd(1), _CMP_LT_OS);
    x = x * (npoints - 1);
    auto lo = _mm512_cvttpd_epi32(x);
    x = x - _mm512_cvtepi32_pd(lo);
    __mmask8 ok = und_ok & ovr_ok;
    auto vlo = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, lo, points, 8);
    auto vhi = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), ok, lo,
                                        (points + 1), 8);
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm512_mask_blend_pd(und_ok, _mm512_set1_pd(points[0]), res);
    res = _mm512_mask_blend_pd(ovr_ok, _mm512_set1_pd(points[npoints - 1]), res);
    return res;
}

__attribute__((target("avx512f,avx512dq"), always_inline, flatten))
static inline __m512d linearInterpolate8(__m512d x, __m512d x0, __m512d dx,
                                         uint32_t npoints, const double *points)
{
    return linearInterpolate8((x - x0) / dx, npoints, points);
}

#  if NACS_OS_WINDOWS
#    if !defined(__GNUC__) || defined(__clang__)
__attribute__((target("sse2"))) NACS_EXPORT()
void linearInterpolate2_sse2(__m128d &res, const __m128d &x, uint32_t npoints,
                             const double *points)
{
    res = linearInterpolate2(x, npoints, points);
}

__attribute__((target("sse2"))) NACS_EXPORT()
void linearInterpolate2_sse2(__m128d &res, const __m128d &x, const __m128d &x0,
                             const __m128d &dx, uint32_t npoints, const double *points)
{
    res = linearInterpolate2((x - x0) / dx, npoints, points);
}

__attribute__((target("avx2,fma"))) NACS_EXPORT()
void linearInterpolate2_avx2(__m128d &res, const __m128d &x, uint32_t npoints,
                             const double *points)
{
    res = _linearInterpolate2_avx2(x, npoints, points);
}

__attribute__((target("avx2,fma"))) NACS_EXPORT()
void linearInterpolate2_avx2(__m128d &res, const __m128d &x, const __m128d &x0,
                             const __m128d &dx, uint32_t npoints, const double *points)
{
    res = _linearInterpolate2_avx2((x - x0) / dx, npoints, points);
}

__attribute__((target("avx"))) NACS_EXPORT()
void linearInterpolate4_avx(__m256d &res, const __m256d &x, uint32_t npoints,
                            const double *points)
{
    res = linearInterpolate4(x, npoints, points);
}

__attribute__((target("avx"))) NACS_EXPORT()
void linearInterpolate4_avx(__m256d &res, const __m256d &x, const __m256d &x0,
                            const __m256d &dx, uint32_t npoints, const double *points)
{
    res = linearInterpolate4(x, x0, dx, npoints, points);
}

__attribute__((target("avx2,fma"))) NACS_EXPORT()
void linearInterpolate4_avx2(__m256d &res, const __m256d &x, uint32_t npoints,
                             const double *points)
{
    res = _linearInterpolate4_avx2(x, npoints, points);
}

__attribute__((target("avx2,fma"))) NACS_EXPORT()
void linearInterpolate4_avx2(__m256d &res, const __m256d &x, const __m256d &x0,
                             const __m256d &dx, uint32_t npoints, const double *points)
{
    res = _linearInterpolate4_avx2((x - x0) / dx, npoints, points);
}
#    endif

__attribute__((target("avx512f,avx512dq"))) NACS_EXPORT()
void linearInterpolate8_avx512f(__m512d &res, const __m512d &x, uint32_t npoints,
                                const double *points)
{
    res = linearInterpolate8(x, npoints, points);
}

__attribute__((target("avx512f,avx512dq"))) NACS_EXPORT()
void linearInterpolate8_avx512f(__m512d &res, const __m512d &x, const __m512d &x0,
                                const __m512d &dx, uint32_t npoints, const double *points)
{
    res = linearInterpolate8(x, x0, dx, npoints, points);
}
#  endif

#  if !NACS_SIMD_USE_REF
NACS_VECTORCALL __attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, uint32_t npoints, const double *points)
{
    return linearInterpolate2(x, npoints, points);
}

NACS_VECTORCALL __attribute__((target("sse2")))
NACS_EXPORT() __m128d linearInterpolate2_sse2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate2((x - x0) / dx, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, uint32_t npoints, const double *points)
{
    return _linearInterpolate2_avx2(x, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx2,fma")))
NACS_EXPORT() __m128d linearInterpolate2_avx2(__m128d x, __m128d x0, __m128d dx,
                                              uint32_t npoints, const double *points)
{
    return _linearInterpolate2_avx2((x - x0) / dx, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx")))
NACS_EXPORT() __m256d linearInterpolate4_avx(__m256d x, uint32_t npoints, const double *points)
{
    return linearInterpolate4(x, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx")))
NACS_EXPORT() __m256d linearInterpolate4_avx(__m256d x, __m256d x0, __m256d dx,
                                             uint32_t npoints, const double *points)
{
    return linearInterpolate4(x, x0, dx, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx2,fma")))
NACS_EXPORT() __m256d linearInterpolate4_avx2(__m256d x, uint32_t npoints, const double *points)
{
    return _linearInterpolate4_avx2(x, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx2,fma")))
NACS_EXPORT() __m256d linearInterpolate4_avx2(__m256d x, __m256d x0, __m256d dx,
                                              uint32_t npoints, const double *points)
{
    return _linearInterpolate4_avx2((x - x0) / dx, npoints, points);
}
#  endif
#  if !NACS_OS_WINDOWS
__attribute__((target("avx512f,avx512dq")))
NACS_EXPORT() __m512d linearInterpolate8_avx512f(__m512d x, uint32_t npoints,
                                                 const double *points)
{
    return linearInterpolate8(x, npoints, points);
}

__attribute__((target("avx512f,avx512dq")))
NACS_EXPORT() __m512d linearInterpolate8_avx512f(__m512d x, __m512d x0, __m512d dx,
                                                 uint32_t npoints, const double *points)
{
    return linearInterpolate8(x, x0, dx, npoints, points);
}
#  endif

#elif NACS_CPU_AARCH64

static NACS_INLINE float64x2_t linearInterpolate2(float64x2_t x, uint32_t npoints,
                                                  const double *points)
{
    auto und_ok = x > 0;
    auto ovr_ok = x < 1;
    auto ok = ovr_ok & und_ok;
    x = x * (npoints - 1);
    x = float64x2_t(ok & uint64x2_t(x));
    auto lo = vcvtq_s64_f64(x);
    x = x - vcvtq_f64_s64(lo);
    auto vlohi = vld2q_f64(&points[lo[0]]);
    vlohi = vld2q_lane_f64(&points[lo[1]], vlohi, 1);
    auto res = x * vlohi.val[1] + (1 - x) * vlohi.val[0];
    res = vbslq_f64(uint64x2_t(ovr_ok), res, vld1q_dup_f64(&points[npoints - 1]));
    return res;
}

NACS_EXPORT() float64x2_t linearInterpolate2_asimd(float64x2_t x, uint32_t npoints,
                                                   const double *points)
{
    return linearInterpolate2(x, npoints, points);
}

NACS_EXPORT() float64x2_t linearInterpolate2_asimd(float64x2_t x, float64x2_t x0,
                                                   float64x2_t dx, uint32_t npoints,
                                                   const double *points)
{
    return linearInterpolate2((x - x0) / dx, npoints, points);
}

// The following versions are not using a wider vector width but should
// remove some duplicated instructions and allow better parallel execution.
// Clang seems to be doing a better job than GCC on this, especially for the x2x4 version.
NACS_EXPORT() float64x2x2_t linearInterpolate4_asimd(float64x2x2_t x, uint32_t npoints,
                                                     const double *points)
{
    return {linearInterpolate2(x.val[0], npoints, points),
            linearInterpolate2(x.val[1], npoints, points)};
}

NACS_EXPORT() float64x2x2_t linearInterpolate4_asimd(float64x2x2_t x, float64x2x2_t x0,
                                                     float64x2x2_t dx, uint32_t npoints,
                                                     const double *points)
{
    return {linearInterpolate2((x.val[0] - x0.val[0]) / dx.val[0], npoints, points),
            linearInterpolate2((x.val[1] - x0.val[1]) / dx.val[1], npoints, points)};
}

NACS_EXPORT() float64x2x4_t linearInterpolate8_asimd(float64x2x4_t x, uint32_t npoints,
                                                     const double *points)
{
    return {linearInterpolate2(x.val[0], npoints, points),
            linearInterpolate2(x.val[1], npoints, points),
            linearInterpolate2(x.val[2], npoints, points),
            linearInterpolate2(x.val[3], npoints, points)};
}

// Pass `dx` by reference since it needs to go through memory according to the ABI anyway
// Give the caller more freedom where that memory should be
NACS_EXPORT() float64x2x4_t linearInterpolate8_asimd(float64x2x4_t x, float64x2x4_t x0,
                                                     const float64x2x4_t &dx,
                                                     uint32_t npoints, const double *points)
{
    return {linearInterpolate2((x.val[0] - x0.val[0]) / dx.val[0], npoints, points),
            linearInterpolate2((x.val[1] - x0.val[1]) / dx.val[1], npoints, points),
            linearInterpolate2((x.val[2] - x0.val[2]) / dx.val[2], npoints, points),
            linearInterpolate2((x.val[3] - x0.val[3]) / dx.val[3], npoints, points)};
}

#endif

}
