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

static NACS_INLINE double _linearInterpolate(double x, uint32_t npoints, const double *points)
{
    if (unlikely(x >= 1))
        return points[npoints - 1];
    x = max(x, 0);
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

static NACS_INLINE __m128d linearInterpolate2(__m128d x, uint32_t npoints, const double *points)
{
    // From benchmark, on a Skylake CPU, this is actually faster than a vectorized version.
    // (This has higher instruction count but a even higher IPC)
    return __m128d{_linearInterpolate(x[0], npoints, points),
            _linearInterpolate(x[1], npoints, points)};
}

__attribute__((target("avx2,fma")))
static NACS_INLINE __m128d _linearInterpolate2_avx2(__m128d x, uint32_t npoints,
                                                    const double *points)
{
    auto end_pt = _mm_set1_pd(points[npoints - 1]);
    auto ovr_ok = (__m128d)(x < 1);
    x = _mm_max_pd(x, _mm_set1_pd(0));
    x = x * (npoints - 1);
    auto lo = _mm_cvttpd_epi32(x);
    x = x - _mm_cvtepi32_pd(lo);
    auto vlo = _mm_mask_i32gather_pd(end_pt, points, lo, ovr_ok, 8);
    auto vhi = _mm_mask_i32gather_pd(end_pt, (points + 1), lo, ovr_ok, 8);
    return x * vhi + (1 - x) * vlo;
}

typedef int v4si __attribute__((__vector_size__(16)));

__attribute__((target("avx"), flatten))
static NACS_INLINE __m256d linearInterpolate4(__m256d x, uint32_t npoints, const double *points)
{
    auto ovr_ok = (__m256)(x < 1);
    x = _mm256_max_pd(x, _mm256_set1_pd(0));
    x = x * (npoints - 1);
    x = (__m256d)_mm256_and_ps((__m256)x, ovr_ok);
    auto lo = (v4si)_mm256_cvttpd_epi32(x);
    x = x - _mm256_cvtepi32_pd((__m128i)lo);
    auto vlo = (__m256d){points[lo[0]], points[lo[1]], points[lo[2]], points[lo[3]]};
    auto vhi = (__m256d){points[lo[0] + 1], points[lo[1] + 1],
                         points[lo[2] + 1], points[lo[3] + 1]};
    auto res = x * vhi + (1 - x) * vlo;
    res = _mm256_blendv_pd(_mm256_broadcast_sd(&points[npoints - 1]), res, (__m256d)ovr_ok);
    return res;
}

__attribute__((target("avx"), flatten))
static NACS_INLINE __m256d linearInterpolate4(__m256d x, __m256d x0, __m256d dx,
                                              uint32_t npoints, const double *points)
{
    return linearInterpolate4((x - x0) / dx, npoints, points);
}

__attribute__((target("avx2,fma")))
static NACS_INLINE __m256d _linearInterpolate4_avx2(__m256d x, uint32_t npoints,
                                                    const double *points)
{
    auto end_pt = _mm256_broadcast_sd(&points[npoints - 1]);
    auto ovr_ok = (__m256d)(x < 1);
    x = _mm256_max_pd(x, _mm256_set1_pd(0));
    x = x * (npoints - 1);
    auto lo = _mm256_cvttpd_epi32(x);
    x = x - _mm256_cvtepi32_pd(lo);
    auto vlo = _mm256_mask_i32gather_pd(end_pt, points, lo, ovr_ok, 8);
    auto vhi = _mm256_mask_i32gather_pd(end_pt, (points + 1), lo, ovr_ok, 8);
    return x * vhi + (1 - x) * vlo;
}

__attribute__((target("avx512f,avx512dq"), flatten))
static NACS_INLINE __m512d linearInterpolate8(__m512d x, uint32_t npoints, const double *points)
{
    auto end_pt = _mm512_set1_pd(points[npoints - 1]);
    auto ovr_ok = _mm512_cmp_pd_mask(x, _mm512_set1_pd(1), _CMP_LT_OS);
    x = _mm512_max_pd(x, _mm512_set1_pd(0));
    x = x * (npoints - 1);
    auto lo = _mm512_cvttpd_epi32(x);
    x = x - _mm512_cvtepi32_pd(lo);
    auto vlo = _mm512_mask_i32gather_pd(end_pt, ovr_ok, lo, points, 8);
    auto vhi = _mm512_mask_i32gather_pd(end_pt, ovr_ok, lo, (points + 1), 8);
    return x * vhi + (1 - x) * vlo;
}

__attribute__((target("avx512f,avx512dq"), flatten))
static NACS_INLINE __m512d linearInterpolate8(__m512d x, __m512d x0, __m512d dx,
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
#    endif
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

NACS_VECTORCALL __attribute__((target("avx512f,avx512dq")))
NACS_EXPORT() __m512d linearInterpolate8_avx512f(__m512d x, uint32_t npoints,
                                                 const double *points)
{
    return linearInterpolate8(x, npoints, points);
}

NACS_VECTORCALL __attribute__((target("avx512f,avx512dq")))
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
    auto ovr_ok = x < 1;
    x = vmaxnmq_f64(x, vdupq_n_f64(0));
    x = x * (npoints - 1);
    x = float64x2_t(ovr_ok & uint64x2_t(x));
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

#ifdef NACS_UTILS_HAS_SVE
static NACS_INLINE __attribute__((target("+sve")))
svfloat64_t _linearInterpolate_sve(svfloat64_t x, uint32_t npoints,
                                   const double *points)
{
    auto ptrue = svptrue_b64();
    auto ovr_ok = svcmplt(ptrue, x, 1);
    x = svmaxnm_x(ptrue, x, 0);
    x = svmul_x(ptrue, x, npoints - 1);
    auto lo = svcvt_s64_x(ptrue, x);
    x = svsub_x(ptrue, x, svcvt_f64_x(ptrue, lo));
    auto vlo = svld1_gather_index(ovr_ok, points, lo);
    auto vhi = svld1_gather_index(ovr_ok, points + 1, lo);
    // Load the overflow value into vlo where needed
    // The mask for svmad_m will make sure
    // these are the values being used in the final result.
    vlo = svdup_f64_m(vlo, svnot_z(ptrue, ovr_ok), points[npoints - 1]);
    return svmad_m(ovr_ok, vlo, svsubr_x(ptrue, x, 1), svmul_x(ptrue, x, vhi));
}

NACS_EXPORT() __attribute__((target("+sve")))
svfloat64_t linearInterpolate_sve(svfloat64_t x, uint32_t npoints, const double *points)
{
    return _linearInterpolate_sve(x, npoints, points);
}

NACS_EXPORT() __attribute__((target("+sve")))
svfloat64_t linearInterpolate_sve(svfloat64_t x, svfloat64_t x0, svfloat64_t dx,
                                  uint32_t npoints, const double *points)
{
    auto ptrue = svptrue_b64();
    x = svdiv_x(ptrue, svsub_x(ptrue, x, x0), dx);
    return _linearInterpolate_sve(x, npoints, points);
}
#endif

#endif

}
