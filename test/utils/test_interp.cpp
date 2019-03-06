/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/utils/number.h"
#include "../../lib/utils/timer.h"
#include <assert.h>
#include <math.h>

using namespace NaCs;

const double points[] = {0, 0.1, 0.2, 0.6};

static bool approx(double a, double b)
{
    double diff = abs(a - b);
    double avg = abs(a + b) / 2;
    return diff < 2e-8 || diff / avg < 2e-8;
}

static void print_avg(Timer &timer, const char *prefix, size_t nele)
{
    auto ttotal = timer.elapsed();
    auto tavg = double(ttotal) / double(nele);
    Log::log("%s: %f ns\n", prefix, tavg);
}

#if NACS_CPU_X86 || NACS_CPU_X86_64
static bool has_avx = false;
static bool has_avx2 = false;
static bool has_avx512 = false;

static void init_cpu()
{
    __builtin_cpu_init();
    assert(__builtin_cpu_supports("sse2"));
    has_avx = __builtin_cpu_supports("avx");
    has_avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    has_avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
}

__attribute__((target("sse2"))) static void test_sse2()
{
    auto res = linearInterpolate2_sse2(__m128d{2.3, 3.4}, __m128d{2, 2},
                                       __m128d{3, 3}, 4, points);
    assert(approx(res[0], 0.03));
    assert(approx(res[1], 0.14));
    res = linearInterpolate2_sse2(__m128d{2, 5}, __m128d{2, 2},
                                  __m128d{3, 3}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.6));

    res = linearInterpolate2_sse2(__m128d{0.3, 1.4} / 3, 4, points);
    assert(approx(res[0], 0.03));
    assert(approx(res[1], 0.14));
    res = linearInterpolate2_sse2(__m128d{-0.1, 1.1}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.6));

    Timer timer;
    for (int i = 0; i < 100000000 / 2; i++)
        linearInterpolate2_sse2(__m128d{0.3, 1.4} / 3, 4, points);
    print_avg(timer, "SSE2<2>", 100000000);
}

__attribute__((target("avx"))) void test_avx()
{
    auto res = linearInterpolate4_avx(__m256d{1.0, 2.3, 3.4, 5.0}, __m256d{2, 2, 2, 2},
                                      __m256d{3, 3, 3, 3}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.03));
    assert(approx(res[2], 0.14));
    assert(approx(res[3], 0.6));
    res = linearInterpolate4_avx(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.03));
    assert(approx(res[2], 0.14));
    assert(approx(res[3], 0.6));

    Timer timer;
    for (int i = 0; i < 100000000 / 4; i++)
        linearInterpolate4_avx(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    print_avg(timer, "AVX<4>", 100000000);
}

__attribute__((target("avx2,fma"))) void test_avx2()
{
    auto res2 = linearInterpolate2_avx2(__m128d{2.3, 3.4}, __m128d{2, 2},
                                        __m128d{3, 3}, 4, points);
    assert(approx(res2[0], 0.03));
    assert(approx(res2[1], 0.14));
    res2 = linearInterpolate2_avx2(__m128d{2, 5}, __m128d{2, 2},
                                   __m128d{3, 3}, 4, points);
    assert(approx(res2[0], 0.0));
    assert(approx(res2[1], 0.6));

    res2 = linearInterpolate2_avx2(__m128d{0.3, 1.4} / 3, 4, points);
    assert(approx(res2[0], 0.03));
    assert(approx(res2[1], 0.14));
    res2 = linearInterpolate2_avx2(__m128d{-0.1, 1.1}, 4, points);
    assert(approx(res2[0], 0.0));
    assert(approx(res2[1], 0.6));

    auto res4 = linearInterpolate4_avx2(__m256d{1.0, 2.3, 3.4, 5.0}, __m256d{2, 2, 2, 2},
                                        __m256d{3, 3, 3, 3}, 4, points);
    assert(approx(res4[0], 0.0));
    assert(approx(res4[1], 0.03));
    assert(approx(res4[2], 0.14));
    assert(approx(res4[3], 0.6));
    res4 = linearInterpolate4_avx2(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    assert(approx(res4[0], 0.0));
    assert(approx(res4[1], 0.03));
    assert(approx(res4[2], 0.14));
    assert(approx(res4[3], 0.6));

    Timer timer;
    for (int i = 0; i < 100000000 / 2; i++)
        linearInterpolate2_avx2(__m128d{0.3, 1.4} / 3, 4, points);
    print_avg(timer, "AVX2<2>", 100000000);

    timer.restart();
    for (int i = 0; i < 100000000 / 4; i++)
        linearInterpolate4_avx2(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    print_avg(timer, "AVX2<4>", 100000000);
}

__attribute__((target("avx512f,avx512dq"))) void test_avx512()
{
    auto res = linearInterpolate8_avx512f(__m512d{1.0, 2.3, 3.4, 5.0, 1.0, 2.3, 3.4, 5.0},
                                          __m512d{2, 2, 2, 2, 2, 2, 2, 2},
                                          __m512d{3, 3, 3, 3, 3, 3, 3, 3}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.03));
    assert(approx(res[2], 0.14));
    assert(approx(res[3], 0.6));
    assert(approx(res[4], 0.0));
    assert(approx(res[5], 0.03));
    assert(approx(res[6], 0.14));
    assert(approx(res[7], 0.6));
    res = linearInterpolate8_avx512f(__m512d{-0.1, 0.3 / 3, 1.4 / 3, 1.1,
                                             -0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    assert(approx(res[0], 0.0));
    assert(approx(res[1], 0.03));
    assert(approx(res[2], 0.14));
    assert(approx(res[3], 0.6));
    assert(approx(res[4], 0.0));
    assert(approx(res[5], 0.03));
    assert(approx(res[6], 0.14));
    assert(approx(res[7], 0.6));

    Timer timer;
    for (int i = 0; i < 100000000 / 8; i++)
        linearInterpolate8_avx512f(__m512d{-0.1, 0.3 / 3, 1.4 / 3, 1.1,
                                           -0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    print_avg(timer, "AVX512<8>", 100000000);
}
#elif NACS_CPU_AARCH64
static void test_asimd()
{
    auto res2 = linearInterpolate2_asimd(float64x2_t{2.3, 3.4}, float64x2_t{2, 2},
                                         float64x2_t{3, 3}, 4, points);
    assert(approx(res2[0], 0.03));
    assert(approx(res2[1], 0.14));
    res2 = linearInterpolate2_asimd(float64x2_t{2, 5}, float64x2_t{2, 2},
                                    float64x2_t{3, 3}, 4, points);
    assert(approx(res2[0], 0.0));
    assert(approx(res2[1], 0.6));

    res2 = linearInterpolate2_asimd(float64x2_t{0.3, 1.4} / 3, 4, points);
    assert(approx(res2[0], 0.03));
    assert(approx(res2[1], 0.14));
    res2 = linearInterpolate2_asimd(float64x2_t{-0.1, 1.1}, 4, points);
    assert(approx(res2[0], 0.0));
    assert(approx(res2[1], 0.6));

    auto res4 = linearInterpolate4_asimd({{float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0}}},
                                         {{float64x2_t{2, 2}, float64x2_t{2, 2}}},
                                         {{float64x2_t{3, 3}, float64x2_t{3, 3}}}, 4, points);
    assert(approx(res4.val[0][0], 0.0));
    assert(approx(res4.val[0][1], 0.03));
    assert(approx(res4.val[1][0], 0.14));
    assert(approx(res4.val[1][1], 0.6));
    res4 = linearInterpolate4_asimd({{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}},
                                    4, points);
    assert(approx(res4.val[0][0], 0.0));
    assert(approx(res4.val[0][1], 0.03));
    assert(approx(res4.val[1][0], 0.14));
    assert(approx(res4.val[1][1], 0.6));

    auto res8 = linearInterpolate8_asimd(
        {{float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0},
          float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0}}},
        {{float64x2_t{2, 2}, float64x2_t{2, 2}, float64x2_t{2, 2}, float64x2_t{2, 2}}},
        {{float64x2_t{3, 3}, float64x2_t{3, 3}, float64x2_t{3, 3}, float64x2_t{3, 3}}},
        4, points);
    assert(approx(res8.val[0][0], 0.0));
    assert(approx(res8.val[0][1], 0.03));
    assert(approx(res8.val[1][0], 0.14));
    assert(approx(res8.val[1][1], 0.6));
    assert(approx(res8.val[2][0], 0.0));
    assert(approx(res8.val[2][1], 0.03));
    assert(approx(res8.val[3][0], 0.14));
    assert(approx(res8.val[3][1], 0.6));
    res8 = linearInterpolate8_asimd(
        {{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1},
          float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}}, 4, points);
    assert(approx(res8.val[0][0], 0.0));
    assert(approx(res8.val[0][1], 0.03));
    assert(approx(res8.val[1][0], 0.14));
    assert(approx(res8.val[1][1], 0.6));
    assert(approx(res8.val[2][0], 0.0));
    assert(approx(res8.val[2][1], 0.03));
    assert(approx(res8.val[3][0], 0.14));
    assert(approx(res8.val[3][1], 0.6));


    Timer timer;
    for (int i = 0; i < 100000000 / 2; i++)
        linearInterpolate2_asimd(float64x2_t{0.3, 1.4} / 3, 4, points);
    print_avg(timer, "ASIMD<2>", 100000000);

    timer.restart();
    for (int i = 0; i < 100000000 / 4; i++)
        linearInterpolate4_asimd({{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}},
                                 4, points);
    print_avg(timer, "ASIMD<4>", 100000000);

    timer.restart();
    for (int i = 0; i < 100000000 / 8; i++)
        linearInterpolate8_asimd(
            {{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1},
              float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}}, 4, points);
    print_avg(timer, "ASIMD<8>", 100000000);
}
#endif

int main()
{
    assert(approx(linearInterpolate(0, 2, 3, 4, points), 0));
    assert(approx(linearInterpolate(2, 2, 3, 4, points), 0));
    assert(approx(linearInterpolate(2.3, 2, 3, 4, points), 0.03));
    assert(approx(linearInterpolate(3.4, 2, 3, 4, points), 0.14));
    assert(approx(linearInterpolate(4.5, 2, 3, 4, points), 0.4));
    assert(approx(linearInterpolate(5, 2, 3, 4, points), 0.6));
    assert(approx(linearInterpolate(7, 2, 3, 4, points), 0.6));

    assert(approx(linearInterpolate(-1, 4, points), 0));
    assert(approx(linearInterpolate(0, 4, points), 0));
    assert(approx(linearInterpolate(0.1, 4, points), 0.03));
    assert(approx(linearInterpolate(1.4 / 3, 4, points), 0.14));
    assert(approx(linearInterpolate(2.5 / 3, 4, points), 0.4));
    assert(approx(linearInterpolate(1, 4, points), 0.6));
    assert(approx(linearInterpolate(2, 4, points), 0.6));

    Timer timer;
    for (int i = 0; i < 100000000; i++)
        linearInterpolate((3.4 - 2) / 3, 4, points);
    print_avg(timer, "Scalar", 100000000);

#if NACS_CPU_X86 || NACS_CPU_X86_64
    init_cpu();
    test_sse2();
    if (has_avx)
        test_avx();
    if (has_avx2)
        test_avx2();
    if (has_avx512)
        test_avx512();
#elif NACS_CPU_AARCH64
    test_asimd();
#endif

    return 0;
}
