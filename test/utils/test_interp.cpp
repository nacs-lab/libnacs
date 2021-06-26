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

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN

#include "../../lib/utils/number.h"
#include <math.h>

#include <catch2/catch.hpp>

using namespace NaCs;

const double points[] = {0, 0.1, 0.2, 0.6};

#if NACS_CPU_X86 || NACS_CPU_X86_64
static bool has_avx = false;
static bool has_avx2 = false;
static bool has_avx512 = false;

#if defined(__clang__) && __clang_major__ < 10
// Clang < 10.0 has problem handling inline asm
// with vector register enabled by target attribute.
#  define ENABLE_VECTOR_BENCH 0
#else
#  define ENABLE_VECTOR_BENCH 1
#endif

static void init_cpu()
{
    __builtin_cpu_init();
    REQUIRE(__builtin_cpu_supports("sse2"));
    has_avx = __builtin_cpu_supports("avx");
    has_avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    has_avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq");
}

__attribute__((target("sse2"))) static void test_sse2()
{
    auto res = linearInterpolate2_sse2(__m128d{2.3, 3.4}, __m128d{2, 2},
                                       __m128d{3, 3}, 4, points);
    REQUIRE(res[0] == Approx(0.03));
    REQUIRE(res[1] == Approx(0.14));
    res = linearInterpolate2_sse2(__m128d{2, 5}, __m128d{2, 2},
                                  __m128d{3, 3}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.6));

    res = linearInterpolate2_sse2(__m128d{0.3, 1.4} / 3, 4, points);
    REQUIRE(res[0] == Approx(0.03));
    REQUIRE(res[1] == Approx(0.14));
    res = linearInterpolate2_sse2(__m128d{-0.1, 1.1}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.6));

    BENCHMARK("SSE2<2>") __attribute__((target("sse2"))) {
        auto res = linearInterpolate2_sse2(__m128d{0.3, 1.4} / 3, 4, points);
        // Do not return the result since the caller may not have the correct ABI declared
        // Instead, use an inline assembly to convince the compiler that the result is used.
        asm volatile ("" :: "v"(res));
    };
}

__attribute__((target("avx"))) void test_avx()
{
    auto res = linearInterpolate4_avx(__m256d{1.0, 2.3, 3.4, 5.0}, __m256d{2, 2, 2, 2},
                                      __m256d{3, 3, 3, 3}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.03));
    REQUIRE(res[2] == Approx(0.14));
    REQUIRE(res[3] == Approx(0.6));
    res = linearInterpolate4_avx(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.03));
    REQUIRE(res[2] == Approx(0.14));
    REQUIRE(res[3] == Approx(0.6));

#if ENABLE_VECTOR_BENCH
    BENCHMARK("AVX<4>") __attribute__((target("avx"))) {
        auto res = linearInterpolate4_avx(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
        // Do not return the result since the caller may not have the correct ABI declared
        // Instead, use an inline assembly to convince the compiler that the result is used.
        asm volatile ("" :: "v"(res));
    };
#endif
}

__attribute__((target("avx2,fma"))) void test_avx2()
{
    auto res2 = linearInterpolate2_avx2(__m128d{2.3, 3.4}, __m128d{2, 2},
                                        __m128d{3, 3}, 4, points);
    REQUIRE(res2[0] == Approx(0.03));
    REQUIRE(res2[1] == Approx(0.14));
    res2 = linearInterpolate2_avx2(__m128d{2, 5}, __m128d{2, 2},
                                   __m128d{3, 3}, 4, points);
    REQUIRE(res2[0] == Approx(0.0));
    REQUIRE(res2[1] == Approx(0.6));

    res2 = linearInterpolate2_avx2(__m128d{0.3, 1.4} / 3, 4, points);
    REQUIRE(res2[0] == Approx(0.03));
    REQUIRE(res2[1] == Approx(0.14));
    res2 = linearInterpolate2_avx2(__m128d{-0.1, 1.1}, 4, points);
    REQUIRE(res2[0] == Approx(0.0));
    REQUIRE(res2[1] == Approx(0.6));

    auto res4 = linearInterpolate4_avx2(__m256d{1.0, 2.3, 3.4, 5.0}, __m256d{2, 2, 2, 2},
                                        __m256d{3, 3, 3, 3}, 4, points);
    REQUIRE(res4[0] == Approx(0.0));
    REQUIRE(res4[1] == Approx(0.03));
    REQUIRE(res4[2] == Approx(0.14));
    REQUIRE(res4[3] == Approx(0.6));
    res4 = linearInterpolate4_avx2(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    REQUIRE(res4[0] == Approx(0.0));
    REQUIRE(res4[1] == Approx(0.03));
    REQUIRE(res4[2] == Approx(0.14));
    REQUIRE(res4[3] == Approx(0.6));

#if ENABLE_VECTOR_BENCH
    BENCHMARK("AVX2<2>") __attribute__((target("avx2,fma"))) {
        auto res = linearInterpolate2_avx2(__m128d{0.3, 1.4} / 3, 4, points);
        // Do not return the result since the caller may not have the correct ABI declared
        // Instead, use an inline assembly to convince the compiler that the result is used.
        asm volatile ("" :: "v"(res));
    };
    BENCHMARK("AVX2<4>") __attribute__((target("avx2,fma"))) {
        auto res = linearInterpolate4_avx2(__m256d{-0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
        // Do not return the result since the caller may not have the correct ABI declared
        // Instead, use an inline assembly to convince the compiler that the result is used.
        asm volatile ("" :: "v"(res));
    };
#endif
}

__attribute__((target("avx512f,avx512dq"))) void test_avx512()
{
    auto res = linearInterpolate8_avx512f(__m512d{1.0, 2.3, 3.4, 5.0, 1.0, 2.3, 3.4, 5.0},
                                          __m512d{2, 2, 2, 2, 2, 2, 2, 2},
                                          __m512d{3, 3, 3, 3, 3, 3, 3, 3}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.03));
    REQUIRE(res[2] == Approx(0.14));
    REQUIRE(res[3] == Approx(0.6));
    REQUIRE(res[4] == Approx(0.0));
    REQUIRE(res[5] == Approx(0.03));
    REQUIRE(res[6] == Approx(0.14));
    REQUIRE(res[7] == Approx(0.6));
    res = linearInterpolate8_avx512f(__m512d{-0.1, 0.3 / 3, 1.4 / 3, 1.1,
                                             -0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
    REQUIRE(res[0] == Approx(0.0));
    REQUIRE(res[1] == Approx(0.03));
    REQUIRE(res[2] == Approx(0.14));
    REQUIRE(res[3] == Approx(0.6));
    REQUIRE(res[4] == Approx(0.0));
    REQUIRE(res[5] == Approx(0.03));
    REQUIRE(res[6] == Approx(0.14));
    REQUIRE(res[7] == Approx(0.6));

#if ENABLE_VECTOR_BENCH
    BENCHMARK("AVX512<8>") __attribute__((target("avx512f,avx512dq"))) {
        auto res = linearInterpolate8_avx512f(__m512d{-0.1, 0.3 / 3, 1.4 / 3, 1.1,
                -0.1, 0.3 / 3, 1.4 / 3, 1.1}, 4, points);
        // Do not return the result since the caller may not have the correct ABI declared
        // Instead, use an inline assembly to convince the compiler that the result is used.
        asm volatile ("" :: "v"(res));
    };
#endif
}
#elif NACS_CPU_AARCH64
static void test_asimd()
{
    auto res2 = linearInterpolate2_asimd(float64x2_t{2.3, 3.4}, float64x2_t{2, 2},
                                         float64x2_t{3, 3}, 4, points);
    REQUIRE(res2[0] == Approx(0.03));
    REQUIRE(res2[1] == Approx(0.14));
    res2 = linearInterpolate2_asimd(float64x2_t{2, 5}, float64x2_t{2, 2},
                                    float64x2_t{3, 3}, 4, points);
    REQUIRE(res2[0] == Approx(0.0));
    REQUIRE(res2[1] == Approx(0.6));

    res2 = linearInterpolate2_asimd(float64x2_t{0.3, 1.4} / 3, 4, points);
    REQUIRE(res2[0] == Approx(0.03));
    REQUIRE(res2[1] == Approx(0.14));
    res2 = linearInterpolate2_asimd(float64x2_t{-0.1, 1.1}, 4, points);
    REQUIRE(res2[0] == Approx(0.0));
    REQUIRE(res2[1] == Approx(0.6));

    auto res4 = linearInterpolate4_asimd({{float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0}}},
                                         {{float64x2_t{2, 2}, float64x2_t{2, 2}}},
                                         {{float64x2_t{3, 3}, float64x2_t{3, 3}}}, 4, points);
    REQUIRE(res4.val[0][0] == Approx(0.0));
    REQUIRE(res4.val[0][1] == Approx(0.03));
    REQUIRE(res4.val[1][0] == Approx(0.14));
    REQUIRE(res4.val[1][1] == Approx(0.6));
    res4 = linearInterpolate4_asimd({{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}},
                                    4, points);
    REQUIRE(res4.val[0][0] == Approx(0.0));
    REQUIRE(res4.val[0][1] == Approx(0.03));
    REQUIRE(res4.val[1][0] == Approx(0.14));
    REQUIRE(res4.val[1][1] == Approx(0.6));

    auto res8 = linearInterpolate8_asimd(
        {{float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0},
          float64x2_t{1.0, 2.3}, float64x2_t{3.4, 5.0}}},
        {{float64x2_t{2, 2}, float64x2_t{2, 2}, float64x2_t{2, 2}, float64x2_t{2, 2}}},
        {{float64x2_t{3, 3}, float64x2_t{3, 3}, float64x2_t{3, 3}, float64x2_t{3, 3}}},
        4, points);
    REQUIRE(res8.val[0][0] == Approx(0.0));
    REQUIRE(res8.val[0][1] == Approx(0.03));
    REQUIRE(res8.val[1][0] == Approx(0.14));
    REQUIRE(res8.val[1][1] == Approx(0.6));
    REQUIRE(res8.val[2][0] == Approx(0.0));
    REQUIRE(res8.val[2][1] == Approx(0.03));
    REQUIRE(res8.val[3][0] == Approx(0.14));
    REQUIRE(res8.val[3][1] == Approx(0.6));
    res8 = linearInterpolate8_asimd(
        {{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1},
          float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}}, 4, points);
    REQUIRE(res8.val[0][0] == Approx(0.0));
    REQUIRE(res8.val[0][1] == Approx(0.03));
    REQUIRE(res8.val[1][0] == Approx(0.14));
    REQUIRE(res8.val[1][1] == Approx(0.6));
    REQUIRE(res8.val[2][0] == Approx(0.0));
    REQUIRE(res8.val[2][1] == Approx(0.03));
    REQUIRE(res8.val[3][0] == Approx(0.14));
    REQUIRE(res8.val[3][1] == Approx(0.6));

    BENCHMARK("ASIMD<2>") {
        return linearInterpolate2_asimd(float64x2_t{0.3, 1.4} / 3, 4, points);
    };
    BENCHMARK("ASIMD<4>") {
        return linearInterpolate4_asimd({{float64x2_t{-0.1, 0.3 / 3},
                    float64x2_t{1.4 / 3, 1.1}}}, 4, points);
    };
    BENCHMARK("ASIMD<8>") {
        return linearInterpolate8_asimd(
            {{float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1},
              float64x2_t{-0.1, 0.3 / 3}, float64x2_t{1.4 / 3, 1.1}}}, 4, points);
    };
}
#endif

TEST_CASE("Scalar") {
    REQUIRE(linearInterpolate(0, 2, 3, 4, points) == Approx(0));
    REQUIRE(linearInterpolate(2, 2, 3, 4, points) == Approx(0));
    REQUIRE(linearInterpolate(2.3, 2, 3, 4, points) == Approx(0.03));
    REQUIRE(linearInterpolate(3.4, 2, 3, 4, points) == Approx(0.14));
    REQUIRE(linearInterpolate(4.5, 2, 3, 4, points) == Approx(0.4));
    REQUIRE(linearInterpolate(5, 2, 3, 4, points) == Approx(0.6));
    REQUIRE(linearInterpolate(7, 2, 3, 4, points) == Approx(0.6));

    REQUIRE(linearInterpolate(-1, 4, points) == Approx(0));
    REQUIRE(linearInterpolate(0, 4, points) == Approx(0));
    REQUIRE(linearInterpolate(0.1, 4, points) == Approx(0.03));
    REQUIRE(linearInterpolate(1.4 / 3, 4, points) == Approx(0.14));
    REQUIRE(linearInterpolate(2.5 / 3, 4, points) == Approx(0.4));
    REQUIRE(linearInterpolate(1, 4, points) == Approx(0.6));
    REQUIRE(linearInterpolate(2, 4, points) == Approx(0.6));

    BENCHMARK("Scalar") {
        return linearInterpolate((3.4 - 2) / 3, 4, points);
    };
}

#if NACS_CPU_X86 || NACS_CPU_X86_64
TEST_CASE("Vector (x86)") {
    init_cpu();
    test_sse2();
    if (has_avx)
        test_avx();
    if (has_avx2)
        test_avx2();
    if (has_avx512)
        test_avx512();
}
#elif NACS_CPU_AARCH64
TEST_CASE("Vector (aarch64)") {
    test_asimd();
}
#endif
