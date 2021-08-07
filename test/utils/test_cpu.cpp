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

#define CATCH_CONFIG_MAIN

#include "../../lib/nacs-utils/processor.h"

#include <catch2/catch.hpp>

using namespace NaCs;

TEST_CASE("Host") {
    auto &host_info = CPUInfo::get_host();
    host_info.dump_llvm();
    REQUIRE(&host_info == &CPUInfo::get_host());
    auto dup_host = CPUInfo::create(host_info.get_arch().c_str(),
                                    std::string(host_info).c_str());
    REQUIRE(dup_host->get_arch() == host_info.get_arch());
    REQUIRE(std::string(*dup_host) == std::string(host_info));
#if NACS_CPU_X86
    REQUIRE(host_info.get_arch() == "i386");
#elif NACS_CPU_X86_64
    REQUIRE(host_info.get_arch() == "x86-64");
#elif NACS_CPU_AARCH32
    REQUIRE(host_info.get_arch() == "arm");
#elif NACS_CPU_AARCH64
    REQUIRE(host_info.get_arch() == "aarch64");
#endif
}

TEST_CASE("x86-64") {
    auto skx = CPUInfo::create("x86-64", "skylake-avx512");
    REQUIRE(skx->get_arch() == "x86-64");
    REQUIRE(skx->get_name() == "skylake-avx512");
    REQUIRE(skx->get_vector_size() == 64);
    REQUIRE(skx->test_feature(X86::Feature::avx512f));
    REQUIRE(skx->test_feature(X86::Feature::avx512dq));

    auto skx_no512 = CPUInfo::create("x86-64", "skylake-avx512,-avx512f");
    REQUIRE(skx_no512->get_arch() == "x86-64");
    REQUIRE(skx_no512->get_name() == "skylake-avx512");
    REQUIRE(skx_no512->get_vector_size() == 32);
    REQUIRE(!skx_no512->test_feature(X86::Feature::avx512f));
    REQUIRE(!skx_no512->test_feature(X86::Feature::avx512dq));

    auto skx_no2 = CPUInfo::create("x86-64", "skylake-avx512,-avx2");
    REQUIRE(skx_no2->get_arch() == "x86-64");
    REQUIRE(skx_no2->get_name() == "skylake-avx512");
    REQUIRE(skx_no2->get_vector_size() == 32);
    REQUIRE(skx_no2->test_feature(X86::Feature::avx));
    REQUIRE(skx_no2->test_feature(X86::Feature::fma));
    REQUIRE(!skx_no2->test_feature(X86::Feature::avx2));
    REQUIRE(!skx_no2->test_feature(X86::Feature::avx512f));
    REQUIRE(!skx_no2->test_feature(X86::Feature::avx512dq));

    auto skx_noavx = CPUInfo::create("x86-64", "skylake-avx512,-avx");
    REQUIRE(skx_noavx->get_arch() == "x86-64");
    REQUIRE(skx_noavx->get_name() == "skylake-avx512");
    REQUIRE(skx_noavx->get_vector_size() == 16);
    REQUIRE(!skx_noavx->test_feature(X86::Feature::avx));
    REQUIRE(!skx_noavx->test_feature(X86::Feature::fma));
    REQUIRE(!skx_noavx->test_feature(X86::Feature::avx2));
    REQUIRE(!skx_noavx->test_feature(X86::Feature::avx512f));
    REQUIRE(!skx_noavx->test_feature(X86::Feature::avx512dq));
}

TEST_CASE("arm") {
    auto cortex_a9 = CPUInfo::create("arm", "cortex-a9");
    REQUIRE(cortex_a9->get_arch() == "arm");
    REQUIRE(cortex_a9->get_name() == "cortex-a9");
    REQUIRE(cortex_a9->get_vector_size() == 8);
    REQUIRE(!cortex_a9->test_feature(AArch32::Feature::neon));

    auto cortex_a7 = CPUInfo::create("arm", "cortex-a7");
    REQUIRE(cortex_a7->get_arch() == "arm");
    REQUIRE(cortex_a7->get_name() == "cortex-a7");
    REQUIRE(cortex_a7->get_vector_size() == 16);
    REQUIRE(cortex_a7->test_feature(AArch32::Feature::neon));
}

TEST_CASE("aarch64") {
    auto cortex_a75 = CPUInfo::create("aarch64", "cortex-a75");
    REQUIRE(cortex_a75->get_arch() == "aarch64");
    REQUIRE(cortex_a75->get_name() == "cortex-a75");
    REQUIRE(cortex_a75->get_vector_size() == 16);
    REQUIRE(cortex_a75->test_feature(AArch64::Feature::v8_2a));
    REQUIRE(cortex_a75->test_feature(AArch64::Feature::crc));
    REQUIRE(!cortex_a75->test_feature(AArch64::Feature::aes));
    REQUIRE(!cortex_a75->test_feature(AArch64::Feature::sha2));
}
