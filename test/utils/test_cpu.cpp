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

#include "../../lib/utils/processor.h"
#include <assert.h>

using namespace NaCs;

int main()
{
    auto &host_info = CPUInfo::get_host();
    assert(&host_info == &CPUInfo::get_host());
    auto dup_host = CPUInfo::create(host_info.get_arch().c_str(),
                                    std::string(host_info).c_str());
    assert(dup_host->get_arch() == host_info.get_arch());
    assert(std::string(*dup_host) == std::string(host_info));
#if NACS_CPU_X86
    assert(host_info.get_arch() == "i386");
#elif NACS_CPU_X86_64
    assert(host_info.get_arch() == "x86-64");
#elif NACS_CPU_AARCH32
    assert(host_info.get_arch() == "arm");
#elif NACS_CPU_AARCH64
    assert(host_info.get_arch() == "aarch64");
#endif

    auto skx = CPUInfo::create("x86-64", "skylake-avx512");
    assert(skx->get_arch() == "x86-64");
    assert(skx->get_name() == "skylake-avx512");
    assert(skx->get_vector_size() == 64);
    assert(skx->test_feature(X86::Feature::avx512f));
    assert(skx->test_feature(X86::Feature::avx512dq));

    auto skx_no512 = CPUInfo::create("x86-64", "skylake-avx512,-avx512f");
    assert(skx_no512->get_arch() == "x86-64");
    assert(skx_no512->get_name() == "skylake-avx512");
    assert(skx_no512->get_vector_size() == 32);
    assert(!skx_no512->test_feature(X86::Feature::avx512f));
    assert(!skx_no512->test_feature(X86::Feature::avx512dq));

    auto skx_no2 = CPUInfo::create("x86-64", "skylake-avx512,-avx2");
    assert(skx_no2->get_arch() == "x86-64");
    assert(skx_no2->get_name() == "skylake-avx512");
    assert(skx_no2->get_vector_size() == 32);
    assert(skx_no2->test_feature(X86::Feature::avx));
    assert(skx_no2->test_feature(X86::Feature::fma));
    assert(!skx_no2->test_feature(X86::Feature::avx2));
    assert(!skx_no2->test_feature(X86::Feature::avx512f));
    assert(!skx_no2->test_feature(X86::Feature::avx512dq));

    auto skx_noavx = CPUInfo::create("x86-64", "skylake-avx512,-avx");
    assert(skx_noavx->get_arch() == "x86-64");
    assert(skx_noavx->get_name() == "skylake-avx512");
    assert(skx_noavx->get_vector_size() == 16);
    assert(!skx_noavx->test_feature(X86::Feature::avx));
    assert(!skx_noavx->test_feature(X86::Feature::fma));
    assert(!skx_noavx->test_feature(X86::Feature::avx2));
    assert(!skx_noavx->test_feature(X86::Feature::avx512f));
    assert(!skx_noavx->test_feature(X86::Feature::avx512dq));

    auto cortex_a9 = CPUInfo::create("arm", "cortex-a9");
    assert(cortex_a9->get_arch() == "arm");
    assert(cortex_a9->get_name() == "cortex-a9");
    assert(cortex_a9->get_vector_size() == 8);
    assert(!cortex_a9->test_feature(AArch32::Feature::neon));

    auto cortex_a7 = CPUInfo::create("arm", "cortex-a7");
    assert(cortex_a7->get_arch() == "arm");
    assert(cortex_a7->get_name() == "cortex-a7");
    assert(cortex_a7->get_vector_size() == 16);
    assert(cortex_a7->test_feature(AArch32::Feature::neon));

    auto cortex_a75 = CPUInfo::create("aarch64", "cortex-a75");
    assert(cortex_a75->get_arch() == "aarch64");
    assert(cortex_a75->get_name() == "cortex-a75");
    assert(cortex_a75->get_vector_size() == 16);
    assert(cortex_a75->test_feature(AArch64::Feature::v8_2a));
    assert(cortex_a75->test_feature(AArch64::Feature::crc));
    assert(cortex_a75->test_feature(AArch64::Feature::aes));
    assert(cortex_a75->test_feature(AArch64::Feature::sha2));

    return 0;
}
