/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_TEST_SEQ_ZYNQ_HELPER_H__
#define __NACS_TEST_SEQ_ZYNQ_HELPER_H__

#include "../../lib/nacs-seq/zynq/bc_gen.h"
#include "../../lib/nacs-seq/zynq/bytecode.h"

#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/processor.h"

#include <sstream>

#include <catch2/catch.hpp>

using namespace NaCs;
using namespace NaCs::Seq;

namespace {

NACS_UNUSED static unsigned vector_size = [] {
#if NACS_CPU_X86 || NACS_CPU_X86_64
    auto &host = CPUInfo::get_host();
    if (host.test_feature(X86::Feature::avx512f) &&
        host.test_feature(X86::Feature::avx512dq)) {
        return 8;
    }
    else if (host.test_feature(X86::Feature::avx)) {
        return 4;
    }
    else {
        return 2;
    }
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    return 2;
#else
    return 1;
#endif
}();

NACS_UNUSED static auto dump_bytecode(const uint8_t *bytecode, size_t size, uint32_t ver)
{
    uint64_t len;
    memcpy(&len, bytecode, sizeof(uint64_t));
    bytecode += sizeof(uint64_t);
    size -= sizeof(uint64_t);
    uint32_t ttl_mask;
    memcpy(&ttl_mask, bytecode, sizeof(uint32_t));
    std::stringstream stm;
    stm << "# len_ns=" << len << std::endl;
    stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
    Zynq::ByteCode::print_raw(stm, bytecode + sizeof(uint32_t),
                              size - sizeof(uint32_t), ver);
    return stm.str();
}

NACS_UNUSED static auto dump_bytecode(const std::vector<uint8_t> &bytecode, uint32_t ver)
{
    return dump_bytecode(bytecode.data(), bytecode.size(), ver);
}

NACS_UNUSED static auto dump_bytecode(const Zynq::BCGen &bc_gen)
{
    return dump_bytecode(bc_gen.bytecode, bc_gen.version());
}

struct Checker {
    Checker(const uint8_t *data, size_t size)
        : data(data),
          size(size)
    {}
    Checker(const std::vector<uint8_t> &data)
        : Checker(data.data(), data.size())
    {}
    template<typename T>
    bool cmp(T v)
    {
        REQUIRE(sizeof(T) <= size);
        bool res = memcmp(data, &v, sizeof(T)) == 0;
        data += sizeof(T);
        size -= sizeof(T);
        return res;
    }
    bool is_end() const
    {
        return size == 0;
    }

    const uint8_t *data;
    size_t size;
};

}

#endif
