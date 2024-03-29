/*************************************************************************
 *   Copyright (c) 2018 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/nacs-seq/zynq/cmdlist.h"
#include "../../lib/nacs-utils/streams.h"
#include "../../lib/nacs-utils/errors.h"
#include "../../lib/nacs-utils/log.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <catch2/catch.hpp>

using namespace NaCs;
using namespace NaCs::Seq::Zynq;

static void test_file_eq(const std::string &fname, const std::string &cmp)
{
    std::ifstream stm(fname);
    REQUIRE(stm.good());
    std::string fstr(std::istreambuf_iterator<char>(stm), {});
    REQUIRE(cmp == fstr);
}

static uint64_t test_cmdlist_eq(const std::string &cmdlist, uint32_t ttl_mask,
                                const std::string &cmp)
{
    uint32_t ver = 1;
    uint64_t len_ns = CmdList::total_time((uint8_t*)cmdlist.data(), cmdlist.size(), ver) * 10;
    auto str_data = (const uint8_t*)cmp.data();
    auto str_sz = cmp.size();
    REQUIRE(memcmp(str_data, &ver, 4) == 0);
    str_data += 4;
    str_sz -= 4;
    REQUIRE(memcmp(str_data, &len_ns, 8) == 0);
    str_data += 8;
    str_sz -= 8;
    REQUIRE(memcmp(str_data, &ttl_mask, 4) == 0);
    str_data += 4;
    str_sz -= 4;
    REQUIRE(str_sz == cmdlist.size());
    REQUIRE(memcmp(str_data, cmdlist.data(), str_sz) == 0);
    return len_ns;
}

static void test(const std::string &dir, const std::string &name)
{
    Log::log("Testing: %s\n", name.c_str());
    auto path = dir + name;
    std::ifstream istm(path);
    REQUIRE(istm.good());
    string_ostream vstm;
    uint32_t ver = 1;
    try {
        uint32_t ttl_mask = CmdList::parse(vstm, istm, ver);
        auto vec = vstm.get_buf();
        std::ifstream bstm(path + ".cmdbin", std::ios::binary);
        REQUIRE(bstm.good());
        std::string binstr(std::istreambuf_iterator<char>(bstm), {});
        uint64_t len_ns = test_cmdlist_eq(vec, ttl_mask, binstr);

        string_ostream tstm;
        tstm << "# " << len_ns << " ns" << std::endl;
        CmdList::print(tstm, (uint8_t*)vec.data(), vec.size(), ttl_mask, ver);
        auto text = tstm.get_buf();
        test_file_eq(path + ".txt", text);

        const_istream tistm(text);
        auto ttl_mask2 = CmdList::parse(vstm, tistm, ver);
        REQUIRE(ttl_mask == ttl_mask2);
        REQUIRE(vec == vstm.get_buf());
    }
    catch (const SyntaxError &err) {
        string_ostream sstr;
        sstr << err;
        test_file_eq(path + ".err", sstr.get_buf());
    }
}

TEST_CASE("CmdList") {
    std::string dir(getenv("TEST_SOURCE_DIR"));
    dir += "/cmdlists/";

    test(dir, "ttl_mask_err1");
    test(dir, "ttl_mask_err2");
    test(dir, "ttl_mask_err3");

    test(dir, "invalid_cmd");
    test(dir, "invalid_cmd2");

    test(dir, "invalid_ttl");
    test(dir, "invalid_ttl1_space");
    test(dir, "invalid_ttlall_space");

    test(dir, "invalid_ttl1_1");
    test(dir, "invalid_ttl1_2");
    test(dir, "invalid_ttl1_3");
    test(dir, "invalid_ttl1_4");

    test(dir, "invalid_ttl_t1");
    test(dir, "invalid_ttl_t2");
    test(dir, "invalid_ttl_t3");

    test(dir, "invalid_wait");
    test(dir, "invalid_wait2");
    test(dir, "invalid_wait3");

    test(dir, "ttl_time");
    test(dir, "ttl_mask_wait");

    test(dir, "dds_freq");
    test(dir, "invalid_dds");
    test(dir, "invalid_freq");
    test(dir, "invalid_freq2");
    test(dir, "invalid_freq_unit");
    test(dir, "invalid_freq_val");
    test(dir, "invalid_freq_val2");

    test(dir, "no_newline_invalid_start");
    test(dir, "no_newline_valid");
}
