/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/seq/cmdlist.h"
#include "../../lib/utils/streams.h"
#include "../../lib/utils/errors.h"
#include "../../lib/utils/log.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace NaCs;

static void test_file_eq(const std::string &fname, const std::string &cmp)
{
    std::ifstream stm(fname);
    assert(stm.good());
    std::string fstr(std::istreambuf_iterator<char>(stm), {});
    assert(cmp == fstr);
}

static uint64_t test_cmdlist_eq(const std::string &cmdlist, uint32_t ttl_mask,
                                const std::string &cmp)
{
    uint64_t len_ns = Seq::CmdList::total_time((uint8_t*)cmdlist.data(), cmdlist.size()) * 10;
    auto str_data = (const uint8_t*)cmp.data();
    auto str_sz = cmp.size();
    uint32_t ver = 1;
    assert(memcmp(str_data, &ver, 4) == 0);
    str_data += 4;
    str_sz -= 4;
    assert(memcmp(str_data, &len_ns, 8) == 0);
    str_data += 8;
    str_sz -= 8;
    assert(memcmp(str_data, &ttl_mask, 4) == 0);
    str_data += 4;
    str_sz -= 4;
    assert(str_sz == cmdlist.size());
    assert(memcmp(str_data, cmdlist.data(), str_sz) == 0);
    return len_ns;
}

static void test(const std::string &dir, const std::string &name)
{
    Log::log("Testing: %s\n", name.c_str());
    auto path = dir + name;
    std::ifstream istm(path);
    assert(istm.good());
    string_ostream vstm;
    try {
        uint32_t ttl_mask = Seq::CmdList::parse(vstm, istm);
        auto vec = vstm.get_buf();
        std::ifstream bstm(path + ".cmdbin");
        assert(bstm.good());
        std::string binstr(std::istreambuf_iterator<char>(bstm), {});
        uint64_t len_ns = test_cmdlist_eq(vec, ttl_mask, binstr);

        string_ostream tstm;
        tstm << "# " << len_ns << " ns" << std::endl;
        Seq::CmdList::print(tstm, (uint8_t*)vec.data(), vec.size(), ttl_mask);
        auto text = tstm.get_buf();
        test_file_eq(path + ".txt", text);

        const_istream tistm(text);
        auto ttl_mask2 = Seq::CmdList::parse(vstm, tistm);
        assert(ttl_mask == ttl_mask2);
        assert(vec == vstm.get_buf());
    }
    catch (const SyntaxError &err) {
        string_ostream sstr;
        sstr << err;
        test_file_eq(path + ".err", sstr.get_buf());
    }
}

int main(int argc, char **argv)
{
    Log::printPID(false);
    if (argc != 2) {
        Log::error("ERROR: wrong number of arguments.\n");
        return 1;
    }

    std::string dir(argv[1]);
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

    return 0;
}
