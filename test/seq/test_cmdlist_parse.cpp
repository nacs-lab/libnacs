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

#include <iostream>
#include <fstream>
#include <sstream>

using namespace NaCs;

void test(std::string name)
{
    std::ifstream istm(name);
    assert(istm.good());
    vector_ostream vstm;
    try {
        uint32_t ttl_mask = Seq::CmdList::parse(vstm, istm);
        auto vec = vstm.get_buf();
        uint64_t len_ns = Seq::CmdList::total_time((uint8_t*)vec.data(), vec.size()) * 10;
        std::ifstream bstm(name + ".bin");
        assert(bstm.good());
        std::string str(std::istreambuf_iterator<char>(bstm), {});
        auto str_data = (const uint8_t*)str.data();
        auto str_sz = str.size();
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
        assert(str_sz == vec.size());
        assert(memcmp(str_data, vec.data(), str_sz));
    }
    catch (const SyntaxError &err) {
        string_ostream sstr;
        sstr << err;
        std::ifstream estm(name + ".err");
        std::string expected(std::istreambuf_iterator<char>(estm), {});
        assert(sstr.get_buf() == expected);
    }
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "ERROR: wrong number of arguments." << std::endl;
        return 1;
    }

    std::string dir(argv[1]);
    dir += "/cmdlists/";

    test(dir + "ttl_mask_err1");
    test(dir + "ttl_mask_err2");
    test(dir + "ttl_mask_err3");

    test(dir + "invalid_cmd");
    test(dir + "invalid_cmd2");

    test(dir + "invalid_ttl");
    test(dir + "invalid_ttl1_space");
    test(dir + "invalid_ttlall_space");

    test(dir + "invalid_ttl1_1");
    test(dir + "invalid_ttl1_2");
    test(dir + "invalid_ttl1_3");
    test(dir + "invalid_ttl1_4");

    test(dir + "invalid_ttl_t1");
    test(dir + "invalid_ttl_t2");
    test(dir + "invalid_ttl_t3");

    test(dir + "invalid_wait");
    test(dir + "invalid_wait2");
    test(dir + "invalid_wait3");

    return 0;
}
