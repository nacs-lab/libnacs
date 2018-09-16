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

#include "../../lib/utils/streams.h"

#include <string>
#include <sstream>

#include <iostream>

#include <assert.h>

using namespace NaCs;

void test_print(std::ostream &stm)
{
    stm << 12345 << " abcde " << 1.23;
}

void test_seek(std::ostream &stm)
{
    stm << 12345 << " abcde " << 1.23 << "asdfjkl;";
    auto n0 = stm.tellp();
    stm << 12345 << " abcde " << 1.23 << "asdfjkl;";
    stm << 12345 << " abcde " << 1.23 << "asdfjkl;";
    auto n = stm.tellp();
    stm.seekp(-5, std::ios_base::cur);
    stm << 1;
    stm.seekp(-30, std::ios_base::end);
    stm << 'o';
    stm.seekp(n0);
    stm << ',';
    stm.seekp(n);
}

template<typename F>
void test_streams(F &&f)
{
    std::stringstream sstm;
    f(sstm);
    auto const res = sstm.str();

    for (int i = 0; i < 2; i++) {
        if (i != 0) {
            std::stringstream sstm;
            f(sstm);
            sstm.seekp(0);
            assert(res == sstm.str());
        }

        malloc_ostream mstm;
        f(mstm);
        if (i != 0)
            mstm.seekp(0);
        size_t sz;
        auto p = mstm.get_buf(sz);
        assert(sz == res.size());
        assert(memcmp(p, &res[0], res.size()) == 0);
        free(p);

        vector_ostream vstm;
        f(vstm);
        if (i != 0)
            vstm.seekp(0);
        auto v = vstm.get_buf();
        assert(v.size() == res.size());
        assert(memcmp(&v[0], &res[0], res.size()) == 0);

        string_ostream sstm2;
        f(sstm2);
        if (i != 0)
            sstm2.seekp(0);
        auto s = sstm2.get_buf();
        assert(s == res);
    }
}

int main()
{
    test_streams(test_print);
    test_streams(test_seek);
    return 0;
}
