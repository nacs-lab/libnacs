/*************************************************************************
 *   Copyright (c) 2015 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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

#include <assert.h>

using namespace NaCs;

void test_print(std::ostream &stm)
{
    stm << 12345 << " abcde " << 1.23;
}

int main()
{
    std::stringstream sstm;
    test_print(sstm);
    auto const res = sstm.str();

    malloc_ostream mstm;
    test_print(mstm);
    size_t sz;
    auto p = mstm.get_buf(sz);
    assert(sz == res.size());
    assert(memcmp(p, &res[0], res.size()) == 0);
    free(p);

    vector_ostream vstm;
    test_print(vstm);
    auto v = vstm.get_buf();
    assert(v.size() == res.size());
    assert(memcmp(&v[0], &res[0], res.size()) == 0);

    return 0;
}
