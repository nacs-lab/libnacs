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

#ifdef NDEBUG
#  undef NDEBUG
#endif

#include "../../lib/utils/timer.h"
#include "../../lib/seq/misc.h"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace NaCs;

int main(int argc, char **argv)
{
    assert(argc >= 2);

    std::ifstream istm(argv[1]);

    WavemeterParser parser;
    size_t sz;

    Timer timer;
    parser.parse(istm, &sz, true);
    timer.print(true);

    auto ptrs = parser.parse(istm, &sz, true);
    timer.print(true);

    parser.parse(istm, &sz, true);
    timer.print(true);

    std::cout << (void*)ptrs.first << " "
              << (void*)ptrs.second << std::endl;

    return 0;
}
