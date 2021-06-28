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

#include "../../lib/nacs-utils/timer.h"
#include "../../lib/nacs-utils/wavemeter.h"

#include <iostream>
#include <fstream>
#include <limits>

#include <assert.h>

using namespace NaCs;

int main(int argc, char **argv)
{
    assert(argc >= 2);

    std::ifstream istm(argv[1]);

    Wavemeter parser(0, std::numeric_limits<double>::max());
    size_t sz;

    Timer timer;
    parser.parse(istm, &sz, 0, std::numeric_limits<double>::max());
    timer.print(true);

    const double *times_ptr;
    const double *datas_ptr;
    const double *heights_ptr;
    std::tie(times_ptr, datas_ptr, heights_ptr) =
        parser.parse(istm, &sz, 0, std::numeric_limits<double>::max());
    timer.print(true);

    parser.parse(istm, &sz, 0, std::numeric_limits<double>::max());
    timer.print(true);

    std::cout << (void*)times_ptr << " "
              << (void*)datas_ptr << std::endl;

    return 0;
}
