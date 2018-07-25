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

#include "../../lib/utils/number.h"
#include "../../lib/utils/timer.h"
#include <assert.h>

#include <iostream>

using namespace NaCs;

template<typename T>
static inline void
test_single_get_bits(T i)
{
    unsigned bits = getBits(i);
    // std::cout << i << ", " << bits << std::endl;
    if (bits < sizeof(i) * 8) {
        assert(i < (T(1) << bits));
    }
    if (bits > 0) {
        assert(i >= (T(1) << (bits - 1)));
    }
}

template<typename T>
static inline void
test_get_bits()
{
    Timer timer;
    T i(-1);
    do {
        test_single_get_bits(i);
    } while (i-- > 0);
    timer.print();
}

int
main()
{
    test_get_bits<uint16_t>();
    test_get_bits<uint32_t>();
    return 0;
}
