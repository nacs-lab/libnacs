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

#define CATCH_CONFIG_MAIN

#include "../../lib/utils/number.h"
#include "../../lib/utils/timer.h"

#include <catch2/catch.hpp>

#include <iostream>

using namespace NaCs;

template<typename T>
static inline void
test_single_get_bits(T i)
{
    unsigned bits = getBits(i);
    if (bits < sizeof(i) * 8) {
        REQUIRE(i < (T(1) << bits));
    }
    if (bits > 0) {
        REQUIRE(i >= (T(1) << (bits - 1)));
    }
}

template<typename T>
static inline void
test_get_bits()
{
    T i(-1);
    do {
        test_single_get_bits(i);
    } while (i-- > 0);
}

TEST_CASE("16 bits") {
    uint16_t i(-1);
    do {
        test_single_get_bits(i);
    } while (i-- > 0);
}

TEST_CASE("32 bits") {
    uint32_t i(-1);
    do {
        test_single_get_bits(i);
        if (i < 17)
            break;
        i -= 17;
    } while (true);
    i = 100;
    do {
        test_single_get_bits(i);
    } while (i-- > 0);
}
