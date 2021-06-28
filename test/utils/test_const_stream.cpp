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

#define CATCH_CONFIG_MAIN

#include "../../lib/nacs-utils/streams.h"

#include <catch2/catch.hpp>

using namespace NaCs;

const char str[] = "123 3.4 0x123";

TEST_CASE("const_istream") {
    int a;
    double b;
    unsigned c;
    const_istream istm(str, str + sizeof(str));
    istm >> a >> b >> std::hex >> c;
    REQUIRE(a == 123);
    REQUIRE(b == 3.4);
    REQUIRE(c == 0x123);
    REQUIRE(istm.good());
}
