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

#include "../../lib/utils/container.h"

#include <catch2/catch.hpp>

using namespace NaCs;

ANON_TEST_CASE() {
    FixedQueue<unsigned, 8> queue;

    REQUIRE(!queue.full());
    REQUIRE(queue.empty());
    for (unsigned i = 0; i < 8; i++) {
        REQUIRE(queue.size() == i);
        queue.push(i);
        REQUIRE(queue.size() == i + 1);
    }
    REQUIRE(queue.full());
    REQUIRE(!queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        REQUIRE(queue.size() == 8 - i);
        REQUIRE(queue.front() == i);
        auto v = queue.pop();
        REQUIRE(v == i);
        REQUIRE(queue.size() == 7 - i);
    }
    REQUIRE(!queue.full());
    REQUIRE(!queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        REQUIRE(queue.size() == 3 + i);
        queue.push(i);
        REQUIRE(queue.size() == 4 + i);
    }
    REQUIRE(queue.full());
    REQUIRE(!queue.empty());
    for (unsigned i = 5; i < 8; i++) {
        REQUIRE(queue.size() == 13 - i);
        REQUIRE(queue.front() == i);
        auto v = queue.pop();
        REQUIRE(v == i);
        REQUIRE(queue.size() == 12 - i);
    }
    REQUIRE(!queue.full());
    REQUIRE(!queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        REQUIRE(queue.size() == 5 - i);
        REQUIRE(queue.front() == i);
        auto v = queue.pop();
        REQUIRE(v == i);
        REQUIRE(queue.size() == 4 - i);
    }
    REQUIRE(!queue.full());
    REQUIRE(queue.empty());
}
