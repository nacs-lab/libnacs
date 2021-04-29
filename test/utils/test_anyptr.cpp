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

#include "test_helpers.h"

#include "../../lib/utils/container.h"

#include <catch2/catch.hpp>

#include <utility>

using namespace NaCs;

TEST_CASE("AnyPtr") {
    int val = 0;
    // Make sure the increment class works.
    {
        Increment inc1(&val);
        REQUIRE(val == 0);
    }
    REQUIRE(val == 1);
    {
        Increment inc1;
    }
    REQUIRE(val == 1);
    {
        Increment inc1;
        {
            Increment inc2(&val);
            inc1 = std::move(inc2);
        }
        REQUIRE(val == 1);
    }
    REQUIRE(val == 2);
    {
        UnmovableIncrement inc1(&val);
        REQUIRE(val == 2);
    }
    REQUIRE(val == 3);
    {
        UnmovableIncrement inc1;
    }
    REQUIRE(val == 3);

    // Actually testing AnyPtr
    {
        AnyPtr ptr(new Increment(&val));
        REQUIRE(val == 3);
    }
    REQUIRE(val == 4);
    {
        // Should work for unmovable classes too
        AnyPtr ptr(new UnmovableIncrement(&val));
        REQUIRE(val == 4);
    }
    REQUIRE(val == 5);
    {
        // Only work for movable class
        AnyPtr ptr{Increment(&val)};
        REQUIRE(val == 5);
    }
    REQUIRE(val == 6);
    {
        // Arbitrary callback
        AnyPtr ptr{&val, [] (void *p) { ++*(int*)p; }};
        REQUIRE(val == 6);
    }
    REQUIRE(val == 7);
    {
        // Make sure we are default constructable.
        AnyPtr ptr;
    }
}
