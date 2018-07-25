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

#include "test_helpers.h"

#include "../../lib/utils/container.h"

#include <assert.h>

#include <utility>

using namespace NaCs;

int main()
{
    int val = 0;
    // Make sure the increment class works.
    {
        Increment inc1(&val);
        assert(val == 0);
    }
    assert(val == 1);
    {
        Increment inc1;
    }
    assert(val == 1);
    {
        Increment inc1;
        {
            Increment inc2(&val);
            inc1 = std::move(inc2);
        }
        assert(val == 1);
    }
    assert(val == 2);
    {
        UnmovableIncrement inc1(&val);
        assert(val == 2);
    }
    assert(val == 3);
    {
        UnmovableIncrement inc1;
    }
    assert(val == 3);

    // Actually testing AnyPtr
    {
        AnyPtr ptr(new Increment(&val));
        assert(val == 3);
    }
    assert(val == 4);
    {
        // Should work for unmovable classes too
        AnyPtr ptr(new UnmovableIncrement(&val));
        assert(val == 4);
    }
    assert(val == 5);
    {
        // Only work for movable class
        AnyPtr ptr{Increment(&val)};
        assert(val == 5);
    }
    assert(val == 6);
    {
        // Arbitrary callback
        AnyPtr ptr{&val, [] (void *p) { ++*(int*)p; }};
        assert(val == 6);
    }
    assert(val == 7);
    return 0;
}
