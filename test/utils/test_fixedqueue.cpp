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

#include "../../lib/utils/container.h"

#include <assert.h>

using namespace NaCs;

int main()
{
    FixedQueue<unsigned, 8> queue;

    assert(!queue.full() && queue.empty());
    for (unsigned i = 0; i < 8; i++) {
        assert(queue.size() == i);
        queue.push(i);
        assert(queue.size() == i + 1);
    }
    assert(queue.full() && !queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        assert(queue.size() == 8 - i);
        assert(queue.front() == i);
        auto v = queue.pop();
        assert(v == i);
        assert(queue.size() == 7 - i);
    }
    assert(!queue.full() && !queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        assert(queue.size() == 3 + i);
        queue.push(i);
        assert(queue.size() == 4 + i);
    }
    assert(queue.full() && !queue.empty());
    for (unsigned i = 5; i < 8; i++) {
        assert(queue.size() == 13 - i);
        assert(queue.front() == i);
        auto v = queue.pop();
        assert(v == i);
        assert(queue.size() == 12 - i);
    }
    assert(!queue.full() && !queue.empty());
    for (unsigned i = 0; i < 5; i++) {
        assert(queue.size() == 5 - i);
        assert(queue.front() == i);
        auto v = queue.pop();
        assert(v == i);
        assert(queue.size() == 4 - i);
    }
    assert(!queue.full() && queue.empty());

    return 0;
}
