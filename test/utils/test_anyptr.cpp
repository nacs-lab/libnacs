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

#include <nacs-utils/container.h>

#include <assert.h>

#include <utility>

using namespace NaCs;

struct Increment {
    Increment(int *ptr, int inc=1)
        : m_ptr(ptr),
          m_inc(inc)
    {
    }
    Increment()
    {}
    Increment &operator=(Increment &&other)
    {
        std::swap(m_ptr, other.m_ptr);
        std::swap(m_inc, other.m_inc);
        return *this;
    }
    Increment &operator=(const Increment &other) = delete;
    Increment(Increment &&other)
    {
        *this = std::move(other);
    }
    Increment(const Increment &other) = delete;
    ~Increment()
    {
        if (m_ptr) {
            *m_ptr += m_inc;
        }
    }

private:
    int *m_ptr = nullptr;
    int m_inc = 0;
};

struct UnmovableIncrement : Increment {
    using Increment::Increment;
    UnmovableIncrement(UnmovableIncrement &&other) = delete;
};

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
