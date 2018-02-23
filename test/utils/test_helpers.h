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

#include <utility>

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
