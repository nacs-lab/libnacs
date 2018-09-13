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

#include "streams.h"

#include <utility>

namespace NaCs {

vector_streambuf::vector_streambuf(std::vector<char> &&buf, size_t offset)
    : m_buf(std::move(buf))
{
    setp(&m_buf[0], &m_buf[m_buf.size()]);
    pbump(offset);
}

vector_streambuf::~vector_streambuf()
{
}

std::streamsize vector_streambuf::xsputn(const char *s, std::streamsize count)
{
    auto p = extend(count);
    memcpy(p, s, count);
    pbump(count);
    return count;
}

auto vector_streambuf::overflow(int_type ch) -> int_type
{
    // The document on this function is really confusing.
    // In terms of what it should do when `ch` is `eof` and
    // what to do about the current pointer
    //
    // [cppreference.com](https://en.cppreference.com/w/cpp/io/basic_streambuf/overflow)
    // says:
    // > Ensures that there is space at the put area for at least one character
    // which seems to imply that the reservation needs to be done even when `ch` is `eof`.
    //
    // [cplusplus.com](http://www.cplusplus.com/reference/streambuf/streambuf/overflow/)
    // says:
    // > put a character into the controlled output sequence
    // > without changing the current position.
    // Mentioning nothing about what needs to be done with `eof` as input and
    // seems to suggest that the current position should not be updated
    // after a character is added.
    //
    // From [a stackoverflow question](https://stackoverflow.com/questions/19849143/what-is-wrong-with-my-implementation-of-overflow)
    // and `libc++`s implementation it seems that `eof` input should just be ignored and
    // the current pointer should be updated after the character is written
    // so that's what we'll do....
    if (traits_type::eq_int_type(ch, traits_type::eof()))
        return traits_type::not_eof(ch);
    *extend(1) = (char)ch;
    pbump(1);
    return traits_type::not_eof(ch);
}

std::vector<char> vector_streambuf::get_buf()
{
    m_buf.resize(pptr() - pbase());
    auto res = std::move(m_buf);
    auto p = &m_buf[m_buf.size()];
    setp(p, p);
    return res;
}

char *vector_streambuf::extend(size_t sz)
{
    auto oldbase = pbase();
    auto oldptr = pptr();
    auto oldsz = oldptr - oldbase;
    // overallocate.
    m_buf.resize((oldsz + sz) * 3 / 2);
    setp(&m_buf[0], &m_buf[m_buf.size()]);
    pbump(oldsz);
    return &m_buf[oldsz];
}

malloc_streambuf::malloc_streambuf()
{
    setp(nullptr, nullptr);
}

malloc_streambuf::~malloc_streambuf()
{
}

std::streamsize malloc_streambuf::xsputn(const char *s, std::streamsize count)
{
    auto p = extend(count);
    memcpy(p, s, count);
    pbump(count);
    return count;
}

auto malloc_streambuf::overflow(int_type ch) -> int_type
{
    // See also `vector_streambuf::overflow`
    if (traits_type::eq_int_type(ch, traits_type::eof()))
        return traits_type::not_eof(ch);
    *extend(1) = (char)ch;
    pbump(1);
    return traits_type::not_eof(ch);
}

char *malloc_streambuf::get_buf(size_t &sz)
{
    sz = pptr() - pbase();
    auto buf = m_buf.release();
    setp(nullptr, nullptr);
    return buf;
}

char *malloc_streambuf::extend(size_t sz)
{
    auto oldbase = pbase();
    auto oldptr = pptr();
    auto oldsz = oldptr - oldbase;
    // overallocate.
    auto new_sz = (oldsz + sz) * 3 / 2;
    auto buf = (char*)realloc(m_buf.release(), new_sz);
    m_buf.reset(buf);
    setp(buf, &buf[new_sz]);
    pbump(oldsz);
    return &buf[oldsz];
}

vector_ostream::vector_ostream(std::vector<char> &&buf, size_t offset)
    : std::ostream(&m_buf), m_buf(std::move(buf), offset)
{
}

vector_ostream::~vector_ostream()
{
}

malloc_ostream::malloc_ostream()
    : std::ostream(&m_buf), m_buf()
{
}

malloc_ostream::~malloc_ostream()
{
}

}
