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

#include <assert.h>

namespace NaCs {

std::streamsize buff_streambuf::xsputn(const char *s, std::streamsize count)
{
    auto p = extend(count);
    memcpy(p, s, count);
    pbump((int)count);
    update_size();
    return count;
}

auto buff_streambuf::overflow(int_type ch) -> int_type
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
    update_size();
    return traits_type::not_eof(ch);
}

auto buff_streambuf::seekpos(pos_type pos, std::ios_base::openmode which) -> pos_type
{
    if (which != std::ios_base::out)
        return pos_type(-1);
    return _seekpos(pos);
}

inline auto buff_streambuf::_seekpos(pos_type pos) -> pos_type
{
    if (pos < 0)
        return pos_type(-1);
    // Update before changing the pointer as well as after changing the pointer
    // so that we can catch seeking back.
    update_size();
    auto base = pbase();
    auto end = epptr();
    if (unlikely(base + (off_type)pos > end)) {
        extend(base + (off_type)pos - end);
        base = pbase();
    }
    pbump(int((off_type)pos - (pptr() - base)));
    update_size();
    return pos;
}

auto buff_streambuf::seekoff(off_type off, std::ios_base::seekdir dir,
                             std::ios_base::openmode which) -> pos_type
{
    if (which != std::ios_base::out)
        return pos_type(-1);
    pos_type pos;
    if (dir == std::ios_base::beg) {
        pos = off;
    }
    else if (dir == std::ios_base::cur) {
        pos = pptr() - pbase();
        if (off == 0)
            return pos;
        pos = pos + off;
    }
    else if (dir == std::ios_base::end) {
        if (off > 0)
            return pos_type(-1);
        pos = m_end + off;
    }
    else {
        return pos_type(-1);
    }
    return _seekpos(pos);
}

inline int buff_streambuf::sync()
{
    update_size();
    return 0;
}

inline void buff_streambuf::update_size()
{
    auto sz = pptr() - pbase();
    if (sz > m_end) {
        m_end = sz;
    }
}

#if !NACS_OS_WINDOWS
template class basic_vector_streambuf<std::vector<char>>;
template class basic_vector_streambuf<std::vector<unsigned char>>;
template class basic_vector_streambuf<std::string>;
#endif

malloc_streambuf::malloc_streambuf()
{
    setp(nullptr, nullptr);
}

malloc_streambuf::~malloc_streambuf()
{
}

char *malloc_streambuf::get_buf(size_t &sz)
{
    sz = m_end;
    auto buf = m_buf.release();
    setp(nullptr, nullptr);
    m_end = 0;
    return buf;
}

char *malloc_streambuf::extend(size_t sz)
{
    auto oldbase = pbase();
    auto oldptr = pptr();
    auto oldsz = oldptr - oldbase;
    // overallocate.
    auto new_sz = (oldsz + sz) * 3 / 2;
    if (oldbase + new_sz <= epptr())
        return &m_buf.get()[oldsz];
    auto buf = (char*)realloc(m_buf.release(), new_sz);
    m_buf.reset(buf);
    setp(buf, &buf[new_sz]);
    pbump((int)oldsz);
    return &buf[oldsz];
}


const_streambuf::const_streambuf(const void *begin, const void *end)
    : m_begin((const char*)begin),
      m_end((const char*)end),
      m_current((const char*)begin)
{
    assert(m_begin <= m_end);
}

const_streambuf::~const_streambuf()
{
}

auto const_streambuf::underflow() -> int_type
{
    if (m_current == m_end)
        return traits_type::eof();

    return traits_type::to_int_type(*m_current);
}

auto const_streambuf::uflow() -> int_type
{
    if (m_current == m_end)
        return traits_type::eof();

    return traits_type::to_int_type(*m_current++);
}

auto const_streambuf::pbackfail(int_type ch) -> int_type
{
    if (m_current == m_begin || (ch != traits_type::eof() && ch != m_current[-1]))
        return traits_type::eof();

    return traits_type::to_int_type(*--m_current);
}

std::streamsize const_streambuf::showmanyc()
{
    assert(m_current <= m_end);
    return m_end - m_current;
}

auto const_streambuf::seekpos(pos_type pos, std::ios_base::openmode which) -> pos_type
{
    if (which != std::ios_base::in)
        return pos_type(-1);
    if (m_begin + std::streamoff(pos) > m_end)
        return pos_type(-1);
    m_current = m_begin + std::streamoff(pos);
    return pos;
}

auto const_streambuf::seekoff(off_type off, std::ios_base::seekdir dir,
                              std::ios_base::openmode which) -> pos_type
{
    if (which != std::ios_base::in)
        return pos_type(-1);
    const char *ptr;;
    if (dir == std::ios_base::beg) {
        ptr = m_begin;
    }
    else if (dir == std::ios_base::cur) {
        ptr = m_current;
    }
    else if (dir == std::ios_base::end) {
        ptr = m_end;
    }
    else {
        return pos_type(-1);
    }
    ptr += off;
    if (ptr < m_begin || ptr > m_end)
        return pos_type(-1);
    m_current = ptr;
    return ptr - m_begin;
}

#if !NACS_OS_WINDOWS
template class basic_vector_ostream<std::vector<char>>;
template class basic_vector_ostream<std::vector<unsigned char>>;
template class basic_vector_ostream<std::string>;
#endif

malloc_ostream::malloc_ostream()
    : buff_ostream(&m_buf), m_buf()
{
}

malloc_ostream::~malloc_ostream()
{
}

const_istream::const_istream(const void *begin, const void *end)
    : std::istream(&m_buf),
    m_buf(begin, end)
{}

const_istream::~const_istream()
{}

}
