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

std::streamsize buff_streambuf::xsputn(const char *s, std::streamsize count)
{
    auto p = extend(count);
    memcpy(p, s, count);
    pbump(count);
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

inline void buff_streambuf::update_size()
{
    auto sz = pptr() - pbase();
    if (sz > m_end) {
        m_end = sz;
    }
}

template<typename charT>
basic_vector_streambuf<charT>::basic_vector_streambuf(std::vector<charT> &&buf, size_t offset)
    : buff_streambuf(offset),
      m_buf(std::move(buf))
{
    setp((char*)&m_buf[0], (char*)&m_buf[m_buf.size()]);
    pbump(offset);
}

template<typename charT>
basic_vector_streambuf<charT>::~basic_vector_streambuf()
{
}

template<typename charT>
std::vector<charT> basic_vector_streambuf<charT>::get_buf()
{
    m_buf.resize(pptr() - pbase());
    auto res = std::move(m_buf);
    setp((char*)&m_buf[0], (char*)&m_buf[m_buf.size()]);
    m_end = m_buf.size();
    return res;
}

template<typename charT>
char *basic_vector_streambuf<charT>::extend(size_t sz)
{
    auto oldbase = pbase();
    auto oldptr = pptr();
    auto oldsz = oldptr - oldbase;
    auto newsz = (oldsz + sz) * 3 / 2;
    if (newsz <= m_buf.size())
        return oldptr;
    // overallocate.
    m_buf.resize(newsz);
    setp((char*)&m_buf[0], (char*)&m_buf[m_buf.size()]);
    pbump(oldsz);
    return (char*)&m_buf[oldsz];
}

template class basic_vector_streambuf<char>;
template class basic_vector_streambuf<unsigned char>;

malloc_streambuf::malloc_streambuf()
{
    setp(nullptr, nullptr);
}

malloc_streambuf::~malloc_streambuf()
{
}

char *malloc_streambuf::get_buf(size_t &sz)
{
    sz = pptr() - pbase();
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
    pbump(oldsz);
    return &buf[oldsz];
}

template<typename charT>
basic_vector_ostream<charT>::basic_vector_ostream(std::vector<charT> &&buf, size_t offset)
    : buff_ostream(&m_buf), m_buf(std::move(buf), offset)
{
}

template<typename charT>
basic_vector_ostream<charT>::~basic_vector_ostream()
{
}

template class basic_vector_ostream<char>;
template class basic_vector_ostream<unsigned char>;

malloc_ostream::malloc_ostream()
    : buff_ostream(&m_buf), m_buf()
{
}

malloc_ostream::~malloc_ostream()
{
}

}
