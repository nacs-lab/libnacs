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

#ifndef __NACS_UTILS_STREAMS_H__
#define __NACS_UTILS_STREAMS_H__

#include "utils.h"

#include <istream>
#include <limits>
#include <memory>
#include <ostream>
#include <streambuf>
#include <vector>

namespace NaCs {

// Similar to `std::stringstream` but is also indexable.
// The subtypes also support multiple different back storages
// and can generally transfer the ownership of the buffer
// which is not possible with `std::stringstream`.
class NACS_EXPORT() buff_streambuf : public std::streambuf {
public:
    buff_streambuf(size_t sz=0)
        : m_end(sz)
    {}
    pos_type tellg() const
    {
        return pptr() - pbase();
    }
    const char &operator[](size_t i) const
    {
        return pbase()[i];
    }
    char &operator[](size_t i)
    {
        return pbase()[i];
    }

private:
    std::streamsize xsputn(const char* s, std::streamsize count) override;
    int_type overflow(int_type ch) override;
    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                     std::ios_base::openmode which) override;
    pos_type seekpos(pos_type pos, std::ios_base::openmode which) override;

    pos_type _seekpos(pos_type pos);
    void update_size();

    // The base class defines most of the interface with the stream framework
    // and subclasses only need to define the `extend` method, which should
    // resize the buffer to fit at least `sz` bytes from the current pointer
    // without loosing any existing content (up to `m_end`).
    virtual char *extend(size_t sz) = 0;

protected:
    // This is the last location accessed on the stream.
    // In another word, this is the length of the file.
    ssize_t m_end;
};

template<typename T>
class NACS_EXPORT() basic_vector_streambuf : public buff_streambuf {
public:
    basic_vector_streambuf(T &&buf, size_t offset);
    basic_vector_streambuf(T &&buf=T())
        : basic_vector_streambuf(std::move(buf), buf.size())
    {}
    ~basic_vector_streambuf() override;

    T get_buf();

private:
    char *extend(size_t sz) override;

    T m_buf;
};

template<typename T>
basic_vector_streambuf<T>::basic_vector_streambuf(T &&buf, size_t offset)
    : buff_streambuf(offset),
      m_buf(std::move(buf))
{
    setp((char*)&m_buf[0], (char*)&m_buf[m_buf.size()]);
    pbump(offset);
}

template<typename T>
basic_vector_streambuf<T>::~basic_vector_streambuf()
{
}

template<typename T>
T basic_vector_streambuf<T>::get_buf()
{
    m_buf.resize(m_end);
    auto res = std::move(m_buf);
    setp((char*)&m_buf[0], (char*)&m_buf[m_buf.size()]);
    m_end = m_buf.size();
    return res;
}

template<typename T>
char *basic_vector_streambuf<T>::extend(size_t sz)
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

extern template class basic_vector_streambuf<std::vector<char>>;
extern template class basic_vector_streambuf<std::vector<unsigned char>>;
extern template class basic_vector_streambuf<std::string>;

using vector_streambuf = basic_vector_streambuf<std::vector<char>>;
using uvector_streambuf = basic_vector_streambuf<std::vector<unsigned char>>;
using string_streambuf = basic_vector_streambuf<std::string>;

class NACS_EXPORT() malloc_streambuf : public buff_streambuf {
public:
    malloc_streambuf();
    ~malloc_streambuf() override;

    char *get_buf(size_t &sz);

private:
    char *extend(size_t sz) override;

    std::unique_ptr<char, CDeleter> m_buf;
};

// This is a streambuf for istream from constant memory.
// Based on http://www.voidcn.com/article/p-vjnlygmc-gy.html
class NACS_EXPORT() const_streambuf : public std::streambuf {
public:
    const_streambuf(const void *begin, const void *end);
    ~const_streambuf() override;

private:
    int_type underflow() override;
    int_type uflow() override;
    int_type pbackfail(int_type ch) override;
    std::streamsize showmanyc() override;

    const char *const m_begin;
    const char *const m_end;
    const char *m_current;
};

// Here we have the (i|o)stream s that uses/wraps the streambufs above
class buff_ostream : public std::ostream {
public:
    buff_ostream(buff_streambuf *buf)
        : std::ostream(buf)
    {}
    pos_type tellg()
    {
        return static_cast<buff_streambuf*>(rdbuf())->tellg();
    }
    const char &operator[](size_t i) const
    {
        return (*static_cast<const buff_streambuf*>(rdbuf()))[i];
    }
    char &operator[](size_t i)
    {
        return (*static_cast<buff_streambuf*>(rdbuf()))[i];
    }
};

template<typename T>
class NACS_EXPORT() basic_vector_ostream : public buff_ostream {
public:
    basic_vector_ostream(T &&buf, size_t offset);
    basic_vector_ostream(T &&buf=T())
        : basic_vector_ostream(std::move(buf), buf.size())
    {}
    ~basic_vector_ostream();

    T get_buf()
    {
        flush();
        return m_buf.get_buf();
    }

private:
    basic_vector_streambuf<T> m_buf;
};

template<typename T>
basic_vector_ostream<T>::basic_vector_ostream(T &&buf, size_t offset)
    : buff_ostream(&m_buf), m_buf(std::move(buf), offset)
{
}

template<typename T>
basic_vector_ostream<T>::~basic_vector_ostream()
{
}

extern template class basic_vector_ostream<std::vector<char>>;
extern template class basic_vector_ostream<std::vector<unsigned char>>;
extern template class basic_vector_ostream<std::string>;

using vector_ostream = basic_vector_ostream<std::vector<char>>;
using uvector_ostream = basic_vector_ostream<std::vector<unsigned char>>;
using string_ostream = basic_vector_ostream<std::string>;

class NACS_EXPORT() malloc_ostream : public buff_ostream {
public:
    malloc_ostream();
    ~malloc_ostream();

    char *get_buf(size_t &sz)
    {
        flush();
        return m_buf.get_buf(sz);
    }

private:
    malloc_streambuf m_buf;
};

class NACS_EXPORT() const_istream : public std::istream {
public:
    const_istream(const void *begin, const void *end);
    const_istream(const std::string &str)
        : const_istream(&str[0], &str[0] + str.size())
    {}
    template<typename T>
    const_istream(const std::vector<T> &vec)
        : const_istream(&vec[0], &vec[0] + vec.size())
    {}
    ~const_istream() override;

private:
    const_streambuf m_buf;
};

constexpr auto eofc = std::istream::traits_type::eof();

static inline std::istream &ignore_line(std::istream &stm, char c)
{
    stm.ignore(std::numeric_limits<std::streamsize>::max(), c);
    return stm;
}

static inline std::istream &ignore_line(std::istream &stm)
{
    // Use a separate signature so that `stm >> ignore_line` could work.
    return ignore_line(stm, '\n');
}

static inline std::istream &ignore_space(std::istream &stm)
{
    while (true) {
        auto c = stm.peek();
        if (c == eofc)
            break;
        if (c != ' ' && c != '\t')
            break;
        stm.get();
    }
    return stm;
}

}

#endif
