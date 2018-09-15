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

    virtual char *extend(size_t sz) = 0;

protected:
    ssize_t m_end;
};

template<typename charT>
class NACS_EXPORT() basic_vector_streambuf : public buff_streambuf {
    static_assert(sizeof(charT) == 1);
public:
    basic_vector_streambuf(std::vector<charT> &&buf, size_t offset);
    basic_vector_streambuf(std::vector<charT> &&buf=std::vector<charT>())
        : basic_vector_streambuf(std::move(buf), buf.size())
    {}
    ~basic_vector_streambuf() override;

    std::vector<charT> get_buf();

private:
    char *extend(size_t sz) override;

    std::vector<charT> m_buf;
};

extern template class basic_vector_streambuf<char>;
extern template class basic_vector_streambuf<unsigned char>;

using vector_streambuf = basic_vector_streambuf<char>;
using uvector_streambuf = basic_vector_streambuf<unsigned char>;

class NACS_EXPORT() malloc_streambuf : public buff_streambuf {
public:
    malloc_streambuf();
    ~malloc_streambuf() override;

    char *get_buf(size_t &sz);

private:
    char *extend(size_t sz) override;

    std::unique_ptr<char, CDeleter> m_buf;
};

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

template<typename charT>
class NACS_EXPORT() basic_vector_ostream : public buff_ostream {
public:
    basic_vector_ostream(std::vector<charT> &&buf, size_t offset);
    basic_vector_ostream(std::vector<charT> &&buf=std::vector<charT>())
        : basic_vector_ostream(std::move(buf), buf.size())
    {}
    ~basic_vector_ostream();

    std::vector<charT> get_buf()
    {
        flush();
        return m_buf.get_buf();
    }

private:
    basic_vector_streambuf<charT> m_buf;
};

extern template class basic_vector_ostream<char>;
extern template class basic_vector_ostream<unsigned char>;

using vector_ostream = basic_vector_ostream<char>;
using uvector_ostream = basic_vector_ostream<unsigned char>;

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