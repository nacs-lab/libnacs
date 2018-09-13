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

class NACS_EXPORT() vector_streambuf : public buff_streambuf {
public:
    vector_streambuf(std::vector<char> &&buf, size_t offset);
    vector_streambuf(std::vector<char> &&buf=std::vector<char>())
        : vector_streambuf(std::move(buf), buf.size())
    {}
    ~vector_streambuf() override;

    std::vector<char> get_buf();

private:
    char *extend(size_t sz) override;

    std::vector<char> m_buf;
};

class NACS_EXPORT() malloc_streambuf : public buff_streambuf {
public:
    malloc_streambuf();
    ~malloc_streambuf() override;

    char *get_buf(size_t &sz);

private:
    char *extend(size_t sz) override;

    std::unique_ptr<char, CDeleter> m_buf;
};

class NACS_EXPORT() vector_ostream : public std::ostream {
public:
    vector_ostream(std::vector<char> &&buf, size_t offset);
    vector_ostream(std::vector<char> &&buf=std::vector<char>())
        : vector_ostream(std::move(buf), buf.size())
    {}
    ~vector_ostream();

    std::vector<char> get_buf()
    {
        flush();
        return m_buf.get_buf();
    }

private:
    vector_streambuf m_buf;
};

class NACS_EXPORT() malloc_ostream : public std::ostream {
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

}

#endif
