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

#ifndef _NACS_UTILS_ERRORS_H_
#define _NACS_UTILS_ERRORS_H_

#include "utils.h"

#include <ostream>
#include <system_error>

namespace NaCs {

// Helper function for C style system calls.
// (return -1 for error and type of error is in errno)
template<typename T, typename... Arg>
static inline T checkErrno(T res, Arg&&... arg)
{
    if (res == -1)
        throw std::system_error(errno, std::system_category(), std::forward<Arg>(arg)...);
    return res;
}

// Class that represent an error raised during parsing.
// It carries information about the reason and location of the error
// which includes the line number, the line of text and optionally column number / range.
// Both line and column numbers / ranges are 1-based.
class NACS_EXPORT(utils) SyntaxError : public std::runtime_error {
public:
    SyntaxError(std::string what, std::string line, int lineno,
                int colnum=-1, int colstart=-1, int colend=-1);
    ~SyntaxError() override;
    const std::string &msg() const
    {
        return m_msg;
    }
    const std::string &line() const
    {
        return m_line;
    }
    int lineno() const
    {
        return m_lineno;
    }
    int columns(int *colstart=nullptr, int *colend=nullptr) const
    {
        if (colstart)
            *colstart = m_colstart;
        if (colend)
            *colend = m_colend;
        return m_colnum;
    }

private:
    std::string m_msg;
    std::string m_line;
    int m_lineno;
    int m_colnum;
    int m_colstart;
    int m_colend;
};

// Pretty print the error
// The error section in the code will be highlighted
// based on the column info if availabe.
std::ostream &operator<<(std::ostream &stm, const SyntaxError &err);

}

#endif
