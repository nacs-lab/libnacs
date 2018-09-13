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

#include "errors.h"

namespace NaCs {

namespace {

static std::string get_message(const std::string &what, const std::string &line,
                               int lineno, int colnum)
{
    std::string msg = "SyntaxError: " + what + '\n';
    if (lineno > 0) {
        msg.push_back('L');
        msg += std::to_string(lineno);
        if (colnum > 0) {
            msg.push_back(':');
            msg += std::to_string(colnum);
        }
        msg.push_back(' ');
    }
    msg += line;
    return msg;
}

}

SyntaxError::SyntaxError(std::string msg, std::string line, int lineno,
                         int colnum, int colstart, int colend)
    : std::runtime_error(get_message(msg, line, lineno, colnum)),
    m_msg(std::move(msg)),
    m_line(std::move(line)),
    m_lineno(lineno),
    m_colnum(colnum),
    m_colstart(colstart),
    m_colend(colend)
{
}

SyntaxError::~SyntaxError()
{
}

NACS_EXPORT() std::ostream &operator<<(std::ostream &stm, const SyntaxError &err)
{
    // This is not unicode aware but is good enough for now........
    stm << "SyntaxError: " << err.msg() << std::endl;
    auto lineno = err.lineno();
    int colstart, colend;
    int colnum = err.columns(&colstart, &colend);
    if (lineno > 0) {
        stm << "Line: " << lineno;
    }
    if (colnum > 0) {
        if (lineno > 0)
            stm << ", ";
        stm << "Column: " << colnum;
    }
    if (colnum > 0 || lineno > 0)
        stm << std::endl;
    const auto &line = err.line();
    stm << line << std::endl;
    if (colstart > 0 && colend > 0 && colstart <= colend &&
        colstart <= line.size() && colend <= line.size()) {
        for (int i = 1; i < colstart; i++)
            stm.put(' ');
        for (int i = colstart; i <= colend; i++)
            stm.put(i == colnum ? '^' : '~');
        stm << std::endl;
    }
    else if (colnum > 0 && colnum <= line.size()) {
        for (int i = 1; i < colnum; i++)
            stm.put(' ');
        stm << '^' << std::endl;
    }
    return stm;
}

}
