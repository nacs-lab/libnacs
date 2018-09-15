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

#include "parser_p.h"
#include "../utils/errors.h"

namespace NaCs {
namespace Seq {

ParserBase::ParserBase(std::istream &istm)
    : istm(istm)
{
    next_line();
}

bool ParserBase::next_line()
{
    std::getline(istm, line);
    lineno++;
    colno = 0;
    return istm.good();
}

void ParserBase::syntax_error(std::string msg, int colnum, int colstart, int colend)
{
    throw SyntaxError(std::move(msg), line, lineno, colnum, colstart, colend);
}

void ParserBase::skip_whitespace()
{
    int linelen = (int)line.size();
    for (; colno < linelen; colno++) {
        if (!std::isspace(line[colno])) {
            return;
        }
    }
}

bool ParserBase::skip_comments()
{
    while (true) {
        skip_whitespace();
        if (colno < (int)line.size() && line[colno] != '#')
            return true;
        if (!next_line()) {
            return false;
        }
    }
}

std::pair<const std::string&,int> ParserBase::read_name(bool allow_num0)
{
    buff.clear();
    skip_whitespace();
    int startcol = colno;
    int linelen = (int)line.size();
    if (colno >= linelen)
        return {buff, -1};
    auto isalphanum = [] (char c, bool num=true) {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_') {
            return true;
        }
        else if (!num) {
            return false;
        }
        return c >= '0' && c <= '9';
    };

    auto c0 = line[colno];
    if (!isalphanum(c0, allow_num0))
        return {buff, -1};
    buff.push_back(c0);
    colno++;
    for (; colno < linelen; colno++) {
        auto c = line[colno];
        if (isalphanum(c)) {
            buff.push_back(c);
        }
        else {
            break;
        }
    }
    return {buff, startcol};
}

std::pair<uint64_t,int> ParserBase::read_hex(uint64_t lo, uint64_t hi)
{
    skip_whitespace();
    int startcol = colno;
    if (peek() != '0' || std::tolower(peek(1)) != 'x')
        return {0, -1};
    auto startptr = &line[startcol];
    char *endptr;
    auto res = strtoull(startptr, &endptr, 16);
    if (endptr <= startptr + 2)
        syntax_error("Empty or invalid hex literal", startcol + 3,
                     startcol + 1, startcol + 3);
    if (res < lo || res > hi || errno == ERANGE)
        syntax_error("Hex literal out of range [" + to_hexstring(lo) + ", " +
                     to_hexstring(hi) + "]", -1,
                     startcol + 1, startcol + int(endptr - startptr));
    colno = startcol + (endptr - startptr);
    return {res, startcol};
}

template std::pair<int,int> ParserBase::read_dec<int>(int lo, int hi);

std::pair<double,int> ParserBase::read_float(double lo, double hi)
{
    skip_whitespace();
    int startcol = colno;
    auto startptr = &line[startcol];
    int hassign = peek() == '+' || peek() == '-';
    // Disallow parsing of hex floating point numbers.
    // This is not a nice way to do it but I couldn't find another function
    // that's easy to use...
    if (peek(hassign) == '0' && std::tolower(peek(hassign + 1)) == 'x')
        return {0, -1};
    char *endptr;
    auto res = strtod(startptr, &endptr);
    if (endptr == startptr)
        return {0, -1};
    if (errno == ERANGE)
        syntax_error("Floating point literal overflow", -1,
                     startcol + 1, startcol + int(endptr - startptr));
    if (res < lo || res > hi)
        syntax_error("Floating point literal out of range [" + std::to_string(lo) + ", " +
                     std::to_string(hi) + "]", -1,
                     startcol + 1, startcol + int(endptr - startptr));
    colno = startcol + (endptr - startptr);
    return {res, startcol};
}

bool ParserBase::checked_next_line()
{
    skip_whitespace();
    auto c = peek();
    if (c && c != '#')
        syntax_error("Unexpected charater at the end of line", colno + 1);
    return next_line();
}

uint8_t ParserBase::read_ddschn(const char *name)
{
    if (peek() != '(')
        syntax_error(std::string("Invalid ") + name + " command: expecting `(`", colno + 1);
    colno++;
    uint8_t chn;
    int chn_start;
    std::tie(chn, chn_start) = read_dec(0, 21);
    if (chn_start == -1)
        syntax_error("Missing DDS channel", colno + 1);
    skip_whitespace();
    if (peek() != ')')
        syntax_error("Expecting `)` after DDS channel", colno + 1);
    colno++;
    skip_whitespace();
    return chn;
}

}
}
