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

#ifndef __NACS_SEQ_PARSER_P_H__
#define __NACS_SEQ_PARSER_P_H__

#include "../utils/utils.h"

#include <cmath>
#include <istream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace NaCs {
namespace Seq {

struct ParserBase {
    std::istream &istm;
    std::string line;
    std::string buff;
    int lineno = 0;
    int colno = 0;

    template<typename T>
    static std::string to_hexstring(T v)
    {
        std::ostringstream stm;
        stm << "0x" << std::hex << v;
        return stm.str();
    }
    char peek(int idx=0)
    {
        if (colno + idx >= (int)line.size())
            return 0;
        return line[colno + idx];
    }

    ParserBase(std::istream &istm);
    bool next_line();
    NACS_NORETURN void syntax_error(std::string msg, int colnum, int colstart=-1, int colend=-1);
    void skip_whitespace();
    bool skip_comments();
    std::pair<const std::string&,int> read_name(bool allow_num0=false);
    std::pair<uint64_t,int> read_hex(uint64_t lo=0, uint64_t hi=UINT64_MAX);
    template<typename T>
    std::pair<T,int> read_dec(T lo=std::numeric_limits<T>::min(),
                              T hi=std::numeric_limits<T>::max());
    std::pair<double,int> read_float(double lo=-INFINITY, double hi=INFINITY);
    bool checked_next_line();

    uint8_t read_ddschn(const char *name);

    uint64_t read_waittime();
    uint64_t read_waitcmd();
    uint64_t read_ttlwait();

    uint32_t read_ttlall();
    std::pair<uint8_t,bool> read_ttl1();

    std::pair<uint8_t,uint32_t> read_freqcmd();
    std::pair<uint8_t,uint16_t> read_ampcmd();
    std::pair<uint8_t,std::pair<bool,uint16_t>> read_phasecmd();

    std::pair<uint8_t,uint16_t> read_daccmd();
    uint8_t read_clockcmd();

    std::pair<bool,uint32_t> read_ttlmask();
};

template<typename T>
std::pair<T,int> ParserBase::read_dec(T lo, T hi)
{
    skip_whitespace();
    int startcol = colno;
    auto startptr = &line[startcol];
    char *endptr;
    T res;
    if (std::is_signed<T>::value) {
        res = (T)strtoll(startptr, &endptr, 10);
    }
    else {
        res = (T)strtoull(startptr, &endptr, 10);
    }
    if (endptr == startptr)
        return {0, -1};
    if (res < lo || res > hi || errno == ERANGE)
        syntax_error("Number literal out of range [" + std::to_string(lo) + ", " +
                     std::to_string(hi) + "]", -1,
                     startcol + 1, startcol + int(endptr - startptr));
    // `strtoull` allows `-`, which may give surprising result for us.
    if (lo >= 0 && peek() == '-')
        syntax_error("Unexpected negative number", startcol + 1,
                     startcol + 1, startcol + int(endptr - startptr));
    colno = startcol + (endptr - startptr);
    return {res, startcol};
}

extern template std::pair<int,int> ParserBase::read_dec<int>(int lo, int hi);

}
}

#endif // __NACS_SEQ_PARSER_P_H__
