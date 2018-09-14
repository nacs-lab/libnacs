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

#include "exehelper_p.h"

#include "cmdlist.h"

#include "../utils/streams.h"
#include "../utils/errors.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <type_traits>

namespace NaCs {
namespace Seq {
namespace CmdList {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len)
{
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b;
        assert(op < 10);
        auto len = cmd_size[op];
        i += len;
        count++;
    }
    return count;
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len)
{
    print(stm, code, code_len, 0);
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask)
{
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
    Printer printer{stm};
    ExeState state;
    state.run(printer, code, code_len);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len)
{
    TimeKeeper keeper;
    ExeState state;
    state.run(keeper, code, code_len);
    return keeper.total_t;
}

namespace {

/**
 * The writer class maintains all the states for cmdlist generation.
 */
class Writer {
    // States
    uint32_t all_ttl_mask;

    uint8_t max_time_left = 0;
    ssize_t last_timed_inst = 0;

    buff_ostream &stm;

    template<typename Inst>
    ssize_t addInst(Inst inst)
    {
        auto len = stm.tellg();
        stm.write((char*)&inst, sizeof(inst));
        max_time_left = 0;
        return (ssize_t)len;
    }

    // Increase the wait time encoded in the last instruction by `t`.
    // The caller should have checked that the time fits in the maximum time possible to be
    // encoded.
    void incLastTime(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t op = stm[last_timed_inst];
        uint8_t *tp;
        uint8_t tmask;
        int tshift;
        if (op == 1) {
            // TTL1
            tmask = 0xc0;
            tp = (uint8_t*)&stm[last_timed_inst + 1];
            tshift = 6;
        }
        else {
            assert(0 && "Invalid command to increase time.");
            abort();
        }
        uint8_t b = *tp;
        t = uint8_t(t + ((b & tmask) >> tshift));
        *tp = uint8_t((b & ~tmask) | ((t << tshift) & tmask));
    }

public:
    // The `add***` functions below provides an abstraction to the cmdlist and hides
    // the detail about the cmdlist encoding.
    // This is the same level as the API of `ExeState`.
    void addTTL(uint32_t ttl)
    {
        all_ttl_mask = uint32_t(-1);
        addInst(Inst::TTLAll{OpCode::TTLAll, ttl});
    }

    void addTTL1(uint8_t chn, bool val)
    {
        chn = chn & 0x1f;
        all_ttl_mask = all_ttl_mask | (1 << chn);
        last_timed_inst = addInst(Inst::TTL1{OpCode::TTL1, 0, val, chn});
        max_time_left = 3;
    }

    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        if (dt <= max_time_left) {
            incLastTime(uint8_t(dt));
            return;
        }
        if (max_time_left) {
            dt -= max_time_left;
            incLastTime(max_time_left);
        }
        addInst(Inst::Wait{OpCode::Wait, dt});
    }

    void addClock(uint8_t period)
    {
        addInst(Inst::Clock{OpCode::Clock, period});
    }

    void addDDSFreq(uint8_t chn, uint32_t freq)
    {
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        addInst(Inst::DDSFreq{OpCode::DDSFreq, chn, freq});
    }

    void addDDSAmp(uint8_t chn, uint16_t amp)
    {
        if (amp > 4095)
            amp = 4095;
        addInst(Inst::DDSAmp{OpCode::DDSAmp, chn, amp});
    }

    void addDDSPhase(uint8_t chn, uint16_t phase)
    {
        addInst(Inst::DDSPhase{OpCode::DDSPhase, chn, phase});
    }

    void addDDSDetPhase(uint8_t chn, uint16_t det_phase)
    {
        addInst(Inst::DDSDetPhase{OpCode::DDSDetPhase, chn, det_phase});
    }

    void addDDSReset(uint8_t chn)
    {
        addInst(Inst::DDSReset{OpCode::DDSReset, chn});
    }

    void addDAC(uint8_t chn, uint16_t amp)
    {
        addInst(Inst::DAC{OpCode::DAC, chn, amp});
    }

    Writer(buff_ostream &stm, uint32_t ttl_mask)
        : all_ttl_mask(ttl_mask),
          stm(stm)
    {}

    uint32_t get_ttl_mask() const
    {
        return all_ttl_mask;
    }
};

struct Parser {
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

    Parser(std::istream &istm)
        : istm(istm)
    {
        next_line();
    }
    bool next_line()
    {
        std::getline(istm, line);
        lineno++;
        colno = 0;
        return istm.good();
    }
    void syntax_error(std::string msg, int colnum, int colstart=-1, int colend=-1)
    {
        throw SyntaxError(std::move(msg), line, lineno, colnum, colstart, colend);
    }
    void skip_whitespace()
    {
        int linelen = (int)line.size();
        for (; colno < linelen; colno++) {
            if (!std::isspace(line[colno])) {
                return;
            }
        }
    }
    bool skip_comments()
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
    char peek(int idx=0)
    {
        if (colno + idx >= (int)line.size())
            return 0;
        return line[colno + idx];
    }
    std::pair<const std::string&,int> read_name()
    {
        buff.clear();

        skip_whitespace();

        int startcol = colno;

        int linelen = (int)line.size();
        if (colno >= linelen)
            return {buff, -1};
        auto c0 = line[colno];
        if (!(c0 >= 'a' && c0 <= 'z') && !(c0 >= 'A' && c0 <= 'Z') && c0 != '_')
            return {buff, -1};
        buff.push_back(c0);
        colno++;
        for (; colno < linelen; colno++) {
            auto c = line[colno];
            if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                (c >= '0' && c <= '9') || c == '_') {
                buff.push_back(c);
            }
            else {
                break;
            }
        }
        return {buff, startcol};
    }
    std::pair<uint64_t,int> read_hex(uint64_t lo=0, uint64_t hi=UINT64_MAX)
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
        if (res < lo || res > hi ||
            (res == std::numeric_limits<decltype(res)>::max() && errno == ERANGE))
            syntax_error("Hex literal out of range [" + to_hexstring(lo) + ", " +
                         to_hexstring(hi) + "]", -1,
                         startcol + 1, startcol + int(endptr - startptr));
        colno = startcol + (endptr - startptr);
        return {res, startcol};
    }
    template<typename T>
    std::pair<T,int> read_dec(T lo=std::numeric_limits<T>::min(),
                              T hi=std::numeric_limits<T>::max())
    {
        skip_whitespace();
        int startcol = colno;
        auto startptr = &line[startcol];
        char *endptr;
        T res;
        if (std::is_signed<T>::value) {
            res = strtoll(startptr, &endptr, 10);
        }
        else {
            res = strtoull(startptr, &endptr, 10);
        }
        if (endptr == startptr)
            return {0, -1};
        if (res < lo || res > hi ||
            (res == std::numeric_limits<T>::max() && errno == ERANGE))
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
    std::pair<double,int> read_float(double lo=-INFINITY, double hi=INFINITY)
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
    bool checked_next_line()
    {
        skip_whitespace();
        auto c = peek();
        if (c && c != '#')
            syntax_error("Unexpected charater at the end of line", colno + 1);
        return next_line();
    }
};

}

NACS_EXPORT() uint32_t parse(buff_ostream &ostm, std::istream &istm)
{
    Parser parser(istm);
    if (!parser.skip_comments())
        return 0;
    auto res = parser.read_name();
    uint32_t ttl_mask = 0;
    if (res.first == "ttl_mask") {
        parser.skip_whitespace();
        if (parser.peek() != '=')
            parser.syntax_error("Expecting `=` after `ttl_mask`", parser.colno + 1);
        parser.colno++;
        parser.skip_whitespace();
        int oldcol;
        std::tie(ttl_mask, oldcol) = parser.read_hex(0, UINT32_MAX);
        if (oldcol < 0)
            parser.syntax_error("Expecting hex literal", parser.colno + 1);
        if (!parser.checked_next_line() || !parser.skip_comments()) {
            return ttl_mask;
        }
    }
    else {
        parser.colno = res.second;
    }

    Writer writer(ostm, ttl_mask);

    return writer.get_ttl_mask();
}

}
}
}
