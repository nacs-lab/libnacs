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

#include <assert.h>

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


uint64_t ParserBase::read_waittime()
{
    auto t_hex = read_hex(3);
    if (t_hex.second != -1) {
        assert(t_hex.first >= 3);
        return t_hex.first;
    }
    auto t_flt = read_float();
    if (t_flt.second == -1)
        syntax_error("Invalid time", colno + 1);
    skip_whitespace();
    auto unit = read_name();
    if (unit.second == -1)
        syntax_error("Missing time unit", colno + 1);
    double t_ns = t_flt.first;
    if (unit.first == "ns") {
    }
    else if (unit.first == "us") {
        t_ns = t_ns * 1000.0;
    }
    else if (unit.first == "ms") {
        t_ns = t_ns * 1000000.0;
    }
    else if (unit.first == "s") {
        t_ns = t_ns * 1000000000.0;
    }
    else {
        syntax_error("Unknown time unit", -1, unit.second + 1, colno);
    }
    if (t_ns < 30)
        syntax_error("Time too short (min 30ns)", -1, t_flt.second + 1, colno);
    return ((uint64_t)t_ns) / 10;
}

uint64_t ParserBase::read_waitcmd()
{
    skip_whitespace();
    if (peek() != '(')
        syntax_error("Invalid wait command: expecting `(`", colno + 1);
    colno++;
    auto res = read_waittime();
    skip_whitespace();
    if (peek() != ')')
        syntax_error("Invalid wait command: expecting `)`", colno + 1);
    colno++;
    return res;
}

uint64_t ParserBase::read_ttlwait()
{
    auto c0 = peek();
    if (!c0 || c0 == '#')
        return 0;
    if (!std::isspace(c0))
        syntax_error("Expecting space after TTL value", colno + 1);
    auto tres = read_name();
    if (tres.first != "t") {
        if (tres.second != -1)
            colno = tres.second;
        return 0;
    }
    skip_whitespace();
    if (peek() != '=')
        syntax_error("Invalid ttl time: expecting `=`", colno + 1);
    colno++;
    return read_waittime() - 3;
}

uint32_t ParserBase::read_ttlall()
{
    assert(peek() == '=');
    colno++;
    auto res = read_hex(0, UINT32_MAX);
    if (res.second == -1)
        syntax_error("Missing TTL value", colno + 1);
    return (uint32_t)res.first;
}

std::pair<uint8_t,bool> ParserBase::read_ttl1()
{
    assert(peek() == '(');
    colno++;
    uint8_t chn;
    int chn_start;
    std::tie(chn, chn_start) = read_dec(0, 31);
    if (chn_start == -1)
        syntax_error("Missing TTL channel", colno + 1);
    skip_whitespace();
    if (peek() != ')')
        syntax_error("Expecting `)` after TTL channel", colno + 1);
    colno++;
    skip_whitespace();
    if (peek() != '=')
        syntax_error("Expecting `=` before TTL value", colno + 1);
    colno++;
    bool val;
    auto rres = read_name(true);
    if (rres.first == "true" || rres.first == "True" || rres.first == "TRUE" ||
        rres.first == "on" || rres.first == "On" || rres.first == "ON" ||
        rres.first == "1") {
        val = true;
    }
    else if (rres.first == "false" || rres.first == "False" || rres.first == "FALSE" ||
             rres.first == "off" || rres.first == "Off" || rres.first == "OFF" ||
             rres.first == "0") {
        val = false;
    }
    else if (rres.second == -1) {
        syntax_error("Expecting TTL value after `=`", colno + 1);
    }
    else {
        syntax_error("Invalid TTL value", -1, rres.second + 1, colno);
    }
    return {chn, val};
}

std::pair<uint8_t,uint32_t> ParserBase::read_freqcmd()
{
    auto chn = read_ddschn("freq");
    if (peek() != '=')
        syntax_error("Expecting `=` before frequency value", colno + 1);
    colno++;
    uint32_t freq;
    int freq_start;
    std::tie(freq, freq_start) = read_hex(0, 0x7fffffff);
    if (freq_start == -1) {
        double freq_hz;
        int freq_hz_start;
        std::tie(freq_hz, freq_hz_start) = read_float(0);
        if (freq_hz_start == -1)
            syntax_error("Invalid frequency", colno + 1);
        auto unit = read_name();
        if (unit.second == -1)
            syntax_error("Missing frequency unit", colno + 1);
        if (unit.first == "Hz") {
        }
        else if (unit.first == "kHz") {
            freq_hz = freq_hz * 1000.0;
        }
        else if (unit.first == "MHz") {
            freq_hz = freq_hz * 1000000.0;
        }
        else if (unit.first == "GHz") {
            freq_hz = freq_hz * 1000000000.0;
        }
        else {
            syntax_error("Unknown frequency unit", -1, unit.second + 1, colno);
        }
        if (freq_hz > 1.75e9)
            syntax_error("Frequency too high (max 1.75GHz)", -1, freq_hz_start + 1, colno);
        constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;
        freq = uint32_t(0.5 + freq_hz * freq_factor);
    }
    return {chn, freq};
}

std::pair<uint8_t,uint16_t> ParserBase::read_ampcmd()
{
    auto chn = read_ddschn("amp");
    if (peek() != '=')
        syntax_error("Expecting `=` before amplitude value", colno + 1);
    colno++;
    uint16_t amp;
    int amp_start;
    std::tie(amp, amp_start) = read_hex(0, 4095);
    if (amp_start == -1) {
        double ampf;
        int ampf_start;
        std::tie(ampf, ampf_start) = read_float(0, 1);
        if (ampf_start == -1)
            syntax_error("Invalid amplitude", colno + 1);
        amp = uint16_t(ampf * 4095.0 + 0.5);
    }
    return {chn, amp};
}

std::pair<uint8_t,std::pair<bool,uint16_t>> ParserBase::read_phasecmd()
{
    auto chn = read_ddschn("phase");
    bool det = false;
    bool neg = false;
    if (peek() == '+') {
        if (peek(1) != '=')
            syntax_error("Expecting `=`, `-=` or `+=` before phase value", colno + 1);
        det = true;
        colno += 1;
    }
    else if (peek() == '-') {
        if (peek(1) != '=')
            syntax_error("Expecting `=`, `-=` or `+=` before phase value", colno + 1);
        det = true;
        neg = true;
        colno += 1;
    }
    else if (peek() != '=') {
        syntax_error("Expecting `=`, `-=` or `+=` before phase value", colno + 1);
    }
    colno++;
    uint16_t phase;
    int phase_start;
    std::tie(phase, phase_start) = read_hex(0, UINT16_MAX);
    if (phase_start == -1) {
        double phase_deg;
        int phase_deg_start;
        std::tie(phase_deg, phase_deg_start) = read_float();
        if (phase_deg_start == -1)
            syntax_error("Invalid phase", colno + 1);
        auto unit = read_name();
        if (unit.second == -1) {
            if (peek() != '%')
                syntax_error("Missing phase unit", colno + 1);
            colno++;
            phase_deg = phase_deg * 3.60;
        }
        else if (unit.first == "deg") {
        }
        else if (unit.first == "pi") {
            phase_deg = phase_deg * 180;
        }
        else if (unit.first == "rad") {
            phase_deg = phase_deg * (180 / M_PI);
        }
        else {
            syntax_error("Unknown phase unit", -1, unit.second + 1, colno);
        }
        if (!(abs(phase_deg) <= 360 * 10))
            syntax_error("Phase too high (max +-1000%)", -1, phase_deg_start + 1, colno);
        phase_deg = fmod(phase_deg, 360);
        constexpr double phase_factor = (1 << 14) / 90.0;
        phase = uint16_t(0.5 + phase_deg * phase_factor);
    }
    if (neg)
        phase = uint16_t(-phase);
    return {chn, {det, phase}};
}

std::pair<uint8_t,uint16_t> ParserBase::read_daccmd()
{
    if (peek() != '(')
        syntax_error("Invalid dac command: expecting `(`", colno + 1);
    colno++;
    uint8_t chn;
    int chn_start;
    std::tie(chn, chn_start) = read_dec(0, 4);
    if (chn_start == -1)
        syntax_error("Missing DAC channel", colno + 1);
    skip_whitespace();
    if (peek() != ')')
        syntax_error("Expecting `)` after DAC channel", colno + 1);
    colno++;
    skip_whitespace();
    if (peek() != '=')
        syntax_error("Expecting `=` before DAC value", colno + 1);
    colno++;
    uint16_t dac;
    int dac_start;
    std::tie(dac, dac_start) = read_hex(0, UINT16_MAX);
    if (dac_start == -1) {
        double dac_mv;
        int dac_mv_start;
        std::tie(dac_mv, dac_mv_start) = read_float();
        if (dac_mv_start == -1)
            syntax_error("Invalid DAC voltage", colno + 1);
        auto unit = read_name();
        if (unit.second == -1)
            syntax_error("Missing voltage unit", colno + 1);
        if (unit.first == "mV") {
        }
        else if (unit.first == "V") {
            dac_mv = dac_mv * 1000;
        }
        else {
            syntax_error("Unknown voltage unit", -1, unit.second + 1, colno);
        }
        if (!(abs(dac_mv) <= 10000))
            syntax_error("DAC voltage too high (max +-10V)", -1, dac_start + 1, colno);
        // this is for the DAC8814 chip in SPI0
        constexpr double scale = 65535 / 20000.0;
        constexpr double offset = 10000.0;
        dac = uint16_t(((offset - dac_mv) * scale) + 0.5);
    }
    return {chn, dac};
}

uint8_t ParserBase::read_clockcmd()
{
    if (peek() != '(')
        syntax_error("Invalid clock command: expecting `(`", colno + 1);
    colno++;
    uint8_t clock;
    auto state = read_name();
    if (state.second != -1) {
        if (state.first != "off")
            syntax_error("Invalid clock state", -1, state.second + 1, colno);
        clock = 255;
    }
    else {
        int clock_start;
        std::tie(clock, clock_start) = read_dec(0, 255);
        if (clock_start == -1) {
            syntax_error("Missing clock state", colno + 1);
        }
    }
    skip_whitespace();
    if (peek() != ')')
        syntax_error("Expecting `)` after clock state", colno + 1);
    colno++;
    return clock;
}

std::pair<bool,uint32_t> ParserBase::read_ttlmask()
{
    if (!skip_comments())
        return {false, 0};
    auto res = read_name();
    uint32_t ttl_mask = 0;
    if (res.first == "ttl_mask") {
        skip_whitespace();
        if (peek() != '=')
            syntax_error("Expecting `=` after `ttl_mask`", colno + 1);
        colno++;
        skip_whitespace();
        int oldcol;
        std::tie(ttl_mask, oldcol) = read_hex(0, UINT32_MAX);
        if (oldcol < 0)
            syntax_error("Expecting hex literal", colno + 1);
        if (!checked_next_line() || !skip_comments()) {
            return {false, ttl_mask};
        }
    }
    else if (res.second != -1) {
        colno = res.second;
    }
    return {true, ttl_mask};
}

}
}
