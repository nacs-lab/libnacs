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
        all_ttl_mask = all_ttl_mask | (1 << chn);
        last_timed_inst = addInst(Inst::TTL1{OpCode::TTL1, 0, val, uint8_t(chn & 0x1f)});
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
    std::pair<const std::string&,int> read_name(bool allow_num0=false)
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
        if (res < lo || res > hi || errno == ERANGE)
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

    void parse_ttlall(Writer &writer)
    {
        assert(peek() == '=');
        colno++;
        writer.addTTL((uint32_t)read_hex(0, UINT32_MAX).first);
    }

    void parse_ttl1(Writer &writer)
    {
        assert(peek() == '(');
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 31).first;
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
        writer.addTTL1(chn, val);
    }

    void parse_waittime(Writer &writer, bool isttl=false)
    {
        auto t_hex = read_hex(3);
        if (t_hex.second != -1) {
            assert(t_hex.first >= 3);
            writer.addWait(t_hex.first - (isttl ? 3 : 0));
            return;
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
        writer.addWait(((uint64_t)t_ns) / 10 - (isttl ? 3 : 0));
    }

    void parse_ttl(Writer &writer)
    {
        skip_whitespace();
        auto c0 = peek();
        if (c0 == '=') {
            parse_ttlall(writer);
        }
        else if (c0 == '(') {
            parse_ttl1(writer);
        }
        else {
            syntax_error("Invalid ttl command: expecting `(` or `=`", colno + 1);
        }
        c0 = peek();
        if (!c0 || c0 == '#')
            return;
        if (!std::isspace(c0))
            syntax_error("Expecting space after TTL value", colno + 1);
        auto tres = read_name();
        if (tres.first != "t") {
            if (tres.second != -1)
                colno = tres.second;
            return;
        }
        skip_whitespace();
        if (peek() != '=')
            syntax_error("Invalid ttl time: expecting `=`", colno + 1);
        colno++;
        parse_waittime(writer, true);
    }

    void parse_wait(Writer &writer)
    {
        skip_whitespace();
        if (peek() != '(')
            syntax_error("Invalid wait command: expecting `(`", colno + 1);
        colno++;
        parse_waittime(writer);
        skip_whitespace();
        if (peek() != ')')
            syntax_error("Invalid wait command: expecting `)`", colno + 1);
        colno++;
    }

    void parse_freq(Writer &writer)
    {
        if (peek() != '(')
            syntax_error("Invalid freq command: expecting `(`", colno + 1);
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 21).first;
        skip_whitespace();
        if (peek() != ')')
            syntax_error("Expecting `)` after DDS channel", colno + 1);
        colno++;
        skip_whitespace();
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
                syntax_error("Frequency too high (max 1.75GHz)", -1, freq_start + 1, colno);
            constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;
            freq = uint32_t(0.5 + freq_hz * freq_factor);
        }
        writer.addDDSFreq(chn, freq);
    }

    void parse_amp(Writer &writer)
    {
        if (peek() != '(')
            syntax_error("Invalid amp command: expecting `(`", colno + 1);
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 21).first;
        skip_whitespace();
        if (peek() != ')')
            syntax_error("Expecting `)` after DDS channel", colno + 1);
        colno++;
        skip_whitespace();
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
        writer.addDDSAmp(chn, amp);
    }

    void parse_phase(Writer &writer)
    {
        if (peek() != '(')
            syntax_error("Invalid phase command: expecting `(`", colno + 1);
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 21).first;
        skip_whitespace();
        if (peek() != ')')
            syntax_error("Expecting `)` after DDS channel", colno + 1);
        colno++;
        skip_whitespace();
        bool det = false;
        if (peek() == '+') {
            if (peek(1) != '=')
                syntax_error("Expecting `=` or `+=` before phase value", colno + 1);
            colno += 1;
        }
        else if (peek() != '=') {
            syntax_error("Expecting `=` or `+=` before phase value", colno + 1);
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
                syntax_error("Phase too high (max +-1000%)", -1, phase_start + 1, colno);
            phase_deg = fmod(phase_deg, 360);
            constexpr double phase_factor = (1 << 14) / 90.0;
            phase = uint16_t(0.5 + phase_deg * phase_factor);
        }
        if (det) {
            writer.addDDSDetPhase(chn, phase);
        }
        else {
            writer.addDDSPhase(chn, phase);
        }
    }

    void parse_reset(Writer &writer)
    {
        if (peek() != '(')
            syntax_error("Invalid reset command: expecting `(`", colno + 1);
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 21).first;
        skip_whitespace();
        if (peek() != ')')
            syntax_error("Expecting `)` after DDS channel", colno + 1);
        colno++;
        writer.addDDSReset(chn);
    }

    void parse_dac(Writer &writer)
    {
        if (peek() != '(')
            syntax_error("Invalid dac command: expecting `(`", colno + 1);
        colno++;
        uint8_t chn = (uint8_t)read_dec(0, 4).first;
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
        writer.addDAC(chn, dac);
    }

    bool parse_cmd(Writer &writer)
    {
        skip_whitespace();
        auto nres = read_name();
        if (nres.second == -1)
            syntax_error("Expecting command name", colno + 1);
        if (nres.first == "ttl") {
            parse_ttl(writer);
        }
        else if (nres.first == "wait") {
            parse_wait(writer);
        }
        else if (nres.first == "freq") {
            parse_freq(writer);
        }
        else if (nres.first == "amp") {
            parse_amp(writer);
        }
        else if (nres.first == "phase") {
            parse_phase(writer);
        }
        else if (nres.first == "reset") {
            parse_reset(writer);
        }
        else if (nres.first == "dac") {
            parse_dac(writer);
        }
        else {
            syntax_error("Unknown command name", -1, nres.second + 1, colno);
        }
        if (!checked_next_line())
            return false;
        return skip_comments();
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
    else if (res.second != -1) {
        parser.colno = res.second;
    }

    Writer writer(ostm, ttl_mask);
    while (parser.parse_cmd(writer)) {
    }
    return writer.get_ttl_mask();
}

}
}
}
