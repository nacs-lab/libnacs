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

// This parses the output format from the `Printer` in `exechelper_p.h`
// The parser read a line of text each time and keeps track of the current position
// inside the line.
// There is at most one command per line and comments are allowed following `#`
// anywhere in a line.
// The line number and column number are useful for error reporting.
// `lineno` is 1-based, `colno` is 0-based.
// For all of the `read_*` functions, white space are insignificant unless otherwise noted
// and are ignored between syntax elements.
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
    // Peek a few characters ahead within the same line.
    char peek(int idx=0)
    {
        if (colno + idx >= (int)line.size())
            return 0;
        return line[colno + idx];
    }

    ParserBase(std::istream &istm);
    // Go to the next line. If EOF is reached, return `false`.
    bool next_line();
    // Throw a `SyntaxError` based on the current `line` and `lineno`.
    NACS_NORETURN void syntax_error(std::string msg, int colnum,
                                    int colstart=-1, int colend=-1);
    // Skip white space within the current line
    // `peek()` should be either `0` or non-space after this returns.
    void skip_whitespace();
    // Skip comments and empty lines.
    // If EOF is reached, return `false`.
    bool skip_comments();
    // Make sure the line is finished and switch to the next line.
    // If EOF is reached, return `false`.
    // If there are still non-whitespace and non-comment characters in the line,
    // throw an error.
    bool checked_next_line();

    /**
     * Basic read functions:
     *
     * The read functions in this section tries to read the respective format
     * but generall do not throw if no match is found.
     * They only throw error if at least a partial input could match
     * or when the values read are out of range.
     */

    // Skip white space and read a string from the current line.
    // Characters allowed in the string are `A-Za-z0-9` and `_`.
    // If `allow_num0` is `true`, the first character can be a number.
    // Otherwise, only letters and `_` are allowed as the first character.
    // Returns the string and the column index where the string starts.
    std::pair<const std::string&,int> read_name(bool allow_num0=false);
    // Skip white space and read a hex number from the current line.
    // Both upper and lower cases are allowed.
    // The low and high limit (including) can be specified.
    // Return the number followed by the column index where the number started.
    //
    // Return `{0, -1}` without error if no `0x` (or `0X`) prefix is detected.
    // Throws an error if no valid digit after the prefix is detected or
    // if the result is outside the valid range.
    // Negative numbers (i.e. `-0x...`) are not supported.
    std::pair<uint64_t,int> read_hex(uint64_t lo=0, uint64_t hi=UINT64_MAX);
    // Skip white space and read a decimal number from the current line.
    // Can be either signed or unsigned.
    // If the lower limit is no less than `0`, negative sign is not allowed.
    // Return the number followed by the column index where the number started.
    //
    // Return `{0, -1}` without error if no valid digit is detected.
    // Throw error if the value read is out of range.
    template<typename T>
    std::pair<T,int> read_dec(T lo=std::numeric_limits<T>::min(),
                              T hi=std::numeric_limits<T>::max());
    // Skip white space and read a decimal floating point number from the current line.
    // Return the number followed by the column index where the number started.
    //
    // Return `{0, -1}` without error if no valid digit is detected.
    // Throw error if the value read is out of range.
    std::pair<double,int> read_float(double lo=-INFINITY, double hi=INFINITY);

    /**
     * Formatted read functions.
     *
     * Different from the ones above,
     * the read functions in this section require a match unless otherwise specified.
     * An error will be throw if the desired information can't be read from the input.
     *
     * Functions below that reads commands generally assume
     * that the command name has already been from the line in the dispatch function.
     * They will start by parsing the arguments of the command
     * which usually begin with `(` or `=`.
     *
     * Some commands accepts both the low level binary representation and phyisical unit
     * to specify the parameter. In this case, the binary representation will always be
     * specified in hex with `0x` (`0X`) prefix and
     * the physical representation will usually require a unit to minimize ambiguity.
     */

    // Read a DDS channel number in the format of `(<n>)` and return it.
    uint8_t read_ddschn(const char *name);

    // Read a (wait)time specification.
    // This can either be specified as a hex number of cycles (10ns per cycle)
    // or as a decimal floating point followed by a unit.
    // The unit can be either `s`, `ms`, `us` or `ns`.
    // The smallest time allowed is `PulseTime::Min` cycles or `30ns`.
    // Return the number of cycles correspond to the time.
    uint64_t read_waittime();
    // Read a wait command in the format of `(<time>)`.
    // Return the number of cycles to wait.
    uint64_t read_waitcmd();
    // Read a time specification in a TTL command.
    // The returned number is the number of cycles to wait in additional the the
    // minimum time of a TTL pulse.
    // The time specification could either be empty, in which case `0` is returned,
    // or `t=<time>`, in which case the cycle number minus PulseTime::Min is returned.
    uint64_t read_ttlwait();

    // Read a TTL all command. Assume (and assert) the first character is `=`.
    // The TTL value returned is parsed from a hex literal following the `=`.
    uint32_t read_ttlall();
    // Read a single TTL setting command. Assume (and assert) the first character is `(`.
    // The format is `(<chn>)=<value>` where `chn` is a decimal integer in `[0, 31]`
    // and `value` can be one of `on`, `On`, `ON`, `true`, `True`, `TRUE` or `1`
    // for true logic, or one of `off`, `Off`, `OFF`, `false`, `False`, `FALSE` or `0`
    // for false logic.
    // Returns `{chn, value}`
    std::pair<uint8_t,bool> read_ttl1();

    // Read a DDS frequency command.
    // Format is `(<chn>)=<freq>`. `chn` is a decimal integer in `[0, 21]`.
    // `freq` can be specified in DDS frequency word in hex
    // or as decimal floating point using units `Hz`, `kHz`, `MHz` or `GHz`.
    // Frequency must be between `[0, 1.75GHz]`.
    // Returns `{chn, freq}`
    std::pair<uint8_t,uint32_t> read_freqcmd();
    // Read a DDS amplitude command.
    // Format is `(<chn>)=<amp>`. `chn` is a decimal integer in `[0, 21]`.
    // `amp` can be specified in DDS amplitude word in hex or as decimal floating point.
    // Amplitude must be between `[0, 1]`, making the decimal point required in most cases.
    // This is enough to avoid ambiguity with the hex representation so the unit isn't
    // required for this command.
    // Returns `{chn, amp}`
    std::pair<uint8_t,uint16_t> read_ampcmd();
    // Read a DDS amplitude command.
    // Format is `(<chn>)(|-|+)=<phase>`. `chn` is a decimal integer in `[0, 21]`.
    // `phase` can be specified in DDS phase word in hex
    // or as decimal floating point using units `%` (`100%` is `360deg`), `deg`,
    // `pi` (`2pi` is `360deg`) or `rad` (`2pi rad` (not an allowed syntax) is `360deg`).
    // Phase must be between `[-3600, 3600]`.
    // If `+` or `-` are present before the `=`,
    // this is a phase shift command rather than a set phase command.
    // Returns `{chn, {is_phase_shift, phase}}`
    std::pair<uint8_t,std::pair<bool,uint16_t>> read_phasecmd();

    // Read a DAC command.
    // Format is `(<chn>)=<V>`. `chn` is a decimal integer in `[0, 3]`.
    // `V` can be specified in DAC number in hex
    // or as decimal floating point using either `mV` or `V` as unit.
    // Amplitude must be between `[-10, 10]`.
    // Returns `{chn, V}`
    std::pair<uint8_t,uint16_t> read_daccmd();
    // Read a clock command.
    // Format is `(<div>)`. `div` is a decimal integer in `[0, 255]` or the string `off`,
    // which is equivalent to `255`.
    // Returns the clock divider.
    uint8_t read_clockcmd();

    // Read a `ttl_mask` specification.
    // This may be present at the beginning of the text representation.
    // The final `ttl_mask` in effect will also include
    // all the TTL channels used in the cmdlist (a `ttl=...` command uses all the channels).
    // Return `0` if no `ttl_mask` spec is found.
    // The spec format is `ttl_mask=<val>` where `val` is a 32bit hex literal.
    // Only throw error if the parsing failed after matching `ttl_mask`.
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
    colno = startcol + int(endptr - startptr);
    return {res, startcol};
}

extern template std::pair<int,int> ParserBase::read_dec<int>(int lo, int hi);

}
}

#endif // __NACS_SEQ_PARSER_P_H__
