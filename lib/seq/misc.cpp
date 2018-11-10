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

#include "misc.h"

#include "../utils/streams.h"

#include <ctime>
#include <iomanip>
#include <limits>
#include <fstream>

namespace NaCs {

bool WavemeterParser::do_parse(std::istream &stm, bool inc)
{
    // Append to `m_time` and `m_data`.
    // Update `m_last_pos` and `m_last_line` if the parsing succeeded.
    bool started = inc;
    pos_type prev_read = 0;
    if (inc)
        prev_read = m_last_pos;
    while (true) {
        auto pos = stm.tellg();
        double tsf;
        double val;
        if (!try_parseline(stm, tsf, val)) {
            if (started) {
                // If we've already started, only allow failure on the last line
                stm.clear();
                stm >> ignore_line;
                if (!stm.eof())
                    return false;
                m_last_pos = prev_read;
                stm.seekg(prev_read);
                if (stm)
                    std::getline(stm, m_last_line);
                if (!stm)
                    m_last_line.clear();
                return true;
            }
            if (stm.eof())
                return false;
            stm.clear();
            stm >> ignore_line;
            continue;
        }
        m_time.push_back(tsf);
        m_data.push_back(val);
        started = true;
        prev_read = pos;
        stm >> ignore_line;
    }
}

NACS_INTERNAL bool WavemeterParser::try_parsetime(std::istream &stm, double &tsf)
{
    std::tm timedate;
    memset(&timedate, 0, sizeof(timedate));
    stm >> ignore_space >> std::get_time(&timedate, "%Y-%m-%dT%H:%M:%S");
    // Read time failed
    if (stm.eof() || !stm)
        return false;
    timedate.tm_isdst = 0;
#if NACS_OS_WINDOWS
    auto ts = _mkgmtime(&timedate);
#else
    auto ts = timegm(&timedate);
#endif
    // Convert time failed
    if (ts == -1)
        return false;
    tsf = (double)ts;
    if (stm.peek() == '.') {
        double f = 0;
        stm >> f;
        // Read second fraction failed
        if (stm.fail())
            return false;
        tsf += f;
    }
    stm >> ignore_space;
    // Extra characters after timestamp
    if (stm.peek() != ',')
        return false;
    stm.get();
    return true;
}

NACS_INTERNAL bool WavemeterParser::try_parsenumber(std::istream &stm, double &val, bool &eol)
{
    stm >> val;
    // Read value failed
    if (stm.fail())
        return false;
    stm >> ignore_space;
    auto nc = stm.peek();
    if (nc == '\n' || nc == eofc) {
        eol = true;
    }
    else if (nc != ',') {
        // Not the end of current record
        return false;
    }
    return true;
}

NACS_INTERNAL bool WavemeterParser::try_parseval_withlim(std::istream &stm, double &val,
                                                         double lo, double hi)
{
    bool found_val = false;
    double max_pos = 0;
    double max_height = std::numeric_limits<double>::lowest();

    bool eol;
    auto parsenum = [&] (double &val) { return try_parsenumber(stm, val, eol); };
    auto getpair = [&] {
        double pos;
        double height;
        if (!parsenum(pos))
            return false;
        if (eol)
            return false;
        stm.get();
        if (!parsenum(height))
            return false;
        if (lo <= pos && pos <= hi) {
            if (!found_val || max_height < height) {
                max_height = height;
                max_pos = pos;
                found_val = true;
            }
        }
        else if (!found_val) {
            // The line is valid.
            // If there's no data in range, make sure we return 0 instead of error.
            found_val = true;
        }
        // height in dB could legally be 0
        if (pos == 0)
            eol = true;
        if (!eol)
            stm.get();
        return true;
    };
    while (getpair()) {
        if (eol) {
            break;
        }
    }
    val = max_pos;
    return found_val;
}

bool WavemeterParser::try_parseline(std::istream &stm, double &tsf, double &val) const
{
    if (!try_parsetime(stm, tsf))
        return false;
    tsf = tsf / 86400 + 719529;
    if (m_lo < m_hi)
        return try_parseval_withlim(stm, val, m_lo, m_hi);
    return try_parsenumber(stm, val);
}

bool WavemeterParser::match_cache(std::istream &stm) const
{
    if (m_last_line.empty())
        return false;
    stm.seekg(m_last_pos);
    if (!stm)
        return false;
    std::string l;
    std::getline(stm, l);
    return stm && m_last_line == l;
}

NACS_EXPORT() std::pair<const double*,const double*>
WavemeterParser::parse(std::istream &stm, size_t *sz, bool allow_cache)
{
    if (allow_cache)
        allow_cache = match_cache(stm);
    if (!allow_cache) {
        m_time.clear();
        m_data.clear();
    }
    auto old_size = m_time.size();
    if (!do_parse(stm, allow_cache)) {
        // `m_last_pos` and `m_last_line` are only modified if the return value is `true`
        m_time.resize(old_size);
        m_data.resize(old_size);
        *sz = 0;
        return {nullptr, nullptr};
    }
    *sz = m_time.size();
    return {m_time.data(), m_data.data()};
}

NACS_EXPORT() WavemeterParser::WavemeterParser()
{
}

NACS_EXPORT() WavemeterParser::WavemeterParser(double lo, double hi)
: m_lo(lo), m_hi(hi)
{
}

}

extern "C" {

using namespace NaCs;

NACS_EXPORT() void *nacs_seq_new_wavemeter_parser(void)
{
    return new WavemeterParser;
}

NACS_EXPORT() void *nacs_seq_new_wavemeter_parser_withlim(double lo, double hi)
{
    return new WavemeterParser(lo, hi);
}

NACS_EXPORT() size_t nacs_seq_wavemeter_parse(void *_parser, const char *name,
                                              const double **ts, const double **data,
                                              int cache)
{
    auto parser = (WavemeterParser*)_parser;
    size_t sz = 0;
    std::ifstream stm(name);
    auto ptrs = parser->parse(stm, &sz, cache);
    *ts = ptrs.first;
    *data = ptrs.second;
    return sz;
}

NACS_EXPORT() void nacs_seq_free_wavemeter_parser(void *parser)
{
    delete (WavemeterParser*)parser;
}

}
