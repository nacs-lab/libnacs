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

#include <ctime>
#include <iomanip>
#include <limits>
#include <fstream>

namespace NaCs {

namespace {

constexpr auto eofc = std::istream::traits_type::eof();

static std::istream &ignore_line(std::istream &stm)
{
    stm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return stm;
}

static std::istream &ignore_space(std::istream &stm)
{
    while (true) {
        auto c = stm.peek();
        if (c == eofc)
            break;
        if (c != ' ' && c != '\t')
            break;
        stm.get();
    }
    return stm;
}

}

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
        if (!try_parseline(stm)) {
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
        started = true;
        prev_read = pos;
        stm >> ignore_line;
    }
}

bool WavemeterParser::try_parseline(std::istream &stm)
{
    std::tm timedate;
    memset(&timedate, 0, sizeof(timedate));
    stm >> ignore_space >> std::get_time(&timedate, "%Y-%m-%dT%H:%M:%S");
    // Read time failed
    if (stm.eof() || !stm)
        return false;
    timedate.tm_isdst = 0;
#ifdef NACS_OS_WINDOWS
    auto ts = _mkgmtime(&timedate);
#else
    auto ts = timegm(&timedate);
#endif
    // Convert time failed
    if (ts == -1)
        return false;
    auto tsf = (double)ts;
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
    double val = 0;
    stm >> val;
    // Read value failed
    if (stm.fail())
        return false;
    stm >> ignore_space;
    auto nc = stm.peek();
    // Not the end of current record
    if (nc != ',' && nc != '\n' && nc != eofc)
        return false;
    m_time.push_back(tsf / 86400 + 719529);
    m_data.push_back(val);
    return true;
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

}

extern "C" {

using namespace NaCs;

NACS_EXPORT() void *nacs_seq_new_wavemeter_parser(void)
{
    return new WavemeterParser;
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
