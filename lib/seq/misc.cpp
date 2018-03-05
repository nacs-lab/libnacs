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
            return stm;
        if (c != ' ' && c != '\t')
            return stm;
        stm.get();
    }
}

}

bool WavemeterParser::do_parse(std::istream &stm, bool inc)
{
    // Append to `m_time` and `m_data`.
    // Update `m_last_pos` and `m_last_line` if the parsing succeeded.
    bool started = inc;
    pos_type prev_read = 0;
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
                std::getline(stm, m_last_line);
                return true;
            }
            if (stm.eof())
                return false;
            stm.clear();
            continue;
        }
        started = true;
        prev_read = pos;
    }
}

bool WavemeterParser::try_parseline(std::istream &stm)
{
    std::tm timedate;
    memset(&timedate, 0, sizeof(timedate));
    stm >> ignore_space >> std::get_time(&timedate, "%Y-%m-%dT%H:%M:%S");
    // Read time failed
    if (stm.fail())
        return false;
    timedate.tm_isdst = -1;
    auto ts = std::mktime(&timedate);
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
    m_time.push_back(tsf);
    m_data.push_back(val);
    return true;
}

}
