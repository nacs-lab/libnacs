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

#include "wavemeter.h"

#include "streams.h"
#include "utils.h"

#include <ctime>
#include <iomanip>

#include <time.h>

namespace NaCs {

NACS_INTERNAL bool Wavemeter::try_parsetime(std::istream &stm, double *tsf)
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
    *tsf = (double)ts;
    if (stm.peek() == '.') {
        double f = 0;
        stm >> f;
        // Read second fraction failed
        if (stm.fail())
            return false;
        *tsf += f;
    }
    stm >> ignore_space;
    // Extra characters after timestamp
    if (stm.peek() != ',')
        return false;
    stm.get();
    return true;
}

NACS_INTERNAL bool Wavemeter::try_parsenumber(std::istream &stm, double *val, bool *eol)
{
    stm >> *val;
    // Read value failed
    if (stm.fail())
        return false;
    stm >> ignore_space;
    auto nc = stm.peek();
    if (nc == '\n' || nc == eofc) {
        *eol = true;
    }
    else if (nc != ',') {
        // Not the end of current record
        return false;
    }
    return true;
}

NACS_INTERNAL bool Wavemeter::try_parseval(std::istream &stm, double *val,
                                           double lo, double hi)
{
    bool found_val = false;
    double max_pos = 0;
    double max_height = std::numeric_limits<double>::lowest();

    bool eol;
    auto parsenum = [&] (double *val) { return try_parsenumber(stm, val, &eol); };
    auto getpair = [&] {
        double pos;
        double height;
        if (!parsenum(&pos))
            return false;
        if (eol)
            return false;
        stm.get();
        if (!parsenum(&height))
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
    *val = max_pos;
    return found_val;
}

bool Wavemeter::try_parseline(std::istream &stm, double *tsf, double *val) const
{
    if (!try_parsetime(stm, tsf))
        return false;
    *tsf = *tsf / 86400 + 719529;
    return try_parsenumber(stm, val);
}

static const std::istream::pos_type pos_error = std::streamoff(-1);

NACS_INTERNAL auto Wavemeter::find_linestart(std::istream &stm, pos_type ub,
                                             pos_type lb) -> pos_type
{
    char buff[150];
    pos_type loc = ub;
    while (loc > lb) {
        size_t sz = loc - lb;
        if (sz > sizeof(buff))
            sz = sizeof(buff);
        loc -= sz;
        stm.seekg(loc);
        stm.read(buff, sizeof(buff));
        if (!stm)
            return pos_error;
        for (int i = sz - 1; i >= 0; i--) {
            if (buff[i] == '\n') {
                return loc + std::streamoff(i);
            }
        }
    }
    return lb;
}

}
