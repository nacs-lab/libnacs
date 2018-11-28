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

NACS_INTERNAL bool Wavemeter::try_parseline(std::istream &stm, double *tsf, double *val) const
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
            throw std::runtime_error("Error finding line start");
        for (int i = sz - 1; i >= 0; i--) {
            if (buff[i] == '\n') {
                return loc + std::streamoff(i);
            }
        }
    }
    return lb;
}

NACS_INTERNAL auto Wavemeter::try_parse_at(std::istream &stm, pos_type pos, pos_type lb,
                                           double *tsf, double *val) const
    -> std::pair<bool,pos_type>
{
    auto ls = find_linestart(stm, pos, lb);
    stm.seekg(ls);
    auto res = try_parseline(stm, tsf, val);
    stm.clear();
    stm >> ignore_line;
    return {res, ls};
}

NACS_INTERNAL auto Wavemeter::find_pos_range(double t) const
    -> std::pair<pos_type,pos_type>
{
    // Empty case
    if (m_pos_cache.empty())
        return {0, pos_error};
    auto it = m_pos_cache.upper_bound(t);
    if (likely(it == m_pos_cache.end())) {
        // No range after us, check the last range.
        auto rit = m_pos_cache.rbegin();
        if (rit->second.tend >= t)
            return {rit->second.pstart, rit->second.pend};
        return {rit->second.pend, pos_error};
    }
    pos_type ub = it->second.pstart;
    --it;
    // The one above us is the first range
    if (it == m_pos_cache.end())
        return {0, ub};
    // We are within the previous range
    if (it->second.tend >= t)
        return {it->second.pstart, it->second.pend};
    // We are in between two ranges.
    return {it->second.pend, ub};
}

NACS_INTERNAL void Wavemeter::add_pos_range(double tstart, double tend,
                                            pos_type pstart, pos_type pend)
{
    auto add_range = [&] {
        return m_pos_cache.emplace(tstart, PosRange{tend, pstart, pend}).first;
    };
    // Empty case
    if (m_pos_cache.empty()) {
        add_range();
        return;
    }
    auto it = m_pos_cache.upper_bound(tstart);
    if (likely(it == m_pos_cache.end())) {
        // No range after us, try to merge with the last range
        auto rit = m_pos_cache.rbegin();
        // Already covered.
        if (unlikely(rit->second.tend >= tend))
            return;
        // Not overlapping with the last range
        if (rit->second.tend < tstart) {
            add_range();
            return;
        }
        // Modify the last range
        rit->second.tend = tend;
        rit->second.pend = pend;
        return;
    }
    // For the following cases, we may need to check if any range(s) after us should be merged.
    --it;
    if (it == m_pos_cache.end() || it->second.tend < tstart) {
        // We aren't overlapping with any (previous) ranges.
        it = add_range();
    }
    else if (unlikely(it->second.tend >= tend)) {
        // No new range.
        return;
    }
    else {
        it->second.tend = tend;
        it->second.pend = pend;
    }
    // `it` is the range that contains the new range.
    // We need to check if any range after `it` have overlap.
    for (auto it2 = it++; it2 != m_pos_cache.end(); it2 = m_pos_cache.erase(it2)) {
        // Not overlapping anymore
        if (it2->first > it->second.tend)
            return;
        if (it2->second.tend > it->second.tend) {
            // The end of the new range is after us so we can stop after update.
            it->second.tend = it2->second.tend;
            it->second.pend = it2->second.pend;
            m_pos_cache.erase(it2);
            return;
        }
    }
}

NACS_EXPORT() Wavemeter::Wavemeter(double lo, double hi)
: m_lo(lo), m_hi(hi)
{
}

}
