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
#include <fstream>
#include <iomanip>
#include <tuple>

#include <assert.h>
#include <time.h>

namespace NaCs {

NACS_INTERNAL bool Wavemeter::parsetime(std::istream &stm, double *tsf)
{
    std::tm timedate;
    memset(&timedate, 0, sizeof(timedate));
    stm >> ignore_space >> std::get_time(&timedate, "%Y-%m-%dT%H:%M:%S");
    // Read time failed
    if (!stm.good())
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
        if (!stm.good())
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

NACS_INTERNAL bool Wavemeter::parsenumber(std::istream &stm, double *val, bool *eol)
{
    stm >> *val;
    // Read value failed
    if (!stm.good())
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

NACS_INTERNAL bool Wavemeter::parseval(std::istream &stm, double *val,
                                       double lo, double hi)
{
    bool found_val = false;
    double max_pos = 0;
    double max_height = std::numeric_limits<double>::lowest();

    bool eol;
    auto parsenum = [&] (double *val) { return parsenumber(stm, val, &eol); };
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
    stm.clear();
    stm >> ignore_line;
    if (!stm.good())
        return false;
    *val = max_pos;
    return found_val;
}

NACS_INTERNAL bool Wavemeter::parseline(std::istream &stm, double *tsf, double *val) const
{
    if (!parsetime(stm, tsf)) {
        stm.clear();
        stm >> ignore_line;
        return false;
    }
    *tsf = *tsf / 86400 + 719529;
    return parseval(stm, val, m_lo, m_hi);
}

static const std::istream::pos_type pos_error = std::streamoff(-1);

NACS_INTERNAL auto Wavemeter::find_linestart(std::istream &stm, pos_type ub,
                                             pos_type lb) -> pos_type
{
    pos_type loc = ub;
    while (loc > lb) {
        char buff[150];
        size_t sz = loc - lb;
        if (sz > sizeof(buff))
            sz = sizeof(buff);
        loc -= sz;
        stm.seekg(loc);
        stm.read(buff, sizeof(buff));
        if (!stm.good())
            throw std::runtime_error("Error finding line start");
        for (int i = sz - 1; i >= 0; i--) {
            if (buff[i] == '\n') {
                return loc + std::streamoff(i);
            }
        }
    }
    return lb;
}

NACS_INTERNAL auto Wavemeter::parse_at(std::istream &stm, pos_type pos, pos_type lb,
                                       double *tsf, double *val) const
    -> std::pair<bool,pos_type>
{
    auto ls = find_linestart(stm, pos, lb);
    stm.seekg(ls);
    return {parseline(stm, tsf, val), ls};
}

NACS_INTERNAL void Wavemeter::parse_until(std::istream &stm, double tmax, pos_type pos_max,
                                          std::vector<double> &times,
                                          std::vector<double> &datas) const
{
    while (true) {
        auto pos = stm.tellg();
        if (pos_max != pos_error && pos >= pos_max)
            return;
        double tsf;
        double val;
        auto res = parseline(stm, &tsf, &val);
        if (!unlikely(res)) {
            // Treat IO error as not possible to read forward
            if (!stm.good()) {
                stm.clear();
                return;
            }
            continue;
        }
        times.push_back(tsf);
        datas.push_back(val);
        if (tsf >= tmax) {
            return;
        }
    }
}

// Find the location that is right before `tstart` and starts parsing there.
// Only look in `[pstart, pend)`.
// The starting time may be after `tstart` if no point before `tstart` exists in the range.
// The returned range may be larger than needed.
// If return value is `true`, `tsf` and `val` contains the number corresponds to the first point,
// and `loc` is the start of the corresponding line. Note that the data is only guaranteed
// to be the first valid line following `loc`, there could be any number of garbage lines
// before it.
NACS_INTERNAL bool Wavemeter::start_parse(std::istream &stm, double tstart,
                                          pos_type pstart, pos_type pend,
                                          double *tsf, double *val, pos_type *loc) const
{
    pos_type lb = pstart;
    pos_type ub = pend;
    bool lb_valid = false;
    std::tuple<double,double,pos_type> lb_res{0, 0, 0};
    while (true) {
        auto sz = ub - lb;
        // Good enough
        if (sz < 10240) {
            stm.seekg(lb);
            std::tie(*tsf, *val, *loc) = lb_res;
            return lb_valid;
        }
        pos_type mid = lb + sz / 2;
        double t, v;
        auto res = parse_at(stm, mid, lb, &t, &v);
        if (unlikely(!res.first)) {
            while (!res.first && stm.good() && stm.tellg() < ub) {
                stm.clear();
                res.first = parseline(stm, &t, &v);
            }
            if (!res.first) {
                // Nothing above us is useful, just update ub
                // This can cause the invalid part to be reparse by the caller but
                // we don't really expect that to happen (in a performance important way)
                // anyway so it's more important to keep the code simpler.
                ub = res.second;
                stm.clear();
                continue;
            }
        }
        assert(res.first);
        if (t > tstart) {
            // mid is ub
            ub = res.second;
        }
        else {
            // mid is lb
            lb = stm.tellg();
            lb_valid = true;
            lb_res = {t, v, res.second};
        }
    }
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

NACS_INTERNAL void Wavemeter::extend_segment(std::istream &stm, seg_ent_t &ent,
                                             double tend, pos_type pend)
{
    if (ent.second.times.back() >= tend)
        return;
    stm.seekg(ent.first + std::streamoff(ent.second.size));
    parse_until(stm, tend, pend, ent.second.times, ent.second.datas);
    ent.second.size = stm.tellg() - ent.first;
}

NACS_INTERNAL auto Wavemeter::new_segment_end(std::istream &stm, double tstart, double tend,
                                              pos_type lb, pos_type ub) -> const Segment*
{
    double tsf, val;
    pos_type loc;
    bool valid = start_parse(stm, tstart, lb, ub, &tsf, &val, &loc);
    if (!stm.good())
        return nullptr;
    std::vector<double> times;
    std::vector<double> datas;
    if (valid) {
        times.push_back(tsf);
        datas.push_back(val);
    }
    else {
        loc = stm.tellg();
    }
    parse_until(stm, tend, ub, times, datas);
    if (times.empty())
        return nullptr;
    auto it = m_segments.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(loc),
                                 std::forward_as_tuple(std::move(times), std::move(datas),
                                                       stm.tellg() - loc)).first;
    return &it->second;
}

NACS_INTERNAL auto Wavemeter::get_segment(std::istream &stm, double tstart,
                                          double tend) -> const Segment*
{
    // TODO
    return nullptr;
}

NACS_EXPORT() std::pair<const double*,const double*>
Wavemeter::parse(std::istream &stm, size_t *sz, double tstart, double tend)
{
    if (unlikely(tend <= tstart)) {
        *sz = 0;
        return {nullptr, nullptr};
    }
    auto seg = get_segment(stm, tstart, tend);
    auto it1 = std::lower_bound(seg->times.begin(), seg->times.end(), tstart);
    if (it1 == seg->times.end())
        return {nullptr, nullptr};
    auto idx1 = it1 - seg->times.begin();
    auto it2 = std::lower_bound(it1, seg->times.end(), tend);
    auto idx2 = it2 == seg->times.end() ? seg->times.size() - 1 : it2 - seg->times.begin();
    *sz = idx2 - idx1;
    return {&*it1, &seg->datas[idx1]};
}

NACS_EXPORT() Wavemeter::Wavemeter(double lo, double hi)
: m_lo(lo), m_hi(hi)
{
}

}

extern "C" {

using namespace NaCs;

NACS_EXPORT() void *nacs_utils_new_wavemeter(double lo, double hi)
{
    return new Wavemeter(lo, hi);
}

NACS_EXPORT() size_t nacs_utils_wavemeter_parse(void *_parser, const char *name,
                                                const double **ts, const double **data,
                                                double tstart, double tend)
{
    auto parser = (Wavemeter*)_parser;
    size_t sz = 0;
    std::ifstream stm(name);
    auto ptrs = parser->parse(stm, &sz, tstart, tend);
    *ts = ptrs.first;
    *data = ptrs.second;
    return sz;
}

NACS_EXPORT() void nacs_utils_free_wavemeter(void *parser)
{
    delete (Wavemeter*)parser;
}

}
