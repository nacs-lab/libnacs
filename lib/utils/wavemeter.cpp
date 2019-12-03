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

#include "number.h"
#include "streams.h"

#include <ctime>
#include <fstream>
#include <iomanip>
#include <tuple>

#include <assert.h>
#include <time.h>

namespace NaCs {

// Read stream from current position
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

// Read stream from current position
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

// Read stream from current position
NACS_INTERNAL bool Wavemeter::parseval(std::istream &stm, double *val, double *height,
                                       double lo, double hi)
{
    bool line_valid = true;
    double max_pos = 0;
    double max_height = std::numeric_limits<double>::lowest();

    auto getpair = [&] {
        bool eol = false;
        auto parsenum = [&] (double *val) { return parsenumber(stm, val, &eol); };
        double pos;
        double height;
        if (!parsenum(&pos))
            return false;
        // Only check pos since height in dB could legally be 0
        if (eol || pos == 0)
            return false;
        stm.get();
        if (lo > pos || pos > hi) {
            while (true) {
                auto c = stm.peek();
                if (c == ',')
                    break;
                if (c == '\n')
                    return false;
                if (unlikely(c == eofc)) {
                    line_valid = false;
                    return false;
                }
                stm.get();
            }
            stm.get();
            return true;
        }
        if (!parsenum(&height))
            return false;
        if (max_height < height) {
            max_height = height;
            max_pos = pos;
        }
        if (!eol)
            stm.get();
        return !eol;
    };
    while (getpair()) {
    }
    stm.clear();
    stm >> ignore_line;
    if (!stm.good())
        return false;
    *val = max_pos;
    *height = max_height;
    return line_valid;
}

// Read stream from current position
NACS_INTERNAL bool Wavemeter::parseline(std::istream &stm, double *tsf,
                                        double *val, double *height) const
{
    if (!parsetime(stm, tsf)) {
        stm.clear();
        stm >> ignore_line;
        return false;
    }
    return parseval(stm, val, height, m_lo, m_hi);
}

static const std::istream::pos_type pos_error = std::streamoff(-1);

// Always seek the stream.
NACS_INTERNAL auto Wavemeter::find_linestart(std::istream &stm, pos_type ub,
                                             pos_type lb) -> pos_type
{
    pos_type loc = ub;
    while (loc > lb) {
        char buff[150];
        std::streamoff sz = loc - lb;
        if (sz > std::streamoff(sizeof(buff)))
            sz = sizeof(buff);
        loc -= sz;
        stm.seekg(loc);
        stm.read(buff, sizeof(buff));
        if (!stm.good())
            throw std::runtime_error("Error finding line start");
        for (std::streamoff i = sz - 1; i >= 0; i--) {
            if (buff[i] == '\n') {
                return loc + std::streamoff(i);
            }
        }
    }
    return lb;
}

// Always seek the stream.
NACS_INTERNAL auto Wavemeter::parse_at(std::istream &stm, pos_type pos, pos_type lb,
                                       double *tsf, double *val, double *height) const
    -> std::pair<bool,pos_type>
{
    auto ls = find_linestart(stm, pos, lb);
    stm.seekg(ls);
    return {parseline(stm, tsf, val, height), ls};
}

// Read stream from current position
NACS_INTERNAL double Wavemeter::parse_until(std::istream &stm, double tmax, pos_type pos_max,
                                            std::vector<double> &times,
                                            std::vector<double> &datas,
                                            std::vector<double> &heights) const
{
    double last_time = 0;
    while (true) {
        auto pos = stm.tellg();
        if (pos_max != pos_error && pos >= pos_max)
            return last_time;
        double tsf;
        double val;
        double height;
        auto res = parseline(stm, &tsf, &val, &height);
        if (!unlikely(res)) {
            // Treat IO error as not possible to read forward
            if (!stm.good()) {
                stm.clear();
                return last_time;
            }
            continue;
        }
        last_time = tsf;
        if (val != 0) {
            times.push_back(tsf);
            datas.push_back(val);
            heights.push_back(height);
        }
        if (tsf >= tmax) {
            return last_time;
        }
    }
}

static constexpr double time_threshold = 120;
static constexpr std::streamoff pos_threshold = 10240;

// Find the location that is right before `tstart` and starts parsing there.
// Only look in `[pstart, pend)`.
// The starting time may be after `tstart` if no point before `tstart` exists in the range.
// The returned range may be larger than needed.
// If return value is `true`, `tsf` and `val` contains the number corresponds to the first point,
// and `loc` is the start of the corresponding line. Note that the data is only guaranteed
// to be the first valid line following `loc`, there could be any number of garbage lines
// before it.
// Always seek the stream.
NACS_INTERNAL bool Wavemeter::start_parse(std::istream &stm, double tstart,
                                          pos_type pstart, pos_type pend,
                                          double *tsf, double *val, double *height,
                                          pos_type *loc) const
{
    pos_type lb = pstart;
    pos_type ub = pend;
    if (ub == pos_error)
        ub = m_file_len;
    bool lb_valid = false;
    std::tuple<double,double,double,pos_type> lb_res{0, 0, 0, 0};
    while (true) {
        auto sz = ub - lb;
        // Good enough
        if (sz < pos_threshold) {
            stm.seekg(lb);
            std::tie(*tsf, *val, *height, *loc) = lb_res;
            return lb_valid;
        }
        pos_type mid = lb + sz / 2;
        double t, v, h;
        auto res = parse_at(stm, mid, lb, &t, &v, &h);
        if (unlikely(!res.first)) {
            while (!res.first && stm.good() && stm.tellg() < ub) {
                stm.clear();
                res.first = parseline(stm, &t, &v, &h);
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
            lb_res = {t, v, h, res.second};
        }
    }
}

// Always seek the stream.
NACS_INTERNAL void Wavemeter::extend_segment(std::istream &stm, seg_iterator seg,
                                             double tend, pos_type pend)
{
    if (seg->tend >= tend)
        return;
    stm.seekg(seg->pend);
    size_t old_size = seg->times.size();
    double last_time = parse_until(stm, tend, pend, seg->times, seg->datas, seg->heights);
    m_cache_size += seg->times.size() - old_size;
    if (last_time > 0)
        seg->tend = last_time;
    seg->pend = stm.tellg();
}

// Always seek the stream.
NACS_INTERNAL auto Wavemeter::new_segment(std::istream &stm, double tstart, double tend,
                                          pos_type lb, pos_type ub, seg_iterator prev)
    -> seg_iterator
{
    double tsf, val, height;
    pos_type loc;
    bool valid;
    try {
        valid = start_parse(stm, tstart, lb, ub, &tsf, &val, &height, &loc);
    }
    catch (...) {
        return m_segments.end();
    }
    if (!stm.good())
        return m_segments.end();
    if (!valid)
        loc = stm.tellg();
    if (loc == lb && prev != m_segments.end()) {
        if (valid && val != 0) {
            prev->times.push_back(tsf);
            prev->datas.push_back(val);
            prev->heights.push_back(height);
        }
        double last_time = parse_until(stm, tend, ub, prev->times, prev->datas, prev->heights);
        if (last_time > 0)
            prev->tend = last_time;
        prev->pend = stm.tellg();
        return prev;
    }
    std::vector<double> times;
    std::vector<double> datas;
    std::vector<double> heights;
    if (valid && val != 0) {
        times.push_back(tsf);
        datas.push_back(val);
        heights.push_back(height);
    }
    double last_time = parse_until(stm, tend, ub, times, datas, heights);
    if (times.empty())
        return m_segments.end();
    m_cache_size += times.size();
    return m_segments.emplace(loc, stm.tellg(), last_time,
                              std::move(times), std::move(datas), std::move(heights)).first;
}

// Always seek the stream.
NACS_INTERNAL auto Wavemeter::get_segment(std::istream &stm, double tstart,
                                          double tend) -> seg_iterator
{
    if (unlikely(m_segments.empty()))
        return new_segment(stm, tstart, tend, 0, pos_error);
    auto it = m_segments.end();
    --it;
    // Handle the most likely case first, i.e. reading from the end of file.
    if (it->times.front() <= tstart) {
        if (it->tend + time_threshold >= tstart) {
            extend_segment(stm, it, tend, pos_error);
            return it;
        }
        return new_segment(stm, tstart, tend, it->pend, pos_error, it);
    }

    // Now the generic case

    // First, determine which block to start.
    // We know that the cache isn't empty already.
    it = m_segments.upper_bound(tstart);
    if (it == m_segments.end()) {
        // This means that `lastit->times.front() > tstart`,
        // which disagrees with what we checked above.
        // Maybe this can happen if the time is actually not sorted
        // (due to daylight saving for example) which we are not handling very well now.
        // Return an error in this case for now instead of crashing...
        // Admittedly, when we are here we've probably already had UB so this is just a best
        // effort...
        return m_segments.end();
    }
    auto it2 = it;
    pos_type lb = 0;
    if (it2 != m_segments.begin()) {
        --it2;
        auto endt2 = it2->tend;
        if (endt2 >= tend)
            return it2; // `front() <= tstart` and `back() >= tend`
        if (endt2 + time_threshold >= tstart) {
            // `front() <= tstart <= back() + time_threshold` and `back() < tend`
            extend_segment(stm, it2, tend, it->pstart);
            it = it2;
            ++it2;
            goto segment_started;
        }
        lb = it2->pend;
    }
    else {
        it2 = m_segments.end();
    }
    // We now know that we aren't overlapping with/closed enough to any previous ones.
    {
        auto tmp = it;
        it = new_segment(stm, tstart, tend, lb, it->pstart, it2);
        it2 = tmp;
        if (it == m_segments.end()) {
            // `new_segment` may return `end()` if it finds no data before
            // the next segment starts, in this case, the next one is what we should use.
            it = tmp;
            ++it2;
        }
    }

segment_started:
    // `it` is the current segment,
    // `it2` is the next segement that might need merging.
    while (true) {
        if (it2 == m_segments.end()) {
            extend_segment(stm, it, tend, pos_error);
            break;
        }
        extend_segment(stm, it, tend, it2->pstart);
        // Gap big enough
        if (it2->pstart - it->pend > pos_threshold &&
            it2->times.front() - it->tend > time_threshold)
            break;
        if (it->pend != it2->pstart)
            extend_segment(stm, it, it2->tend, it2->pstart);
        // Merge the two segments
        it->pend = it2->pend;
        it->tend = it2->tend;
        it->times.insert(it->times.end(), it2->times.begin(), it2->times.end());
        it->datas.insert(it->datas.end(), it2->datas.begin(), it2->datas.end());
        auto tmp = it2;
        ++it2;
        m_segments.erase(tmp);
        if (it->tend >= tend) {
            break;
        }
    }

    return it;
}

// Always seek the stream.
NACS_INTERNAL void Wavemeter::check_cache(std::istream &stm)
{
    static constexpr int max_len = 128;
    stm.seekg(0, std::ios::end);
    if (!stm.good())
        throw std::runtime_error("Unable to seek stream.");
    pos_type len = stm.tellg();
    auto update = [&] {
        m_file_len = len;
        if (len > max_len) {
            stm.seekg(-max_len, std::ios::end);
            m_file_end.resize(max_len);
        }
        else {
            stm.seekg(0);
            m_file_end.resize((size_t)len);
        }
        stm.read(&m_file_end[0], m_file_end.size());
        if (!stm.good()) {
            m_file_end.resize(0);
            throw std::runtime_error("Unable to read from file end.");
        }
    };
    if (len < m_file_len) {
        clear();
    }
    else {
        stm.seekg(m_file_len - (off_t)m_file_end.size());
        char buff[max_len];
        assert(max_len >= m_file_end.size());
        stm.read(buff, m_file_end.size());
        if (!stm.good())
            throw std::runtime_error("Unable to read from file end.");
        if (memcmp(buff, m_file_end.data(), m_file_end.size()) != 0) {
            clear();
        }
    }
    update();
}

// Trigger GC based on the size and the time range of the cache
NACS_INTERNAL void Wavemeter::check_gc(double tstart, double tend)
{
    if (unlikely(m_segments.empty()))
        return;
    // Trigger cache if
    // 1. The cache size is larger than `8M` entries, and
    // 2. The time range of the cache is larger than 100 days and 10 x the read time.
    size_t sz_thresh = 8 * 1024 * 1024;
    if (likely(m_cache_size <= sz_thresh))
        return;
    double cache_begin_t = m_segments.begin()->times.front();
    double cache_end_t = m_segments.rbegin()->times.back();
    double t_thresh = 100 * max(86400, tstart - tend);
    if (likely(cache_end_t - cache_begin_t <= t_thresh))
        return;
    // Now drop everything that's too far from the current read range
    t_thresh /= 10;
    sz_thresh /= 10;
    size_t new_sz = 0;
    for (auto next = m_segments.begin(), end = m_segments.end(); next != end;) {
        auto it = next++;
        double t0 = it->times.front();
        double t1 = it->times.back();
        if (t1 >= tstart - t_thresh && t0 <= tend + t_thresh) {
            if (new_sz == 0)
                cache_begin_t = t0;
            cache_end_t = t1;
            // We could calculate the delta as well.
            // However, we are dropping more time period than keeping and recounting
            // can fix any miscalculation we have had (just in case I'm missing something)
            new_sz += it->times.size();
            continue;
        }
        m_segments.erase(it);
    }
    if (new_sz <= sz_thresh || cache_end_t - cache_begin_t <= t_thresh) {
        // We've dropped enough entries.
        m_cache_size = new_sz;
        return;
    }
    assert(!m_segments.empty());
    for (auto it = m_segments.begin(); it->times.back() < tstart;) {
        new_sz -= it->times.size();
        m_segments.erase(it);
        if (m_segments.empty()) {
            m_cache_size = 0;
            return;
        }
        it = m_segments.begin();
        cache_begin_t = it->times.front();
    }
    if (new_sz <= sz_thresh || cache_end_t - cache_begin_t <= t_thresh) {
        m_cache_size = new_sz;
        return;
    }
    // We use `--m_segments.end()` instead of `m_segments.rbegin()` since
    // `erase` does not accept reverse iterator.
    for (auto it = --m_segments.end(); it->times.front() > tend;) {
        new_sz -= it->times.size();
        m_segments.erase(it);
        if (m_segments.empty()) {
            m_cache_size = 0;
            return;
        }
        it = --m_segments.end();
        cache_end_t = it->times.back();
    }
    // The check for whether we've dropped enough for this one is merged into the loop below.

    // Currently the location where each line is located isn't cached so we can't
    // drop partial segment and have to drop a whole segment at a time.
    // We prefer dropping from the beginning since the end is more likely to be used later
    // so we'll just drop segments from the beginning until we've dropped enough.
    for (auto next = m_segments.begin(), end = m_segments.end(); next != end;) {
        auto it = next++;
        cache_begin_t = it->times.front();
        if (new_sz <= sz_thresh || cache_end_t - cache_begin_t <= t_thresh) {
            // We've dropped enough entries.
            m_cache_size = new_sz;
            return;
        }
        new_sz -= it->times.size();
        m_segments.erase(it);
    }
    // All dropped.
    m_cache_size = 0;
}

// Always seek the stream.
NACS_EXPORT() std::tuple<const double*,const double*,const double*>
Wavemeter::parse(std::istream &stm, size_t *sz, double tstart, double tend)
{
    check_cache(stm);
    if (unlikely(tend <= tstart)) {
        *sz = 0;
        return {nullptr, nullptr, nullptr};
    }
    check_gc(tstart, tend);
    auto seg = get_segment(stm, tstart, tend);
    if (seg == m_segments.end())
        return {nullptr, nullptr, nullptr};
    auto it1 = std::lower_bound(seg->times.begin(), seg->times.end(), tstart);
    if (it1 == seg->times.end())
        return {nullptr, nullptr, nullptr};
    auto idx1 = it1 - seg->times.begin();
    auto it2 = std::lower_bound(it1, seg->times.end(), tend);
    auto idx2 = it2 == seg->times.end() ? seg->times.size() - 1 : it2 - seg->times.begin();
    *sz = idx2 - idx1;
    return {&*it1, &seg->datas[idx1], &seg->heights[idx1]};
}

NACS_EXPORT_ Wavemeter::Wavemeter(double lo, double hi)
    : m_lo(lo), m_hi(hi)
{
}

NACS_EXPORT() void Wavemeter::clear()
{
    m_segments.clear();
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
                                                const double **height,
                                                double tstart, double tend)
{
    auto parser = (Wavemeter*)_parser;
    size_t sz = 0;
    std::ifstream stm(name);
    std::tie(*ts, *data, *height) = parser->parse(stm, &sz, tstart, tend);
    return sz;
}

NACS_EXPORT() void nacs_utils_wavemeter_clear(void *_parser)
{
    auto parser = (Wavemeter*)_parser;
    parser->clear();
}

NACS_EXPORT() void nacs_utils_free_wavemeter(void *parser)
{
    delete (Wavemeter*)parser;
}

}
