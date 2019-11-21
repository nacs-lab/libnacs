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

#ifndef __NACS_UTILS_WAVEMETER_H__
#define __NACS_UTILS_WAVEMETER_H__

#include "utils.h"

#include <istream>
#include <set>
#include <vector>

#include <assert.h>

namespace NaCs {

class Wavemeter {
    using pos_type = std::istream::pos_type;
    struct Segment {
        pos_type pstart;
        mutable pos_type pend;
        // These two should never be empty.
        mutable std::vector<double> times;
        mutable std::vector<double> datas;
        Segment(pos_type pstart, pos_type pend,
                std::vector<double> times, std::vector<double> datas)
            : pstart(pstart),
              pend(pend),
              times(std::move(times)),
              datas(std::move(datas))
        {
            assert(!this->times.empty());
        }
    };
    struct SegComp {
        using is_transparent = void;
        bool operator()(const pos_type &pos1, const Segment &seg2) const
        {
            return pos1 < seg2.pstart;
        }
        bool operator()(const Segment &seg1, const pos_type &pos2) const
        {
            return seg1.pstart < pos2;
        }
        bool operator()(const double &t1, const Segment &seg2) const
        {
            return t1 < seg2.times.front();
        }
        bool operator()(const Segment &seg1, const double &t2) const
        {
            return seg1.times.front() < t2;
        }
        bool operator()(const Segment &seg1, const Segment &seg2) const
        {
            return seg1.pstart < seg2.pstart;
        }
    };
    using seg_map_t = std::set<Segment,SegComp>;
    using seg_iterator = seg_map_t::iterator;

    // Stateless parsing functions
    // Parse the time stamp
    static bool parsetime(std::istream &stm, double *tsf);
    // Parse a single number
    static bool parsenumber(std::istream &stm, double *val, bool *eol);
    static inline bool parsenumber(std::istream &stm, double *val)
    {
        bool eol;
        return parsenumber(stm, val, &eol);
    }
    // Parse the data section of the log file and find the best match
    // given the upper and lower bound.
    // The line must be terminated by a `\n`.
    static bool parseval(std::istream &stm, double *val, double lo, double hi);
    // Parse both the time and the data from a complete line.
    // A `good()` `stm` after the function return guarantees forward progress.
    bool parseline(std::istream &stm, double *tsf, double *val) const;
    // Find the beginning of the line that includes `ub`
    // Do not look back more than `lb`.
    static pos_type find_linestart(std::istream &stm, pos_type ub, pos_type lb);
    // Parse the line that includes `pos`. Do not look back more than `lb`
    std::pair<bool,pos_type> parse_at(std::istream &stm, pos_type pos, pos_type lb,
                                      double *tsf, double *val) const;
    // Parse until the time and position limit.
    void parse_until(std::istream &stm, double tmax, pos_type pos_max,
                     std::vector<double> &times, std::vector<double> &datas) const;
    bool start_parse(std::istream &stm, double tstart, pos_type pstart, pos_type pend,
                     double *tsf, double *val, pos_type *loc) const;

    void extend_segment(std::istream &stm, seg_iterator seg, double tend, pos_type pend);
    // If `prev` is not NULL, it's a segment that ends at `lb`.
    seg_iterator new_segment(std::istream &stm, double tstart, double tend,
                             pos_type lb, pos_type ub, seg_iterator prev);
    seg_iterator new_segment(std::istream &stm, double tstart, double tend,
                             pos_type lb, pos_type ub)
    {
        return new_segment(stm, tstart, tend, lb, ub, m_segments.end());
    }
    // Parse and cache the result for a block.
    seg_iterator get_segment(std::istream &stm, double tstart, double tend);

    // Check if the cache is still valid and clear it if not.
    // Also update the `m_file_len` and `m_file_end` fields.
    void check_cache(std::istream &stm);

    // TODO: GC of cache

public:
    NACS_EXPORT(utils) Wavemeter(double lo, double hi);
    NACS_EXPORT(utils) std::pair<const double*,const double*>
    parse(std::istream &stm, size_t *sz, double tstart, double tend);
    NACS_EXPORT(utils) void clear();

private:
    seg_map_t m_segments;

    const double m_lo = 0;
    const double m_hi = 0;

    // The file length we got last time and the corresponding (up to) last 100 characters.
    // If the file is shorter than the length or if the content doesn't match
    // we'll assume the cache is invalide and clear it.
    pos_type m_file_len = 0;
    std::string m_file_end;
};

}

#endif // __NACS_UTILS_WAVEMETER_H__
