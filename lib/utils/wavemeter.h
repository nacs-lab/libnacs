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

#include <istream>
#include <map>
#include <vector>

namespace NaCs {

class Wavemeter {
    using pos_type = std::istream::pos_type;
    struct PosRange {
        double tend;
        pos_type pstart;
        pos_type pend;
    };
    struct Segment {
        std::vector<double> times{};
        std::vector<double> datas{};
        size_t size = 0;
    };

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
    // given the upper and lower bound
    static bool parseval(std::istream &stm, double *val, double lo, double hi);
    // Parse both the time and the data from a complete line.
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
    bool start_parse(std::istream &stm, double tstart,
                     pos_type pstart, pos_type pend, double *tsf, double *val) const;

    // Time -> position
    std::pair<pos_type,pos_type> find_pos_range(double t) const;
    void add_pos_range(double tstart, double tend, pos_type pstart, pos_type pend);

    // Parse and cache the result for a block.
    const Segment *get_segment(std::istream &stm, double tstart, double tend);

public:
    Wavemeter(double lo, double hi);
    std::pair<const double*,const double*>
    parse(std::istream &stm, size_t *sz, double tstart, double tend);

private:
    std::map<double,PosRange> m_pos_cache;
    std::map<pos_type,Segment> m_segments;

    const double m_lo = 0;
    const double m_hi = 0;
};

}

#endif // __NACS_UTILS_WAVEMETER_H__
