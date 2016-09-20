/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include <stdint.h>
#include <stdlib.h>
#include <vector>

#ifndef __NACS_SEQ_TIMING_H__
#define __NACS_SEQ_TIMING_H__

namespace NaCs {
namespace Seq {
namespace Time {

struct Constraints {
    // preferred pulse spacing when outputting contineous updates
    uint64_t prefer_dt;
    // minimum average pulse spacing and the average window
    uint64_t avg_dt;
    size_t avg_window;
};

struct Keeper {
public:
    Keeper(const Constraints &);
    void addPulse(uint64_t dt);
    uint64_t minDt(uint64_t min_dt) const;
    uint64_t minDt(void) const
    {
        return minDt(m_cons.prefer_dt);
    }
    void reset(void);
private:
    Constraints m_cons;
    size_t m_npulses;
    uint64_t m_total_len;
    std::vector<uint64_t> m_plens;
};

}
}
}

#endif
