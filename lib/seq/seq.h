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

#include "timing.h"

#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <stdint.h>

#ifndef __NACS_SEQ_SEQ_H__
#define __NACS_SEQ_SEQ_H__

namespace NaCs {
namespace Seq {

template<typename Cid, typename Cb>
struct BasePulse {
    uint64_t t;
    uint64_t len;
    Cid chn;
    // Called with `cb(uint64_t t, ValT start, uint64_t len) -> ValT`
    // where `t` is the time within the pulse,
    // `start` is the value of the channel at the begining of the pulse.
    // `len` is the length pulse.
    Cb cb;
};

template<typename Cid, typename Cb>
static void sort(std::vector<BasePulse<Cid,Cb>> &seq)
{
    std::stable_sort(seq.begin(), seq.end(),
                     [] (auto &p1, auto &p2) { return p1.t < p2.t; });
}

enum class Event {
    start,
    end
};

struct Clock {
    // Time index of the first clock edge. The clock pulse should happen `div` time points
    // before this time.
    uint64_t t;
    // Length of the pulse in time index. The clock end pulse should be `len` time points after
    // the start pulse.
    uint64_t len;
    // Clock division. Each clock cycle consists of `div` number of time points at high
    // and `div` number of time points at low.
    uint32_t div;
};

}
}

#endif
