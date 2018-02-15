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

#include <nacs-utils/ir.h>

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

struct Channel {
    enum Type {
        TTL = 1,
        DDS_FREQ = 2,
        DDS_AMP = 3,
        DAC = 4,
        CLOCK = 5
    };
    Type typ;
    int id;
};

static inline bool operator<(const Channel &id1, const Channel id2)
{
    if (id1.typ < id2.typ) {
        return true;
    } else if (id1.typ > id2.typ) {
        return false;
    }
    return id1.id < id2.id;
}

struct Val {
    IR::GenVal val;
    Val(IR::GenVal v=getGenVal<double>(0)) : val(v)
    {}
    template<typename T, typename T2> static inline Val get(T2 val=0)
    {
        Val v;
        v.val = getGenVal<T>(val);
        return v;
    }
    template<typename T> static inline Val get()
    {
        return get<T>(0);
    }
private:
    template<typename T, typename T2>
    static inline std::enable_if_t<std::is_same<T,double>::value, IR::GenVal>
    getGenVal(T2 val)
    {
        IR::GenVal gv;
        gv.f64 = val;
        return gv;
    }
    template<typename T, typename T2>
    static inline std::enable_if_t<std::is_same<T,uint32_t>::value, IR::GenVal>
    getGenVal(T2 val)
    {
        IR::GenVal gv;
        gv.i32 = val;
        return gv;
    }
};

struct PulseData {
    typedef IR::ExeContext::Func<double(double,double)> func_t;
    PulseData(PulseData&&) = default;
    PulseData &operator=(PulseData&&) = default;
    PulseData()
        : m_val(),
          m_cb()
    {}
    PulseData(const Val &val)
        : m_val(val),
          m_cb()
    {}
    PulseData(Val &&val)
        : m_val(val),
          m_cb()
    {}
    PulseData(func_t func)
        : m_val(),
          m_cb(std::move(func))
    {}
    template<typename T>
    PulseData(T &&v)
        : PulseData(func_t(std::forward<T>(v)))
    {
    }
    Val operator()(uint64_t t, Val start) const
    {
        if (m_cb)
            return IR::TagVal(m_cb(double(t) * 10e-9, start.val.f64)).val;
        return m_val;
    }
private:
    Val m_val;
    mutable func_t m_cb;
};

struct Pulse {
    uint64_t t;
    uint64_t len;
    Channel chn;
    // Called with `cb(uint64_t t, ValT start, uint64_t len) -> ValT`
    // where `t` is the time within the pulse,
    // `start` is the value of the channel at the begining of the pulse.
    // `len` is the length pulse.
    PulseData cb;
};

static inline void sort(std::vector<Pulse> &seq)
{
    std::stable_sort(seq.begin(), seq.end(),
                     [] (auto &p1, auto &p2) { return p1.t < p2.t; });
}

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

struct Sequence {
    std::vector<Pulse> pulses;
    std::map<Channel,Val> defaults;
    std::vector<Clock> clocks;
    std::unique_ptr<IR::ExeContext> exectx;
    Sequence(std::vector<Pulse> &&_pulses, std::map<Channel,Val> &&_defaults,
             std::vector<Clock> &&_clocks={},
             std::unique_ptr<IR::ExeContext> _exectx=IR::ExeContext::get())
        : pulses(std::move(_pulses)),
          defaults(std::move(_defaults)),
          clocks(std::move(_clocks)),
          exectx(std::move(_exectx))
    {
    }
    static Sequence fromBase64(const uint8_t *data, size_t len);
    static Sequence fromBinary(const uint32_t *data, size_t len);
    std::vector<uint8_t> toByteCode(uint32_t *ttl_mask);
    uint8_t *toByteCode(size_t *sz, uint32_t *ttl_mask);
};

}
}

#endif
