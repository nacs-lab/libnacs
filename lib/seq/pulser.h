/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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
#include <nacs-seq/seq.h>

#include <type_traits>
#include <functional>
#include <memory>

#include <math.h>

#ifndef __NACS_SEQ_PULSER_H__
#define __NACS_SEQ_PULSER_H__

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

typedef BasePulse<Channel, PulseData> Pulse;

struct Filter {
    bool operator()(Channel)
    {
        return true;
    }
    bool operator()(Channel cid, Val val1, Val val2)
    {
        switch (cid.typ) {
        case Channel::TTL:
            return val1.val.i32 != val2.val.i32;
        case Channel::DDS_FREQ:
            return fabs(val1.val.f64 - val2.val.f64) > 0.2;
        case Channel::DDS_AMP:
            return fabs(val1.val.f64 - val2.val.f64) > 0.00005;
        case Channel::DAC:
            return fabs(val1.val.f64 - val2.val.f64) > 0.0002;
        default:
            return val1.val.f64 != val2.val.f64;
        }
    }
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

struct PulsesBuilder {
    typedef std::function<uint64_t(Channel,Val,uint64_t,uint64_t)> cb_t;
    typedef std::function<uint64_t(PulsesBuilder&,uint64_t,Event)> seq_cb_t;
    template<typename T>
    PulsesBuilder(T &&_cb)
        : cb(std::forward<T>(_cb))
    {}
    uint64_t operator()(Channel chn, Val val, uint64_t t, uint64_t tlim)
    {
        return cb(chn, val, t, tlim);
    }
    __attribute__((deprecated)) static Sequence fromBase64(const uint8_t *data, size_t len);
    // In `ttl_mask_out`, when not NULL, will return the mask of all TTL channels
    // used by the sequence. It is guaranteed that the slot is filled before
    // actual scheduling starts or the start callback is invoked.
    void schedule(Sequence&, seq_cb_t seq_cb, uint32_t *ttl_mask_out=nullptr,
                  Time::Constraints t_cons={100, 40, 4096});
    void schedule(Sequence &&seq, seq_cb_t seq_cb, uint32_t *ttl_mask_out=nullptr,
                  Time::Constraints t_cons={100, 40, 4096})
    {
        schedule(seq, seq_cb, ttl_mask_out, t_cons);
    }
private:
    cb_t cb;
};

}
}

#endif
