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
#include <nacs-seq/seq.h>

#include <type_traits>
#include <functional>
#include <memory>

#ifndef __NACS_SEQ_PULSER_H__
#define __NACS_SEQ_PULSER_H__

namespace NaCs {
namespace Seq {

struct Channel {
    enum Type {
        TTL = 1,
        DDS_FREQ = 2,
        DDS_AMP = 3,
        DAC = 4
    };
    Type typ;
    int id;
};

bool operator<(const Channel &id1, const Channel id2)
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
    typedef std::function<Val(uint64_t,Val,uint64_t)> func_t;
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
    template<typename T>
    PulseData(T &&v)
        : m_val(),
          m_cb(new func_t(std::forward<T>(v)))
    {}
    Val operator()(uint64_t t, Val start, uint64_t len) const
    {
        if (m_cb)
            return (*m_cb)(t, start, len);
        return m_val;
    }
private:
    Val m_val;
    mutable std::unique_ptr<func_t> m_cb;
};

typedef BasePulse<Channel, PulseData> Pulse;

struct Filter {
    bool operator()(Channel)
    {
        return true;
    }
    bool operator()(Channel cid, Val val1, Val val2)
    {
        if (cid.typ == Channel::TTL)
            return val1.val.i32 != val2.val.i32;
        return val1.val.f64 != val2.val.f64;
    }
};

struct PulsesBuilder {
    typedef std::function<uint64_t(Channel,Val,uint64_t)> cb_t;
    typedef std::function<uint64_t(PulsesBuilder&,uint64_t,Event)> seq_cb_t;
    template<typename T>
    PulsesBuilder(T &&_cb)
        : cb(std::forward<T>(_cb))
    {}
    uint64_t operator()(Channel chn, Val val, uint64_t t)
    {
        return cb(chn, val, t);
    }
    static std::pair<std::vector<Pulse>,std::map<Channel,Val>>
    fromBase64(const uint8_t *data, size_t len);
    void schedule(std::vector<Pulse> &seq,
                  const std::map<Channel,Val> &defaults,
                  seq_cb_t seq_cb,
                  Time::Constraints t_cons={100, 40, 4096});
private:
    cb_t cb;
};

}
}

#endif
