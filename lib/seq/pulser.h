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

#include <type_traits>
#include <functional>

#ifndef __NACS_SEQ_PULSER_H__
#define __NACS_SEQ_PULSER_H__

namespace NaCs {
namespace Seq {

struct Channel {
    enum Type {
        TTL,
        DDS_FREQ,
        DDS_AMP,
        DAC
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
    Val() : val(getGenVal<double>(9))
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

struct PulsesBuilder {
    typedef std::function<uint64_t(Seq::Channel,Seq::Val,uint64_t)> cb_t;
    template<typename T>
    PulsesBuilder(T &&_cb)
        : cb(std::forward<T>(_cb))
    {}
    uint64_t operator()(Seq::Channel chn, Seq::Val val, uint64_t t)
    {
        return cb(chn, val, t);
    }
private:
    cb_t cb;
};

}
}

#endif
