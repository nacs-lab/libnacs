/*************************************************************************
 *   Copyright (c) 2016 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <stdint.h>

#ifndef __NACS_SEQ_ZYNQ_LEGACY_SEQ_H__
#define __NACS_SEQ_ZYNQ_LEGACY_SEQ_H__

namespace NaCs::Seq::Zynq {

struct LegacySeq {
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
        friend inline bool operator<(const Channel &id1, const Channel id2)
        {
            if (id1.typ < id2.typ) {
                return true;
            } else if (id1.typ > id2.typ) {
                return false;
            }
            return id1.id < id2.id;
        }
    };
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
        static inline std::enable_if_t<std::is_same_v<T,double>, IR::GenVal>
        getGenVal(T2 val)
        {
            IR::GenVal gv;
            gv.f64 = val;
            return gv;
        }
        template<typename T, typename T2>
        static inline std::enable_if_t<std::is_same_v<T,uint32_t>, IR::GenVal>
        getGenVal(T2 val)
        {
            IR::GenVal gv;
            gv.i32 = val;
            return gv;
        }
    };
    struct Pulse {
        typedef IR::ExeContext::Func<double(double,double)> func_t;
        uint64_t t;
        uint64_t len;
        Channel chn;
        Pulse(Pulse&&) = default;
        Pulse &operator=(Pulse&&) = default;
        Pulse()
        {}
        Pulse(uint64_t t, uint64_t len, Channel chn, Val val)
            : t(t), len(len), chn(chn), m_val(val), m_cb()
        {
        }
        Pulse(uint64_t t, uint64_t len, Channel chn, func_t func)
            : t(t), len(len), chn(chn), m_val(), m_cb(std::move(func))
        {
        }
        Val operator()(uint64_t t) const
        {
            // `t` is the time within the pulse,
            // `start` is the value of the channel at the begining of the pulse.
            if (m_cb)
                return IR::TagVal(m_cb(double(t) * 10e-9, m_val.val.f64)).val;
            return m_val;
        }
        void set(Val v)
        {
            m_val = v;
            m_cb = func_t();
        }
        void set_start(Val start)
        {
            if (m_cb) {
                m_val = start;
            }
        }
    private:
        Val m_val; // this is the start value if m_cb is not null.
        func_t m_cb;
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
    std::vector<Pulse> pulses;
    std::map<Channel,Val> defaults;
    std::vector<Clock> clocks;
    std::unique_ptr<IR::ExeContext> exectx;
    LegacySeq(std::vector<Pulse> &&_pulses, std::map<Channel,Val> &&_defaults,
              std::vector<Clock> &&_clocks={},
              std::unique_ptr<IR::ExeContext> _exectx=IR::ExeContext::get())
        : pulses(std::move(_pulses)),
          defaults(std::move(_defaults)),
          clocks(std::move(_clocks)),
          exectx(std::move(_exectx))
    {
    }
    static LegacySeq fromBinary(const uint32_t *data, size_t len);
    static void dumpBinary(std::ostream &stm, const uint32_t *data, size_t len);
    std::vector<uint8_t> toByteCode(uint32_t *ttl_mask);
    uint8_t *toByteCode(size_t *sz, uint32_t *ttl_mask);
};

}

extern "C" {

uint8_t *nacs_seq_bin_to_bytecode(const uint32_t *data, size_t data_len,
                                  size_t *code_len, uint32_t *ttl_mask);
uint64_t nacs_seq_bytecode_total_time(const uint8_t *code, size_t code_len);

}

#endif
