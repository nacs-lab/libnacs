/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_SEQ_PULSE_H__
#define __NACS_SEQ_PULSE_H__

#include "event_time.h"

namespace NaCs::Seq {

class BasicSeq;

class Pulse {
public:
    // id must be unique in the whole sequence
    // A pulse with the same time but a larger ID will be treated as happens after
    // one with a smaller ID.
    Pulse(uint32_t id, EventTime &start, Var *len, Var *val, Var *cond, bool is_measure);
    uint32_t id() const
    {
        return m_id;
    }
    bool is_measure() const
    {
        return m_is_measure;
    }
    // Start time of the pulse.
    // For pulses of a single value (instead of a ramp) or a measurement,
    // this is the time of the output/measure.
    const EventTime &start() const
    {
        return *m_start;
    }
    // May be null if the pulse contains no ramp or is a measurement.
    Var *len() const
    {
        return m_len.get();
    }
    Var *cond() const
    {
        return m_cond.get();
    }
    // Note that `val()` and `endval()` are only for when the pulse is enabled.
    // If `cond()` evaluates to false, their values must be ignored.

    // `val()` is a function of time and old value for output,
    // and the slot that represent the result of the measure for measurement pulses.
    // When the optimizater is able to determine the exact expression for a pulse,
    // the value should be replaced with the corresponding expression
    // and the measurement pulse should be removed from the sequence.
    Var *val() const
    {
        return m_val.get();
    }
    Var *endval() const
    {
        assert(!m_is_measure);
        if (!m_endval)
            m_endval.reset(compute_endval());
        return m_endval.get();
    }

    bool optimize();
    bool needs_oldval() const;
    bool clear_unused_args();
    // Assumes `needs_oldval()`
    bool set_oldval(Var *oldval);
    void check_start() const;

    void print(std::ostream &stm, Var *startval=nullptr, bool newline=false) const;

private:
    Pulse(Pulse&&) = delete;
    Var *compute_endval() const;

    const uint32_t m_id;
    const bool m_is_measure;
public:
    // Property set by BasicSeq. Only valid for output (i.e. non-measure) pulse.
    // A ordered pulse means all other pulses on the same channel
    // are either known before or after this one.
    bool is_ordered = false;
private:
    EventTime::Ref m_start;
    Var::Ref m_len;
    Var::Ref m_val;
    Var::Ref m_cond;

    mutable Var::Ref m_endval;
    friend class BasicSeq;
};

static inline std::ostream &operator<<(std::ostream &stm, const Pulse &p)
{
    p.print(stm);
    return stm;
}

}

#endif
