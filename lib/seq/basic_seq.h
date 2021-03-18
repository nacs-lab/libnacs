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

#ifndef __NACS_SEQ_BASIC_SEQ_H__
#define __NACS_SEQ_BASIC_SEQ_H__

#include "pulse.h"

#include <list>
#include <map>

namespace NaCs::Seq {

class Seq;

class BasicSeq {
    struct ChannelInfo {
        std::list<Pulse> pulses;
        Var::Ref startval;
        Var::Ref endval;
        std::map<const Pulse*,Var::Ref> pulse_start_vars;
    };
public:
    enum Errno : uint8_t {
        ExternMeasure,
        ExternMeasureLength,
        MeasureOrder,
    };
    struct Assignment {
        Var::Ref val;
        // ID is for error reporting only
        uint32_t id;
    };
    struct Assumption {
        EventTime::Sign sign;
        // ID is for error reporting only
        // This ID is in the same namespace as EventTime and will be reported as so.
        uint32_t id = 0;
    };
    struct Branch {
        Var::Ref cond;
        BasicSeq *target;
        // ID is for error reporting only
        uint32_t id;
    };
    BasicSeq(uint32_t id);
    uint32_t id() const
    {
        return m_id;
    }
    Pulse *add_pulse(uint32_t chn, uint32_t id, EventTime &start, Var *len, Var *val);
    Pulse *add_measure(uint32_t chn, uint32_t id, EventTime &start, Var *val);
    Var *new_measure(Env &env, uint32_t measure_id) const;
    // NULL target means termination of sequence
    void add_branch(Var *cond, BasicSeq *target, uint32_t id);
    // NULL target means termination of sequence
    void set_default_branch(BasicSeq *target);
    Var *endval(uint32_t chn) const
    {
        auto it = m_channels.find(chn);
        if (it == m_channels.end())
            return nullptr;
        return it->second.endval.get();
    }
    Var *startval(uint32_t chn) const
    {
        auto it = m_channels.find(chn);
        if (it == m_channels.end())
            return nullptr;
        return it->second.startval.get();
    }
    const std::list<Pulse> *get_pulses(uint32_t chn) const
    {
        auto it = m_channels.find(chn);
        if (it == m_channels.end())
            return nullptr;
        return &it->second.pulses;
    }
    const std::map<const Pulse*,Var::Ref> *get_pulse_start_vars(uint32_t chn) const
    {
        auto it = m_channels.find(chn);
        if (it == m_channels.end())
            return nullptr;
        return &it->second.pulse_start_vars;
    }
    const std::map<uint32_t,Assignment> &get_assigns() const
    {
        return m_assign;
    }
    // Whether the assumption triggers or now for an unused basic seq is implementation defined.
    const auto &get_assumes() const
    {
        return m_assume;
    }
    llvm::ArrayRef<Branch> get_branches() const
    {
        return m_branches;
    }
    BasicSeq *get_default_branch() const
    {
        return m_default_branch;
    }
    const std::list<EventTime::Ref> &get_endtimes() const
    {
        return m_endtimes;
    }
    const std::list<EventTime> &get_eventtimes() const
    {
        return m_eventtimes;
    }
    const std::map<uint32_t,ChannelInfo> &get_channel_infos() const
    {
        return m_channels;
    }

    bool has_output(uint32_t chn) const;
    void assign_global(uint32_t global_id, Var *val, uint32_t assignment_id);
    void add_assume(EventTime::Sign sign, Var *val, uint32_t assume_id);
    void add_endtime(EventTime &t);
    EventTime &track_time(const EventTime &t);

    void prepare(Seq &seq);
    void print(std::ostream &stm) const;

    // For testing only
    bool needs_oldval(Pulse *pulse) const;

private:
    BasicSeq(BasicSeq&&) = delete;
    void mark_recursive();
    bool preoptimize_pulse(uint32_t chn, Env &env);
    bool optimize_order(uint32_t chn);
    bool optimize_endtimes();
    bool optimize_vars(Seq &seq);
    bool optimize_branch();
    bool preoptimize_eventtimes();
    bool postoptimize_eventtimes();
    bool optimize_final(Seq &seq);

    Var *alloc_startval(Env &env)
    {
        uint64_t start_id = (uint64_t(m_id) << 32) | uint32_t(--m_startval_count);
        assert((start_id >> 32) == m_id);
        return env.new_extern({IR::Type::Float64, start_id});
    };

    const uint32_t m_id;
    bool m_used = false; // Used by GC
    std::map<uint32_t,ChannelInfo> m_channels;
    std::map<uint32_t,Assignment> m_assign;
    std::map<Var::Ref,Assumption> m_assume;
    // We use a different list to keep track of the sequence times
    // since in the real sequence there are likely many duplicates
    // between times on different channels
    // and most of the times will be known to be smaller than other ones.
    // Using a different list allows us to not reprocess all of those
    // in each optimization cycle.
    std::list<EventTime::Ref> m_endtimes;
    // All the times used in the sequence
    std::list<EventTime> m_eventtimes;
    std::vector<Var::Ref> m_global_escapes;

    std::vector<Branch> m_branches;
    BasicSeq *m_default_branch = nullptr;

    int32_t m_startval_count = 0;

    friend class Seq;
};

static inline std::ostream &operator<<(std::ostream &stm, const BasicSeq &seq)
{
    seq.print(stm);
    return stm;
}

}

#endif
