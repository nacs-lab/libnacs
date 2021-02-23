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

class BasicSeq {
public:
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
    Pulse *add_pulse(uint32_t chn, uint32_t id, EventTime &&start, Var *len, Var *val);
    Pulse *add_measure(uint32_t chn, uint32_t id, EventTime &&start, Var *val);
    Var *new_measure(Env &env, uint32_t measure_id) const;
    // NULL target means termination of sequence
    void add_branch(Var *cond, BasicSeq *target, uint32_t id);
    // NULL target means termination of sequence
    void set_default_branch(BasicSeq *target);
    Var *endval(uint32_t chn) const
    {
        auto it = m_endval.find(chn);
        if (it == m_endval.end())
            return nullptr;
        return it->second.get();
    }
    Var *startval(uint32_t chn) const
    {
        auto it = m_startval.find(chn);
        if (it == m_startval.end())
            return nullptr;
        return it->second.get();
    }
    const std::list<Pulse> *get_pulses(uint32_t chn) const
    {
        auto it = m_pulses.find(chn);
        if (it == m_pulses.end())
            return nullptr;
        return &it->second;
    }
    llvm::ArrayRef<Branch> get_branches() const
    {
        return m_branches;
    }
    BasicSeq *get_default_branch() const
    {
        return m_default_branch;
    }
    bool has_output(uint32_t chn) const;

private:
    BasicSeq(BasicSeq&&) = delete;
    void mark_recursive();
    bool optimize_pulse(uint32_t chn);
    bool optimize_order(uint32_t chn);
    void optimize_vars();
    bool optimize_branch();

    const uint32_t m_id;
    bool m_used = false; // Used by GC
    std::map<uint32_t,std::list<Pulse>> m_pulses;
    // The maps below is guaranteed to not have empty Var::Ref
    std::map<uint32_t,Var::Ref> m_startval;
    std::map<uint32_t,Var::Ref> m_endval;

    std::vector<Branch> m_branches;
    BasicSeq *m_default_branch = nullptr;
};

}

#endif
