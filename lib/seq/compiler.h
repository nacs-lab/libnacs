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

#ifndef __NACS_SEQ_COMPILER_H__
#define __NACS_SEQ_COMPILER_H__

#include "basic_seq.h"

#include <functional>
#include <map>
#include <vector>

#include <assert.h>

namespace NaCs::LLVM::Exe {
class Engine;
}

namespace NaCs::Seq {

class Seq;
class HostSeq;

/**
 * Compile a `Seq` to a `HostSeq`
 * Keep some information about the mapping for the backends to access.
 */
class Compiler {
    enum class VarType : uint8_t {
        None,
        Const, // Includes ramps
        Global,
        GlobalValue,
        Channel,
        BSeqMeasure,
        BSeqDirect,
        BSeqNeedOrder,
    };
    struct TimeVal;
    using val_eval_t = void (*)(void*);
    using ramp_eval_t = double (*)(double, void*);
    using assume_entry_t = std::pair<const Var::Ref,BasicSeq::Assumption>;

    template<typename Map, typename Key>
    static inline decltype(auto) find_val(const Map &map, Key &&key)
    {
        auto it = map.find(std::forward<Key>(key));
        assert(it != map.end());
        return it->second;
    }

public:
    Compiler(HostSeq &host_seq, Seq &seq, LLVM::Exe::Engine &engine,
             std::function<uintptr_t(const std::string&)> resolver);
    ~Compiler();

    uint64_t compile();
    uint32_t var_slot(const Var *var) const
    {
        return m_var_slots[var->varid()];
    }
    uint32_t time_slot(const EventTime *time) const
    {
        return find_val(m_time_slots, time);
    }

private:
    void prescan_vars();
    void prescan_times();
    void classify_times();
    void assign_shared_vars();
    void assign_shared_times();
    void assign_defaults();
    void create_bseqs();
    void assign_bseq_vars();
    void cache_time_slots();
    void populate_bseqs();
    val_eval_t *find_fptr_slot(uint32_t id, bool is_var) const;
    void generate_fptrs();

    HostSeq &m_host_seq;
    Seq &m_seq;

    uint32_t nglobal_vals_var;
    uint32_t nglobal_vals_time;

    uint32_t global_offset;
    uint32_t global_value_offset;
    uint32_t channels_offset;

    LLVM::Exe::Engine &m_engine;
    std::function<uintptr_t(const std::string&)> m_resolver;

    std::vector<uint32_t> m_var_slots;
    std::map<const EventTime*,uint32_t> m_time_slots;
    std::vector<Var*> m_all_vars;
    std::vector<VarType> m_var_types;
    std::vector<uint32_t> m_var_owners;
    // In principle a single ramp function could be used in multiple basic sequences and pulses
    std::map<uint32_t,std::vector<std::pair<uint32_t,uint32_t>>> m_var_bseqidx_pulseid;
    std::vector<bool> m_var_isramps;
    // Dependency of global_vals on each other.
    // The runtime doesn't care about this but we need this to compute
    // the complete set of global_vals used by each basic sequence.
    std::vector<std::vector<uint32_t>> m_global_vals_inter_deps;

    std::vector<const EventTime*> m_event_times;
    std::vector<const assume_entry_t*> m_assumptions;
    std::map<const EventTime*,uint32_t> m_event_time_ids;
    std::map<const assume_entry_t*,uint32_t> m_assumption_ids;

    std::vector<TimeVal> m_timevals;
    // Mapping from the index of `EventTime` and `Assumption`s to index in `m_timevals`.
    std::vector<uint32_t> m_time_timevals;
    std::vector<VarType> m_timeval_types;
    std::vector<uint32_t> m_timeval_owners;
    std::vector<uint32_t> m_timeval_slots;

    std::map<uint32_t,uint32_t> m_bseq_idxs;
    std::vector<const BasicSeq*> m_bseqs;
    uint64_t m_obj_id = 0;
};

}

#endif
