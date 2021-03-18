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

#include "basic_seq.h"

#include "error.h"
#include "seq.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/number.h"

#include <algorithm>
#include <set>
#include <utility>

namespace NaCs::Seq {

struct BasicSeq::TimeInfo {
    EventTime *time;
    int32_t ordered_pos = -1;
    int32_t ordered_nonneg = -1;
};

struct BasicSeq::TimeOrderIterator : FieldIterator<EventTime*> {
    TimeOrderIterator(EventTime *time, BasicSeq &bseq)
        : m_time(time)
    {
        auto it = bseq.m_time_order.find(time);
        if (it != bseq.m_time_order.end()) {
            m_prev_time = &it->second;
            m_it = m_prev_time->begin();
        }
    }
    EventTime *get() const
    {
        if (!m_prev_time)
            return nullptr;
        return const_cast<EventTime*>(m_it->first);
    }
    EventTime *parent() const
    {
        return m_time;
    }
    bool is_end() const
    {
        return !m_prev_time || m_it == m_prev_time->end();
    }
    TimeOrderIterator &operator++()
    {
        ++m_it;
        return *this;
    }

private:
    EventTime *m_time;
    std::map<const EventTime*,bool> *m_prev_time = nullptr;
    std::map<const EventTime*,bool>::iterator m_it;
};

struct BasicSeq::TimeOrderVisitor: DFSVisitor<EventTime*,TimeOrderIterator,StackVector> {
    TimeOrderVisitor(BasicSeq &bseq)
        : bseq(bseq)
    {
    }
    TimeOrderIterator begin(EventTime *time)
    {
        return TimeOrderIterator(time, bseq);
    }

    BasicSeq &bseq;
};

// Used only by `Seq::add_basicseq`
BasicSeq::BasicSeq(uint32_t id)
    : m_id(id)
{
    assert(id);
}

BasicSeq::~BasicSeq()
{
}

NACS_EXPORT() Pulse *BasicSeq::add_pulse(uint32_t chn, uint32_t id,
                                         EventTime &start, Var *len, Var *val, Var *cond)
{
    return &m_channels[chn].pulses.emplace_back(id, start, len, val, cond, false);
}

NACS_EXPORT() Pulse *BasicSeq::add_measure(uint32_t chn, uint32_t id,
                                           EventTime &start, Var *val)
{
    assert(val->is_extern());
    assert((val->get_extern().second >> 32) == m_id);
    return &m_channels[chn].pulses.emplace_back(id, start, nullptr, val, nullptr, true);
}

NACS_EXPORT() Var *BasicSeq::new_measure(Env &env, uint32_t _measure_id) const
{
    uint64_t measure_id = (uint64_t(m_id) << 32) | _measure_id;
    return env.new_extern({IR::Type::Float64, measure_id});
}

NACS_EXPORT() void BasicSeq::add_branch(Var *cond, BasicSeq *target, uint32_t id)
{
    assert(cond);
    m_branches.push_back({cond->ref(), target, id});
}

NACS_EXPORT() void BasicSeq::set_default_branch(BasicSeq *target)
{
    m_default_branch = target;
}

NACS_EXPORT() bool BasicSeq::has_output(uint32_t chn) const
{
    auto it = m_channels.find(chn);
    if (it == m_channels.end())
        return false;
    for (auto &pulse: it->second.pulses) {
        if (!pulse.is_measure()) {
            return true;
        }
    }
    return false;
}

NACS_EXPORT() void BasicSeq::assign_global(uint32_t global_id, Var *val, uint32_t assignment_id)
{
    assert(val);
    auto &assign = m_assign[global_id];
    assign.val.reset(val);
    assign.id = assignment_id;
}

NACS_EXPORT() void BasicSeq::add_assume(Sign sign, Var *val, uint32_t assume_id)
{
    assert(sign == Sign::Pos || sign == Sign::NonNeg);
    m_assume.emplace(val->ref(), Assumption{sign, assume_id});
}

NACS_EXPORT() void BasicSeq::add_endtime(EventTime &t)
{
    m_endtimes.push_back(t.ref());
}

NACS_EXPORT() void BasicSeq::add_time_order(const EventTime *front, const EventTime *back,
                                            bool may_equal)
{
    auto &map = m_time_order[back];
    if (may_equal) {
        map.try_emplace(front, true);
    }
    else {
        map.insert_or_assign(front, false);
    }
}

NACS_EXPORT() EventTime &BasicSeq::track_time(const EventTime &t)
{
    return m_eventtimes.emplace_back(t);
}

NACS_EXPORT() Sign BasicSeq::compare_time(const EventTime *t1,
                                                     const EventTime *t2) const
{
    int32_t tid1 = t1->m_order_id;
    int32_t tid2 = t2->m_order_id;
    if (tid1 > tid2)
        return Sign::Unknown;
    if (tid1 == tid2)
        return Sign::NonNeg; // equal
    auto tinfo2 = m_time_sorted[tid2];
    if (tid1 <= tinfo2.ordered_pos)
        return Sign::Pos;
    auto res = Sign::Unknown;
    if (tid1 <= tinfo2.ordered_nonneg)
        res = Sign::NonNeg;
    auto it = m_time_order.find(tinfo2.time);
    if (it == m_time_order.end())
        return res;
    auto norm_t1 = m_time_sorted[tid1].time;
    auto it2 = it->second.find(norm_t1);
    if (it2 == it->second.end())
        return res;
    if (!it2->second)
        return Sign::Pos;
    return res == Sign::Pos ? res : Sign::NonNeg;
}

NACS_EXPORT() Sign BasicSeq::compare_time(const Pulse *p1, const Pulse *p2) const
{
    auto sign = compare_time(&p1->start(), &p2->start());
    if (sign == Sign::NonNeg)
        sign = p1->id() < p2->id() ? Sign::Pos : Sign::Unknown;
    return sign;
}

NACS_EXPORT() void BasicSeq::prepare(Seq &seq)
{
    preoptimize_eventtimes();
    discover_time_order();
    sort_times();
    for (auto &[chn, info]: m_channels) {
        for (auto &pulse: info.pulses) {
            if (auto norm_t = normalize(&pulse.start())) {
                pulse.m_start.reset(norm_t);
            }
        }
    }
    for (auto &time: m_endtimes) {
        if (auto norm_t = normalize(time.get())) {
            time.reset(norm_t);
        }
    }
    // Pre-compile check to make sure the sequence is computable.
    // 1. Measurements are not accessed from a different basic sequence.
    // 2. Measurements are used only by pulses after themselves.
    // Even if this sequence contains no measure, we still need to check since it might
    // reference illegal measures.
    struct MeasureUseVisitor : Env::Visitor {
        MeasureUseVisitor(const BasicSeq &seq)
            : seq(seq)
        {
        }

        void merge_depend(Var *var, std::set<Var*> &dep)
        {
            if (var->is_extern()) {
                // We do time check on each measure.
                // Since the time check is transitive,
                // we don't have to check the dependency of the measure.
                // Global variables (hi 32bit == 0) also doesn't have any dependencies
                // so we can ignore those.
                if ((var->get_extern().second >> 32) != 0)
                    dep.insert(var);
                return;
            }
            if (var->is_const())
                return;
            auto it = depends.find(var);
            if (it == depends.end())
                return;
            dep.insert(it->second.begin(), it->second.end());
        }

        bool previsit(Var *var)
        {
            if (var->is_const())
                return false;
            if (var->is_extern()) {
                auto id = uint32_t(var->get_extern().second >> 32);
                if (id != 0 && id != seq.id())
                    throw Error(Error::Type::BasicSeq, Error::BasicSeq::ExternMeasure,
                                Error::Type::Measure, var->get_extern().second,
                                Error::Type::BasicSeq, seq.id(),
                                "Measurement value can only be used "
                                "in the same basic sequence");
                return false;
            }
            assert(var->is_call());
            if (var->ref_none())
                return false;
            // We unconditionally initialize the depends entry for the variable
            // in postvisit. Since we are doing DFS on a DAG,
            // we should never visit a node twice before finishing it
            // so this should be enough to check if this node has been scanned.
            if (depends.find(var) != depends.end())
                return false;
            return true;
        }
        void postvisit(Var *var)
        {
            assert(var->is_call());
            int nargs = var->args().size();
            // If the node has been scanned, it should have been finished
            // before we start this one and rejected by previsit already.
            assert(depends.find(var) == depends.end());
            auto &dep = depends[var];
            for (int i = -1; i < nargs; i++) {
                auto ref = var->get_ref(i);
                if (!ref)
                    continue;
                merge_depend(ref, dep);
            }
        }

        void scan_deps(Var *var, std::set<Var*> &dep)
        {
            if (var->is_const())
                return;
            if (previsit(var))
                Env::scan_dfs(var, *this);
            merge_depend(var, dep);
        }
        void scan_deps(const EventTime &time, std::set<Var*> &dep)
        {
            for (auto &term: time.terms) {
                try {
                    scan_deps(term.var.get(), dep);
                }
                catch (Error &error) {
                    if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                        assert(error.type == Error::Type::BasicSeq);
                        error.type2 = Error::Type::EventTime;
                        error.id2 = term.id;
                    }
                    throw error;
                }
            }
        }
        const BasicSeq &seq;
        std::map<Var*,std::set<Var*>> depends;
    };
    MeasureUseVisitor visitor(*this);
    std::map<const Pulse*,std::set<Var*>> pulse_depends;
    std::map<Var*,const Pulse*> measures;
    for (auto &[chn, info]: m_channels) {
        auto &pulses = info.pulses;
        for (auto &pulse: pulses) {
            auto val = pulse.val();
            assert(val);
            auto &dep = pulse_depends[&pulse];
            visitor.scan_deps(pulse.start(), dep);
            if (pulse.is_measure()) {
                assert(val->is_extern());
                assert((val->get_extern().second >> 32) == m_id);
                assert(!pulse.len());
                measures[val] = &pulse;
                continue;
            }
            try {
                if (auto len = pulse.len()) {
                    visitor.scan_deps(len, dep);
                }
            }
            catch (Error &error) {
                if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                    assert(error.type == Error::Type::BasicSeq);
                    error.code = uint16_t(Error::BasicSeq::ExternMeasureLength);
                    error.type2 = Error::Type::Pulse;
                    error.id2 = pulse.id();
                }
                throw error;
            }
            try {
                if (auto cond = pulse.cond()) {
                    visitor.scan_deps(cond, dep);
                }
            }
            catch (Error &error) {
                if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                    assert(error.type == Error::Type::BasicSeq);
                    error.code = uint16_t(Error::BasicSeq::ExternMeasureCond);
                    error.type2 = Error::Type::Pulse;
                    error.id2 = pulse.id();
                }
                throw error;
            }
            try {
                visitor.scan_deps(val, dep);
            }
            catch (Error &error) {
                if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                    assert(error.type == Error::Type::BasicSeq);
                    error.type2 = Error::Type::Pulse;
                    error.id2 = pulse.id();
                }
                throw error;
            }
        }
    }
    for (auto &v: pulse_depends) {
        auto pulse = v.first;
        for (auto var: v.second) {
            assert(var->is_extern());
            assert((var->get_extern().second >> 32) == m_id);
            auto it = measures.find(var);
            assert(it != measures.end());
            // Require the measure to happen before the use
            // and also require the ID to be ordered the same way.
            // With the additional requirement on the ID ordering,
            // we can provide a consistent requirement for the user
            // when writing a subsequence.
            // This is because when a subsequence is disabled by a condition,
            // all of the pulses will happen at the same time.
            // By requiring the ID to be ordered, we also guarantee that
            // the pulses themselves are ordered so as long as the subsequence
            // is valid when enabled, it is also valid when disabled.
            // Additionally, assuming the pulse IDs are assigned
            // by their creation order in the code,
            // this requirement should be satisfied trivially in most cases.
            // This does mean that a measure
            // cannot be used to position a floating step/sequence
            // that was created before it, but I don't think that's a big problem.
            if (it->second->id() >= pulse->id() ||
                compare_time(&*it->second, pulse) != Sign::Pos) {
                if (pulse->is_measure())
                    throw Error(Error::Type::BasicSeq, Error::BasicSeq::MeasureOrder,
                                Error::Type::Measure, var->get_extern().second,
                                Error::Type::Measure, pulse->val()->get_extern().second,
                                "Use of measure must be later than the measure.");
                throw Error(Error::Type::BasicSeq, Error::BasicSeq::MeasureOrder,
                            Error::Type::Measure, var->get_extern().second,
                            Error::Type::Pulse, pulse->id(),
                            "Use of measure must be later than the measure.");
            }
        }
    }
    // Make sure assignment values are valid.
    for (auto &assign: m_assign) {
        auto var = assign.second.val.get();
        assert(var);
        try {
            if (!visitor.previsit(var))
                continue;
            Env::scan_dfs(var, visitor);
        }
        catch (Error &error) {
            if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                assert(error.type == Error::Type::BasicSeq);
                error.type2 = Error::Type::Assignment;
                error.id2 = assign.second.id;
            }
            throw error;
        }
    }
    // Make sure assumption values are valid.
    for (auto &[varref, as]: m_assume) {
        auto var = varref.get();
        assert(var);
        try {
            if (!visitor.previsit(var))
                continue;
            Env::scan_dfs(var, visitor);
        }
        catch (Error &error) {
            if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                assert(error.type == Error::Type::BasicSeq);
                error.type2 = Error::Type::EventTime;
                error.id2 = as.id;
            }
            throw error;
        }
    }
    struct CondVisitor : Env::Visitor {
        CondVisitor(Seq &seq, BasicSeq &bseq)
            : seq(seq),
              bseq(bseq)
        {
        }

        bool previsit(Var *var)
        {
            if (var->is_const())
                return false;
            if (var->is_extern()) {
                auto id = uint32_t(var->get_extern().second >> 32);
                if (id == 0) {
                    ref_types[var] = {true, false};
                }
                else {
                    ref_types[var] = {false, true};
                }
                return false;
            }
            assert(var->is_call());
            if (var->ref_none())
                return false;
            return ref_types.find(var) == ref_types.end();
        }
        void postvisit(Var *var)
        {
            assert(var->is_call());
            int nargs = var->args().size();
            assert(ref_types.find(var) == ref_types.end());
            auto &ref_type = ref_types[var];
            ref_type = {false, false};
            for (int i = -1; i < nargs; i++) {
                auto ref = var->get_ref(i);
                if (!ref)
                    continue;
                auto it = ref_types.find(ref);
                if (it == ref_types.end())
                    continue;
                ref_type.first |= it->second.first;
                ref_type.second |= it->second.second;
            }
            if (!ref_type.first || !ref_type.second)
                return;
            auto callee = var->get_callee();
            assert(var->nfreeargs() == 0);
            auto _args = var->args();
            llvm::SmallVector<Arg,4> args(_args.begin(), _args.end());
            for (int i = -1; i < nargs; i++) {
                auto ref = var->get_ref(i);
                if (!ref)
                    continue;
                auto it = ref_types.find(ref);
                if (it == ref_types.end() || !it->second.second)
                    continue;
                ref = escape_var(ref);
                if (i < 0) {
                    callee.var = ref;
                }
                else {
                    args[i] = Arg::create_var(ref);
                }
            }
            auto newvar = callee.is_llvm ? seq.env().new_call(callee.llvm, args, 0) :
                seq.env().new_call(callee.var, args, 0);
            global_escapes.emplace(var, newvar);
        }
        Var *escape_var(Var *var)
        {
            auto &res = global_escapes[var];
            if (!res) {
                uint32_t offset = bseq.m_global_escapes.size();
                assert(offset + 1 <= global_escapes.size());
                uint32_t global_id = seq.npublic_globals() + offset;
                bseq.assign_global(global_id, var, 0);
                // Create global `Var` with a sequence specific pointer identity
                // so that we can replace it with a normal variable later.
                // We will allocate the actual global variables
                // after all the optimizations are done in `optimize_final`.
                res = seq.env().new_extern({var->type(), global_id});
                // We need to keep a reference to it
                // so that `Env` will honer the pointer identity
                bseq.m_global_escapes.push_back(res->ref());
                assert(bseq.m_global_escapes.size() <= global_escapes.size());
            }
            return res;
        }
        Var *check_var(Var *var)
        {
            if (var->is_const())
                return var;
            if (previsit(var))
                Env::scan_dfs(var, *this);
            auto it = global_escapes.find(var);
            if (it != global_escapes.end())
                return it->second;
            auto it2 = ref_types.find(var);
            if (it2 == ref_types.end() || !it2->second.second)
                return var;
            return escape_var(var);
        }

        Seq &seq;
        BasicSeq &bseq;
        std::map<Var*,std::pair<bool,bool>> ref_types; // ref_global, ref_measure
        std::map<Var*,Var*> global_escapes;
    };
    CondVisitor cond_visitor(seq, *this);

    // Make sure branch conditions are valid.
    for (auto &br: m_branches) {
        auto var = br.cond.get();
        assert(var);
        try {
            if (!visitor.previsit(var))
                continue;
            Env::scan_dfs(var, visitor);
        }
        catch (Error &error) {
            if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                assert(error.type == Error::Type::BasicSeq);
                error.type2 = Error::Type::Branch;
                error.id2 = br.id;
            }
            throw error;
        }
        // Unlike everything else, the use of the branch condition
        // happens **AFTER** user callback, which may change global variables
        // and we would like to reflect such change in the condition.
        // However, if the condition depends on a measure in the sequence
        // that ends up being a function of a global variable, we do not want to update those.
        // Therefore, create an assignment of global variable from the measure results.
        // and make the branch condition use those instead.
        // This way, the result of the measure would be frozen after the sequence.
        // Due to the explicit escape through global variable assignment
        // we won't automatically optimize measured results in the conditions.
        // However, we explicitly check if the assignment is constant
        // and replace the use of the global variable if that is the case.
        br.cond.reset(cond_visitor.check_var(var));
    }
    // Make sure endtimes are valid.
    for (auto &time: m_endtimes) {
        for (auto &term: time->terms) {
            try {
                auto var = term.var.get();
                if (var->is_const() || !visitor.previsit(var))
                    continue;
                Env::scan_dfs(var, visitor);
            }
            catch (Error &error) {
                if (error.code == uint16_t(Error::BasicSeq::ExternMeasure)) {
                    assert(error.type == Error::Type::BasicSeq);
                    error.type2 = Error::Type::EventTime;
                    error.id2 = term.id;
                }
                throw error;
            }
        }
    }
}

void BasicSeq::discover_time_order()
{
    m_eventtimes.sort([] (auto &t1, auto &t2) {
        if (t1.tconst < t2.tconst)
            return true;
        if (t1.tconst > t2.tconst)
            return false;
        auto nt1 = t1.terms.size();
        auto nt2 = t2.terms.size();
        if (nt1 < nt2)
            return true;
        return false;
    });
    // These variables stores the current status of the ordering among the ones we've checked.
    // This is a vector of times that forms a total order.
    llvm::SmallVector<EventTime*,32> ordered;
    // The latest element on each compariable branch.
    // None of the elements in `ordered` are in `latest`.
    // All elements in `latest` are known to be after all elements in `ordered`.
    std::set<EventTime*> latest;
    // For each new time, if it is known to be after everything we've seen before,
    // then we can simply add it to the list and assume it has a unique next time.
    // Otherwise, we need to find the last time it is definitely after
    // and remove all times after that time (without adding this one to the list).
    auto known_after = [&] (auto *t1, auto *t2) {
        // assume `t1` was sorted before `t2`;
        assert(t1->tconst <= t2->tconst);
        auto isless = t1->isless_terms(*t2);
        if (isless == Sign::Unknown)
            return false;
        if (isless == Sign::NonNeg && t1->tconst < t2->tconst)
            isless = Sign::Pos;
        add_time_order(t1, t2, isless != Sign::Pos);
        if (isless != Sign::Pos) {
            auto isless2 = t2->isless(*t1);
            if (isless2 != Sign::Unknown) {
                assert(isless2 == Sign::NonNeg);
                add_time_order(t2, t1, true);
            }
        }
        return isless == Sign::Pos;
    };
    auto num_ordered = [&] (EventTime *time) {
        assert(!ordered.empty());
        assert(!known_after(ordered.back(), time));
        if (ordered.size() == 1 || !known_after(ordered.front(), time))
            return 0u;
        unsigned lo = 0;
        unsigned hi = ordered.size() - 1;
        // The index of the last element that is before `time` is `[lo, hi)`.
        // `lo` is known before `time` whereas `hi` isn't.
        while (hi > lo + 1) {
            unsigned mid = (lo + hi) / 2;
            if (known_after(ordered[mid], time)) {
                lo = mid;
            }
            else {
                hi = mid;
            }
        }
        // `lo` is the index of the last element that is known before `time`.
        return lo + 1; // Total number of times known before `time`.
    };
    auto add_time = [&] (EventTime *time) {
        auto had_latest = !latest.empty();
        for (auto it = latest.begin(); it != latest.end();) {
            if (known_after(*it, time)) {
                it = latest.erase(it);
            }
            else {
                ++it;
            }
        }
        if (ordered.empty() || known_after(ordered.back(), time)) {
            if (latest.empty()) {
                // Time is after everything we've seen
                ordered.push_back(time);
            }
            else {
                latest.insert(time);
            }
            return;
        }
        // If `latest` wasn't empty but is empty now,
        // then time is after everything that was in `latest`
        // which is after the last element in `ordered` so it cannot happen in this branch.
        if (had_latest)
            assert(!latest.empty());
        latest.insert(time);
        // If `latest` wasn't empty we know that either `time` or another one in `latest`
        // is already after everything in `ordered` so we don't need to add anything
        // from `ordered` to `latest`.
        if (!had_latest)
            latest.insert(ordered.back());
        ordered.resize(num_ordered(time));
    };
    for (auto &time: m_eventtimes) {
        add_time(&time);
    }
}

void BasicSeq::sort_times()
{
    // Do a topological sort of the event times
    // and use Tarjan's algorithm to detect times that are equivalent.
    struct Visitor : TimeOrderVisitor {
        using TimeOrderVisitor::TimeOrderVisitor;

        bool previsit(EventTime *time)
        {
            if (time->m_order_id >= 0 || time->m_on_stack)
                return false;
            time->m_on_stack = true;
            time->m_scc_index = ++scc_count;
            time->m_scc_index_lowest = time->m_scc_index;
            scc_stack.push_back(time);
            return true;
        }

        void set_scc_lowest(EventTime *time)
        {
            auto it = bseq.m_time_order.find(time);
            if (it == bseq.m_time_order.end())
                return;
            auto res = time->m_scc_index_lowest;
            for (auto [earlier, may_equal]: it->second) {
                if (!earlier->m_on_stack)
                    continue;
                if (!may_equal)
                    throw std::runtime_error("Time loop detected.");
                res = min(res, earlier->m_scc_index_lowest);
            }
            time->m_scc_index_lowest = res;
        }

        void postvisit(EventTime *time)
        {
            set_scc_lowest(time);
            // Not a root of the SCC
            if (time->m_scc_index_lowest != time->m_scc_index)
                return;
            time->m_order_id = time_order_id++;
            if (scc_stack.back() == time) {
                scc_stack.pop_back();
                bseq.m_time_sorted.push_back(TimeInfo{time});
                time->m_on_stack = false;
            }
            else {
                auto best_time = time;
                uint32_t stack_sz = scc_stack.size();
                assert(stack_sz > 1);
                uint32_t diff;
                for (diff = 1; ; diff++) {
                    assert(diff <= stack_sz);
                    auto t = scc_stack[stack_sz - diff];
                    assert(t->m_on_stack);
                    t->m_on_stack = false;
                    if (t == time)
                        break;
                    assert(t->m_scc_index_lowest == time->m_scc_index_lowest);
                    t->m_order_id = time->m_order_id;
                    if (t->tconst < best_time->tconst ||
                        (t->tconst == best_time->tconst &&
                         t->terms.size() == best_time->terms.size()))
                        continue;
                    best_time = t;
                }
                // Merge all ordering to the normalized time
                for (uint32_t i = stack_sz - diff; i < stack_sz; i++) {
                    auto t = scc_stack[stack_sz - diff];
                    if (t == best_time)
                        continue;
                    auto it = bseq.m_time_order.find(t);
                    assert(it != bseq.m_time_order.end());
                    for (auto [earlier, may_equal]: it->second) {
                        if (earlier == best_time)
                            continue;
                        bseq.add_time_order(best_time, earlier, may_equal);
                    }
                    bseq.m_time_order.erase(it);
                }
                scc_stack.resize(stack_sz - diff);
                bseq.m_time_sorted.push_back(TimeInfo{best_time});
            }
            assert(bseq.m_time_sorted.size() == uint32_t(time->m_order_id + 1));
        }

        uint32_t time_order_id = 0;
        uint32_t scc_count = 0;
        llvm::SmallVector<EventTime*,8> scc_stack;
    };
    Visitor visitor(*this);
    for (auto &time: m_eventtimes) {
        assert(visitor.scc_stack.empty());
        if (!visitor.previsit(&time))
            continue;
        visit_dfs(&time, visitor);
        assert(visitor.scc_stack.empty());
    }

    assert(visitor.time_order_id == m_time_sorted.size());
    assert(visitor.scc_count == m_eventtimes.size());
    assert(visitor.time_order_id <= visitor.scc_count);

    // keys in m_time_order should all be normalized now.
    // However, the values aren't necessarily so.
    // Normalize all the times in the values now.
    for (auto &[time, order_map]: m_time_order) {
        for (auto it = order_map.begin(), end = order_map.end(); it != end;) {
            auto id = it->first->m_order_id;
            auto norm_t = m_time_sorted[id].time;
            if (it->first == norm_t) {
                ++it;
                continue;
            }
            auto may_equal = it->second;
            it = order_map.erase(it);
            if (may_equal) {
                order_map.try_emplace(norm_t, true);
            }
            else {
                order_map.insert_or_assign(norm_t, false);
            }
        }
    }

    // Now `m_time_sorted` contains a topological sort of the times
    // We can scan the tree to optimize the datastructure for future lookup.
    uint32_t ntimes = m_time_sorted.size();
    for (uint32_t scan_id = 0; scan_id < ntimes; scan_id++) {
        auto &info = m_time_sorted[scan_id];
        auto order_it = m_time_order.find(info.time);
        if (order_it == m_time_order.end())
            continue;
        auto &order_map = order_it->second;
        auto update_ordered_idx = [&] {
            for (uint32_t prev_idx = info.ordered_pos + 1; prev_idx < scan_id; prev_idx++) {
                auto pt = m_time_sorted[prev_idx].time;
                auto it = order_map.find(pt);
                if (it == order_map.end())
                    return;
                if (!it->second)
                    info.ordered_pos = prev_idx;
                info.ordered_nonneg = max(int32_t(prev_idx), info.ordered_nonneg);
                order_map.erase(it);
            }
        };
        update_ordered_idx();
        for (uint32_t diff = 1; diff <= scan_id; diff++) {
            int32_t idx = scan_id - diff;
            if (idx <= info.ordered_pos)
                break;
            auto &prev_info = m_time_sorted[idx];
            auto order = Sign::Unknown;
            auto it = order_map.find(prev_info.time);
            if (it != order_map.end()) {
                order = it->second ? Sign::NonNeg : Sign::Pos;
            }
            else if (idx <= info.ordered_nonneg) {
                order = Sign::NonNeg;
            }
            else {
                continue;
            }
            auto may_equal = order == Sign::NonNeg;
            info.ordered_nonneg = max(info.ordered_nonneg,
                                      may_equal ? prev_info.ordered_nonneg :
                                      prev_info.ordered_pos);
            info.ordered_pos = max(info.ordered_pos, prev_info.ordered_pos);
            assert(info.ordered_pos <= info.ordered_nonneg);
            auto it2 = m_time_order.find(prev_info.time);
            if (it2 != m_time_order.end()) {
                for (auto [earlier, _may_equal]: it2->second) {
                    _may_equal = _may_equal && may_equal;
                    if (_may_equal) {
                        if (earlier->m_order_id <= info.ordered_nonneg)
                            continue;
                        order_map.emplace(earlier, true);
                    }
                    else {
                        if (earlier->m_order_id <= info.ordered_pos)
                            continue;
                        order_map.emplace(earlier, false);
                    }
                }
            }
            update_ordered_idx();
        }
        // Remove unneeded entries from the map.
        for (auto it = order_map.begin(), end = order_map.end(); it != end;) {
            int32_t id = it->first->m_order_id;
            if (id > info.ordered_nonneg || (id > info.ordered_pos && it->second)) {
                ++it;
                continue;
            }
            it = order_map.erase(it);
        }
    }
}

void BasicSeq::mark_recursive()
{
    if (m_used)
        return;
    m_used = true;
    // Simply use the native stack since I don't expect there to be too many basic sequence
    // for stack overflow to be a problem...
    if (m_default_branch) {
        m_default_branch->m_has_branchin = true;
        m_default_branch->mark_recursive();
    }
    for (auto &br: m_branches) {
        if (br.target) {
            br.target->m_has_branchin = true;
            br.target->mark_recursive();
        }
    }
}

bool BasicSeq::preoptimize_pulse(uint32_t chn, Env &env)
{
    // This makes sure a channel info is allocated for all channels.
    bool changed = false;
    auto &info = m_channels[chn];
    auto &pulses = info.pulses;
    auto &pulse_start_vars = info.pulse_start_vars;
    info.startval = alloc_startval(env)->ref();
    for (auto it = pulses.begin(), end = pulses.end(); it != end;) {
        auto &pulse = *it;
        pulse.check_start();
        if (!pulse.is_measure()) {
            changed |= pulse.clear_unused_args();
            changed |= pulse.optimize();
            if (auto cond = pulse.cond()) {
                if (cond->is_const() && !cond->get_const().get<bool>()) {
                    it = pulses.erase(it);
                    changed = true;
                    continue;
                }
            }
            if (pulse.needs_oldval()) {
                auto startval = alloc_startval(env);
                pulse_start_vars.try_emplace(&pulse, startval->ref());
                pulse.set_oldval(startval);
            }
            assert(!pulse.needs_oldval());
        }
        else {
            auto m = pulse.val();
            assert(m->extern_used() >= 1);
            if (m->extern_used() == 1 && !m->used(false)) {
                it = pulses.erase(it);
                changed = true;
                continue;
            }
        }
        ++it;
    }
    return changed;
}

bool BasicSeq::optimize_order(uint32_t chn)
{
    bool changed = false;
    auto &info = m_channels[chn];
    auto &pulses = info.pulses;
    if (pulses.empty()) {
        if (info.endval)
            return false;
        info.endval.reset(info.startval.get());
        return true;
    }
    auto &pulse_start_vars = info.pulse_start_vars;
    bool should_continue = false;
    for (auto it = pulses.begin(), end = pulses.end(); it != end;) {
        auto &pulse = *it;
        pulse.check_start();
        if (pulse.is_measure()) {
            auto m = pulse.val();
            assert(m->extern_used() >= 1);
            if (m->extern_used() == 1 && !m->used(false)) {
                it = pulses.erase(it);
                changed = true;
                continue;
            }
            should_continue = true;
        }
        else {
            // Requires preoptimize_pulse
            assert(!pulse.needs_oldval());
            changed |= pulse.clear_unused_args();
            changed |= pulse.optimize();
            if (auto cond = pulse.cond()) {
                // We do not forward the value of a disabled pulse
                // to anything else so it should be safe to delete the pulse
                if (cond->is_const() && !cond->get_const().get<bool>()) {
                    auto start_it = pulse_start_vars.find(&pulse);
                    it = pulses.erase(it);
                    if (start_it != pulse_start_vars.end()) {
                        assert(start_it->second->extern_used() == 1 &&
                               !start_it->second->used(false));
                        pulse_start_vars.erase(start_it);
                    }
                    changed = true;
                    continue;
                }
            }
        }
        ++it;
    }
    for (auto it = pulse_start_vars.begin(), end = pulse_start_vars.end(); it != end;) {
        auto startval = it->second.get();
        assert(startval->extern_used() >= 1);
        if (startval->extern_used() == 1 && !startval->used(false)) {
            it = pulse_start_vars.erase(it);
            changed = true;
            continue;
        }
        should_continue = true;
        ++it;
    }
    // Nothing more to optimize here...
    if (!should_continue && endval(chn))
        return changed;
    pulses.sort([] (const Pulse &p1, const Pulse &p2) {
        auto &t1 = p1.start();
        auto &t2 = p2.start();
        if (t1.m_order_id != t2.m_order_id)
            return t1.m_order_id < t2.m_order_id;
        return p1.id() < p2.id();
    });
    // These variables stores the current status of the ordering among the ones we've checked.
    // This is a vector of pulses that forms a total order.
    llvm::SmallVector<Pulse*,32> ordered;
    // Whether there's a unique first pulse
    // Whether each pulse have a unique next pulse,
    // there is one more element than `orderred`.
    // The first element record if the first pulse in the sequence is unique.
    llvm::SmallVector<bool,32> next_unique(1, true);
    // Measures that are interested in the particular pulse
    // Also has one more element than `ordered`.
    llvm::SmallVector<llvm::SmallVector<Pulse*,1>,32> measures(1);
    // The latest element on each compariable branch.
    // None of the elements in `ordered` are in `latest`.
    // All elements in `latest` are known to be after all elements in `ordered`.
    std::set<Pulse*> latest;
    // For each new pulse, if it is known to be after everything we've seen before,
    // then we can simply add it to the list and assume it has a unique next pulse.
    // Otherwise, we need to find the last pulse it is definitely after
    // and remove all pulses after that pulse (without adding this one to the list).
    auto known_after_real = [&] (Pulse *p1, Pulse *p2) {
        // assume `p1` was sorted before `p2`
        auto &t1 = p1->start();
        auto &t2 = p2->start();
        auto isless = compare_time(&t1, &t2);
        if (isless == Sign::NonNeg)
            return p1->id() < p2->id();
        return isless == Sign::Pos;
    };
    auto known_after = [&] (Pulse *p1, Pulse *p2) {
        // assume `p1` was sorted before `p2`;
        if ((p1->is_ordered && !p2->is_measure()) || (p2->is_ordered && !p1->is_measure())) {
            // For ordered pulse everything has a defined known before/after relation.
            // Since `p1` was sorted before `p2` we must have `p1` is known before `p2`
            assert(known_after_real(p1, p2));
            return true;
        }
        return known_after_real(p1, p2);
    };
    auto num_ordered = [&] (Pulse *pulse) {
        assert(!ordered.empty());
        assert(!known_after(ordered.back(), pulse));
        if (ordered.size() == 1 || !known_after(ordered.front(), pulse))
            return 0u;
        unsigned lo = 0;
        unsigned hi = ordered.size() - 1;
        // The index of the last element that is before `pulse` is `[lo, hi)`.
        // `lo` is known before `pulse` whereas `hi` isn't.
        while (hi > lo + 1) {
            unsigned mid = (lo + hi) / 2;
            if (known_after(ordered[mid], pulse)) {
                lo = mid;
            }
            else {
                hi = mid;
            }
        }
        // `lo` is the index of the last element that is known before `pulse`.
        return lo + 1; // Total number of pulses known before `pulse`.
    };
    auto add_measure = [&] (Pulse *pulse) {
        // If a measure pulse isn't after the last orderred one,
        // or the last one isn't ordered, simply ignore it.
        // We already know that we can't decide it's value
        // and it won't affect any time ordering of the normal pulses.
        if (!latest.empty() || (!ordered.empty() && !known_after(ordered.back(), pulse)))
            return;
        if (next_unique.back()) {
            measures.back().push_back(pulse);
        }
    };
    auto check_measures = [&] (Pulse *pulse) {
        assert(ordered.empty() || known_after(ordered.back(), pulse));
        assert(!measures.empty());
        assert(measures.size() == ordered.size() + 1);
        auto &measure = measures.back();
        if (measure.empty())
            return;
        measure.erase(std::remove_if(measure.begin(), measure.end(), [&] (const auto &m) {
            return !known_after(m, pulse);
        }), measure.end());
    };
    auto add_pulse = [&] (Pulse *pulse) {
        if (pulse->is_measure()) {
            add_measure(pulse);
            return;
        }
        if (pulse->is_ordered) {
#ifndef NDEBUG
            for (auto p: latest)
                assert(known_after(p, pulse));
#endif
            latest.clear();
            assert(ordered.empty() || known_after(ordered.back(), pulse));
            check_measures(pulse);
            // Pulse is after everything we've seen
            ordered.push_back(pulse);
            next_unique.push_back(true);
            measures.emplace_back();
            return;
        }
        auto had_latest = !latest.empty();
        for (auto it = latest.begin(); it != latest.end();) {
            if (known_after(*it, pulse)) {
                it = latest.erase(it);
            }
            else {
                ++it;
            }
        }
        if (ordered.empty() || known_after(ordered.back(), pulse)) {
            check_measures(pulse);
            if (latest.empty()) {
                // Pulse is after everything we've seen
                ordered.push_back(pulse);
                next_unique.push_back(true);
                measures.emplace_back();
            }
            else {
                latest.insert(pulse);
            }
            return;
        }
        // If `latest` wasn't empty but is empty now,
        // then pulse is after everything that was in `latest`
        // which is after the last element in `ordered` so it cannot happen in this branch.
        if (had_latest)
            assert(!latest.empty());
        latest.insert(pulse);
        // If `latest` wasn't empty we know that either `pulse` or another one in `latest`
        // is already after everything in `ordered` so we don't need to add anything
        // from `ordered` to `latest`.
        if (!had_latest)
            latest.insert(ordered.back());
        auto no = num_ordered(pulse);
        ordered.resize(no);
        next_unique.resize(no + 1);
        next_unique[no] = false;
        measures.resize(no + 1);
        check_measures(pulse);
    };
    for (auto &pulse: pulses)
        add_pulse(&pulse);
    latest.clear(); // Don't need anymore

    Var *curval = info.startval.get();
    assert(curval);
    bool measure_replaced = false;
    for (auto m: measures.front()) {
        measure_replaced = true;
        m->val()->assign_var(curval);
        assert(!m->val()->is_extern());
    }
    if (!next_unique.front())
        curval = nullptr;
    auto use_val = [&] (Var *val, Pulse *measure) {
        measure_replaced = true;
        measure->val()->assign_var(val);
        assert(!measure->val()->is_extern());
    };
    auto use_pulse_end = [&] (Pulse *pulse, Pulse *measure) {
        assert(!pulse->cond());
        auto val = pulse->endval();
        assert(val); // since we know `!pulse->needs_oldval()`
        use_val(val, measure);
    };
    auto replace_measure = [&] (Pulse *pulse, Pulse *measure) {
        assert(!pulse->cond());
        assert(!pulse->needs_oldval());
        auto pulse_len = pulse->len();
        if (!pulse_len) {
            use_pulse_end(pulse, measure);
            return;
        }
        // Caller should have optimized each pulse already (clear_unused_args)
        // so a pulse that doesn't need time should already have a NULL length.
        assert(pulse->val()->is_call());
        assert(pulse->val()->nfreeargs() == 1);
        assert(!pulse->val()->argument_unused(0));
        auto &tp = pulse->start();
        auto &mp = measure->start();
        assert(tp.tconst <= mp.tconst);
        assert(known_after(pulse, measure));
        auto tdiff = mp - tp;
        bool needs_branch = true;
        if (pulse_len->is_const()) {
            int64_t t = round<int64_t>(pulse_len->get_const().get<double>());
            if (t <= tdiff.tconst) {
                use_pulse_end(pulse, measure);
                return;
            }
            needs_branch = !tdiff.terms.empty();
        }
        else {
            for (auto &term: tdiff.terms) {
                if (term.var.get() == pulse_len) {
                    use_pulse_end(pulse, measure);
                    return;
                }
            }
        }
        auto &env = measure->val()->env();
        auto toffset = tdiff.to_var(env);
        if (needs_branch) {
            auto mod = env.llvm_module();
            auto &ctx = mod->getContext();
            auto cgctx = env.cg_context();
            llvm::Function *f;
            {
                llvm::FunctionType *ftype;
                llvm::Type *fsig[] = { cgctx->T_f64, cgctx->T_f64 };
                ftype = llvm::FunctionType::get(cgctx->T_f64, fsig, false);
                f = llvm::Function::Create(ftype, llvm::GlobalValue::ExternalLinkage,
                                           "m", mod);
            }
            f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
            f->addFnAttr(llvm::Attribute::AlwaysInline);
            f->addFnAttr(llvm::Attribute::Speculatable);
            f->addFnAttr(llvm::Attribute::NoRecurse);
            f->addFnAttr(llvm::Attribute::NoUnwind);
            f->addFnAttr(llvm::Attribute::ReadOnly);

            {
                auto b0 = llvm::BasicBlock::Create(ctx, "top", f);
                llvm::IRBuilder<> builder(b0);
                llvm::FastMathFlags fmf;
                // Do not set noinf since we may use divide by 0
                fmf.setNoNaNs();
                fmf.setNoSignedZeros();
                fmf.setAllowReciprocal();
#if LLVM_VERSION_MAJOR >= 7
                fmf.setAllowContract();
#endif
                builder.setFastMathFlags(fmf);

                auto intrin = llvm::Intrinsic::getDeclaration(mod, llvm::Intrinsic::minnum,
                                                              { cgctx->T_f64 });
                builder.CreateRet(builder.CreateCall(intrin, { f->getArg(0), f->getArg(1) }));
            }
            toffset = env.new_call(f, { Arg::create_var(toffset),
                    Arg::create_var(pulse_len) });
        }
        use_val(env.new_call(pulse->val(), { Arg::create_var(toffset) }), measure);
    };
    auto nordered = ordered.size();
    for (unsigned i = 0; i < nordered; i++) {
        auto pulse = ordered[i];
        pulse->is_ordered = true;
        assert(!pulse->is_measure());
        if (curval) {
            auto it = pulse_start_vars.find(pulse);
            if (it != pulse_start_vars.end()) {
                it->second->assign_var(curval);
                pulse_start_vars.erase(it);
                changed = true;
            }
        }
        // If a pulse may be disabled, we cannot use it's value for measures
        // since they may not be valid
        // (their values will be the end value of the previous pulse
        //  which is wrong if the pulse is disabled
        //  but starts in the middle of the previous ramp pulse)
        // In principle, using it as the start value of the next pulse is fine
        // since we want the end value anyway.
        // However, I'd like to be more conservative here
        // and this also allow us to make sure the pulse start value of this pulse
        // isn't used anywhere else so we can safely delete the pulse
        // when we know it's disabled.
        if (pulse->cond()) {
            curval = nullptr;
            continue;
        }
        for (auto m: measures[i + 1])
            replace_measure(pulse, m);
        curval = next_unique[i + 1] ? pulse->endval() : nullptr;
    }
    if (curval) {
        auto &old_endval = info.endval;
        if (!old_endval) {
            old_endval.reset(curval);
            changed = true;
        }
    }

    // remove measurements that are filled in
    if (measure_replaced) {
        changed = true;
        pulses.remove_if([&] (auto &pulse) {
            return pulse.is_measure() && !pulse.val()->is_extern();
        });
    }
    return changed;
}

bool BasicSeq::optimize_endtimes()
{
    bool changed = false;
    // Sort in reverse order so that we'll get the best bound quickest
    m_endtimes.sort([] (const auto &t1, const auto &t2) {
        return t1->m_order_id > t2->m_order_id;
    });
    int32_t ordered_nonneg = -1;
    std::set<const EventTime*> known_earlier;
    for (auto it = m_endtimes.begin(), end = m_endtimes.end(); it != end;) {
        auto time = &**it;
        if (time->m_order_id <= ordered_nonneg ||
            (time->terms.size() == 0 && time->tconst <= 0) ||
            known_earlier.find(time) != known_earlier.end()) {
            changed = true;
            it = m_endtimes.erase(it);
            continue;
        }
        ++it;
        auto info = m_time_sorted[time->m_order_id];
        ordered_nonneg = max(ordered_nonneg, info.ordered_nonneg);
        auto order_it = m_time_order.find(info.time);
        if (order_it == m_time_order.end())
            continue;
        for (auto [earlier, may_equal]: order_it->second) {
            if (earlier->m_order_id > ordered_nonneg) {
                known_earlier.insert(earlier);
            }
        }
    }
    return changed;
}

bool BasicSeq::preoptimize_eventtimes()
{
    bool changed = false;
    for (auto &time: m_eventtimes)
        changed |= time.normalize();
    return changed;
}

bool BasicSeq::postoptimize_eventtimes()
{
    bool changed = false;
    for (auto it = m_eventtimes.begin(), end = m_eventtimes.end(); it != end;) {
#ifndef NDEBUG
        if (m_time_sorted[it->m_order_id].time != &*it)
            assert(it->m_ref_count == 0);
#endif
        if (it->m_ref_count == 0) {
            changed = true;
            // We never need to traverse more than 1 level in `m_time_order`
            // so there's no need to forward the value in `m_time_order`
            // for the time to be deleted.
            // However, we do need to delete it from all the values...
            m_time_order.erase(&*it);
            for (auto &[time, order_map]: m_time_order)
                order_map.erase(&*it);
            if (m_time_sorted[it->m_order_id].time == &*it)
                m_time_sorted[it->m_order_id] = TimeInfo{nullptr};
            it = m_eventtimes.erase(it);
        }
        else {
            ++it;
        }
    }
    return changed;
}

EventTime *BasicSeq::normalize(const EventTime *time) const
{
    auto norm_t = m_time_sorted[time->m_order_id].time;
    if (norm_t == time)
        return nullptr;
    return const_cast<EventTime*>(norm_t);
}

bool BasicSeq::optimize_vars(Seq &seq)
{
    bool changed = false;
    for (auto &[chn, info]: m_channels) {
        if (auto v = info.startval->get_assigned_var())
            info.startval.reset(v);
        if (info.endval) {
            if (auto v = info.endval->get_assigned_var()) {
                info.endval.reset(v);
            }
        }
        for (auto &pulse: info.pulses) {
            if (auto norm_t = normalize(&pulse.start())) {
                pulse.m_start.reset(norm_t);
                changed = true;
            }
        }
    }
    for (auto &time: m_endtimes) {
        if (auto norm_t = normalize(time.get())) {
            time.reset(norm_t);
            changed = true;
        }
    }
    auto npublic_globals = seq.npublic_globals();
    for (auto it = m_assign.begin(), end = m_assign.end(); it != end;) {
        auto &[global_id, assign] = *it;
        assert(assign.val);
        if (auto v = assign.val->get_assigned_var())
            assign.val.reset(v);
        if (global_id < npublic_globals) {
            ++it;
            continue;
        }
        auto local_id = global_id - npublic_globals;
        auto &escape = m_global_escapes[local_id];
        assert(escape);
        assert(escape->extern_used() >= 1);
        if (escape->extern_used() == 1 && !escape->used(false)) {
            // No other of this value anymore
            changed = true;
            escape.reset(nullptr);
            it = m_assign.erase(it);
            continue;
        }
        else if (!assign.val->is_const()) {
            ++it;
            continue;
        }
        changed = true;
        // If the assignment to a global escape value is constant,
        // we can replace it with the constant.
        escape->assign_var(assign.val.get());
        escape.reset(nullptr);
        it = m_assign.erase(it);
    }
    auto check_assumption = [&] (Var *var, auto &as) {
        assert(var->is_const());
        auto v = var->get_const().get<double>();
        if (as.sign == Sign::Pos && !(round<int64_t>(v) > 0)) {
            throw Error(Error::Type::EventTime, Error::EventTime::NonPosTime,
                        Error::Type::EventTime, as.id, "Positive time expected.");
        }
        else if (as.sign == Sign::NonNeg && !(v >= 0)) {
            throw Error(Error::Type::EventTime, Error::EventTime::NegTime,
                        Error::Type::EventTime, as.id, "Non-negative time expected.");
        }
    };
    for (auto it = m_assume.begin(), end = m_assume.end(); it != end;) {
        auto &[var, as] = *it;
        if (auto v = var->get_assigned_var()) {
            if (v->is_const()) {
                check_assumption(v, as);
            }
            else {
                auto [it2, inserted] = m_assume.emplace(v->ref(), as);
                if (!inserted && as.sign == Sign::Pos &&
                    it2->second.sign != Sign::Pos) {
                    it2->second = as;
                }
            }
            it = m_assume.erase(it);
            continue;
        }
        if (!var->is_const()) {
            ++it;
            continue;
        }
        check_assumption(var.get(), as);
        it = m_assume.erase(it);
    }

    for (auto &br: m_branches) {
        if (auto v = br.cond->get_assigned_var()) {
            br.cond.reset(v);
        }
    }
    changed |= preoptimize_eventtimes();
    return changed;
}

bool BasicSeq::optimize_branch()
{
    bool changed = false;
    bool dead = false;
    for (auto &br: m_branches) {
        auto &cond = br.cond;
        if (dead) {
            cond = nullptr;
            continue;
        }
        if (auto v = cond->get_assigned_var())
            cond = v->ref();
        if (cond->is_const()) {
            auto c = cond->get_const();
            changed = true;
            if (c.get<bool>()) {
                dead = true;
                m_default_branch = br.target;
            }
            cond = nullptr;
            continue;
        }
    }
    size_t nbranches = m_branches.size();
    size_t num_non_default = 0;
    for (size_t i = 0; i < nbranches; i++) {
        auto &br = m_branches[i];
        if (!br.cond)
            continue;
        if (br.target != m_default_branch) {
            num_non_default = i + 1;
        }
    }
    if (num_non_default < nbranches) {
        m_branches.resize(num_non_default);
        changed = true;
    }
    m_branches.erase(std::remove_if(m_branches.begin(), m_branches.end(), [&] (const auto &br) {
        return !br.cond;
    }), m_branches.end());
    return changed;
}

bool BasicSeq::optimize_final(Seq &seq)
{
    bool changed = false;
    uint32_t npublic_globals = seq.npublic_globals();
    uint32_t nescape = m_global_escapes.size();
    int32_t cond_use = 0;
    for (uint32_t old_id = 0, new_id = 0; old_id < nescape; old_id++) {
        auto &escape = m_global_escapes[old_id];
        if (!escape) {
            assert(m_assign.find(npublic_globals + old_id) == m_assign.end());
            continue;
        }
        changed = true;
        if (old_id != new_id) {
            assert(new_id < old_id);
            auto node = m_assign.extract(npublic_globals + old_id);
            assert(node);
            node.key() = npublic_globals + new_id;
            auto res = m_assign.insert(std::move(node));
            assert(res.inserted);
            (void)res;
        }
        else {
            assert(m_assign.find(npublic_globals + new_id) != m_assign.end());
        }
        assert(escape->is_extern());
        escape->assign_var(seq.get_slot(escape->get_extern().first, npublic_globals + new_id));
        cond_use += escape->extern_used() - 1;
        escape.reset(nullptr);
        new_id++;
    }
    if (cond_use > 0) {
        assert(changed);
        // If we've created variables assignments, we need to optimize those out.
        // They can either be used by other values (which are then used in conditions)
        // or in a condition directly.
        // The former will be handled by a variable optimization done by the caller,
        // but we need to handle the latter here.
        for (auto &br: m_branches) {
            auto &cond = br.cond;
            if (auto v = cond->get_assigned_var()) {
                cond = v->ref();
                assert(cond->is_extern());
                cond_use--;
            }
            assert(!cond->is_const());
        }
        assert(cond_use == 0);
    }
    return changed;
}

NACS_EXPORT() void BasicSeq::print(std::ostream &stm) const
{
    stm << "BS(" << id() << "):" << std::endl;
    for (auto &[chn, info]: m_channels) {
        auto &pulses = info.pulses;
        stm << "  CH(" << chn << "):" << std::endl;
        if (info.startval) {
            stm << "   -init=";
            info.startval->print(stm, false, true);
        }
        stm << std::endl;
        for (auto &pulse: pulses) {
            auto it = info.pulse_start_vars.find(&pulse);
            auto startval = it != info.pulse_start_vars.end() ? it->second.get() : nullptr;
            stm << "    ";
            pulse.print(stm, startval, true);
        }
        if (info.endval) {
            stm << "   -final=";
            info.endval->print(stm, false, true);
            stm << std::endl;
        }
    }
    if (!m_assign.empty()) {
        stm << "  Assignment:" << std::endl;
        for (auto &assign: m_assign) {
            stm << "    " << assign.first << " <- ";
            assign.second.val->print(stm, true, true);
        }
    }
    if (!m_assume.empty()) {
        stm << "  Assumptions:" << std::endl;
        for (auto &[var, as]: m_assume) {
            stm << "    ";
            var->print(stm, false, true);
            stm << "/";
            if (as.sign == Sign::Pos) {
                stm << "p";
            } else if (as.sign == Sign::NonNeg) {
                stm << "nn";
            } else {
                stm << "u";
            }
            stm << std::endl;
        }
    }
    stm << "  Branch:" << std::endl;
    for (auto &br: m_branches) {
        stm << "    ";
        br.cond->print(stm, false, true);
        stm << ": ";
        if (br.target) {
            stm << "BS(" << br.target->id() << ")";
        }
        else {
            stm << "end";
        }
        stm << std::endl;
    }
    stm << "    default: ";
    if (m_default_branch) {
        stm << "BS(" << m_default_branch->id() << ")";
    }
    else {
        stm << "end";
    }
    stm << std::endl;
    if (!m_endtimes.empty()) {
        stm << "  End Times:" << std::endl;
        for (auto &time: m_endtimes) {
            stm << "    " << *time << std::endl;
        }
    }
}

NACS_EXPORT() bool BasicSeq::needs_oldval(Pulse *p) const
{
    bool found = false;
    for (auto &[chn, info]: m_channels) {
        bool found_in_info = false;
        for (auto &pulse: info.pulses) {
            if (&pulse == p) {
                assert(!found);
                assert(!found_in_info);
                found_in_info = true;
            }
        }
        if (found_in_info)
            found = true;
        for (auto &[pulse, val]: info.pulse_start_vars) {
            assert(val);
            if (pulse == p) {
                assert(found_in_info);
                return true;
            }
        }
    }
    assert(found);
    (void)found;
    return false;
}

}
