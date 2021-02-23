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

#include "../utils/llvm/codegen.h"
#include "../utils/number.h"

#include <algorithm>
#include <set>
#include <utility>

namespace NaCs::Seq {

// Used only by `Seq::add_basicseq`
BasicSeq::BasicSeq(uint32_t id)
    : m_id(id)
{
    assert(id);
}

NACS_EXPORT() Pulse *BasicSeq::add_pulse(uint32_t chn, uint32_t id,
                                         EventTime &start, Var *len, Var *val)
{
    return &m_channels[chn].pulses.emplace_back(id, start, len, val, false);
}

NACS_EXPORT() Pulse *BasicSeq::add_measure(uint32_t chn, uint32_t id,
                                           EventTime &start, Var *val)
{
    assert(val->is_extern());
    assert((val->get_extern().second >> 32) == m_id);
    return &m_channels[chn].pulses.emplace_back(id, start, nullptr, val, true);
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

NACS_EXPORT() void BasicSeq::add_assume(EventTime::Sign sign, Var *val, uint32_t assume_id)
{
    assert(sign == EventTime::Pos || sign == EventTime::NonNeg);
    m_assume.emplace(val->ref(), Assumption{sign, assume_id});
}

NACS_EXPORT() void BasicSeq::add_endtime(EventTime &t)
{
    m_endtimes.push_back(t.ref());
}

NACS_EXPORT() EventTime &BasicSeq::track_time(const EventTime &t)
{
    return m_eventtimes.emplace_back(t);
}

NACS_EXPORT() void BasicSeq::check() const
{
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
                    throw Error(uint32_t(Error::BasicSeq) << 16 | ExternMeasure,
                                Error::Measure, var->get_extern().second,
                                Error::BasicSeq, seq.id(),
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
                    if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                        error.type2 = Error::EventTime;
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
                if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                    error.code = uint32_t(Error::BasicSeq) << 16 | ExternMeasureLength;
                    error.type2 = Error::Pulse;
                    error.id2 = pulse.id();
                }
                throw error;
            }
            try {
                visitor.scan_deps(val, dep);
            }
            catch (Error &error) {
                if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                    error.type2 = Error::Pulse;
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
            if (it->second->known_before(*pulse) != EventTime::Pos) {
                if (pulse->is_measure())
                    throw Error(uint32_t(Error::BasicSeq) << 16 | MeasureOrder,
                                Error::Measure, var->get_extern().second,
                                Error::Measure, pulse->val()->get_extern().second,
                                "Use of measure must be later than the measure.");
                throw Error(uint32_t(Error::BasicSeq) << 16 | MeasureOrder,
                            Error::Measure, var->get_extern().second,
                            Error::Pulse, pulse->id(),
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
            if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                error.type2 = Error::Assignment;
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
            if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                error.type2 = Error::EventTime;
                error.id2 = as.id;
            }
            throw error;
        }
    }
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
            if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                error.type2 = Error::Branch;
                error.id2 = br.id;
            }
            throw error;
        }
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
                if (error.code == (uint32_t(Error::BasicSeq) << 16 | ExternMeasure)) {
                    error.type2 = Error::EventTime;
                    error.id2 = term.id;
                }
                throw error;
            }
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
    if (m_default_branch)
        m_default_branch->mark_recursive();
    for (auto &br: m_branches) {
        if (br.target) {
            br.target->mark_recursive();
        }
    }
}

bool BasicSeq::preoptimize_pulse(uint32_t chn, Env &env)
{
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
        }
        ++it;
    }
    auto &pulse_start_vars = info.pulse_start_vars;
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
        if (t1.tconst < t2.tconst)
            return true;
        if (t1.tconst > t2.tconst)
            return false;
        auto nt1 = t1.terms.size();
        auto nt2 = t2.terms.size();
        if (nt1 < nt2)
            return true;
        if (nt1 > nt2)
            return false;
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
        // assume `p1` was sorted before `p2`;
        auto &t1 = p1->start();
        auto &t2 = p2->start();
        assert(t1.tconst <= t2.tconst);
        auto isless = t1.isless_terms(t2);
        if (isless == EventTime::NonNeg) {
            if (t1.tconst < t2.tconst)
                return true;
            return p1->id() < p2->id();
        }
        return isless == EventTime::Pos;
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
        auto val = pulse->endval();
        assert(val); // since we know `!pulse->needs_oldval()`
        use_val(val, measure);
    };
    auto replace_measure = [&] (Pulse *pulse, Pulse *measure) {
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
                fmf.setNoNaNs();
                fmf.setNoInfs();
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
    // size has constant complexity
    bool changed = false;
    if (m_endtimes.size() > 1)
        m_endtimes.sort([] (const auto &t1, const auto &t2) {
            if (t1->tconst < t2->tconst)
                return true;
            if (t1->tconst > t2->tconst)
                return false;
            auto nt1 = t1->terms.size();
            auto nt2 = t2->terms.size();
            if (nt1 < nt2)
                return true;
            if (nt1 > nt2)
                return false;
            return false;
        });
    auto check_time = [&] (auto tit) {
        for (auto it = m_endtimes.begin(); it != tit;) {
            assert((*it)->tconst <= (*tit)->tconst);
            if ((*it)->isless_terms(**tit) != EventTime::Unknown) {
                changed = true;
                it = m_endtimes.erase(it);
            }
            else {
                ++it;
            }
        }
    };
    for (auto it = m_endtimes.begin(), end = m_endtimes.end(); it != end;) {
        if ((*it)->terms.size() == 0 && (*it)->tconst <= 0) {
            it = m_endtimes.erase(it);
        }
        else {
            check_time(it);
            ++it;
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
        if (it->m_ref_count == 0) {
            changed = true;
            it = m_eventtimes.erase(it);
        }
        else {
            ++it;
        }
    }
    return changed;
}

bool BasicSeq::optimize_vars()
{
    for (auto &[chn, info]: m_channels) {
        if (auto v = info.startval->get_assigned_var())
            info.startval.reset(v);
        if (info.endval) {
            if (auto v = info.endval->get_assigned_var()) {
                info.endval.reset(v);
            }
        }
    }
    for (auto &kv: m_assign) {
        assert(kv.second.val);
        if (auto v = kv.second.val->get_assigned_var()) {
            kv.second.val.reset(v);
        }
    }
    auto check_assumption = [&] (Var *var, auto &as) {
        assert(var->is_const());
        auto v = var->get_const().get<double>();
        if (as.sign == EventTime::Pos && !(round<int64_t>(v) > 0)) {
            throw Error(uint32_t(Error::EventTime) << 16 | EventTime::NonPosTime,
                        Error::EventTime, as.id, "Positive time expected.");
        }
        else if (as.sign == EventTime::NonNeg && !(v >= 0)) {
            throw Error(uint32_t(Error::EventTime) << 16 | EventTime::NegTime,
                        Error::EventTime, as.id, "Non-negative time expected.");
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
                if (!inserted && as.sign == EventTime::Pos &&
                    it2->second.sign != EventTime::Pos) {
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
    return preoptimize_eventtimes();
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
            if (as.sign == EventTime::Pos) {
                stm << "p";
            } else if (as.sign == EventTime::NonNeg) {
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
