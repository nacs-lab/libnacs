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

#include "compiler.h"

#include "host_seq.h"
#include "seq.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/llvm/compile.h"
#include "../nacs-utils/llvm/execute.h"
#include "../nacs-utils/llvm/passes.h"
#include "../nacs-utils/llvm/utils.h"
#include "../nacs-utils/number.h"

#include <algorithm>
#include <numeric>

namespace NaCs::Seq {

namespace {

template<typename V>
struct ValCache {
    std::pair<V&,bool> get(IR::Type type, IR::GenVal value)
    {
        switch (type) {
        case IR::Type::Bool:
            return load_bool(value.b);
        case IR::Type::Int32:
            return load(m_i32_cache, value.i32);
        case IR::Type::Float64:
            return load(m_f64_cache, value.f64);
        default:
            return {m_dummy, false};
        }
    }

    std::pair<V&,bool> get(IR::TagVal value)
    {
        return get(value.typ, value.val);
    }

private:
    template<typename K>
    static std::pair<V&,bool> load(std::map<K,V> &cache, K key)
    {
        auto [it, inserted] = cache.try_emplace(key);
        return {it->second, !inserted};
    }
    std::pair<V&,bool> load_bool(bool _key)
    {
        uint32_t key = _key;
        auto cached = m_b_cached[key];
        m_b_cached[key] = true;
        return {m_b_cache[key], cached};
    }

    V m_dummy;
    bool m_b_cached[2] = {false, false};
    V m_b_cache[2];
    std::map<int32_t,V> m_i32_cache;
    std::map<double,V> m_f64_cache;
};

static inline HostSeq::Type ir2hostseq(IR::Type type)
{
    switch (type) {
    case IR::Type::Bool:
        return HostSeq::Type::Bool;
    case IR::Type::Int32:
        return HostSeq::Type::Int32;
    case IR::Type::Float64:
        return HostSeq::Type::Float64;
    default:
        return HostSeq::Type::Bool;
    }
}

static inline std::pair<HostSeq::Type,HostSeq::Value> ir2hostseq(IR::TagVal value)
{
    std::pair<HostSeq::Type,HostSeq::Value> res;
    switch (value.typ) {
    case IR::Type::Bool:
        res.first = HostSeq::Type::Bool;
        res.second.b = value.val.b;
        break;
    case IR::Type::Int32:
        res.first = HostSeq::Type::Int32;
        res.second.i32 = value.val.i32;
        break;
    case IR::Type::Float64:
        res.first = HostSeq::Type::Float64;
        res.second.f64 = value.val.f64;
        break;
    default:
        res.first = HostSeq::Type::Bool;
        res.second.b = false;
        break;
    }
    return res;
}

}

struct Compiler::TimeVal {
    llvm::SmallVector<uint32_t,4> vars; // varid
    llvm::SmallVector<uint32_t,4> timevals; // timeval id
    int64_t offset;
};

NACS_EXPORT_ Compiler::Compiler(HostSeq &host_seq, Seq &seq, LLVM::Exe::Engine &engine,
                                std::function<uintptr_t(const std::string&)> resolver)
    : m_host_seq(host_seq),
      m_seq(seq),
      m_engine(engine),
      m_resolver(resolver)
{
}

NACS_EXPORT() Compiler::~Compiler()
{
}

NACS_EXPORT() llvm::Function *Compiler::get_ramp_func(Var *val, llvm::Function *f,
                                                      llvm::Type *ret_ty,
                                                      llvm::Type *time_ty, bool _export)
{
    auto args = val->args();
    uint32_t nargs = args.size();
    auto check_func_type = [&] {
        // The optimizer makes sure that all the arguments types
        // are consistent with the input values.
        // However, it doesn't make sure the types are exactly what we want.
        if (f->getReturnType() != ret_ty)
            return false;
        for (uint32_t argi = 0; argi < nargs; argi++) {
            auto arg = args[argi];
            if (arg.is_arg()) {
                assert(arg.get_arg() == 0);
                if (f->getArg(argi)->getType() != time_ty) {
                    return false;
                }
            }
        }
        return true;
    };
    if (check_func_type())
        return f;
    llvm::SmallVector<llvm::Type*,8> arg_types(nargs);
    auto oldfty = f->getFunctionType();
    for (uint32_t i = 0; i < nargs; i++) {
        auto arg = args[i];
        if (arg.is_arg()) {
            arg_types[i] = time_ty;
        }
        else {
            arg_types[i] = oldfty->getParamType(i);
        }
    }
    auto fty = llvm::FunctionType::get(ret_ty, arg_types, false);
    auto newf = llvm::Function::Create(fty,  _export ? llvm::GlobalValue::ExternalLinkage :
                                       llvm::GlobalValue::PrivateLinkage,
                                       f->getName() + ".r", f->getParent());
    if (_export)
        newf->setVisibility(llvm::GlobalValue::ProtectedVisibility);
    newf->setAttributes(llvm::AttributeList::get(
                            f->getContext(), LLVM::getFnAttrs(*f), {}, {}));
    auto *b0 = llvm::BasicBlock::Create(f->getContext(), "top", newf);
    llvm::IRBuilder<> builder(b0);
    llvm::SmallVector<llvm::Value*,8> call_args(nargs);
    for (uint32_t i = 0; i < nargs; i++) {
        auto arg = args[i];
        llvm::Value *larg = newf->getArg(i);
        if (arg.is_arg())
            larg = LLVM::convert_scalar(builder,
                                        f->getArg(i)->getType(), larg);
        call_args[i] = larg;
    }
    auto call = builder.CreateCall(f, call_args);
    builder.CreateRet(LLVM::convert_scalar(builder, ret_ty, call));
    return newf;
}

void Compiler::prescan_vars()
{
    Env &env = m_seq.env();
    uint32_t nvars = env.num_vars();

    m_all_vars.resize(nvars);
    m_var_types.resize(nvars, VarType::None);
    m_var_owners.resize(nvars, 0);
    m_var_slots.resize(nvars, 0);
    m_var_isramps.resize(nvars, false);

    ValCache<uint32_t> const_ids;
    uint32_t nglobal_vals = 0;

    // The only way to distinguish between a measure and a sequence start value
    // is by looking at where it is referred to in the sequence.
    // Therefore, we need to mark all the `Channel` type before we scan the rest of the `Var`s.
    for (auto &bseq: m_seq.get_basicseqs()) {
        m_bseqs.push_back(&bseq);
        auto bseq_id = bseq.id();
        for (auto &[chn, info]: bseq.get_channel_infos()) {
            auto startval = info.startval.get();
            if (startval->is_const())
                continue;
            assert(startval->is_extern());
            assert((startval->get_extern().second >> 32) == bseq_id);
            uint32_t varid = startval->varid();
            assert(varid < nvars);
            m_var_types[varid] = VarType::Channel;
            m_var_owners[varid] = bseq_id;
        }
    }

    // Store the vars in topological order
    for (auto var: env) {
        assert((uint32_t)var->varid() < nvars);
        m_all_vars[var->varid()] = var;
    }

    for (uint32_t varid = 0; varid < nvars; varid++) {
        auto var = m_all_vars[varid];
        if (m_var_types[varid] != VarType::None) {
            assert(m_var_types[varid] == VarType::Channel);
            continue;
        }
        else if (var->is_const()) {
            m_var_types[varid] = VarType::Const;
            auto c = var->get_const();
            auto [slot, cached] = const_ids.get(c);
            if (!cached) {
                auto [typ, val] = ir2hostseq(c);
                slot = m_host_seq.types.size();
                m_host_seq.types.push_back(typ);
                m_host_seq.values.push_back(val);
            }
            m_var_slots[varid] = slot;
            continue;
        }
        else if (var->is_extern()) {
            auto ext_id = var->get_extern().second;
            auto bseq_id = uint32_t(ext_id >> 32);
            if (bseq_id) {
                // basic sequence start values already marked as `Channel`
                m_var_types[varid] = VarType::BSeqMeasure;
                m_var_owners[varid] = bseq_id;
            }
            else {
                m_var_types[varid] = VarType::Global;
            }
            continue;
        }

        assert(var->is_call());
        assert(var->get_callee().is_llvm);
        uint32_t nargs = var->args().size();

        bool refglobal = false;
        bool refchannel = false;
        bool refmeasure = false;
        uint32_t owner = 0;

        for (uint32_t i = 0; i < nargs; i++) {
            auto arg = var->args()[i];
            if (arg.is_arg()) {
                m_var_isramps[varid] = true;
                continue;
            }
            else if (!arg.is_var()) {
                continue;
            }
            auto ref = arg.get_var();
            uint32_t refid = ref->varid();
            assert(refid < varid);
            auto type = m_var_types[refid];
            switch (type) {
            case VarType::Global:
            case VarType::GlobalValue:
                refglobal = true;
                break;
            case VarType::Channel:
            case VarType::BSeqDirect:
                refchannel = true;
                assert(owner == m_var_owners[refid] || owner == 0);
                owner = m_var_owners[refid];
                break;
            case VarType::BSeqMeasure:
            case VarType::BSeqNeedOrder:
                refmeasure = true;
                assert(owner == m_var_owners[refid] || owner == 0);
                owner = m_var_owners[refid];
                break;
            default:
                assert(false && "Invalid type");
            }
            if (refmeasure) {
                break;
            }
        }
        if (refmeasure) {
            m_var_owners[varid] = owner;
            m_var_types[varid] = VarType::BSeqNeedOrder;
        }
        else if (refchannel) {
            m_var_owners[varid] = owner;
            m_var_types[varid] = VarType::BSeqDirect;
        }
        else if (refglobal) {
            nglobal_vals++;
            m_var_types[varid] = VarType::GlobalValue;
        }
        else {
            assert(m_var_isramps[varid]);
            m_var_types[varid] = VarType::Const;
            uint32_t slot = m_host_seq.types.size();
            m_host_seq.types.push_back(ir2hostseq(var->type()));
            HostSeq::Value val;
            val.i64 = 0;
            m_host_seq.values.push_back(val);
            m_var_slots[varid] = slot;
        }
    }

    assert(m_host_seq.values.size() == m_host_seq.types.size());
    m_host_seq.nglobals = m_seq.get_slots().size();
    nglobal_vals_var = nglobal_vals;
    m_host_seq.nchannels = m_seq.get_chn_names().size();
}

namespace {

static uint32_t find_repeat(const std::vector<uint32_t> &ary)
{
    assert(!ary.empty());
    uint32_t nele = ary.size();
    if (nele == 1)
        return 1;
    uint32_t res = 0;
    uint32_t cur_val = ary[0];
    uint32_t cur_rep = 1;
    for (uint32_t i = 1; i < nele; i++) {
        auto val = ary[i];
        if (val == cur_val) {
            cur_rep += 1;
            continue;
        }
        res = std::gcd(res, cur_rep);
        if (res == 1)
            return 1;
        cur_val = val;
        cur_rep = 1;
    }
    return std::gcd(res, cur_rep);
}

struct TermKeyCompare {
    bool operator()(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b) const
    {
        // Sort ones with more use in front
        if (a.size() > b.size())
            return true;
        if (a.size() < b.size())
            return false;
        // Fallback to default order for ones with the same size
        return a < b;
    }
};

}

void Compiler::prescan_times()
{
    std::map<int64_t,uint32_t> cache;
    std::map<const Var*,std::vector<uint32_t>> varuse;

    // First assign all users of the i64 values an ID
    // Sort all the times before all the assumptions.
    // In the mean time, collect all the users of each i64 terms.
    // For event time that is purely constant, translate to a slot directly
    for (auto &bseq: m_seq.get_basicseqs()) {
        for (auto &et: bseq.get_eventtimes()) {
            if (et.terms.empty()) {
                auto t0 = et.tconst;
                auto [it, inserted] = cache.try_emplace(t0);
                if (inserted) {
                    it->second = m_host_seq.types.size();
                    HostSeq::Value v;
                    v.i64 = t0;
                    m_host_seq.types.push_back(HostSeq::Type::Int64);
                    m_host_seq.values.push_back(v);
                }
                m_time_slots.emplace(&et, it->second);
                continue;
            }
            uint32_t timeid = m_event_times.size();
            m_event_times.push_back(&et);
            m_event_time_ids.emplace(&et, timeid);
            for (auto &term: et.terms) {
                varuse[term.var.get()].push_back(timeid);
            }
        }
    }

    uint32_t nevent_times = m_event_times.size();
    for (auto &bseq: m_seq.get_basicseqs()) {
        for (auto &assume: bseq.get_assumes()) {
            uint32_t id = m_assumptions.size() + nevent_times;
            m_assumptions.push_back(&assume);
            m_assumption_ids.emplace(&assume, id);
            varuse[assume.first.get()].push_back(id);
        }
    }

    uint32_t nassumes = m_assumptions.size();
    uint32_t nresults = nassumes + nevent_times;

    // Now we can group the terms into part that will always be evaluated together
    // based on their use.
    // Sort the terms with more use to the front
    // so that they'll also appear in the front when we decompose
    // each final result into terms.
    // We use this as a heuristic to increase the chance of subexpression reuse.
    struct TermInfo {
        llvm::SmallVector<const Var*,2> vars;
        // How many IntermediateResult refers to this term.
        uint32_t refcount = 0;
    };
    std::map<std::vector<uint32_t>,TermInfo,TermKeyCompare> terms;
    for (auto &[var, use]: varuse) {
        // Reduce the use array if possible
        // and turn that into repetition of the variable itself so that we can find
        // other similar varibles that should be in the same group with different repetitions.
        // i.e. we'd like to turn
        //
        //   v1: [1, 1, 2, 2]
        //   v2: [1, 1, 1, 2, 2, 2]
        //
        // in `varuse` into
        //
        //   [1, 2]: [v1, v1, v2, v2, v2]
        //
        // in `terms`
        auto rep = find_repeat(use);
        if (rep == 1) {
            terms[use].vars.push_back(var);
            continue;
        }
        assert(use.size() % rep == 0);
        uint32_t nuse = use.size() / rep;
        std::vector<uint32_t> reduce_use(nuse);
        for (uint32_t i = 0; i < nuse; i++)
            reduce_use[i] = use[i * rep];
        terms[reduce_use].vars.append(rep, var);
    }

    // From this point on we'll use the keys in `terms` to identify the terms
    using term_key_t = const std::vector<uint32_t>*;
    struct ResultInfo {
        uint32_t id;
        std::vector<term_key_t> terms;
        void add_term(term_key_t term)
        {
            terms.push_back(term);
        }
    };
    std::vector<ResultInfo> results(nresults);
    for (uint32_t i = 0; i < nresults; i++)
        results[i].id = i;
    for (auto &[uses, term]: terms) {
        auto id = &uses;
        // The number of repetition in `use` is the number of time each result uses this term.
        for (auto use: uses) {
            results[use].add_term(id);
        }
    }
    // Sort results with fewer terms to the front so that
    // Stable sort here is just to get more predictable result.
    std::stable_sort(results.begin(), results.end(), [&] (const auto &v1, const auto &v2) {
        return v1.terms.size() < v2.terms.size();
    });

    // Now we've identified all the terms, next we need to figure out a good way to evaluate
    // all the results while reusing intermediate results as much as we can.
    // I'm not sure how to find the optimal way, but using a heuristic-based approach
    // (similar to what LLVM does with `reassociate` followed by `cse` passes)
    // we'll sort the terms in each result in the same order and with the most used terms
    // in the front.
    // This should improve the chance of finding common subexpressions and
    // in particular this should cover the "time chain" case that we care about
    // (i.e. a series of `t(n) = t(n-1) + term`)
    // We'll also limit our reduction to have each intermediate result being able to use
    // at most one other intermediate result for simplicity.

    struct IntermediateResult {
        IntermediateResult *parent;
        llvm::SmallVector<term_key_t,2> terms;
        uint32_t refcount = 1;
        uint32_t extrefcount = 0;
        int64_t offset = 0;
        bool offset_set = false;
    };
    std::map<std::vector<term_key_t>,IntermediateResult,std::less<>> intermediates;
    using lookup_func_t = std::function<IntermediateResult*(term_key_t*,size_t)>;
    lookup_func_t create_intermediate = [&] (term_key_t *term_keys, size_t nterms) ->
        IntermediateResult* {
        if (nterms == 0)
            return nullptr;
        ArrayKey<term_key_t> key{term_keys, nterms};
        auto it = intermediates.find(key);
        if (it != intermediates.end()) {
            auto res = &it->second;
            res->refcount++;
            return res;
        }
        auto parent = create_intermediate(term_keys, nterms - 1);
        auto term_key = term_keys[nterms - 1];
        it = intermediates.emplace(std::vector<term_key_t>{term_keys, term_keys + nterms},
                                   IntermediateResult{parent, {term_key}}).first;
        auto term_it = terms.find(*term_key);
        assert(term_it != terms.end());
        term_it->second.refcount += 1;
        return &it->second;
    };

    std::vector<IntermediateResult*> result_intermediates(results.size());
    // First trying to identify all the common subexressions by breaking down each computation
    // to subexpressions and reuse the subexpressions as we go.
    for (auto &result: results) {
        // Assumptions must have only one term.
        if (result.id >= nevent_times)
            assert(result.terms.size() == 1);
        auto inter = create_intermediate(result.terms.data(), result.terms.size());
        assert(result.id < result_intermediates.size());
        result_intermediates[result.id] = inter;
        inter->extrefcount++;
    }
    // Now we got all the terms right, it's time to clean up the unneeded intermediate results
    // We know each terms are used more than once or has external use.
    // However, this doesn't mean each subexpression are used more than once.
    // The subexpressions that has a single user can be inlined into parent.
    // Since each subexpression can only have at most one use of other subexpression,
    // we know that diamond shape dependencies like
    //   v3 - v2 - v1
    //   /    /
    // *v5 - v4
    // cannot exist so the elimination of intermediate can be done with only local information.
    // The lexicographical order guarantees that the parent is sorted before the child
    // so if we iterate the map normally it is guaranteed that the parent is already optimized
    // before we reach the child.
    for (auto &[_term, intermediate]: intermediates) {
        auto parent = intermediate.parent;
        if (intermediate.refcount == 0 || !parent)
            continue;
        if (parent->extrefcount)
            assert(parent->refcount > 1);
        if (parent->refcount > 1)
            continue;
        // Take the memory of parent terms which is likely larger.
        parent->terms.append(intermediate.terms.begin(), intermediate.terms.end());
        intermediate.terms = std::move(parent->terms);
        intermediate.parent = parent->parent;
        // Now I want to delete `parent` but I don't have the iterator to do that.
        // Instead of looking it up here, let's just mark it for deletion and
        // it'll be deleted below when we iterate through it next time.
        parent->refcount = 0;
    }

    // Now we've finalized all the variable terms, we can now assign the constant offset.
    // For other ones, we simply pick the smallest offset of all the users.
    // Note that unlike the time offset in `TimeVal`,
    // the offset here is not cumulative but is rather the total offset for this subexpression.
    for (uint32_t i = 0; i < nevent_times; i++) {
        assert(i < result_intermediates.size());
        auto *intermediate = result_intermediates[i];
        if (intermediate->terms.size() == 1 && !intermediate->parent) {
            auto term_key = intermediate->terms[0];
            auto term_it = terms.find(*term_key);
            assert(term_it != terms.end());
            if (term_it->second.refcount > 1) {
                // A term with refcount > 1 must be used in an increment term
                // in an intermediate result, we should assign offset 0.
                intermediate->offset_set = true;
                intermediate->offset = 0;
                continue;
            }
        }
        assert(i < m_event_times.size());
        if (!intermediate->offset_set) {
            intermediate->offset_set = true;
            intermediate->offset = m_event_times[i]->tconst;
        }
        else {
            intermediate->offset = min(intermediate->offset, m_event_times[i]->tconst);
        }
    }
    for (uint32_t i = 0; i < nassumes; i++) {
        assert(i + nevent_times < result_intermediates.size());
        auto *intermediate = result_intermediates[i + nevent_times];
        assert(intermediate->terms.size() == 1 && !intermediate->parent);
        intermediate->offset_set = true;
        intermediate->offset = 0;
    }

    // Now create the final variables.
    std::map<term_key_t,uint32_t> term_slots;
    // Assign each terms first.
    // The subexpressions in `intermediates` are sorted in their dependencies
    // on each other but not with respect to the terms.
    for (auto &[uses, term]: terms) {
        // Do not create the term if it has a single use
        // Delay it to the actual user so that we can inline it and also
        // know what the final offset is.
        if (term.refcount == 1)
            continue;
        assert(term.refcount > 1);
        assert(uses.size() > 1);
        auto id = &uses;
        auto slot = m_timevals.size();
        auto &timeval = m_timevals.emplace_back(TimeVal{{}, {}, 0});
        for (auto var: term.vars)
            timeval.vars.push_back(var->varid());
        term_slots.emplace(id, slot);
    }
    std::map<IntermediateResult*,uint32_t> intermediate_slots;
    // Now assign the intermediates
    for (auto _it = intermediates.begin(), _end = intermediates.end(); _it != _end;) {
        auto &intermediate = _it->second;
        if (intermediate.refcount == 0) {
            // Marked for deletion from above. Do it now.
            _it = intermediates.erase(_it);
            continue;
        }
        ++_it;
        if (intermediate.terms.size() == 1 && !intermediate.parent) {
            auto term_key = intermediate.terms[0];
            auto term_it = terms.find(*term_key);
            assert(term_it != terms.end());
            if (term_it->second.refcount == 1) {
                // This is a single used term which we skipped above,
                // fall through to create it below.
                assert(term_slots.find(term_key) == term_slots.end());
                if (!intermediate.offset_set) {
                    assert(intermediate.offset == 0);
                    intermediate.offset_set = true;
                }
            }
            else {
                assert(intermediate.offset == 0);
                intermediate.offset_set = true;
                // This is a single term, use the slot for that term
                intermediate_slots.emplace(&intermediate, find_val(term_slots, term_key));
                continue;
            }
        }
        int64_t parent_offset = intermediate.parent ? intermediate.parent->offset : 0;
        assert(!intermediate.parent || intermediate.parent->offset_set);
        int64_t rel_offset;
        if (!intermediate.offset_set) {
            intermediate.offset = parent_offset;
            intermediate.offset_set = true;
            rel_offset = 0;
        }
        else {
            rel_offset = intermediate.offset - parent_offset;
        }
        auto slot = m_timevals.size();
        auto &timeval = m_timevals.emplace_back(TimeVal{{}, {}, rel_offset});
        if (intermediate.parent)
            timeval.timevals.push_back(find_val(intermediate_slots, intermediate.parent));
        for (auto term: intermediate.terms) {
            auto it = term_slots.find(term);
            if (it != term_slots.end()) {
                timeval.timevals.push_back(find_val(term_slots, term));
                continue;
            }
            auto term_it = terms.find(*term);
            assert(term_it != terms.end());
            assert(term_it->second.refcount == 1);
            for (auto var: term_it->second.vars) {
                timeval.vars.push_back(var->varid());
            }
        }
        intermediate_slots.emplace(&intermediate, slot);
    }
    m_time_timevals.resize(nresults);
    // All the intermediate variables are created.
    // Find the corresponding one for each users
    // and create new ones with a different offset if necessary.
    for (uint32_t i = 0; i < nevent_times; i++) {
        assert(i < result_intermediates.size());
        assert(i < m_event_times.size());
        assert(i < m_time_timevals.size());
        auto intermediate = result_intermediates[i];
        auto et = m_event_times[i];
        assert(intermediate->offset_set);
        auto intermediate_slot = find_val(intermediate_slots, intermediate);
        auto rel_offset = et->tconst - intermediate->offset;
        if (rel_offset == 0) {
            m_time_timevals[i] = intermediate_slot;
            continue;
        }
        auto slot = m_timevals.size();
        m_timevals.emplace_back(TimeVal{{}, {intermediate_slot}, rel_offset});
        m_time_timevals[i] = slot;
    }
    for (uint32_t i = 0; i < nassumes; i++) {
        assert(i + nevent_times < result_intermediates.size());
        assert(i + nevent_times < m_time_timevals.size());
        auto intermediate = result_intermediates[i + nevent_times];
        assert(intermediate->offset == 0);
        auto intermediate_slot = find_val(intermediate_slots, intermediate);
        m_time_timevals[i + nevent_times] = intermediate_slot;
    }
    assert(m_host_seq.values.size() == m_host_seq.types.size());
}

void Compiler::classify_times()
{
    uint32_t ntimevals = m_timevals.size();
    uint32_t nglobal_vals = 0;
    m_timeval_types.resize(ntimevals, VarType::None);
    m_timeval_owners.resize(ntimevals, 0);
    m_timeval_slots.resize(ntimevals, 0);

    for (uint32_t timevalid = 0; timevalid < ntimevals; timevalid++) {
        auto &timeval = m_timevals[timevalid];
        bool refglobal = false;
        bool refchannel = false;
        bool refmeasure = false;
        uint32_t owner = 0;
        for (auto refid: timeval.vars) {
            assert(refid < m_var_types.size());
            auto type = m_var_types[refid];
            switch (type) {
            case VarType::Global:
            case VarType::GlobalValue:
                refglobal = true;
                break;
            case VarType::Channel:
            case VarType::BSeqDirect:
                refchannel = true;
                assert(owner == m_var_owners[refid] || owner == 0);
                owner = m_var_owners[refid];
                break;
            case VarType::BSeqMeasure:
            case VarType::BSeqNeedOrder:
                refmeasure = true;
                assert(owner == m_var_owners[refid] || owner == 0);
                owner = m_var_owners[refid];
                break;
            default:
                assert(false && "Invalid type");
            }
            if (refmeasure) {
                break;
            }
        }
        if (!refmeasure) {
            for (auto refid: timeval.timevals) {
                assert(refid < m_timeval_types.size());
                auto type = m_timeval_types[refid];
                switch (type) {
                case VarType::Global:
                case VarType::GlobalValue:
                    refglobal = true;
                    break;
                case VarType::Channel:
                case VarType::BSeqDirect:
                    refchannel = true;
                    assert(owner == m_timeval_owners[refid] || owner == 0);
                    owner = m_timeval_owners[refid];
                    break;
                case VarType::BSeqMeasure:
                case VarType::BSeqNeedOrder:
                    refmeasure = true;
                    assert(owner == m_timeval_owners[refid] || owner == 0);
                    owner = m_timeval_owners[refid];
                    break;
                default:
                    assert(false && "Invalid type");
                }
                if (refmeasure) {
                    break;
                }
            }
        }
        if (refmeasure) {
            m_timeval_owners[timevalid] = owner;
            m_timeval_types[timevalid] = VarType::BSeqNeedOrder;
        }
        else if (refchannel) {
            m_timeval_owners[timevalid] = owner;
            m_timeval_types[timevalid] = VarType::BSeqDirect;
        }
        else {
            assert(refglobal);
            (void)refglobal;
            nglobal_vals++;
            m_timeval_types[timevalid] = VarType::GlobalValue;
        }
    }
    nglobal_vals_time = nglobal_vals;
}

void Compiler::assign_shared_vars()
{
    uint32_t nvars = m_all_vars.size();
    global_offset = m_host_seq.nconsts;
    global_value_offset = global_offset + m_host_seq.nglobals;
    channels_offset = global_value_offset + m_host_seq.nglobal_vals;
    uint32_t global_value_id = 0;
    m_host_seq.types.resize(channels_offset);
    m_host_seq.depends.resize(m_host_seq.nglobal_vals);
    m_global_vals_inter_deps.resize(m_host_seq.nglobal_vals);
    std::set<uint32_t> global_val_deps;
    std::set<uint32_t> global_val_inter_deps;

    for (uint32_t varid = 0; varid < nvars; varid++) {
        auto var = m_all_vars[varid];
        switch (m_var_types[varid]) {
        case VarType::Global: {
            auto slot = global_offset + uint32_t(var->get_extern().second);
            assert(var->is_extern());
            assert(slot < m_host_seq.types.size());
            m_host_seq.types[slot] = ir2hostseq(var->get_extern().first);
            m_var_slots[varid] = slot;
            break;
        }
        case VarType::GlobalValue: {
            auto slot = global_value_offset + global_value_id;
            assert(slot < m_host_seq.types.size());
            m_host_seq.types[slot] = ir2hostseq(var->type());
            m_var_slots[varid] = slot;
            assert(var->is_call());
            assert(var->get_callee().is_llvm);
            uint32_t nargs = var->args().size();

            for (uint32_t i = 0; i < nargs; i++) {
                auto ref = var->get_ref(i);
                if (!ref)
                    continue;
                uint32_t refid = ref->varid();
                assert(refid < varid);
                auto refslot = m_var_slots[refid];
                auto type = m_var_types[refid];
                assert(type == VarType::Global || type == VarType::GlobalValue);
                if (type == VarType::Global) {
                    assert(refslot >= global_offset);
                    global_val_deps.insert(refslot - global_offset);
                }
                else {
                    assert(refslot >= global_value_offset);
                    auto ref_val_id = refslot - global_value_offset;
                    assert(ref_val_id < m_host_seq.depends.size());
                    auto &refdeps = m_host_seq.depends[ref_val_id];
                    global_val_deps.insert(refdeps.begin(), refdeps.end());
                    assert(ref_val_id < m_global_vals_inter_deps.size());
                    const auto &ref_interdeps =
                        m_global_vals_inter_deps[ref_val_id];
                    global_val_inter_deps.insert(ref_val_id);
                    global_val_inter_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
                }
            }
            assert(!global_val_deps.empty());
            assert(global_value_id < m_host_seq.depends.size());
            auto &deps = m_host_seq.depends[global_value_id];
            deps.insert(deps.end(), global_val_deps.begin(), global_val_deps.end());
            global_val_deps.clear();
            assert(global_value_id < m_global_vals_inter_deps.size());
            auto &interdeps = m_global_vals_inter_deps[global_value_id];
            interdeps.insert(interdeps.end(), global_val_inter_deps.begin(),
                             global_val_inter_deps.end());
            global_val_inter_deps.clear();
            global_value_id++;
            break;
        }
        default:
            break;
        }
    }

    for (auto &bseq: m_seq.get_basicseqs()) {
        auto bseq_id = bseq.id();
        for (auto &[chn, info]: bseq.get_channel_infos()) {
            auto startval = info.startval.get();
            if (startval->is_const())
                continue;
            assert(startval->is_extern());
            assert((startval->get_extern().second >> 32) == bseq_id);
            (void)bseq_id;
            m_var_slots[startval->varid()] = channels_offset + chn - 1;
        }
    }
}

void Compiler::assign_shared_times()
{
    uint32_t ntimevals = m_timevals.size();
    uint32_t global_value_id = nglobal_vals_var;
    std::set<uint32_t> global_val_deps;
    std::set<uint32_t> global_val_inter_deps;

    for (uint32_t timevalid = 0; timevalid < ntimevals; timevalid++) {
        auto &timeval = m_timevals[timevalid];
        if (m_timeval_types[timevalid] == VarType::GlobalValue) {
            auto slot = global_value_offset + global_value_id;
            assert(slot < m_host_seq.types.size());
            m_host_seq.types[slot] = HostSeq::Type::Int64;
            m_timeval_slots[timevalid] = slot;
            for (auto refid: timeval.vars) {
                assert(refid < m_var_slots.size());
                auto refslot = m_var_slots[refid];
                assert(refslot < slot);
                auto reftype = m_var_types[refid];
                assert(reftype == VarType::Global || reftype == VarType::GlobalValue);
                if (reftype == VarType::Global) {
                    assert(refslot >= global_offset);
                    global_val_deps.insert(refslot - global_offset);
                }
                else {
                    assert(refslot >= global_value_offset);
                    auto ref_val_id = refslot - global_value_offset;
                    assert(ref_val_id < m_host_seq.depends.size());
                    auto &refdeps = m_host_seq.depends[ref_val_id];
                    global_val_deps.insert(refdeps.begin(), refdeps.end());
                    assert(ref_val_id < m_global_vals_inter_deps.size());
                    const auto &ref_interdeps = m_global_vals_inter_deps[ref_val_id];
                    global_val_inter_deps.insert(ref_val_id);
                    global_val_inter_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
                }
            }
            for (auto refid: timeval.timevals) {
                assert(refid < m_timeval_slots.size());
                auto refslot = m_timeval_slots[refid];
                assert(refslot < slot);
                auto reftype = m_timeval_types[refid];
                assert(reftype == VarType::Global || reftype == VarType::GlobalValue);
                if (reftype == VarType::Global) {
                    assert(refslot >= global_offset);
                    global_val_deps.insert(refslot - global_offset);
                }
                else {
                    assert(refslot >= global_value_offset);
                    auto ref_val_id = refslot - global_value_offset;
                    assert(ref_val_id < m_host_seq.depends.size());
                    auto &refdeps = m_host_seq.depends[ref_val_id];
                    global_val_deps.insert(refdeps.begin(), refdeps.end());
                    assert(ref_val_id < m_global_vals_inter_deps.size());
                    const auto &ref_interdeps = m_global_vals_inter_deps[ref_val_id];
                    global_val_inter_deps.insert(ref_val_id);
                    global_val_inter_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
                }
            }
            assert(!global_val_deps.empty());
            assert(global_value_id < m_host_seq.depends.size());
            auto &deps = m_host_seq.depends[global_value_id];
            deps.insert(deps.end(), global_val_deps.begin(), global_val_deps.end());
            global_val_deps.clear();
            assert(global_value_id < m_global_vals_inter_deps.size());
            auto &interdeps = m_global_vals_inter_deps[global_value_id];
            interdeps.insert(interdeps.end(), global_val_inter_deps.begin(),
                             global_val_inter_deps.end());
            global_val_inter_deps.clear();
            global_value_id++;
        }
    }

}

void Compiler::assign_defaults()
{
    m_host_seq.default_values.resize(m_host_seq.nchannels);
    memset(m_host_seq.default_values.data(), 0,
           m_host_seq.nchannels * sizeof(HostSeq::Value));
    for (auto [chn, val]: m_seq.get_defvals()) {
        assert(chn > 0);
        assert(chn <= m_host_seq.default_values.size());
        m_host_seq.default_values[chn - 1].f64 = val;
    }
}

void Compiler::create_bseqs()
{
    uint32_t nbseqs = m_bseqs.size();
    assert(nbseqs == m_seq.get_basicseqs().size());
    m_host_seq.seqs.resize(nbseqs);
    for (uint32_t bseq_idx = 0; bseq_idx < nbseqs; bseq_idx++) {
        auto bseq = m_bseqs[bseq_idx];
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        auto [it, inserted] = m_bseq_idxs.emplace(bseq->id(), bseq_idx);
        assert(inserted); // ID should have no duplicate
        (void)inserted;
        host_bseq.id = bseq->id();
        host_bseq.branches.resize(bseq->get_branches().size());
    }
    // Branches
    auto get_target_idx = [&] (const BasicSeq *bseq) {
        if (!bseq)
            return uint32_t(-1);
        return find_val(m_bseq_idxs, bseq->id());
    };
    for (uint32_t bseq_idx = 0; bseq_idx < nbseqs; bseq_idx++) {
        auto bseq = m_bseqs[bseq_idx];
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        auto branches = bseq->get_branches();
        auto nbranches = branches.size();
        assert(host_bseq.branches.size() == nbranches);
        for (uint32_t i = 0; i < nbranches; i++)
            host_bseq.branches[i].target = get_target_idx(branches[i].target);
        host_bseq.default_branch = get_target_idx(bseq->get_default_branch());
    }
}

void Compiler::assign_bseq_vars()
{
    uint32_t nvars = m_all_vars.size();
    uint32_t ntimevals = m_timevals.size();
    std::set<uint32_t> global_val_deps;
    std::set<uint32_t> cond_global_val_deps;

    uint32_t nbseq_max = 0;
    uint32_t nbseqs = m_bseqs.size();
    // Count and scan dependencies on global values
    for (uint32_t bseq_idx = 0; bseq_idx < nbseqs; bseq_idx++) {
        auto bseq = m_bseqs[bseq_idx];
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        auto seq_id = host_bseq.id;
        uint32_t nmeasure = 0;
        uint32_t ndirect = 0;
        uint32_t nneed_order = 0;
        for (uint32_t varid = 0; varid < nvars; varid++) {
            auto owner = m_var_owners[varid];
            if (owner != seq_id)
                continue;
            switch (m_var_types[varid]) {
            case VarType::BSeqDirect:
                ndirect++;
                break;
            case VarType::BSeqMeasure:
                nmeasure++;
                continue;
            case VarType::BSeqNeedOrder:
                nneed_order++;
                break;
            default:
                continue;
            }
            auto var = m_all_vars[varid];
            assert(var->is_call());
            assert(var->get_callee().is_llvm);
            uint32_t nargs = var->args().size();
            for (uint32_t i = 0; i < nargs; i++) {
                auto ref = var->get_ref(i);
                if (!ref)
                    continue;
                uint32_t refid = ref->varid();
                assert(refid < varid);
                auto reftype = m_var_types[refid];
                if (reftype != VarType::GlobalValue)
                    continue;
                auto refslot = m_var_slots[refid];
                assert(refslot >= global_value_offset);
                auto ref_val_id = refslot - global_value_offset;
                assert(ref_val_id < m_global_vals_inter_deps.size());
                const auto &ref_interdeps = m_global_vals_inter_deps[ref_val_id];
                global_val_deps.insert(ref_val_id);
                global_val_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
            }
        }

        for (uint32_t timevalid = 0; timevalid < ntimevals; timevalid++) {
            auto owner = m_timeval_owners[timevalid];
            if (owner != seq_id)
                continue;
            switch (m_timeval_types[timevalid]) {
            case VarType::BSeqDirect:
                ndirect++;
                break;
            case VarType::BSeqNeedOrder:
                nneed_order++;
                break;
            default:
                continue;
            }
            auto &timeval = m_timevals[timevalid];
            for (auto refid: timeval.vars) {
                assert(refid < m_var_types.size());
                auto reftype = m_var_types[refid];
                if (reftype != VarType::GlobalValue)
                    continue;
                auto refslot = m_var_slots[refid];
                assert(refslot >= global_value_offset);
                auto ref_val_id = refslot - global_value_offset;
                assert(ref_val_id < m_global_vals_inter_deps.size());
                const auto &ref_interdeps = m_global_vals_inter_deps[ref_val_id];
                global_val_deps.insert(ref_val_id);
                global_val_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
            }
            for (auto refid: timeval.timevals) {
                assert(refid < m_timeval_types.size());
                auto reftype = m_timeval_types[refid];
                if (reftype != VarType::GlobalValue)
                    continue;
                auto refslot = m_timeval_slots[refid];
                assert(refslot >= global_value_offset);
                auto ref_val_id = refslot - global_value_offset;
                assert(ref_val_id < m_global_vals_inter_deps.size());
                const auto &ref_interdeps = m_global_vals_inter_deps[ref_val_id];
                global_val_deps.insert(ref_val_id);
                global_val_deps.insert(ref_interdeps.begin(), ref_interdeps.end());
            }
        }
        host_bseq.nmeasure = nmeasure;
        host_bseq.ndirect = ndirect;
        host_bseq.nneed_order = nneed_order;
        host_bseq.types.resize(ndirect + nneed_order);
        host_bseq.evals.resize(ndirect + nneed_order, nullptr);
        host_bseq.deps_count.resize(nneed_order);
        host_bseq.reverse_depends.resize(nmeasure);
        nbseq_max = max(nbseq_max, nmeasure + ndirect + nneed_order);

        auto add_dependency_global_val_id = [&] (uint32_t global_val_id) {
            const auto &interdeps = m_global_vals_inter_deps[global_val_id];
            global_val_deps.insert(global_val_id);
            global_val_deps.insert(interdeps.begin(), interdeps.end());
        };
        auto add_dependency = [&] (auto &&var) {
            auto varid = var->varid();
            auto type = m_var_types[varid];
            if (type != VarType::GlobalValue)
                return;
            add_dependency_global_val_id(m_var_slots[varid] - global_value_offset);
        };
        auto add_dependency_timeid = [&] (uint32_t timeid) {
            auto timevalid = m_time_timevals[timeid];
            if (m_timeval_types[timevalid] != VarType::GlobalValue)
                return;
            add_dependency_global_val_id(m_timeval_slots[timevalid] - global_value_offset);
        };
        auto add_dependency_time = [&] (const EventTime *et) {
            if (et->terms.empty())
                return;
            add_dependency_timeid(find_val(m_event_time_ids, et));
        };

        // We've collected implicit dependency of the basic sequence
        // on global values through basic sequence specific values.
        // We should also collect the direct dependencies
        // from the components of the basic sequence directly.
        for (auto &[global_id, assign]: bseq->get_assigns())
            add_dependency(assign.val);
        for (auto &assume: bseq->get_assumes())
            add_dependency_timeid(find_val(m_assumption_ids, &assume));
        for (auto &endtime: bseq->get_endtimes())
            add_dependency_time(endtime.get());
        for (auto &[chn, info]: bseq->get_channel_infos()) {
            for (auto &pulse: info.pulses) {
                add_dependency_time(&pulse.start());
                if (pulse.is_measure())
                    continue;
                auto len = pulse.len();
                if (len)
                    add_dependency(len);
                auto cond = pulse.cond();
                if (cond)
                    add_dependency(cond);
                add_dependency(pulse.val());
                add_dependency(pulse.endval());
            }
        }

        for (auto &br: bseq->get_branches()) {
            auto varid = br.cond->varid();
            auto type = m_var_types[varid];
            if (type != VarType::GlobalValue) {
                assert(type == VarType::Global);
                continue;
            }
            assert(m_var_slots[varid] >= global_value_offset);
            uint32_t global_val_id = m_var_slots[varid] - global_value_offset;
            assert(global_val_id < m_global_vals_inter_deps.size());
            const auto &interdeps = m_global_vals_inter_deps[global_val_id];
            cond_global_val_deps.insert(global_val_id);
            cond_global_val_deps.insert(interdeps.begin(), interdeps.end());
        }

        auto &deps = host_bseq.global_refs;
        deps.insert(deps.end(), global_val_deps.begin(), global_val_deps.end());
        global_val_deps.clear();
        auto &cond_deps = host_bseq.cond_global_refs;
        cond_deps.insert(cond_deps.end(), cond_global_val_deps.begin(),
                         cond_global_val_deps.end());
        cond_global_val_deps.clear();
    }
    m_host_seq.values.resize(m_host_seq.nshared + nbseq_max);

    // Assign slot, store types, scan dependency on measures
    for (uint32_t bseq_idx = 0; bseq_idx < nbseqs; bseq_idx++) {
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        auto seq_id = host_bseq.id;
        uint32_t measure_id = 0;
        uint32_t direct_id = 0;
        uint32_t need_order_id = 0;
        uint32_t need_order_offset = host_bseq.nmeasure + host_bseq.ndirect;
        std::vector<std::set<uint32_t>> measure_deps(host_bseq.nneed_order);
        for (uint32_t varid = 0; varid < nvars; varid++) {
            auto owner = m_var_owners[varid];
            if (owner != seq_id)
                continue;
            auto type = m_var_types[varid];
            auto var = m_all_vars[varid];
            uint32_t bseq_slot;
            if (type == VarType::BSeqMeasure) {
                m_var_slots[varid] = m_host_seq.nshared + measure_id++;
                continue;
            }
            else if (type == VarType::BSeqDirect) {
                bseq_slot = host_bseq.nmeasure + direct_id++;
            }
            else if (type == VarType::BSeqNeedOrder) {
                bseq_slot = need_order_offset + need_order_id;
                assert(need_order_id < measure_deps.size());
                auto &deps = measure_deps[need_order_id];
                assert(var->is_call());
                assert(var->get_callee().is_llvm);
                uint32_t nargs = var->args().size();
                for (uint32_t i = 0; i < nargs; i++) {
                    auto ref = var->get_ref(i);
                    if (!ref)
                        continue;
                    uint32_t refid = ref->varid();
                    assert(refid < varid);
                    auto refslot = m_var_slots[refid];
                    auto ref_bseq_slot = refslot - m_host_seq.nshared;
                    auto type = m_var_types[refid];
                    if (type == VarType::BSeqMeasure) {
                        assert(refslot >= m_host_seq.nshared);
                        deps.insert(ref_bseq_slot);
                    }
                    else if (type == VarType::BSeqNeedOrder) {
                        assert(refslot >= m_host_seq.nshared);
                        assert(ref_bseq_slot >= need_order_offset);
                        assert(ref_bseq_slot - need_order_offset < measure_deps.size());
                        auto &ref_deps = measure_deps[ref_bseq_slot - need_order_offset];
                        deps.insert(ref_deps.begin(), ref_deps.end());
                    }
                }
                assert(need_order_id < host_bseq.deps_count.size());
                host_bseq.deps_count[need_order_id] = deps.size();
                for (auto dep: deps) {
                    assert(dep < host_bseq.reverse_depends.size());
                    host_bseq.reverse_depends[dep].push_back(need_order_id);
                }
                need_order_id++;
            }
            else {
                continue;
            }
            auto slot = m_host_seq.nshared + bseq_slot;
            assert(bseq_slot >= host_bseq.nmeasure);
            assert(bseq_slot - host_bseq.nmeasure < host_bseq.types.size());
            host_bseq.types[bseq_slot - host_bseq.nmeasure] = ir2hostseq(var->type());
            m_var_slots[varid] = slot;
        }

        for (uint32_t timevalid = 0; timevalid < ntimevals; timevalid++) {
            auto owner = m_timeval_owners[timevalid];
            if (owner != seq_id)
                continue;
            auto type = m_timeval_types[timevalid];
            auto &timeval = m_timevals[timevalid];
            uint32_t bseq_slot;
            if (type == VarType::BSeqDirect) {
                bseq_slot = host_bseq.nmeasure + direct_id++;
            }
            else if (type == VarType::BSeqNeedOrder) {
                bseq_slot = need_order_offset + need_order_id;
                assert(need_order_id < measure_deps.size());
                auto &deps = measure_deps[need_order_id];
                for (auto refid: timeval.vars) {
                    assert(refid < m_var_slots.size());
                    auto refslot = m_var_slots[refid];
                    assert(refslot < bseq_slot + m_host_seq.nshared);
                    auto ref_bseq_slot = refslot - m_host_seq.nshared;
                    auto type = m_var_types[refid];
                    if (type == VarType::BSeqMeasure) {
                        assert(refslot >= m_host_seq.nshared);
                        deps.insert(ref_bseq_slot);
                    }
                    else if (type == VarType::BSeqNeedOrder) {
                        assert(refslot >= m_host_seq.nshared);
                        assert(ref_bseq_slot >= need_order_offset);
                        assert(ref_bseq_slot - need_order_offset < measure_deps.size());
                        auto &ref_deps = measure_deps[ref_bseq_slot - need_order_offset];
                        deps.insert(ref_deps.begin(), ref_deps.end());
                    }
                }
                for (auto refid: timeval.timevals) {
                    assert(refid < m_timeval_slots.size());
                    auto refslot = m_timeval_slots[refid];
                    assert(refslot < bseq_slot + m_host_seq.nshared);
                    auto ref_bseq_slot = refslot - m_host_seq.nshared;
                    auto type = m_timeval_types[refid];
                    if (type == VarType::BSeqMeasure) {
                        assert(refslot >= m_host_seq.nshared);
                        deps.insert(ref_bseq_slot);
                    }
                    else if (type == VarType::BSeqNeedOrder) {
                        assert(refslot >= m_host_seq.nshared);
                        assert(ref_bseq_slot >= need_order_offset);
                        assert(ref_bseq_slot - need_order_offset < measure_deps.size());
                        auto &ref_deps = measure_deps[ref_bseq_slot - need_order_offset];
                        deps.insert(ref_deps.begin(), ref_deps.end());
                    }
                }
                assert(need_order_id < host_bseq.deps_count.size());
                host_bseq.deps_count[need_order_id] = deps.size();
                for (auto dep: deps) {
                    assert(dep < host_bseq.reverse_depends.size());
                    host_bseq.reverse_depends[dep].push_back(need_order_id);
                }
                need_order_id++;
            }
            else {
                continue;
            }
            auto slot = m_host_seq.nshared + bseq_slot;
            assert(bseq_slot >= host_bseq.nmeasure);
            assert(bseq_slot - host_bseq.nmeasure < host_bseq.types.size());
            host_bseq.types[bseq_slot - host_bseq.nmeasure] = HostSeq::Type::Int64;
            m_timeval_slots[timevalid] = slot;
        }
    }
}

void Compiler::cache_time_slots()
{
    uint32_t nevent_times = m_event_times.size();
    for (uint32_t timeid = 0; timeid < nevent_times; timeid++) {
        auto et = m_event_times[timeid];
        auto timevalid = m_time_timevals[timeid];
        auto slot = m_timeval_slots[timevalid];
        m_time_slots.emplace(et, slot);
    }
}

void Compiler::populate_bseqs()
{
    uint32_t nbseqs = m_bseqs.size();
    for (uint32_t bseq_idx = 0; bseq_idx < nbseqs; bseq_idx++) {
        auto bseq = m_bseqs[bseq_idx];
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        for (auto &[global_id, assign]: bseq->get_assigns())
            host_bseq.assignments.push_back({global_id, m_var_slots[assign.val->varid()]});
        // Put assumptions that can be evaluated directly first.
        for (auto &assume: bseq->get_assumes()) {
            auto timeid = find_val(m_assumption_ids, &assume);
            auto timevalid = m_time_timevals[timeid];
            assert(timevalid < m_timeval_types.size());
            if (m_timeval_types[timevalid] == VarType::BSeqNeedOrder)
                continue;
            auto slot = m_timeval_slots[timevalid];
            host_bseq.assumptions.push_back({assume.second.sign, slot, assume.second.id});
        }
        host_bseq.ndirect_assumes = host_bseq.assumptions.size();
        host_bseq.assumptions_idx.resize(host_bseq.nneed_order, uint32_t(-1));
        uint32_t need_order_slot_offset = (m_host_seq.nshared + host_bseq.nmeasure +
                                           host_bseq.ndirect);
        for (auto &assume: bseq->get_assumes()) {
            auto timeid = find_val(m_assumption_ids, &assume);
            auto timevalid = m_time_timevals[timeid];
            assert(timevalid < m_timeval_types.size());
            if (m_timeval_types[timevalid] != VarType::BSeqNeedOrder)
                continue;
            auto slot = m_timeval_slots[timevalid];
            assert(slot >= need_order_slot_offset);
            auto need_order_id = slot - need_order_slot_offset;
            host_bseq.assumptions_idx[need_order_id] = host_bseq.assumptions.size();
            host_bseq.assumptions.push_back({assume.second.sign, slot, assume.second.id});
        }
        auto branches = bseq->get_branches();
        auto nbranches = branches.size();
        for (uint32_t i = 0; i < nbranches; i++)
            host_bseq.branches[i].cond = m_var_slots[branches[i].cond->varid()];
        for (auto &endtime: bseq->get_endtimes())
            host_bseq.endtimes.push_back(find_val(m_time_slots, endtime.get()));
        for (auto &[chn, info]: bseq->get_channel_infos()) {
            for (auto &pulse: info.pulses) {
                HostSeq::Pulse host_pulse;
                host_pulse.id = pulse.id();
                host_pulse.time = find_val(m_time_slots, &pulse.start());
                host_pulse.chn = chn;
                host_pulse.ramp_func = nullptr;
                if (pulse.is_measure()) {
                    host_pulse.len = -2;
                    host_pulse.value = -1;
                    host_pulse.measure = m_var_slots[pulse.val()->varid()] - m_host_seq.nshared;
                    host_pulse.endvalue = -1;
                    host_pulse.cond = -1;
                }
                else {
                    auto len = pulse.len();
                    auto value_id = pulse.val()->varid();
                    if (len) {
                        host_pulse.len = m_var_slots[len->varid()];
                        m_var_bseqidx_pulseid[value_id].emplace_back(
                            bseq_idx, host_bseq.pulses.size());
                    }
                    else {
                        host_pulse.len = uint32_t(-1);
                    }
                    auto cond = pulse.cond();
                    if (cond) {
                        host_pulse.cond = m_var_slots[cond->varid()];
                    }
                    else {
                        host_pulse.cond = uint32_t(-1);
                    }
                    host_pulse.value = m_var_slots[value_id];
                    auto it = info.pulse_start_vars.find(&pulse);
                    if (it == info.pulse_start_vars.end()) {
                        host_pulse.measure = -1;
                    }
                    else {
                        host_pulse.measure = m_var_slots[it->second->varid()] - m_host_seq.nshared;
                    }
                    host_pulse.endvalue = m_var_slots[pulse.endval()->varid()];
                }
                host_bseq.pulses.push_back(host_pulse);
            }
        }
    }
}

auto Compiler::find_fptr_slot(uint32_t id, bool is_var) const -> val_eval_t*
{
    auto &types = is_var ? m_var_types : m_timeval_types;
    auto &slots = is_var ? m_var_slots : m_timeval_slots;
    auto &owners = is_var ? m_var_owners : m_timeval_owners;
    assert(types.size() == slots.size());
    assert(types.size() == owners.size());

    auto type = types[id];
    auto slot = slots[id];
    switch (type) {
    case VarType::GlobalValue:
        assert(slot >= global_value_offset);
        assert(slot - global_value_offset < m_host_seq.global_evals.size());
        return &m_host_seq.global_evals[slot - global_value_offset];
    case VarType::BSeqDirect:
    case VarType::BSeqNeedOrder: {
        auto owner = owners[id];
        assert(owner);
        auto bseq_idx = find_val(m_bseq_idxs, owner);
        auto &host_bseq = m_host_seq.seqs[bseq_idx];
        assert(slot >= m_host_seq.nshared + host_bseq.nmeasure);
        assert(slot - m_host_seq.nshared - host_bseq.nmeasure < host_bseq.evals.size());
        return &host_bseq.evals[slot - m_host_seq.nshared - host_bseq.nmeasure];
    }
    default:
        return nullptr;
    }
}

namespace {
static void dummy_eval_func(void*)
{
}
}

void Compiler::generate_fptrs()
{
    // LLVM setup...
    auto &env = m_seq.env();
    auto env_mod = env.llvm_module();
    auto &llvm_ctx = env_mod->getContext();
    // Create a new module that includes all the code that we want to compile
    // This probably includes most of the functions in the old module
    // but we do not want to add new wrapper functions to the old one.
    // Also most of the old functions will be inlined and won't be emitted directly.
    llvm::Module mod("", llvm_ctx);
    // For copying code into the new module
    LLVM::FunctionMover mover(&mod);
    // A new cgctx that emits code in the new module
    LLVM::Codegen::Context cgctx(&mod);

    struct ValFunc {
        val_eval_t *slot;
        llvm::Function *f;
        std::string name{};
    };
    struct RampFunc {
        ramp_eval_t *slot;
        llvm::Function *f;
        std::string name{};
    };
    std::vector<ValFunc> val_funcs;
    std::vector<RampFunc> ramp_funcs;

    // First, copy the functions we need over and generate wrapper for each slot
    uint32_t nvars = m_all_vars.size();
    for (uint32_t varid = 0; varid < nvars; varid++) {
        auto fptr_slot = find_fptr_slot(varid, true);
        bool isramp = m_var_isramps[varid];
        if (!fptr_slot && !isramp)
            continue;
        assert(!fptr_slot || !*fptr_slot);
        if (isramp && fptr_slot)
            *fptr_slot = dummy_eval_func;
        auto var = m_all_vars[varid];
        assert(var->is_call());
        assert(var->get_callee().is_llvm);
        auto f = mover.clone_function(var->get_callee().llvm);
        LLVM::Codegen::Wrapper wrap{true};
        if (!isramp) {
            assert(var->llvm_type() == f->getReturnType());
            wrap.add_closure(-1, m_var_slots[varid]);
        }
        else {
            f = get_ramp_func(var, f, cgctx.T_f64, cgctx.T_f64);
        }
        bool found_arg = false;
        auto args = var->args();
        uint32_t nargs = args.size();
        for (uint32_t argi = 0; argi < nargs; argi++) {
            auto arg = args[argi];
            if (arg.is_arg()) {
                assert(isramp);
                assert(!found_arg);
                (void)found_arg;
                found_arg = true;
                assert(arg.get_arg() == 0);
                continue;
            }
            assert(arg.is_var());
            assert(arg.get_var()->llvm_type() == f->getArg(argi)->getType());
            wrap.add_closure(argi, m_var_slots[arg.get_var()->varid()]);
        }
        auto wrapf = cgctx.emit_wrapper(f, "0", wrap);
        if (isramp) {
            for (auto [bseq_idx, pulse_id]: find_val(m_var_bseqidx_pulseid, varid)) {
                auto &host_bseq = m_host_seq.seqs[bseq_idx];
                ramp_funcs.push_back(RampFunc{&(host_bseq.pulses[pulse_id].ramp_func), wrapf});
            }
        }
        else {
            assert(fptr_slot);
            val_funcs.push_back(ValFunc{fptr_slot, wrapf});
        }
    }

    // Then do the same for the times
    llvm::FunctionType *val_ftype = llvm::FunctionType::get(
        llvm::Type::getVoidTy(llvm_ctx), {LLVM::get_pointer_type(cgctx.T_i8)}, false);
    auto T_i64 = llvm::Type::getInt64Ty(llvm_ctx);
    // lrint can be lowered directly into SSE instructions
    auto intrin = LLVM::get_intrinsic(&mod, llvm::Intrinsic::lrint, {T_i64, cgctx.T_f64});
    uint32_t ntimevals = m_timevals.size();
    for (uint32_t timevalid = 0; timevalid < ntimevals; timevalid++) {
        auto fptr_slot = find_fptr_slot(timevalid, false);
        if (!fptr_slot)
            continue;
        assert(!*fptr_slot);
        auto &timeval = m_timevals[timevalid];
        auto slot = m_timeval_slots[timevalid];
        llvm::Function *f =
            llvm::Function::Create(val_ftype, llvm::GlobalValue::ExternalLinkage, "+", mod);
        f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        f->addFnAttr(llvm::Attribute::Speculatable);
        f->addFnAttr(llvm::Attribute::NoRecurse);
        f->addFnAttr(llvm::Attribute::NoUnwind);
        auto b0 = llvm::BasicBlock::Create(llvm_ctx, "top", f);
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

        llvm::Value *sumv = nullptr;
        if (timeval.offset != 0)
            sumv = llvm::ConstantInt::get(T_i64, timeval.offset);
        auto clarg = &*f->arg_begin();
        for (auto refid: timeval.vars) {
            auto refslot = m_var_slots[refid];
            auto ptr = builder.CreateConstGEP1_32(cgctx.T_i8, clarg, refslot * 8);
            llvm::LoadInst *load;
            llvm::Value *val;
            auto ref_type = m_all_vars[refid]->type();
            auto ref_ltype = cgctx.llvm_ty(ref_type);
            ptr = builder.CreateBitCast(ptr, LLVM::get_pointer_type(ref_ltype));
            load = builder.CreateLoad(ref_ltype, ptr);
            switch (ref_type) {
            case IR::Type::Bool:
                val = builder.CreateZExt(load, T_i64);
                break;
            case IR::Type::Int32:
                val = builder.CreateSExt(load, T_i64);
                break;
            case IR::Type::Float64:
                val = builder.CreateCall(intrin, {load});
                break;
            default:
                val = nullptr;
            }
            assert(val);
#if LLVM_VERSION_MAJOR >= 11
            load->setAlignment(llvm::Align(alignof(double)));
#elif LLVM_VERSION_MAJOR >= 10
            load->setAlignment(llvm::MaybeAlign(alignof(double)));
#else
            load->setAlignment(alignof(double));
#endif
            load->setMetadata(llvm::LLVMContext::MD_tbaa, cgctx.tbaa_const);
            if (!sumv) {
                sumv = val;
                continue;
            }
            sumv = builder.CreateAdd(sumv, val);
        }
        for (auto refid: timeval.timevals) {
            auto refslot = m_timeval_slots[refid];
            auto ptr = builder.CreateConstGEP1_32(cgctx.T_i8, clarg, refslot * 8);
            llvm::LoadInst *load;
            llvm::Value *val;
            ptr = builder.CreateBitCast(ptr, LLVM::get_pointer_type(T_i64));
            load = builder.CreateLoad(T_i64, ptr);
            val = load;
#if LLVM_VERSION_MAJOR >= 11
            load->setAlignment(llvm::Align(alignof(double)));
#elif LLVM_VERSION_MAJOR >= 10
            load->setAlignment(llvm::MaybeAlign(alignof(double)));
#else
            load->setAlignment(alignof(double));
#endif
            load->setMetadata(llvm::LLVMContext::MD_tbaa, cgctx.tbaa_const);
            if (!sumv) {
                sumv = val;
                continue;
            }
            sumv = builder.CreateAdd(sumv, val);
        }
        auto resptr = builder.CreateConstGEP1_32(cgctx.T_i8, clarg, slot * 8);
        resptr = builder.CreateBitCast(resptr, LLVM::get_pointer_type(T_i64));
        auto store = builder.CreateStore(sumv, resptr);
#if LLVM_VERSION_MAJOR >= 11
        store->setAlignment(llvm::Align(alignof(double)));
#elif LLVM_VERSION_MAJOR >= 10
        store->setAlignment(llvm::MaybeAlign(alignof(double)));
#else
        store->setAlignment(alignof(double));
#endif
        builder.CreateRetVoid();
        val_funcs.push_back(ValFunc{fptr_slot, f});
    }

    // Clean up function exports
    for (auto &go: mod.global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    for (auto &val_func: val_funcs) {
        val_func.f->setLinkage(llvm::GlobalValue::ExternalLinkage);
        val_func.f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        val_func.f->setUnnamedAddr(llvm::GlobalValue::UnnamedAddr::Global);
    }
    for (auto &ramp_func: ramp_funcs) {
        ramp_func.f->setLinkage(llvm::GlobalValue::ExternalLinkage);
        ramp_func.f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        ramp_func.f->setUnnamedAddr(llvm::GlobalValue::UnnamedAddr::Global);
    }

    LLVM::runGlobalRenamePasses(mod);

    // The optimization passes below (in `emit_objfile`) may recreate the functions
    // so the function handle may not be valid anymore.
    // We are done with the function renaming so we can simply get the names now.
    for (auto &val_func: val_funcs) {
        val_func.name = val_func.f->getName();
        assert(!val_func.name.empty());
        val_func.f = nullptr;
    }
    for (auto &ramp_func: ramp_funcs) {
        ramp_func.name = ramp_func.f->getName();
        assert(!ramp_func.name.empty());
        ramp_func.f = nullptr;
    }

    llvm::SmallVector<char,0> vec;
    auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(), &mod, true);
    if (!res)
        throw std::runtime_error("Sequence function compilation failed.");
    m_obj_id = m_engine.load(&vec[0], vec.size(), m_resolver);
    if (!m_obj_id)
        throw std::runtime_error("Sequence function compilation failed.");
    for (auto &val_func: val_funcs) {
        auto sym = m_engine.get_symbol(val_func.name);
        assert(sym);
        *val_func.slot = val_eval_t(sym);
    }
    for (auto &ramp_func: ramp_funcs) {
        auto sym = m_engine.get_symbol(ramp_func.name);
        assert(sym);
        *ramp_func.slot = ramp_eval_t(sym);
    }
    m_engine.reset_dyld();
}

NACS_EXPORT() uint64_t Compiler::compile()
{
    m_host_seq.npublic_globals = m_seq.npublic_globals();
    // Classify and count variables
    // Assign constant variable ID
    prescan_vars();
    // Create slot for constant time, create time variables for event time and assumptions
    prescan_times();
    // All the constants have been created
    m_host_seq.nconsts = m_host_seq.values.size();
    classify_times();
    m_host_seq.nglobal_vals = nglobal_vals_var + nglobal_vals_time;
    m_host_seq.global_evals.resize(m_host_seq.nglobal_vals, nullptr);
    m_host_seq.nshared = (m_host_seq.nconsts + m_host_seq.nglobals +
                          m_host_seq.nglobal_vals + m_host_seq.nchannels);
    // Assign slots to global variables and values
    assign_shared_vars();
    // Assign slots to global time values
    assign_shared_times();
    // Fill in default values
    assign_defaults();
    // Create `BasicSeq`s and assign branch targets
    create_bseqs();
    // Assign slots to `BasicSeq`s variables and values (including times)
    assign_bseq_vars();
    // Populate `m_time_slots` so that we can directly find the slots
    cache_time_slots();
    // Now with most of the dirty work done,
    // create/populate the actual basic sequences...
    populate_bseqs();
    // Finally fill in all the function pointers
    generate_fptrs();
    return m_obj_id;
}

}
