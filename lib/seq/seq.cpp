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

#include "seq.h"

#include "../utils/llvm/codegen.h"

#include <algorithm>
#include <utility>

namespace NaCs::Seq {

NACS_EXPORT_ Seq::Seq(std::unique_ptr<LLVM::Codegen::Context> cgctx)
    : m_env(std::move(cgctx))
{
}

NACS_EXPORT_ Seq::Seq(llvm::LLVMContext &llvm_ctx)
    : m_env(llvm_ctx)
{
}

NACS_EXPORT() Seq::~Seq()
{
}

NACS_EXPORT() BasicSeq *Seq::add_basicseq(uint32_t id)
{
    assert(id);
    return &m_seqs.emplace_back(id);
}

NACS_EXPORT() uint32_t Seq::get_chn_id(llvm::StringRef name, bool create)
{
    if (!create) {
        auto it = m_chnids.find(name);
        if (it == m_chnids.end())
            return 0;
        return it->second;
    }
    auto res = m_chnids.try_emplace(name);
    if (res.second) {
        m_chnnames.push_back(name.str());
        res.first->second = m_chnnames.size();
    }
    return res.first->second;
}

NACS_EXPORT() const std::string &Seq::get_chn_name(uint32_t id) const
{
    return m_chnnames[id - 1];
}

NACS_EXPORT() void Seq::set_defval(uint32_t chn, Var *val)
{
    assert(val->is_const());
    m_defval[chn].reset(val);
}

NACS_EXPORT() Var *Seq::defval(uint32_t chn)
{
    auto &val = m_defval[chn];
    if (!val)
        val.reset(m_env.new_const(false));
    return val.get();
}

NACS_EXPORT() Var *Seq::get_const(IR::TagVal c)
{
    if (c.typ != IR::Type::Float64)
        return m_env.new_const(c);
    auto &cached = m_float_cache[c.val.f64];
    if (!cached)
        cached = m_env.new_const(c);
    return cached;
}

NACS_EXPORT() Var *Seq::get_call(const IR::Function &func, llvm::ArrayRef<Arg> args,
                                 int nfreeargs)
{
    return m_env.new_call(func, args, nfreeargs);
}

NACS_EXPORT() Var *Seq::get_slot(IR::Type typ, uint32_t id)
{
    if (id >= m_slot_vars.size())
        m_slot_vars.resize(id + 1);
    auto &var = m_slot_vars[id];
    if (!var)
        var = m_env.new_extern({typ, id})->ref();
    return var.get();
}

NACS_INTERNAL bool Seq::optimize_chn(uint32_t chn)
{
    bool changed = false;
    std::map<BasicSeq*,Var*> startval;
    auto same_val = [] (Var *v1, Var *v2) {
        if (v1 == v2)
            return true;
        if (!v1 || !v2 || !v1->is_const() || !v2->is_const())
            return false;
        // The conversion to double is guaranteed to be lossless.
        return v1->get_const().get<double>() == v2->get_const().get<double>();
    };
    // Return whether `seq` should be scanned again.
    auto add_startval = [&] (BasicSeq *seq, Var *val) {
        if (!seq)
            return false;
        if (!val) {
            auto &v = startval[seq];
            if (!v)
                return false;
            // We previously believe this sequence has a unique start value.
            // Scan again to make sure we clear that assumption.
            v = nullptr;
            return true;
        }
        auto it = startval.find(seq);
        if (it == startval.end()) {
            startval[seq] = val;
            return true;
        }
        // Same value, good
        if (same_val(it->second, val))
            return false;
        if (!it->second)
            return false;
        it->second = nullptr;
        return true;
    };
    for (auto &seq: m_seqs) {
        changed |= seq.optimize_order(chn);
        seq.m_used = false;
    }
    std::function<void(BasicSeq*, Var*)> scan_seq = [&] (BasicSeq *seq, Var *val) {
        bool should_rescan = add_startval(seq, val);
        if (!should_rescan && (!seq || seq->m_used))
            return;
        seq->m_used = true;
        auto endval = seq->endval(chn);
        if (!endval && val && !seq->has_output(chn))
            endval = val;
        // Only allow constants to propagate.
        if (endval && !endval->is_const())
            endval = nullptr;
        for (auto &br: seq->m_branches)
            scan_seq(br.target, endval);
        scan_seq(seq->m_default_branch, endval);
    };
    scan_seq(&m_seqs.front(), defval(chn));
    for (auto &v: startval) {
        if (!v.second)
            continue;
        auto &sv = v.first->m_startval[chn];
        if (sv) // Already has a value
            continue;
        assert(v.second->is_const());
        sv.reset(v.second);
        changed = true;
    }
    return changed;
}

NACS_INTERNAL bool Seq::optimize_cfg()
{
    bool changed = false;
    for (auto &seq: m_seqs) {
        changed |= seq.optimize_branch();
        seq.m_used = false;
    }
    m_seqs.front().mark_recursive();
    m_seqs.remove_if([&changed] (auto &seq) {
        if (seq.m_used)
            return false;
        changed = true;
        return true;
    });
    return changed;
}

NACS_EXPORT() void Seq::optimize()
{
    m_float_cache.clear();
    // Do some clean up before doing the actual optimization of the variables
    unsigned nchns = m_chnnames.size();
    for (auto &seq: m_seqs) {
        for (unsigned i = 1; i <= nchns; i++) {
            seq.optimize_pulse(i);
        }
    }
    optimize_cfg();
    while (true) {
        m_env.optimize();
        for (auto &seq: m_seqs)
            seq.optimize_vars();
        bool changed = optimize_cfg();
        for (unsigned i = 1; i <= nchns; i++)
            changed |= optimize_chn(i);
        if (!changed) {
            break;
        }
    }
    m_env.gc();
}

}
