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

#include "../nacs-utils/llvm/codegen.h"

#include <algorithm>
#include <utility>

#include <math.h>

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

NACS_EXPORT() void Seq::set_defval(uint32_t chn, double val)
{
    m_defval[chn] = val;
}

NACS_EXPORT() double Seq::defval(uint32_t chn)
{
    auto [it, inserted] = m_defval.emplace(chn, 0);
    return it->second;
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
    auto nan = std::numeric_limits<double>::quiet_NaN();
    bool changed = false;
    std::map<BasicSeq*,double> startval;
    // Return whether `seq` should be scanned again.
    auto add_startval = [&] (BasicSeq *seq, double val) {
        if (!seq)
            return false;
        if (!std::isfinite(val)) {
            auto &v = startval[seq];
            if (!std::isfinite(v))
                return false;
            // We previously believe this sequence has a unique start value.
            // Scan again to make sure we clear that assumption.
            v = nan;
            return true;
        }
        auto [it, inserted] = startval.emplace(seq, val);
        if (inserted)
            return true;
        // Same value, good
        if (it->second == val)
            return false;
        if (!std::isfinite(it->second))
            return false;
        it->second = nan;
        return true;
    };
    for (auto &seq: m_seqs) {
        changed |= seq.optimize_order(chn);
        seq.m_used = false;
    }
    std::function<void(BasicSeq*, double)> scan_seq = [&] (BasicSeq *seq, double val) {
        bool should_rescan = add_startval(seq, val);
        if (!should_rescan && (!seq || seq->m_used))
            return;
        seq->m_used = true;
        auto &info = seq->m_channels.find(chn)->second;
        double endval = nan;
        auto endval_var = info.endval.get();
        if (endval_var == info.startval.get()) { // Endval was a simple forward of startval.
            endval = val;
        }
        else if (endval_var && endval_var->is_const()) {
            // Only allow constants to propagate.
            endval = endval_var->get_const().get<double>();
        }
        for (auto &br: seq->m_branches)
            scan_seq(br.target, endval);
        scan_seq(seq->m_default_branch, endval);
    };
    scan_seq(&m_seqs.front(), defval(chn));
    for (auto [bseq, val]: startval) {
        if (!std::isfinite(val))
            continue;
        auto &sv = bseq->m_channels[chn].startval;
        if (sv->is_const()) // Already has a value
            continue;
        sv->assign_const(IR::TagVal(val));
        changed = true;
    }
    return changed;
}

NACS_INTERNAL bool Seq::optimize_cfg()
{
    bool changed = false;
    for (auto &seq: m_seqs) {
        changed |= seq.optimize_branch();
        seq.reset_used();
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
        for (unsigned i = 1; i <= nchns; i++)
            seq.preoptimize_pulse(i, env());
        seq.preoptimize_eventtimes();
    }
    optimize_cfg();
    while (true) {
        m_env.optimize();
        bool changed = false;
        for (auto &seq: m_seqs)
            changed |= seq.optimize_vars(*this);
        changed |= optimize_cfg();
        for (auto &seq: m_seqs)
            changed |= seq.optimize_endtimes();
        for (unsigned i = 1; i <= nchns; i++)
            changed |= optimize_chn(i);
        for (auto &seq: m_seqs)
            changed |= seq.postoptimize_eventtimes();
        if (!changed) {
            break;
        }
    }
    bool changed = false;
    for (auto &seq: m_seqs)
        changed |= seq.optimize_final(*this);
    if (changed) {
        m_env.optimize();
    }
    else {
        m_env.gc();
    }
}

NACS_EXPORT() void Seq::prepare()
{
    assert(!m_seqs.empty());
    m_npublic_globals = m_slot_vars.size();
    // Do an optimization to normalize the values
    // so that we don't need to deal with call of var's anywhere.
    m_env.optimize();
    for (auto &seq: m_seqs) {
        seq.prepare(*this);
    }
}

NACS_EXPORT() void Seq::print(std::ostream &stm) const
{
    unsigned nchns = m_chnnames.size();
    stm << "Channels:" << std::endl;
    for (unsigned i = 1; i <= nchns; i++)
        stm << "  " << i << ": " << m_chnnames[i - 1] << std::endl;
    for (auto &seq: m_seqs)
        seq.print(stm);
    stm << "Default Value:" << std::endl;
    for (auto &dv: m_defval)
        stm << "  " << dv.first << " = Float64 " << dv.second << std::endl;
    if (!m_slot_vars.empty()) {
        stm << "Global Variable:" << std::endl;
        for (auto &slot: m_slot_vars) {
            if (!slot)
                continue;
            stm << "  ";
            slot->print(stm, true);
        }
    }
    m_env.print(stm);
}

}
