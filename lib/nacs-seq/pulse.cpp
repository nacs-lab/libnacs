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

#include "pulse.h"

#include "error.h"

#include "../nacs-utils/llvm/codegen.h"

#include <algorithm>

namespace NaCs::Seq {

NACS_EXPORT_ Pulse::Pulse(uint32_t id, EventTime &start,
                          Var *len, Var *val, Var *cond, bool is_measure)
    : m_id(id),
      m_is_measure(is_measure),
      m_start(start.ref()),
      m_len(len ? len->ref() : nullptr),
      m_val(val->ref()),
      m_cond(cond ? cond->ref() : nullptr)
{
    assert(val);
    // Caller should check this
    assert(val->nfreeargs() <= 2);
    // Measurement value must be extern with zero length and no condition
    assert(!is_measure || (val->is_extern() && !len && !cond));
    // zero-length cannot be ramp
    assert(len || is_measure || val->nfreeargs() == 0 || val->argument_unused(0));
}

NACS_EXPORT() bool Pulse::needs_oldval() const
{
    return (m_val->nfreeargs() >= 2 && !m_val->argument_unused(1));
}

NACS_EXPORT() bool Pulse::clear_unused_args()
{
    assert(!m_is_measure);
    bool changed = false;
    if (m_cond) {
        if (auto v = m_cond->get_assigned_var())
            m_cond = v->ref();
        if (m_cond->is_const()) {
            if (m_cond->get_const().get<bool>()) {
                m_cond.reset(nullptr);
            }
            else {
                // Clear the value, wait for basic sequence to remove us.
                m_val->assign_const(IR::TagVal(false));
                if (m_endval)
                    m_endval->assign_const(IR::TagVal(false));
                return true;
            }
            changed = true;
        }
    }
    if (m_val->nfreeargs() == 0) {
        if (m_len) {
            m_len.reset(nullptr);
            return true;
        }
        return changed;
    }
    auto &env = m_val->env();
    if (m_val->nfreeargs() == 1) {
        if (!m_val->argument_unused(0))
            return changed;
        Arg args[] = {Arg::create_const(false)};
        m_val.reset(env.new_call(m_val.get(), args, 0));
        m_len.reset(nullptr);
        return true;
    }
    assert(m_val->nfreeargs() == 2);
    bool used0 = !m_val->argument_unused(0);
    if (!m_val->argument_unused(1)) {
        if (!used0 && m_len) {
            m_len.reset(nullptr);
            return true;
        }
        return changed;
    }
    Arg args[2];
    args[0] = used0 ? Arg::create_arg(0) : Arg::create_const(false);
    args[1] = Arg::create_const(false);
    m_val.reset(env.new_call(m_val.get(), args, used0 ? 1 : 0));
    if (!used0)
        m_len.reset(nullptr);
    return true;
}

NACS_EXPORT() bool Pulse::set_oldval(Var *oldval)
{
    assert(!m_is_measure);
    assert(needs_oldval());
    assert(oldval);
    auto &env = m_val->env();
    assert(m_val->nfreeargs() == 2);
    bool used0 = !m_val->argument_unused(0);
    Arg args[2];
    args[0] = used0 ? Arg::create_arg(0) : Arg::create_const(false);
    args[1] = Arg::create_var(oldval);
    m_val = env.new_call(m_val.get(), args, used0 ? 1 : 0)->ref();
    if (!used0)
        m_len.reset(nullptr);
    return true;
}

NACS_EXPORT() void Pulse::print(std::ostream &stm, Var *startval, bool newline) const
{
    if (is_measure()) {
        assert(!startval);
        stm << "M(" << id() << "@";
        start().print(stm);
        stm << ")=";
        val()->print(stm, false, true);
    }
    else {
        stm << "O(" << id() << "@";
        start().print(stm);
        if (m_len) {
            stm << "; +[";
            m_len->print(stm, false, true);
            stm << "]";
        }
        stm << ")=";
        val()->print(stm, false, true);
        if (m_cond) {
            stm << " if ";
            m_cond->print(stm, false, true);
        }
        if (startval) {
            stm << " # I=";
            startval->print(stm, false, true);
        }
        if (m_endval) {
            stm << (startval ? ", E=" : " # E=");
            m_endval->print(stm, false, true);
        }
    }
    if (newline) {
        stm << std::endl;
    }
}

NACS_EXPORT() Var *Pulse::compute_endval() const
{
    assert(!m_is_measure);
    if (m_val->nfreeargs() < 1)
        return m_val.get();
    auto &env = m_val->env();
    auto arg0 = Arg::create_var(m_len.get());
    if (m_val->nfreeargs() == 1) {
        Arg args = {arg0};
        return env.new_call(m_val.get(), args, 0);
    }
    assert(m_val->nfreeargs() == 2);
    Arg args[2];
    args[0] = arg0;
    if (!m_val->argument_unused(1)) {
        return nullptr;
    }
    else {
        args[1] = Arg::create_const(false);
    }
    return env.new_call(m_val.get(), args, 0);
}

NACS_EXPORT() bool Pulse::optimize()
{
    bool changed = false;
    assert(!m_is_measure);
    // Simple assignment can happen since the variable optimization maintains pointer
    // identity for external variables. Update our external references to allow
    // such assgignments to be optimized out.
    if (m_len) {
        if (auto v = m_len->get_assigned_var()) {
            // Here we'll rely on the parent to check if the length is greater than 0.
            // Since pulse length can be optimized out,
            // this isn't a good/reliable place to do the check anyway.
            m_len = v->ref();
        }
    }
    if (auto v = m_val->get_assigned_var())
        m_val = v->ref();
    if (m_endval) {
        if (auto v = m_endval->get_assigned_var()) {
            m_endval = v->ref();
        }
    }
    else {
        // eagerly compute endval so that it can be optimized
        changed |= endval() != nullptr;
    }
    return changed;
}

NACS_EXPORT() void Pulse::check_start() const
{
    if (m_start->terms.empty() && m_start->tconst < 0) {
        throw Error(Error::Type::Pulse, Error::Pulse::NegTime,
                    Error::Type::Pulse, m_id, "Pulse time must be positive.");
    }
}

}
