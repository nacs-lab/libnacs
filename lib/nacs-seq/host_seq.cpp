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

#include "host_seq.h"

#include "error.h"
#include "../nacs-utils/number.h"

#include <algorithm>

#include <assert.h>

namespace NaCs::Seq {

inline void HostSeq::BasicSeq::init() const
{
    m_measure_filled.resize(nmeasure);
    m_deps_count_copy.resize(deps_count.size());
}

inline void HostSeq::BasicSeq::init_run() const
{
    std::fill(m_measure_filled.begin(), m_measure_filled.end(), false);
    std::copy(deps_count.begin(), deps_count.end(), m_deps_count_copy.begin());
}

NACS_EXPORT() void HostSeq::init()
{
    if (m_initialized)
        throw std::runtime_error("Sequence already initialized.");
    for (auto &bseq: seqs)
        bseq.init();
    m_cur_age = 0;
    m_ages.resize(nglobals + nglobal_vals, 0);
    for (uint32_t i = 0; i < nglobals; i++)
        m_ages[i] = 1;
    m_global_dirty = true;
    m_chn_vals.resize(nchannels);
    m_chn_pulses.resize(nchannels);
    start_values.resize(nchannels);
    memset(values.data() + nconsts, 0, nglobals * sizeof(Value));
    m_initialized = true;
}

NACS_EXPORT() void HostSeq::init_run()
{
    if (!m_initialized)
        throw std::runtime_error("Sequence must be initialized first.");
    if (m_cur_seq_idx != uint32_t(-1) || m_first_bseq)
        throw std::runtime_error("Unfinished sequence cannot be restarted again.");
    m_cur_seq_idx = 0;
    first_bseq = true;
    m_first_bseq = true;
    m_running = false;

    static_assert(sizeof(Value) == sizeof(double));
    memcpy(m_chn_vals.data(), default_values.data(), nchannels * sizeof(Value));
}

NACS_EXPORT() void HostSeq::pre_run()
{
    if (!m_initialized)
        throw std::runtime_error("Sequence must be initialized first.");
    if (m_cur_seq_idx == uint32_t(-1))
        throw std::runtime_error("Sequence must be started before running.");
    if (m_running)
        throw std::runtime_error("Sequence already running.");
    m_running = true;
    assert(m_cur_seq_idx < seqs.size());
    auto &bseq = seqs[m_cur_seq_idx];
    update_values(bseq);
    first_bseq = m_first_bseq;
    m_first_bseq = false;
}

NACS_EXPORT() uint32_t HostSeq::post_run()
{
    if (!m_initialized)
        throw std::runtime_error("Sequence must be initialized first.");
    if (m_cur_seq_idx == uint32_t(-1))
        throw std::runtime_error("Sequence must be started before running.");
    if (!m_running)
        throw std::runtime_error("Sequence not running.");
    m_running = false;
    assert(m_cur_seq_idx < seqs.size());
    auto br = eval_branch(seqs[m_cur_seq_idx]);
    m_cur_seq_idx = br;
    if (br == uint32_t(-1))
         return 0;
    assert(br < seqs.size());
    return seqs[m_cur_seq_idx].id;
}

NACS_EXPORT() double HostSeq::get_global(uint32_t i) const
{
    if (!m_initialized)
        throw std::runtime_error("Sequence must be initialized first.");
    if (i >= npublic_globals)
        return 0;
    i += nconsts;
    return convert_value(types[i], values[i]);
}

void HostSeq::_set_global(uint32_t i, double val)
{
    if (i >= nglobals)
        return;
    auto &age = m_ages[i];
    i += nconsts;
    if (assign_value(types[i], values[i], val)) {
        m_global_dirty = true;
        age = m_cur_age + 1;
    }
}

NACS_EXPORT() void HostSeq::set_global(uint32_t i, double val)
{
    if (!m_initialized)
        throw std::runtime_error("Sequence must be initialized first.");
    if (i >= npublic_globals)
        return;
    _set_global(i, val);
}

inline double HostSeq::convert_value(Type type, Value value)
{
    switch (type) {
    case Type::Bool:
        return value.b;
    case Type::Int32:
        return value.i32;
    case Type::Float64:
        return value.f64;
    case Type::Int64:
        return (double)value.i64;
    default:
        return 0;
    }
}

inline bool HostSeq::assign_value(Type type, Value &slot, double val)
{
    switch (type) {
    case Type::Bool: {
        auto v = val != 0;
        if (slot.b != v) {
            slot.b = v;
            return true;
        }
        return false;
    }
    case Type::Int32: {
        auto v = int32_t(val);
        if (slot.i32 != v) {
            slot.i32 = v;
            return true;
        }
        return false;
    }
    case Type::Float64:
        if (slot.f64 != val) {
            slot.f64 = val;
            return true;
        }
        return false;
    case Type::Int64: {
        auto v = int64_t(val);
        if (slot.i64 != v) {
            slot.i64 = v;
            return true;
        }
        return false;
    }
    default:
        return false;
    }
}

bool HostSeq::value_computed(const BasicSeq &bseq, uint32_t i) const
{
    if (i < nshared)
        return true;
    i -= nshared;
    if (i < bseq.nmeasure)
        return bseq.m_measure_filled[i];
    i -= bseq.nmeasure;
    if (i < bseq.ndirect)
        return true;
    i -= bseq.ndirect;
    return bseq.m_deps_count_copy[i] == 0;
}

NACS_EXPORT() HostSeq::Type HostSeq::get_type(uint32_t i, uint32_t seq_idx) const
{
    auto &bseq = seqs[seq_idx];
    if (i < nshared) {
        if (i < types.size())
            return types[i];
        return Type::Float64;
    }
    i -= nshared;
    if (i < bseq.nmeasure)
        return Type::Float64;
    i -= bseq.nmeasure;
    return bseq.types[i];
}

NACS_EXPORT() double HostSeq::get_value(const BasicSeq &bseq, uint32_t i) const
{
    assert(value_computed(bseq, i));
    uint32_t ifull = i;
    if (i < nshared) {
        if (i < types.size())
            return convert_value(types[i], values[i]);
        return values[i].f64;
    }
    i -= nshared;
    if (i < bseq.nmeasure)
        return values[ifull].f64;
    i -= bseq.nmeasure;
    return convert_value(bseq.types[i], values[ifull]);
}

void HostSeq::update_global_value(uint32_t i)
{
    assert(i < nglobal_vals);
    auto age = m_ages[nglobals + i];
    for (auto dep_idx: depends[i]) {
        assert(dep_idx < nglobal_vals);
        if (m_ages[dep_idx] > age) {
            global_evals[i](values.data());
            m_ages[nglobals + i] = m_cur_age;
            return;
        }
    }
}

bool HostSeq::propagate_measure(const BasicSeq &bseq, uint32_t i)
{
    bool has_time = false;
    for (auto idx: bseq.reverse_depends[i]) {
        assert(idx < bseq.nneed_order);
        assert(bseq.m_deps_count_copy[idx] > 0);
        if (--bseq.m_deps_count_copy[idx])
            continue;
        has_time |= bseq.types[bseq.ndirect + idx] == Type::Int64;
        bseq.evals[bseq.ndirect + idx](values.data());
        auto assume_idx = bseq.assumptions_idx[idx];
        if (assume_idx != uint32_t(-1)) {
            check_assumption(bseq, bseq.assumptions[bseq.ndirect_assumes + assume_idx]);
        }
    }
    bseq.m_measure_filled[i] = true;
    return has_time;
}

void HostSeq::check_assumption(const BasicSeq &bseq, const Assumption &as) const
{
    assert(value_computed(bseq, as.value));
#ifndef NDEBUG
    {
        uint32_t i = as.value;
        if (i < nshared) {
            assert(i < types.size());
            assert(types[i] == Type::Int64);
        }
        else {
            i -= nshared;
            assert(i >= bseq.nmeasure);
            i -= bseq.nmeasure;
            assert(bseq.types[i] == Type::Int64);
        }
    }
#endif
    auto t = get_time(as.value);
    if (as.sign == Sign::Pos && t <= 0) {
        throw Error(Error::Type::EventTime, Error::EventTime::NonPosTime,
                    Error::Type::EventTime, as.id, "Positive time expected.");
    }
    else if (as.sign == Sign::NonNeg && t < 0) {
        throw Error(Error::Type::EventTime, Error::EventTime::NegTime,
                    Error::Type::EventTime, as.id, "Non-negative time expected.");
    }
}

void HostSeq::update_values(const BasicSeq &bseq)
{
    // Globals
    if (m_global_dirty) {
        // Flush global age change.
        m_global_dirty = false;
        m_cur_age += 1;
    }
    static_assert(sizeof(Value) == sizeof(double));
    memcpy(values.data() + nconsts + nglobals + nglobal_vals, m_chn_vals.data(),
           nchannels * sizeof(Value));
    memcpy(start_values.data(), m_chn_vals.data(), nchannels * sizeof(Value));
    for (auto gidx: bseq.global_refs)
        update_global_value(gidx);
    // Initialize tracking info
    bseq.init_run();
    // Initialize channel values tracking info
    std::fill(m_chn_pulses.begin(), m_chn_pulses.end(), nullptr);
    auto get_chn_val = [&] (uint32_t chn, int64_t time) {
        auto p = m_chn_pulses[chn - 1];
        if (!p || time == int64_t(-1))
            return m_chn_vals[chn - 1];
        assert(!p->is_measure());
        assert(p->len != uint32_t(-1));
        auto plen = get_value(bseq, p->len);
        auto ptime = get_time(p->time);
        auto delta = time - ptime;
        if ((double)delta < plen)
            return p->ramp_func((double)delta, values.data());
        return m_chn_vals[chn - 1];
    };

    // Values that doesn't depend on time order of things
    auto data = values.data();
    for (uint32_t i = 0; i < bseq.ndirect; i++)
        bseq.evals[i](data);

    for (uint32_t i = 0; i < bseq.ndirect_assumes; i++)
        check_assumption(bseq, bseq.assumptions[i]);

    auto scan_start = bseq.pulses.begin();
    auto pulses_end = bseq.pulses.end();

    int64_t last_time = 0;
    int64_t last_id = -1;

    while (scan_start != pulses_end) {
        // Sort pulses with known time
        std::sort(scan_start, pulses_end, [&] (const Pulse &p1, const Pulse &p2) {
            if (!value_computed(bseq, p1.time))
                return false;
            if (!value_computed(bseq, p2.time))
                return true;
            auto t1 = get_time(p1.time);
            auto t2 = get_time(p2.time);
            if (t1 < t2)
                return true;
            if (t1 > t2)
                return false;
            return p1.id < p2.id;
        });
        bool changed = false;
        while (scan_start != pulses_end) {
            if (!value_computed(bseq, scan_start->time))
                break;
            auto pulse_time = get_time(scan_start->time);
            if (pulse_time < last_time ||
                (pulse_time == last_time && scan_start->id < last_id)) {
                if (pulse_time < 0)
                    throw Error(Error::Type::Pulse, Error::Pulse::NegTime,
                                Error::Type::Pulse, scan_start->id,
                                "Pulse time must be positive.");
                // This is an internal and not a user error since
                // it should be guaranteed by the checks.
                // It could be an assertion but it's pretty cheap
                // to throw an error here since we need the `t < 0` check anyway...
                throw std::runtime_error("Going back in time not allowed.");
            }
            last_time = pulse_time;
            last_id = scan_start->id;
            auto measure_id = scan_start->measure;
            bool time_updated = false;
            if (measure_id != uint32_t(-1)) {
                assert(measure_id < bseq.nmeasure);
                auto val = get_chn_val(scan_start->chn,
                                       scan_start->is_measure() ? pulse_time : -1);
                values[nshared + measure_id].f64 = val;
                time_updated = propagate_measure(bseq, measure_id);
                if (scan_start->is_measure()) {
                    changed = true;
                    ++scan_start;
                    if (time_updated)
                        break;
                    continue;
                }
            }
            assert(!scan_start->is_measure());
            if (scan_start->cond != uint32_t(-1)) {
                if (!value_computed(bseq, scan_start->cond))
                    break;
                if (get_value(bseq, scan_start->cond) == 0) {
                    // Skip disabled pulse
                    changed = true;
                    ++scan_start;
                    continue;
                }
            }
            if (scan_start->len != uint32_t(-1)) {
                if (!value_computed(bseq, scan_start->len))
                    break;
                if (get_value(bseq, scan_start->len) <= 0) {
                    throw Error(Error::Type::Pulse, Error::Pulse::NonPosLen,
                                Error::Type::Pulse, scan_start->id,
                                "Pulse length must be positive");
                }
            }
            assert(scan_start->value != uint32_t(-1));
            if (!value_computed(bseq, scan_start->value))
                break;
            changed = true;
            auto chn = scan_start->chn;
            double endval;
            if (scan_start->len == uint32_t(-1)) {
                m_chn_pulses[chn - 1] = nullptr;
                endval = get_value(bseq, scan_start->value);
            }
            else {
                m_chn_pulses[chn - 1] = &*scan_start;
                endval = scan_start->ramp_func(get_value(bseq, scan_start->len),
                                               values.data());
            }
            m_chn_vals[chn - 1] = endval;
            assert(value_computed(bseq, scan_start->endvalue));
            ++scan_start;
            if (time_updated) {
                // Even though a pulse's initial value isn't user accessible in the sequence
                // and therefore won't be used to determine a pulse time directly,
                // it can still affect time due to measurements that are
                // handled during optimization.
                break;
            }
        }
        // Shouldn't happen
        if (!changed) {
            throw std::runtime_error("Cannot compute sequence.");
        }
    }

#ifndef NDEBUG
    for (uint32_t i = 0; i < bseq.nmeasure + bseq.ndirect + bseq.nneed_order; i++)
        assert(value_computed(bseq, i + nshared));
#endif

    int64_t length = 0;
    for (auto id: bseq.endtimes)
        length = max(length, get_time(id));
    bseq.length = length;

    for (auto &assign: bseq.assignments) {
        assert(assign.global_id < nglobals);
        auto val = get_value(bseq, assign.value);
        _set_global(assign.global_id, val);
    }

    for (auto &pulse: bseq.pulses) {
        assert(value_computed(bseq, pulse.time));
        auto t = get_time(pulse.time);
        if (t != 0)
            break;
        if (pulse.cond != uint32_t(-1) && get_value(bseq, pulse.cond) == 0)
            continue;
        double startval = pulse.len == uint32_t(-1) ? get_value(bseq, pulse.value) :
            pulse.ramp_func(0, values.data());
        start_values[pulse.chn - 1].f64 = startval;
    }
}

uint32_t HostSeq::eval_branch(const BasicSeq &bseq)
{
    // Globals
    if (m_global_dirty) {
        // Flush global age change.
        m_global_dirty = false;
        m_cur_age += 1;
    }
    for (auto gidx: bseq.cond_global_refs)
        update_global_value(gidx);
    for (auto &br: bseq.branches) {
        if (get_value(bseq, br.cond) == 0)
            continue;
        return br.target;
    }
    return bseq.default_branch;
}

}
