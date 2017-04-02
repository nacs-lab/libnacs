/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include "timing.h"

#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <stdint.h>

#ifndef __NACS_SEQ_SEQ_H__
#define __NACS_SEQ_SEQ_H__

namespace NaCs {
namespace Seq {

template<typename Cid, typename Cb>
struct BasePulse {
    uint64_t t;
    uint64_t len;
    Cid chn;
    // Called with `cb(uint64_t t, ValT start, uint64_t len) -> ValT`
    // where `t` is the time within the pulse,
    // `start` is the value of the channel at the begining of the pulse.
    // `len` is the length pulse.
    Cb cb;
};

template<typename Cid, typename Cb>
static void sort(std::vector<BasePulse<Cid,Cb>> &seq)
{
    std::stable_sort(seq.begin(), seq.end(),
                     [] (auto &p1, auto &p2) { return p1.t < p2.t; });
}

enum class Event {
    start,
    end
};

struct Clock {
    // Time index of the first clock edge. The clock pulse should happen `div` time points
    // before this time.
    uint64_t t;
    // Length of the pulse in time index. The clock end pulse should be `len` time points after
    // the start pulse.
    uint64_t len;
    // Clock division. Each clock cycle consists of `div` number of time points at high
    // and `div` number of time points at low.
    uint32_t div;
};

static constexpr int default_clock_div = 100;

// TTL pulses should be folded first.
template<typename Accum, typename Cid, typename Cb, typename Filter,
         typename SeqCB, typename ValT, typename ConvertClock>
static void schedule(Accum &accum, const std::vector<BasePulse<Cid,Cb>> &seq,
                     const Time::Constraints &t_cons,
                     const std::map<Cid,ValT> &defaults,
                     const std::vector<Clock> &clocks,
                     Filter &&filter, SeqCB &&seq_cb, ConvertClock &&convert_clock, Cid clock_chn)
{
    // Complexity O(nchannel * npulse)
    // Assume pulses are sorted according to start time.

    // `accum` is called with
    // * `accum(Cid cid, ValT val, uint64_t t, uint64_t tlim)`
    //     Add a pulse at `t` to output `val` on channel `cid`.
    //     The pulse should finish before `tlim`. Otherwise, no pulse should be added.
    //     Returns the minimum time the pulse will take, `0` if no pulse was added.

    // `filter` is called with
    // * `filter(Cid cid) -> bool`:
    //     Returns if the specified channel should be handled.
    // * `filter(Cid cid, ValT orig_val, ValT new_val) -> bool`
    //     Returns if the new value is significantly different
    //     from the original value for the channel.

    // `seq_cb` is called with
    // * `seq_cb(Accum accum, uint64_t cur_t, Event evt) -> uint64_t`
    //     Add necessary pulses to start or finish the sequence
    //     (turn on/off start trigger/clock output). Returns the new
    //     current time. (t=0 for start)

    Time::Keeper keeper(t_cons);
    std::map<Cid,ValT> start_vals;
    std::map<Cid,ValT> cur_vals;
    // The earliest time we can schedule the next pulse.
    // When we are not hitting the time constraint, this is the time we finish
    // the previous pulse.
    uint64_t next_t = 0;
    // The time of the last pulse.
    uint64_t prev_t = 0;
    // Time offset of start, only the output should care about this.
    uint64_t start_t = 0;
    // The time when we should output the next clock pulse.
    // Set to `UINT64_MAX` when there isn't any left.
    // The time is shifted by `-clock_div` such that it is the time to output the pulse
    // and not the time the clock edge should be delivered.
    uint64_t next_clock_time = UINT64_MAX;
    // The index of the next clock pulse in `clocks` after the current clock period is finished.
    // If one is currently running (started and haven't end) this is the index to the next one
    // that is not running.
    size_t next_clock_idx = 0;
    // `next_clock_div` == 256 means we are finishing a clock period.
    int next_clock_div = 256;
    int clock_neg_offset = 0;

    // Helper functions
    auto forward_clock = [&] () {
        if (next_clock_div != 256) {
            next_clock_div = 256;
            next_clock_time = next_clock_time + clocks[next_clock_idx - 1].len;
            return;
        }
        if (next_clock_idx >= clocks.size()) {
            next_clock_time = UINT64_MAX;
            return;
        }
        auto next_clock = clocks[next_clock_idx];
        next_clock_idx++;
        next_clock_div = next_clock.div;
        next_clock_time = next_clock.t - clock_neg_offset;
    };
    // Get the time to output the next pulse at least `dt` after the previous
    // pulse.
    auto get_next_time = [&] (uint64_t dt) {
        return std::max(prev_t + dt, next_t);
    };
    // Add a pulse at `t`. The caller is expected to check the time limit
    // before calling this function.
    auto output_pulse = [&] (Cid cid, ValT val, uint64_t t) {
        if (!filter(cid, cur_vals[cid], val))
            return;
        while (true) {
            uint64_t tlim = next_clock_time;
            if (tlim != UINT64_MAX)
                tlim += start_t;
            uint64_t min_dt = accum(cid, val, t + start_t, tlim);
            if (min_dt == 0) {
                keeper.addPulse(next_clock_time - prev_t);
                uint64_t min_dt = accum(clock_chn, convert_clock(next_clock_div),
                                        next_clock_time + start_t, UINT64_MAX);
                prev_t = next_clock_time;
                next_t = prev_t + min_dt;
                forward_clock();
                if (t < next_t)
                    t = next_t;
                continue;
            }
            cur_vals[cid] = val;
            uint64_t dt = t - prev_t;
            keeper.addPulse(dt);
            prev_t = t;
            next_t = t + keeper.minDt(min_dt);
            return;
        }
    };

    // Initialize channels
    for (auto &pulse: seq) {
        auto cid = pulse.chn;
        if (filter(cid) && start_vals.find(cid) == start_vals.end()) {
            auto it = defaults.find(cid);
            ValT def_val = it == defaults.end() ? ValT() : it->second;
            start_vals[cid] = def_val;
            cur_vals[cid] = def_val;
            prev_t = get_next_time(t_cons.prefer_dt * 2);
            uint64_t min_dt = accum(cid, def_val, prev_t, UINT64_MAX);
            next_t = prev_t + std::max(t_cons.prefer_dt, min_dt);
        }
    }

    // Start the sequence and restart timer.
    start_t = seq_cb(accum, next_t, Event::start);
    if (clocks.empty()) {
        // Start continuous clock with default divider
        accum(clock_chn, convert_clock(default_clock_div), start_t, UINT64_MAX);
        start_t += default_clock_div;
    }
    else {
        // Check if the first clock period needs to be started and setup book keeping vars.
        next_clock_idx = 1;
        auto first_clock = clocks[0];
        if (first_clock.t <= first_clock.div) {
            clock_neg_offset = int(first_clock.div - first_clock.t);
            accum(clock_chn, convert_clock(first_clock.div), start_t, UINT64_MAX);
            next_clock_time = first_clock.len - clock_neg_offset;
            start_t += clock_neg_offset;
            next_clock_div = 256;
        }
        else {
            clock_neg_offset = first_clock.div;
            next_clock_time = first_clock.t - first_clock.div;
            next_clock_div = first_clock.div;
        }
    }
    keeper.reset();
    prev_t = next_t = 0;

    auto calc_pulse = [&] (size_t id, uint64_t t) -> ValT {
        auto &pulse = seq[id];
        uint64_t rel_t = t < pulse.t ? 0 : t - pulse.t;
        if (rel_t > pulse.len)
            rel_t = pulse.len;
        auto start = start_vals[pulse.chn];
        return pulse.cb(rel_t, start, pulse.len);
    };

    std::map<Cid,size_t> cur_pulses;
    std::set<size_t> finalized;
    size_t cursor = 0;
    size_t npulse = seq.size();

    auto record_pulse = [&] (size_t id, uint64_t t) {
        auto &pulse = seq[id];
        if (pulse.t + pulse.len <= t)
            return;
        cur_pulses[pulse.chn] = id;
    };

    // * Update `start_val` to the pulse's final value
    // * Remove the pulse from `cur_pulses`
    // * Add the pulse to `finalized`
    // * Does **NOT** add any output. If the user want to output the final value
    //   of the pulse, it can read the value from the updated `start_val`.
    //   Similarly, it doesn't update `cur_vals` either since that holds the
    //   value of the latest pulse.
    auto finalize_chn = [&] (Cid cid) {
        auto it = cur_pulses.find(cid);
        if (it == cur_pulses.end())
            return it;
        size_t pid = it->second;
        auto &pulse = seq[pid];
        auto fin_val = calc_pulse(pid, pulse.t + pulse.len);
        start_vals[cid] = fin_val;
        it = cur_pulses.erase(it);
        finalized.insert(pid);
        return it;
    };

    // Check if any channel has overdue changes.
    // This includes new pulses or finishing of pulses that should happen
    // before the next preferred time point.
    std::set<Cid> to_flush;
    std::vector<size_t> pending;
    auto handle_overdue = [&] () {
        uint64_t tlim = next_t;
        // First collect and finalize pulses that should be finished.
        for (auto it = cur_pulses.begin();it != cur_pulses.end();) {
            auto pid = it->second;
            auto &pulse = seq[pid];
            if (pulse.t + pulse.len > tlim) {
                ++it;
                continue;
            }
            to_flush.insert(it->first);
            it = finalize_chn(pulse.chn);
        }

        // Now see if there's any new pulses that needs handling.
        // Remove the corresponding pulse from `to_flush` as we encounter it.
        for (;cursor < npulse;cursor++) {
            auto &pulse = seq[cursor];
            if (pulse.t > tlim)
                break;
            if (!filter(pulse.chn))
                continue;
            auto flush_it = to_flush.find(pulse.chn);
            if (flush_it != to_flush.end())
                to_flush.erase(flush_it);
            finalize_chn(pulse.chn);
            pending.push_back(cursor);
        }
        if (pending.empty() && to_flush.empty())
            return false;
        // Flush the ones to finish.
        if (!to_flush.empty()) {
            for (auto cid: to_flush)
                output_pulse(cid, start_vals[cid], next_t);
            to_flush.clear();
        }
        // Output queued pulses in the queued order, after the finishing ones.
        if (!pending.empty()) {
            for (auto pid: pending) {
                if (finalized.find(pid) != finalized.end())
                    continue;
                auto val = calc_pulse(pid, next_t);
                record_pulse(pid, next_t);
                output_pulse(seq[pid].chn, val, next_t);
            }
            pending.clear();
        }
        return true;
    };

    auto handle_update = [&] () {
        size_t num_updates = cur_pulses.size();
        // Compute a deadline after which we will ignore any new pulses and
        // focuses on updating the on-going updates.
        uint64_t deadline = get_next_time(t_cons.prefer_dt * num_updates * 2);
        // Check repeatedly if there's any channels to be finalized
        bool has_finalize;
        do {
            has_finalize = false;
            for (auto it = cur_pulses.begin();it != cur_pulses.end();) {
                auto pid = it->second;
                auto cid = it->first;
                auto &pulse = seq[pid];
                if (pulse.t + pulse.len > next_t) {
                    ++it;
                    continue;
                }
                has_finalize = true;
                it = finalize_chn(cid);
                output_pulse(cid, start_vals[cid], next_t);
            }
        } while (has_finalize);
        auto cur_copy = cur_pulses;
        while (next_t <= deadline && !cur_copy.empty() && cursor < npulse) {
            uint64_t next_seq_t = get_next_time(t_cons.prefer_dt);
            auto &new_pulse = seq[cursor];
            if (!filter(new_pulse.chn)) {
                cursor++;
                continue;
            }
            if (new_pulse.t <= next_seq_t) {
                uint64_t dt = (new_pulse.t > prev_t ? new_pulse.t - prev_t : 0);
                uint64_t t = get_next_time(dt);
                auto val = calc_pulse(cursor, t);
                record_pulse(cursor, t);
                output_pulse(new_pulse.chn, val, t);
                cursor++;
                continue;
            }
            auto it = cur_copy.begin();
            auto pid = it->second;
            auto cid = it->first;
            cur_copy.erase(it);
            auto &pulse = seq[pid];
            if (pulse.t + pulse.len <= next_seq_t) {
                finalize_chn(cid);
                output_pulse(cid, start_vals[cid], next_seq_t);
            } else {
                output_pulse(cid, calc_pulse(pid, next_seq_t), next_seq_t);
            }
        }
        // If we reached the deadline with pending updates, flush all of them.
        for (auto &it: cur_copy) {
            uint64_t next_seq_t = get_next_time(t_cons.prefer_dt);
            auto pid = it.second;
            auto cid = it.first;
            auto &pulse = seq[pid];
            if (pulse.t + pulse.len <= next_seq_t) {
                finalize_chn(cid);
                output_pulse(cid, start_vals[cid], next_seq_t);
            } else {
                output_pulse(cid, calc_pulse(pid, next_seq_t), next_seq_t);
            }
        }
    };

    bool prev_overdue = false;
    while (cursor < npulse || !cur_pulses.empty()) {
        uint64_t old_next_t = next_t;
        // 1. Check and handle any overdue changes.
        //
        //     If there's output in this step restart the loop, unless the
        //     previous loop is also aborted here.
        bool has_overdue = handle_overdue();
        if (has_overdue && !prev_overdue) {
            prev_overdue = true;
            continue;
        }
        prev_overdue = false;

        // 2. Output pulses
        handle_update();

        // 3. Forward time if there's no on-going pulses.
        if (cur_pulses.empty() && cursor < npulse) {
            auto &pulse = seq[cursor];
            if (!filter(pulse.chn)) {
                cursor++;
                continue;
            }
            // These should be no-op but just to be safe
            uint64_t dt = (pulse.t > prev_t ? pulse.t - prev_t : 0);
            uint64_t t = get_next_time(dt);
            auto val = calc_pulse(cursor, t);
            record_pulse(cursor, t);
            output_pulse(pulse.chn, val, t);
            cursor++;
            continue;
        }
        if (next_t == old_next_t) {
            next_t = next_t + t_cons.prefer_dt;
        }
    }
    while (next_clock_time != UINT64_MAX) {
        uint64_t mindt = accum(clock_chn, convert_clock(next_clock_div),
                               next_clock_time + start_t, UINT64_MAX);
        next_t = next_clock_time + mindt;
        forward_clock();
    }
    seq_cb(accum, next_t + start_t, Event::end);
}

}
}

#endif
