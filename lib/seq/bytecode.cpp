/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "exehelper_p.h"

#include "bytecode.h"

namespace NaCs {
namespace Seq {

namespace ByteCode {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len)
{
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b & 0xf;
        assert(op < 14);
        auto len = inst_size[op];
        i += len;
        count++;
    }
    return count;
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len)
{
    Printer printer{stm};
    ExeState state;
    state.run(printer, code, code_len);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len)
{
    TimeKeeper keeper;
    ExeState state;
    state.run(keeper, code, code_len);
    return keeper.total_t;
}

namespace {

struct Filter {
    bool operator()(Channel)
    {
        return true;
    }
    bool operator()(Channel cid, Val val1, Val val2)
    {
        switch (cid.typ) {
        case Channel::TTL:
            return val1.val.i32 != val2.val.i32;
        case Channel::DDS_FREQ:
            return fabs(val1.val.f64 - val2.val.f64) > 0.2;
        case Channel::DDS_AMP:
            return fabs(val1.val.f64 - val2.val.f64) > 0.00005;
        case Channel::DAC:
            return fabs(val1.val.f64 - val2.val.f64) > 0.0002;
        default:
            return val1.val.f64 != val2.val.f64;
        }
    }
};

template<typename Vec>
struct ScheduleState {
    struct DDS {
        uint32_t freq = 0;
        uint16_t amp = 0;
        bool freq_set = false;
        bool amp_set = false;
    };
    struct DAC {
        uint16_t V = 0;
        bool set = false;
    };

    ScheduleState(Sequence &seq)
        : seq(seq)
    {
        auto &pulses = seq.pulses;
        auto &defaults = seq.defaults;
        auto ttl_it = defaults.find({Channel::TTL, 0});
        if (ttl_it != defaults.end())
            cur_ttl = ttl_it->second.val.i32;

        sort(pulses);

        // Merge TTL pulses that happen at the same time into a single TTL pulse
        uint32_t ttl_val = cur_ttl;
        uint64_t prev_ttl_t = 0;
        ssize_t prev_ttl_idx = -1;
        size_t to = 0, from = 0;
        for (;from < pulses.size();from++, to++) {
            auto &pulse = pulses[from];
            if (pulse.chn.typ != Channel::TTL) {
                if (from != to)
                    pulses[to] = std::move(pulse);
                continue;
            }
            assert(pulse.len == 0);
            assert(pulse.chn.id < 32);
            uint32_t mask = uint32_t(1) << pulse.chn.id;
            all_ttl_mask = all_ttl_mask | mask;
            bool val = pulse.cb(pulse.t, Val::get<double>((ttl_val & mask) != 0)).val.f64 != 0;
            uint32_t new_ttl_val;
            if (val) {
                new_ttl_val = ttl_val | mask;
            } else {
                new_ttl_val = ttl_val & ~mask;
            }
            if (new_ttl_val == ttl_val && prev_ttl_idx != -1) {
                to--;
                continue;
            }
            ttl_val = new_ttl_val;
            if (prev_ttl_idx != -1 && prev_ttl_t == pulse.t) {
                pulses[prev_ttl_idx].cb = PulseData(Val::get<uint32_t>(new_ttl_val));
                to--;
            } else {
                prev_ttl_idx = to;
                prev_ttl_t = pulse.t;
                pulse.chn = {Channel::Type::TTL, 0};
                pulse.cb = PulseData(Val::get<uint32_t>(new_ttl_val));
                if (from != to) {
                    pulses[to] = std::move(pulse);
                }
            }
        }
        pulses.resize(to);
    }

    uint64_t addPulse(Channel chn, Val val, uint64_t t, uint64_t tlim)
    {
        uint64_t mint = 50;
        if (chn.typ == Channel::TTL) {
            mint = 3;
        }
        else if (chn.typ == Channel::CLOCK) {
            mint = 5;
        }
        else if (chn.typ == Channel::DAC) {
            mint = 45;
        }
        if (t + mint > tlim)
            return 0;
        switch (chn.typ) {
        case Channel::TTL:
            return addTTL(t, val.val.i32);
        case Channel::DDS_FREQ:
            return addDDSFreq(t, uint8_t(chn.id), val.val.f64);
        case Channel::DDS_AMP:
            return addDDSAmp(t, uint8_t(chn.id), val.val.f64);
        case Channel::DAC:
            return addDAC(t, uint8_t(chn.id), val.val.f64);
        case Channel::CLOCK:
            return addClock(t, uint8_t(val.val.i32 - 1));
        default:
            throw std::runtime_error("Invalid Pulse.");
        }
        return mint;
    }

    Sequence &seq;
    uint64_t cur_t = 0;
    uint32_t cur_ttl = 0;
    uint32_t all_ttl_mask = start_ttl_mask;
    bool ttl_set = false;
    DDS dds[22] = {};
    DAC dac[4] = {};
    size_t last_timed_inst = 0;
    uint8_t max_time_left = 0;
    Vec code;

    static constexpr int start_ttl = 0;
    static constexpr int start_ttl_mask = (1 << start_ttl);
    static constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;

    template<typename Inst>
    size_t addInst(Inst inst)
    {
        auto len = code.size();
        code.resize(len + sizeof(inst));
        memcpy(&code[len], &inst, sizeof(inst));
        max_time_left = 0;
        return len;
    }

    void incLastTime(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t b = code[last_timed_inst];
        uint8_t op = b & 0xf;
        uint8_t tmask = 0xf0;
        switch (op) {
        case 0: // TTLAll
            break;
        case 1: // TTL2
            tmask = 0x30;
            break;
        case 3: // TTL5
            tmask = 0x70;
            break;
        default:
            abort();
        }
        t = uint8_t(t + ((b & tmask) >> 4));
        code[last_timed_inst] = uint8_t((b & ~tmask) | ((t << 4) & tmask));
    }

    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        cur_t += dt;
        if (dt <= max_time_left) {
            incLastTime(uint8_t(dt));
            return;
        }
        uint16_t times[17];
        memset(times, 0, sizeof(times));
        for (int i = 15; i >= 0; i--) {
            // The first number not representable by the next one is
            // 0x10000 * 2^(3 * i - 3)
            // or
            // 0x2000 << (3 * i)
            if (i > 0) {
                uint64_t tnext_max = uint64_t(0x2000) << (3 * i);
                if (dt < tnext_max) {
                    continue;
                }
            }
            else if (dt <= uint16_t(2047 + max_time_left)) {
                if (max_time_left) {
                    dt -= max_time_left;
                    incLastTime(max_time_left);
                }
                times[16] = uint16_t(dt);
                dt = 0;
                break;
            }
            times[i] = uint16_t(dt >> (3 * i));
            dt -= uint64_t(times[i]) << (3 * i);
            if (dt == 0) {
                break;
            }
            else if (dt <= max_time_left) {
                incLastTime(uint8_t(dt));
                break;
            }
        }
        for (int i = 15; i >= 0; i--) {
            if (!times[i])
                continue;
            addInst(Inst::Wait{OpCode::Wait, uint8_t(i & 0xf), times[i]});
        }
        if (times[16]) {
            addInst(Inst::Wait2{OpCode::Wait2, 1, uint16_t(times[16] & 2047)});
        }
    }

    uint64_t addTTLReal(uint64_t t, uint32_t ttl)
    {
        ttl = ttl & all_ttl_mask;
        cur_ttl = cur_ttl & all_ttl_mask;
        if (ttl == cur_ttl && ttl_set)
            return 1;
        addWait(t - cur_t);
        cur_t += 3;
        auto changes = ttl ^ cur_ttl;
        if (!ttl_set) {
            ttl_set = true;
            changes = ttl;
        }
        cur_ttl = ttl;
        auto nchgs = __builtin_popcountll(changes);
        assert(nchgs > 0);
        if (nchgs == 1) {
            auto bit = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit & 0x1f),
                        uint8_t(bit & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 2) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 3) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            addInst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit3 & 0x1f)});
        }
        else if (nchgs == 4) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit3);
            auto bit4 = __builtin_ffs(changes) - 1;
            addInst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit4 & 0x1f)});
        }
        else if (nchgs == 5) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit3);
            auto bit4 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit4);
            auto bit5 = __builtin_ffs(changes) - 1;
            last_timed_inst = addInst(Inst::TTL5{OpCode::TTL5, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit5 & 0x1f)});
            max_time_left = 7;
        }
        else {
            last_timed_inst = addInst(Inst::TTLAll{OpCode::TTLAll, 0, ttl});
            max_time_left = 15;
        }
        return 3;
    }

    uint64_t addTTLSingle(uint64_t t, uint8_t chn, bool val)
    {
        if (val) {
            return addTTLReal(t, cur_ttl | (1 << chn));
        }
        else {
            return addTTLReal(t, cur_ttl & ~uint32_t(1 << chn));
        }
    }

    uint64_t addTTL(uint64_t t, uint32_t ttl)
    {
        return addTTLReal(t, ttl & ~start_ttl_mask);
    }

    uint64_t addClock(uint64_t t, uint8_t period)
    {
        addWait(t - cur_t);
        cur_t += 5;
        addInst(Inst::Clock{OpCode::Clock, 0, period});
        return 5;
    }

    uint64_t addDDSFreq(uint64_t t, uint8_t chn, double freqf)
    {
        uint32_t freq = uint32_t(0.5 + freqf * freq_factor);
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        if (dds[chn].freq_set && freq == dds[chn].freq)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint32_t dfreq = freq - dds[chn].freq;
        dds[chn].freq = freq;
        dds[chn].freq_set = true;
        if (dfreq <= 0x3f || dfreq >= 0xffffffc0) {
            addInst(Inst::DDSDetFreq2{OpCode::DDSDetFreq2, uint8_t(chn & 0x1f),
                        uint8_t(dfreq & 0x7f)});
        }
        else if (dfreq <= 0x3fff || dfreq >= 0xffffc000) {
            addInst(Inst::DDSDetFreq3{OpCode::DDSDetFreq3, uint8_t(chn & 0x1f),
                        uint16_t(dfreq & 0x7fff)});
        }
        else if (dfreq <= 0x003fffff || dfreq >= 0xffc00000) {
            addInst(Inst::DDSDetFreq4{OpCode::DDSDetFreq4, uint8_t(chn & 0x1f),
                        uint32_t(dfreq & 0x7fffff)});
        }
        else {
            addInst(Inst::DDSFreq{OpCode::DDSFreq, uint8_t(chn & 0x1f),
                        uint32_t(freq & 0x7fffffff)});
        }
        return 50;
    }

    uint64_t addDDSAmp(uint64_t t, uint8_t chn, double ampf)
    {
        uint16_t amp = uint16_t(ampf * 4095.0 + 0.5);
        if (amp > 4095)
            amp = 4095;
        if (dds[chn].amp_set && amp == dds[chn].amp)
            return 1;
        addWait(t - cur_t);
        cur_t += 50;
        uint16_t damp = uint16_t(amp - dds[chn].amp);
        dds[chn].amp = amp;
        dds[chn].amp_set = true;
        if (damp <= 0x3f || damp >= 0xffc0) {
            addInst(Inst::DDSDetAmp{OpCode::DDSDetAmp, uint8_t(chn & 0x1f),
                        uint8_t(damp & 0x7f)});
        }
        else {
            addInst(Inst::DDSAmp{OpCode::DDSAmp, 0, uint8_t(chn & 0x1f),
                        uint16_t(amp & 0xfff)});
        }
        return 50;
    }

    uint64_t addDAC(uint64_t t, uint8_t chn, double Vf)
    {
        uint16_t V;
        if (Vf >= 10) {
            V = 0;
        }
        else if (Vf <= -10) {
            V = 0xffff;
        }
        else {
            // this is for the DAC8814 chip in SPI0
            double scale = 65535 / 20.0;
            double offset = 10.0;
            V = uint16_t(((offset - Vf) * scale) + 0.5);
        }
        if (dac[chn].set && V == dac[chn].V)
            return 1;
        addWait(t - cur_t);
        cur_t += 45;
        uint16_t dV = uint16_t(V - dac[chn].V);
        dac[chn].V = V;
        dac[chn].set = true;
        if (dV <= 0x1ff || dV >= 0xfe00) {
            addInst(Inst::DACDet{OpCode::DACDet, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            addInst(Inst::DAC{OpCode::DAC, 0, uint8_t(chn & 0x3), V});
        }
        return 45;
    }

    uint64_t start(uint64_t t)
    {
        // wait 20us
        t += 2000;
        addTTLSingle(t, start_ttl, true);
        // 1us
        t += 100;
        addTTLSingle(t, start_ttl, false);
        // 5us
        return t + 500;
    }

    void end(uint64_t t)
    {
        // 1us
        t += 100;
        addWait(t - cur_t);
    }

    uint32_t schedule(Time::Constraints t_cons)
    {
        static constexpr int default_clock_div = 100;
        Filter filter{};
        auto &pulses = seq.pulses;
        auto &clocks = seq.clocks;
        auto &defaults = seq.defaults;

        // Complexity O(nchannel * npulse)
        // Assume pulses are sorted according to start time.
        Channel clock_chn{Channel::Type::CLOCK, 0};

        // `filter` is called with
        // * `filter(Channel cid) -> bool`:
        //     Returns if the specified channel should be handled.
        // * `filter(Channel cid, Val orig_val, Val new_val) -> bool`
        //     Returns if the new value is significantly different
        //     from the original value for the channel.

        Time::Keeper keeper(t_cons);
        std::map<Channel,Val> start_vals;
        std::map<Channel,Val> cur_vals;
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
        auto output_pulse = [&] (Channel cid, Val val, uint64_t t) {
            if (!filter(cid, cur_vals[cid], val))
                return;
            while (true) {
                uint64_t tlim = next_clock_time;
                if (tlim != UINT64_MAX)
                    tlim += start_t;
                uint64_t min_dt = addPulse(cid, val, t + start_t, tlim);
                if (min_dt == 0) {
                    keeper.addPulse(next_clock_time - prev_t);
                    uint64_t min_dt = addPulse(clock_chn, Val::get<uint32_t>(next_clock_div),
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
        for (auto &pulse: pulses) {
            auto cid = pulse.chn;
            if (filter(cid) && start_vals.find(cid) == start_vals.end()) {
                auto it = defaults.find(cid);
                Val def_val = it == defaults.end() ? Val() : it->second;
                start_vals[cid] = def_val;
                cur_vals[cid] = def_val;
                prev_t = get_next_time(t_cons.prefer_dt * 2);
                uint64_t min_dt = addPulse(cid, def_val, prev_t, UINT64_MAX);
                next_t = prev_t + std::max(t_cons.prefer_dt, min_dt);
            }
        }

        // Start the sequence and restart timer.
        start_t = start(next_t);
        if (clocks.empty()) {
            // Start continuous clock with default divider
            addPulse(clock_chn, Val::get<uint32_t>(default_clock_div), start_t, UINT64_MAX);
            start_t += default_clock_div;
        }
        else {
            // Check if the first clock period needs to be started and setup book keeping vars.
            next_clock_idx = 1;
            auto first_clock = clocks[0];
            if (first_clock.t <= first_clock.div) {
                clock_neg_offset = int(first_clock.div - first_clock.t);
                addPulse(clock_chn, Val::get<uint32_t>(first_clock.div), start_t, UINT64_MAX);
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

        auto calc_pulse = [&] (size_t id, uint64_t t) -> Val {
            auto &pulse = pulses[id];
            uint64_t rel_t = t < pulse.t ? 0 : t - pulse.t;
            if (rel_t > pulse.len)
                rel_t = pulse.len;
            auto start = start_vals[pulse.chn];
            return pulse.cb(rel_t, start);
        };

        std::map<Channel,size_t> cur_pulses;
        std::set<size_t> finalized;
        size_t cursor = 0;
        size_t npulse = pulses.size();

        auto record_pulse = [&] (size_t id, uint64_t t) {
            auto &pulse = pulses[id];
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
        auto finalize_chn = [&] (Channel cid) {
            auto it = cur_pulses.find(cid);
            if (it == cur_pulses.end())
                return it;
            size_t pid = it->second;
            auto &pulse = pulses[pid];
            auto fin_val = calc_pulse(pid, pulse.t + pulse.len);
            start_vals[cid] = fin_val;
            it = cur_pulses.erase(it);
            finalized.insert(pid);
            return it;
        };

        // Check if any channel has overdue changes.
        // This includes new pulses or finishing of pulses that should happen
        // before the next preferred time point.
        std::set<Channel> to_flush;
        std::vector<size_t> pending;
        auto handle_overdue = [&] () {
            uint64_t tlim = next_t;
            // First collect and finalize pulses that should be finished.
            for (auto it = cur_pulses.begin();it != cur_pulses.end();) {
                auto pid = it->second;
                auto &pulse = pulses[pid];
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
                auto &pulse = pulses[cursor];
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
                    output_pulse(pulses[pid].chn, val, next_t);
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
                    auto &pulse = pulses[pid];
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
                auto &new_pulse = pulses[cursor];
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
                auto &pulse = pulses[pid];
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
                auto &pulse = pulses[pid];
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
                auto &pulse = pulses[cursor];
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
            uint64_t mindt = addPulse(clock_chn, Val::get<uint32_t>(next_clock_div),
                                      next_clock_time + start_t, UINT64_MAX);
            next_t = next_clock_time + mindt;
            forward_clock();
        }
        end(next_t + start_t);
        return all_ttl_mask;
    }
};

} // anonymous

template<typename Vec>
Vec SeqToByteCode(Sequence &seq, uint32_t *ttl_mask)
{
    // The bytecode is guaranteed to not enable any channel that is not present in the mask.
    ScheduleState<Vec> state(seq);
    auto ttl_mask_v = state.schedule({50, 40, 4096});
    if (ttl_mask)
        *ttl_mask = ttl_mask_v;
    return state.code;
}

} // ByteCode

NACS_EXPORT() std::vector<uint8_t> Sequence::toByteCode(uint32_t *ttl_mask)
{
    return ByteCode::SeqToByteCode<std::vector<uint8_t>>(*this, ttl_mask);
}

namespace {

// Simple wrapper around malloc for C interop.
// Implement an API close enough to `std::vector` for our purpose.
// Note that the destructor of this type does **NOT** free the memory.
template<typename ET>
struct MallocVector {
    MallocVector()
        : m_data(nullptr),
          m_len(0)
    {}
    static_assert(std::is_trivial<ET>::value, "");
    ET &operator[](size_t i)
    {
        return m_data[i];
    }
    const ET &operator[](size_t i) const
    {
        return m_data[i];
    }
    size_t size() const
    {
        return m_len;
    }
    void resize(size_t sz)
    {
        m_data = (ET*)realloc(m_data, sizeof(ET) * sz);
        m_len = sz;
    }
private:
    ET *m_data;
    size_t m_len;
};

}

NACS_EXPORT() uint8_t *Sequence::toByteCode(size_t *sz, uint32_t *ttl_mask)
{
    auto vec = ByteCode::SeqToByteCode<MallocVector<uint8_t>>(*this, ttl_mask);
    *sz = vec.size();
    return &vec[0];
}

}
}

using namespace NaCs::Seq;

extern "C" NACS_EXPORT() uint8_t *nacs_seq_bin_to_bytecode(const uint32_t *data,
                                                           size_t data_len,
                                                           size_t *code_len,
                                                           uint32_t *ttl_mask)
{
    return Sequence::fromBinary(data, data_len).toByteCode(code_len, ttl_mask);
}

extern "C" NACS_EXPORT() uint64_t nacs_seq_bytecode_total_time(const uint8_t *code,
                                                               size_t code_len)
{
    return ByteCode::total_time(code, code_len);
}
