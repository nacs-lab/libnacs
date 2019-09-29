/*************************************************************************
 *   Copyright (c) 2018 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../utils/streams.h"
#include "../../utils/number.h"

#include <math.h>

namespace NaCs::Seq::Zynq {

namespace ByteCode {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len)
{
    // This cannot use the `ExeState` helper to avoid fusing of instructions.
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

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask)
{
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
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

namespace Time {

struct Constraints {
    // preferred pulse spacing when outputting contineous updates
    uint64_t prefer_dt;
    // minimum average pulse spacing and the average window
    uint64_t avg_dt;
    size_t avg_window;
};

struct Keeper {
public:
    Keeper(const Constraints &);
    void addPulse(uint64_t dt);
    uint64_t minDt(uint64_t min_dt) const;
    void reset(void);
private:
    Constraints m_cons;
    size_t m_npulses;
    uint64_t m_total_len;
    std::vector<uint64_t> m_plens;
};

Keeper::Keeper(const Constraints &cons)
    : m_cons(cons),
      m_plens(cons.avg_window - 1)
{
    reset();
}

void Keeper::addPulse(uint64_t dt)
{
    auto max_dt = m_cons.avg_dt * max(m_cons.avg_window / 100, (unsigned)10);
    if (dt > max_dt)
        dt = max_dt;
    size_t pidx = m_npulses % (m_cons.avg_window - 1);
    m_total_len = m_total_len + dt - m_plens[pidx];
    m_plens[pidx] = dt;
    m_npulses++;
}

uint64_t Keeper::minDt(uint64_t min_dt) const
{
    uint64_t window_sz;
    if (m_npulses < m_cons.avg_window - 1) {
        window_sz = m_cons.avg_dt * m_npulses * m_npulses / m_cons.avg_window;
    }
    else {
        window_sz = m_cons.avg_dt * m_cons.avg_window;
    }
    min_dt = std::max<uint64_t>(min_dt, 1);
    if (m_total_len >= window_sz)
        return min_dt;
    return std::max(window_sz - m_total_len, min_dt);
}

void Keeper::reset(void)
{
    m_npulses = 0;
    m_total_len = 0;
    memset(m_plens.data(), 0, m_plens.size() * sizeof(uint64_t));
}

}

/**
 * The writer class maintains all the states for bytecode generation.
 */
class Writer {
    // Hard coded start trigger TTL channel.
    static constexpr int start_ttl = 0;
    static constexpr int start_ttl_mask = (1 << start_ttl);

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

    // States
    uint32_t cur_ttl = 0;
    uint32_t all_ttl_mask = start_ttl_mask;
    bool ttl_set = false;
    DDS dds[22] = {};
    DAC dac[4] = {};

    uint64_t cur_t = 0;
    uint8_t max_time_left = 0;
    ssize_t last_timed_inst = 0;

    buff_ostream &stm;

    template<typename Inst>
    ssize_t addInst(Inst inst)
    {
        auto len = stm.tellg();
        stm.write((char*)&inst, sizeof(inst));
        max_time_left = 0;
        return (ssize_t)len;
    }

    // Increase the wait time encoded in the last instruction by `t`.
    // The caller should have checked that the time fits in the maximum time possible to be
    // encoded.
    void incLastTime(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t b = stm[last_timed_inst];
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
        stm[last_timed_inst] = uint8_t((b & ~tmask) | ((t << 4) & tmask));
    }

    // The `add***` functions below provides an abstraction to the bytecode and hides
    // the detail about the bytecode encoding. The code tries to use the most compact
    // encoding when available and can fallback to generic encoding if not.
    // This is the same level as the API of `ExeState`.
    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        cur_t += dt;
        if (dt <= max_time_left) {
            incLastTime(uint8_t(dt));
            return;
        }
        keeper.addPulse(dt);
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

    int addTTLReal(uint64_t t, uint32_t ttl)
    {
        ttl = ttl & all_ttl_mask;
        cur_ttl = cur_ttl & all_ttl_mask;
        if (ttl == cur_ttl && ttl_set)
            return 0;
        addWait(t - cur_t);
        keeper.addPulse(PulseTime::Min);
        cur_t += PulseTime::Min;
        auto changes = ttl ^ cur_ttl;
        if (!ttl_set) {
            ttl_set = true;
            changes = ttl;
        }
        else {
            assert(changes);
        }
        cur_ttl = ttl;
        auto nchgs = __builtin_popcountll(changes);
        // The case where nchgs == 0, i.e. initial setting to 0, will be handled by
        // the generic case below.
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
        return PulseTime::Min;
    }

    int addTTLSingle(uint64_t t, uint8_t chn, bool val)
    {
        if (val) {
            return addTTLReal(t, cur_ttl | (1 << chn));
        }
        else {
            return addTTLReal(t, cur_ttl & ~uint32_t(1 << chn));
        }
    }

    int addTTL(uint64_t t, uint32_t ttl)
    {
        return addTTLReal(t, ttl & ~start_ttl_mask);
    }

    int addClock(uint64_t t, uint8_t period)
    {
        addWait(t - cur_t);
        cur_t += PulseTime::Clock;
        addInst(Inst::Clock{OpCode::Clock, 0, period});
        keeper.addPulse(PulseTime::Clock);
        return PulseTime::Clock;
    }

    static constexpr double freq_factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;
    int addDDSFreq(uint64_t t, uint8_t chn, double freqf)
    {
        uint32_t freq = uint32_t(0.5 + freqf * freq_factor);
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        if (dds[chn].freq_set && freq == dds[chn].freq)
            return 0;
        addWait(t - cur_t);
        keeper.addPulse(PulseTime::DDSFreq);
        cur_t += PulseTime::DDSFreq;
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
        return PulseTime::DDSFreq;
    }

    int addDDSAmp(uint64_t t, uint8_t chn, double ampf)
    {
        uint16_t amp = uint16_t(ampf * 4095.0 + 0.5);
        if (amp > 4095)
            amp = 4095;
        if (dds[chn].amp_set && amp == dds[chn].amp)
            return 0;
        addWait(t - cur_t);
        keeper.addPulse(PulseTime::DDSAmp);
        cur_t += PulseTime::DDSAmp;
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
        return PulseTime::DDSAmp;
    }

    int addDAC(uint64_t t, uint8_t chn, double Vf)
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
            return 0;
        addWait(t - cur_t);
        keeper.addPulse(PulseTime::DAC);
        cur_t += PulseTime::DAC;
        uint16_t dV = uint16_t(V - dac[chn].V);
        dac[chn].V = V;
        dac[chn].set = true;
        if (dV <= 0x1ff || dV >= 0xfe00) {
            addInst(Inst::DACDet{OpCode::DACDet, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            addInst(Inst::DAC{OpCode::DAC, 0, uint8_t(chn & 0x3), V});
        }
        return PulseTime::DAC;
    }

public:
    Writer(buff_ostream &stm, Time::Keeper keeper)
        : stm(stm),
          keeper(keeper)
    {}

    uint32_t get_ttl_mask() const
    {
        return all_ttl_mask;
    }

    void init_ttl(uint32_t ttl, uint32_t mask)
    {
        cur_ttl = ttl;
        all_ttl_mask = start_ttl_mask | mask;
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

    int addPulse(Channel chn, Val val, uint64_t t, uint64_t tlim)
    {
        int mint;
        if (chn.typ == Channel::TTL) {
            mint = PulseTime::Min;
        }
        else if (chn.typ == Channel::CLOCK) {
            mint = PulseTime::Clock;
        }
        else if (chn.typ == Channel::DAC) {
            mint = PulseTime::DAC;
        }
        else if (chn.typ == Channel::DDS_FREQ) {
            mint = PulseTime::DDSFreq;
        }
        else if (chn.typ == Channel::DDS_AMP) {
            mint = PulseTime::DDSAmp;
        }
        else {
            return -1;
        }
        if (t + mint > tlim)
            return -1;
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

    Time::Keeper keeper;
};

// Easier to deal with than static members.
static constexpr Channel clock_chn{Channel::Type::CLOCK, 0};
static constexpr int default_clock_div = 100;

class Scheduler {
    std::vector<ExpSeq::Pulse> &pulses;
    size_t n_pulses;
    std::map<Channel,Val> &defaults;
    std::vector<ExpSeq::Clock> &clocks;
    Writer writer;

    Time::Constraints t_cons;
    Time::Keeper &keeper;

    /**
     * States for value calculations.
     */
    // Final value for each pulses.
    // This value can be computed from `start_vals` and `pulses` but is cached because
    // it is already computed when initializing `start_vals`. See `init_pulse_vals`
    std::vector<Val> end_vals;
    // The value we should initialize the channel to at the beginning of the sequence.
    std::map<Channel,Val> init_vals;

    /**
     * States for timing.
     */
    // The earliest time we can schedule the next pulse.
    // When we are not hitting the time constraint, this is the time we finish
    // the previous pulse.
    uint64_t next_t = 0;
    // The time of the last pulse.
    uint64_t prev_t = 0;
    // Time offset of start, only the output should care about this.
    uint64_t start_t = 0;

    /**
     * States for clock output.
     */
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

    /**
     * States for keeping track of pulses.
     */
    // The pulses that have been started but not finished yet.
    // This is where we look for pulses to output when we've scheduled more urgent
    // stuff (i.e. clocks, new pulses, etc).
    // The only code that mutate this is `record_pulse` and `finalize_chn`.
    // * `record_pulse` adds a new pulse to `cur_pulses` if it's not finished already.
    // * `finalize_chn` removes pulses from `cur_pulses` after their end time.
    std::map<Channel,size_t> cur_pulses;
    size_t cursor = 0;

    /**
     * Helper functions
     */
    void forward_clock()
    {
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
    }

    // Get the time to output the next pulse at least `dt` after the previous
    // pulse.
    uint64_t get_next_time(uint64_t dt)
    {
        // Hardcoded for now.
        // The minimum wait time we support is `PulseTime::Min`
        // so if the pulse does not happen exactly after the last one
        // it need to happen at least `PulseTime::Min` cycles later...
        auto req_t = prev_t + dt;
        if (next_t >= req_t)
            return next_t;
        if (req_t >= next_t + PulseTime::Min)
            return req_t;
        return next_t + PulseTime::Min;
    }

    // Add a pulse at `t`. The caller is expected to check the time limit
    // before calling this function.
    void output_pulse(Channel cid, Val val, uint64_t t)
    {
        while (true) {
            uint64_t tlim = next_clock_time;
            if (tlim != UINT64_MAX)
                tlim += start_t;
            int min_dt = writer.addPulse(cid, val, t + start_t, tlim);
            if (min_dt < 0) {
                int min_dt = writer.addPulse(clock_chn, Val::get<uint32_t>(next_clock_div),
                                             next_clock_time + start_t, UINT64_MAX);
                prev_t = next_clock_time;
                next_t = prev_t + min_dt;
                forward_clock();
                if (t < next_t)
                    t = next_t;
                continue;
            }
            else if (min_dt != 0) {
                prev_t = t;
                next_t = t + keeper.minDt(min_dt);
            }
            return;
        }
    }

    Val calc_pulse(size_t id, uint64_t t)
    {
        auto &pulse = pulses[id];
        uint64_t rel_t = t < pulse.t ? 0 : t - pulse.t;
        if (rel_t > pulse.len)
            return end_vals[id];
        return pulse(rel_t);
    }

    // If `cleanup` is `false`,
    // the caller must make sure everything before the new pulse is finalized on the channel.
    void record_pulse(size_t id, uint64_t t, bool cleanup=false)
    {
        // id: pulse index
        // t: current time
        auto &pulse = pulses[id];
        assert(cleanup || cur_pulses.find(pulse.chn) == cur_pulses.end());
        if (pulse.t + pulse.len <= t) {
            if (cleanup)
                finalize_chn(pulse.chn);
            return;
        }
        cur_pulses[pulse.chn] = id;
    }

    // * Remove the pulse from `cur_pulses`
    // * Does **NOT** add any output. If the user want to output the final value
    //   of the pulse, the value can be read from the `end_vals`.
    auto finalize_chn(Channel cid) -> decltype(cur_pulses.find(cid))
    {
        auto it = cur_pulses.find(cid);
        if (it == cur_pulses.end())
            return it;
        it = cur_pulses.erase(it);
        return it;
    };

    /* These variables are actually local variables for `handle_overdue`.
     * They are declared here to have their allocations cached.
     */
    // Channels that are finished and needs their final value outputted.
    // This will be overwritten if a new pulse is taking over on the channel.
    std::map<Channel,Val> to_flush;
    // New pulse on the channel. We don't output immediately in case
    // the next pulse starts before that.
    // The one we actually need to output is recorded for the channel in
    // `latest_pending`.
    std::vector<size_t> pending;
    std::map<Channel,size_t> latest_pending;
    // Check if any channel has overdue changes.
    // This includes new pulses or finishing of pulses that should happen
    // before the next preferred time point.
    bool handle_overdue()
    {
        uint64_t tlim = next_t;
        // First collect and finalize pulses that should be finished.
        for (auto it = cur_pulses.begin();it != cur_pulses.end();) {
            auto pid = it->second;
            auto &pulse = pulses[pid];
            if (pulse.t + pulse.len > tlim) {
                ++it;
                continue;
            }
            to_flush.insert({it->first, end_vals[pid]});
            it = finalize_chn(pulse.chn);
        }

        // Now see if there's any new pulses that needs handling.
        // Remove the corresponding pulse from `to_flush` as we encounter it.
        for (;cursor < n_pulses;cursor++) {
            auto &pulse = pulses[cursor];
            if (pulse.t > tlim)
                break;
            auto flush_it = to_flush.find(pulse.chn);
            if (flush_it != to_flush.end())
                to_flush.erase(flush_it);
            // We'll finalize the channel automatically in `record_pulse`
            // when handling the pending pulses below.
            pending.push_back(cursor);
            latest_pending[pulse.chn] = cursor;
        }
        if (pending.empty() && to_flush.empty())
            return false;
        // Flush the ones to finish.
        for (auto it: to_flush)
            output_pulse(it.first, it.second, next_t);
        to_flush.clear();
        // Output queued pulses in the queued order, after the finishing ones.
        for (auto pid: pending) {
            auto &pulse = pulses[pid];
            // This is replaced by a new one.
            if (latest_pending[pulse.chn] != pid)
                continue;
            auto val = calc_pulse(pid, next_t);
            record_pulse(pid, next_t, true);
            output_pulse(pulses[pid].chn, val, next_t);
        }
        pending.clear();
        latest_pending.clear();
        return true;
    }

    void handle_update()
    {
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
                output_pulse(cid, end_vals[pid], next_t);
            }
        } while (has_finalize);
        auto cur_copy = cur_pulses;
        while (next_t <= deadline && !cur_copy.empty() && cursor < n_pulses) {
            uint64_t next_seq_t = get_next_time(t_cons.prefer_dt);
            auto &new_pulse = pulses[cursor];
            if (new_pulse.t <= next_seq_t) {
                uint64_t dt = (new_pulse.t > prev_t ? new_pulse.t - prev_t : 0);
                uint64_t t = get_next_time(dt);
                // Found a new pulse to output.
                auto it = cur_copy.find(new_pulse.chn);
                if (it != cur_copy.end())
                    cur_copy.erase(it);
                auto val = calc_pulse(cursor, t);
                record_pulse(cursor, t, true);
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
                output_pulse(cid, end_vals[pid], next_seq_t);
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
                output_pulse(cid, end_vals[pid], next_seq_t);
            } else {
                output_pulse(cid, calc_pulse(pid, next_seq_t), next_seq_t);
            }
        }
    }

    void init_pulse_vals()
    {
        // Populate `start_vals` and `end_vals` for all pulses,
        // and `init_vals` for all channels.
        std::map<Channel,Val> cur_vals;
        for (size_t i = 0; i < n_pulses; i++) {
            auto &pulse = pulses[i];
            auto cid = pulse.chn;
            auto it = cur_vals.find(cid);
            Val val;
            if (it == cur_vals.end()) {
                // This is a new channel, find it's default value and populate `init_vals`
                auto dit = defaults.find(cid);
                val = dit == defaults.end() ? Val() : dit->second;
                it = cur_vals.insert({cid, val}).first;
                init_vals[cid] = val;
            }
            else {
                val = it->second;
            }
            pulse.set_start(val);
            val = pulse(pulse.len);
            end_vals[i] = val;
            it->second = val;
        }
    }

public:
    Scheduler(ExpSeq &seq, buff_ostream &stm, Time::Constraints t_cons)
        : pulses(seq.pulses),
          n_pulses(pulses.size()),
          defaults(seq.defaults),
          clocks(seq.clocks),
          writer(stm, t_cons),
          t_cons(t_cons),
          keeper(writer.keeper),
          end_vals(n_pulses)
    {
    }

    void init()
    {
        auto ttl_it = defaults.find({Channel::TTL, 0});
        uint32_t cur_ttl = 0;
        uint32_t all_ttl_mask = 0;
        if (ttl_it != defaults.end())
            cur_ttl = ttl_it->second.val.i32;

        std::stable_sort(pulses.begin(), pulses.end(),
                         [] (auto &p1, auto &p2) { return p1.t < p2.t; });

        // Merge TTL pulses that happen at the same time into a single TTL pulse
        uint32_t ttl_val = cur_ttl;
        uint64_t prev_ttl_t = 0;
        ssize_t prev_ttl_idx = -1;
        size_t to = 0, from = 0;
        for (;from < n_pulses;from++, to++) {
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
            pulse.set_start(Val::get<double>((ttl_val & mask) != 0));
            bool val = pulse(pulse.t).val.f64 != 0;
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
                pulses[prev_ttl_idx].set(Val::get<uint32_t>(new_ttl_val));
                to--;
            } else {
                prev_ttl_idx = to;
                prev_ttl_t = pulse.t;
                pulse.chn = {Channel::Type::TTL, 0};
                pulse.set(Val::get<uint32_t>(new_ttl_val));
                if (from != to) {
                    pulses[to] = std::move(pulse);
                }
            }
        }
        pulses.resize(to);
        n_pulses = to;
        writer.init_ttl(cur_ttl, all_ttl_mask);
        init_pulse_vals();
    }

    // Complexity O(nchannel * npulse)
    // Assume pulses are sorted according to start time.
    uint32_t schedule()
    {
        // Initialize channels
        for (auto it: init_vals) {
            prev_t = get_next_time(t_cons.prefer_dt * 2);
            int min_dt = writer.addPulse(it.first, it.second, prev_t, UINT64_MAX);
            next_t = prev_t;
            if (min_dt > 0) {
                next_t += std::max(t_cons.prefer_dt, uint64_t(min_dt));
            }
        }

        // Start the sequence and restart timer.
        start_t = writer.start(next_t);
        if (clocks.empty()) {
            // Start continuous clock with default divider
            writer.addPulse(clock_chn, Val::get<uint32_t>(default_clock_div),
                            start_t, UINT64_MAX);
            start_t += default_clock_div;
        }
        else {
            // Check if the first clock period needs to be started and setup book keeping vars.
            next_clock_idx = 1;
            auto first_clock = clocks[0];
            if (first_clock.t <= first_clock.div) {
                clock_neg_offset = int(first_clock.div - first_clock.t);
                writer.addPulse(clock_chn, Val::get<uint32_t>(first_clock.div),
                                start_t, UINT64_MAX);
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

        bool prev_overdue = false;
        while (cursor < n_pulses || !cur_pulses.empty()) {
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
            if (cur_pulses.empty() && cursor < n_pulses) {
                auto &pulse = pulses[cursor];
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
            int mindt = writer.addPulse(clock_chn, Val::get<uint32_t>(next_clock_div),
                                        next_clock_time + start_t, UINT64_MAX);
            next_t = next_clock_time + mindt;
            forward_clock();
        }
        writer.end(next_t + start_t);
        return writer.get_ttl_mask();
    }
};

} // anonymous

static uint32_t SeqToByteCode(buff_ostream &stm, ExpSeq &seq)
{
    // The bytecode is guaranteed to not enable any channel that is not present in the mask.
    Scheduler state(seq, stm, {50, 40, 4096});
    state.init();
    return state.schedule();
}

} // ByteCode

NACS_EXPORT() std::vector<uint8_t> ExpSeq::toByteCode(uint32_t *ttl_mask)
{
    basic_vector_ostream<std::vector<uint8_t>> stm;
    auto tm = ByteCode::SeqToByteCode(stm, *this);
    if (ttl_mask)
        *ttl_mask = tm;
    return stm.get_buf();
}

NACS_EXPORT() uint8_t *ExpSeq::toByteCode(size_t *sz, uint32_t *ttl_mask)
{
    malloc_ostream stm;
    auto tm = ByteCode::SeqToByteCode(stm, *this);
    if (ttl_mask)
        *ttl_mask = tm;
    return (uint8_t*)stm.get_buf(*sz);
}

}

using namespace NaCs::Seq::Zynq;

extern "C" NACS_EXPORT() uint8_t *nacs_seq_bin_to_bytecode(const uint32_t *data,
                                                           size_t data_len,
                                                           size_t *code_len,
                                                           uint32_t *ttl_mask)
{
    return ExpSeq::fromBinary(data, data_len).toByteCode(code_len, ttl_mask);
}

extern "C" NACS_EXPORT() uint64_t nacs_seq_bytecode_total_time(const uint8_t *code,
                                                               size_t code_len)
{
    return ByteCode::total_time(code, code_len);
}
