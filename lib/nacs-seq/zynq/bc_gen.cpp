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

#include "bc_gen.h"

#include "bytecode.h"
#include "dma.h"
#include "utils.h"
#include "../../nacs-utils/processor.h"

#include <algorithm>

#include <string.h>

namespace NaCs::Seq::Zynq {

struct BCGen::Pulse {
    PulseType pulse_type;
    ChnType chn_type;
    uint8_t chn;
    bool need_endvalue;
    uint32_t id;
    uint32_t endvalue;
    int64_t time; // in sequence unit before sorting and 10ns unit after sorting
    // 0 for non-ramps or pulses that are too short
    int64_t len; // in 10ns unit
    void (*ramp_func)(void);
};

namespace {

// This needs to be consistent with the compiler of the ramp functions.
struct StreamBuffer {
    StreamBuffer(uint32_t fpga_clock_div)
        : fpga_clock_div(fpga_clock_div)
    {
#if NACS_CPU_X86 || NACS_CPU_X86_64
        for (int i = 0; i < 8; i++) {
            static_assert(sizeof(time) == 8 * sizeof(double));
            time[i] = i * fpga_clock_div;
        }
#elif NACS_CPU_AARCH64
        for (int i = 0; i < 2; i++) {
            static_assert(sizeof(time) == 2 * sizeof(double));
            time[i] = i * fpga_clock_div;
        }
#endif
    }
#if NACS_CPU_X86 || NACS_CPU_X86_64
    double time[8] __attribute__((aligned(64)));
    static constexpr double zeros[8] __attribute__((aligned(64))) = {};
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    double time[2] __attribute__((aligned(16)));
    static constexpr double zeros[2] __attribute__((aligned(16))) = {};
#endif
    uint32_t fpga_clock_div;
};

}

struct BCGen::DataStream {
    // TODO: use actual pulse spacing to refine the step size prediction.
    // This can happen when the ramp are slow and do not result in updates.
    bool is_vector;
    uint32_t endvalue;
    int64_t time; // in 10ns unit
    void (*ramp_func)(void) = nullptr;

    double rel_time(int64_t t, const StreamBuffer &streambuffer)
    {
        return double((t - time) * streambuffer.fpga_clock_div);
    }

    double compute_scalar(int64_t t, const void *data,
                          const StreamBuffer &streambuffer)
    {
        auto func = (double (*)(double, const void*))ramp_func;
        return func(rel_time(t, streambuffer), data);
    }

    static double compute_start_val(PulseType pulse_type, void (*ramp_func)(void),
                                    const void *data)
    {
        assert(pulse_type == PulseType::Vector || pulse_type == PulseType::Scalar);
        if (pulse_type == PulseType::Scalar) {
            auto func = (double (*)(double, const void*))ramp_func;
            return func(0, data);
        }
#if NACS_CPU_X86 || NACS_CPU_X86_64
        double buffer[8] __attribute__((aligned(64)));
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&StreamBuffer::zeros, data);
        return buffer[0];
#elif NACS_CPU_AARCH64
        double buffer[2] __attribute__((aligned(16)));
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&StreamBuffer::zeros, data);
        return buffer[0];
#else
        assert(pulse_type == PulseType::Scalar);
#endif
    }

    int64_t time_base;
    int64_t time_step;
    int64_t nsteps;

    void clear()
    {
        ramp_func = nullptr;
    }
    void add_pulse(const Pulse &pulse)
    {
        assert(pulse.pulse_type == PulseType::Vector ||
               pulse.pulse_type == PulseType::Scalar);
        is_vector = pulse.pulse_type == PulseType::Vector;
        endvalue = pulse.endvalue;
        time = pulse.time;
        ramp_func = pulse.ramp_func;
        assert(ramp_func);
        nsteps = 0;
    }

    // This needs to be consistent with the compiler of the ramp functions.
#if NACS_CPU_X86 || NACS_CPU_X86_64
    static double (DataStream::*compute_vector)(int64_t, const void*, int64_t,
                                                const StreamBuffer&);
    double buffer[8] __attribute__((aligned(64)));
    static void init()
    {
        // We could use a PLT style callback to resolve this
        // but since pointer to method may not be pointers we may not be able to use
        // a single relaxed atomic for load and store.
        // Using a local static variable has the additional check on a flag.
        auto &host = CPUInfo::get_host();
        if (host.test_feature(X86::Feature::avx512f) &&
            host.test_feature(X86::Feature::avx512dq)) {
            compute_vector = &DataStream::compute_vector_avx512;
        }
        else if (host.test_feature(X86::Feature::avx2) &&
                 host.test_feature(X86::Feature::fma)) {
            compute_vector = &DataStream::compute_vector_avx2;
        }
        else if (host.test_feature(X86::Feature::avx)) {
            compute_vector = &DataStream::compute_vector_avx;
        }
        else {
            compute_vector = &DataStream::compute_vector_sse;
        }
    }
    double compute(int64_t t, const void *data, int64_t step, const StreamBuffer &streambuffer)
    {
        if (!is_vector)
            return compute_scalar(t, data, streambuffer);
        return (this->*compute_vector)(t, data, step, streambuffer);
    }
    double compute_vector_sse(int64_t t, const void *data, int64_t step,
                              const StreamBuffer &streambuffer)
    {
        if (nsteps > 0) {
            // We compute two numbers at the same time for SSE and have at most one extra...
            assert(nsteps == 1);
            if (time_base == t) {
                nsteps = 0;
                return buffer[1];
            }
        }
        auto ts __attribute__((aligned(16))) =
            rel_time(t, streambuffer) + _mm_load_pd(streambuffer.time) * (double)step;
        time_base = t + step;
        nsteps = 1;
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&ts, data);
        return buffer[0];
    }
    __attribute__((target("avx")))
    double compute_vector_avx(int64_t t, const void *data, int64_t step,
                              const StreamBuffer &streambuffer)
    {
        if (nsteps > 0) {
            if (time_base == t) {
                nsteps--;
                time_base += time_step;
                return buffer[3 - nsteps];
            }
            else if (time_base + time_step == t && nsteps > 1) {
                nsteps -= 2;
                time_base += time_step * 2;
                return buffer[3 - nsteps];
            }
        }
        auto ts __attribute__((aligned(32))) =
            rel_time(t, streambuffer) + _mm256_load_pd(streambuffer.time) * (double)step;
        time_base = t + step;
        nsteps = 3;
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&ts, data);
        return buffer[0];
    }
    __attribute__((target("avx2,fma")))
    double compute_vector_avx2(int64_t t, const void *data, int64_t step,
                               const StreamBuffer &streambuffer)
    {
        // Same as the AVX version but with FMA
        if (nsteps > 0) {
            if (time_base == t) {
                nsteps--;
                time_base += time_step;
                return buffer[3 - nsteps];
            }
            else if (time_base + time_step == t && nsteps > 1) {
                nsteps -= 2;
                time_base += time_step * 2;
                return buffer[3 - nsteps];
            }
        }
        auto ts __attribute__((aligned(32))) =
            rel_time(t, streambuffer) + _mm256_load_pd(streambuffer.time) * (double)step;
        time_base = t + step;
        nsteps = 3;
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&ts, data);
        return buffer[0];
    }
    __attribute__((target("avx512f,avx512dq")))
    double compute_vector_avx512(int64_t t, const void *data, int64_t step,
                                 const StreamBuffer &streambuffer)
    {
        if (nsteps > 0) {
            if (time_base == t) {
                nsteps--;
                time_base += time_step;
                return buffer[7 - nsteps];
            }
            else if (time_base + time_step == t && nsteps > 1) {
                nsteps -= 2;
                time_base += time_step * 2;
                return buffer[7 - nsteps];
            }
        }
        auto ts __attribute__((aligned(64))) =
            rel_time(t, streambuffer) + _mm512_load_pd(streambuffer.time) * (double)step;
        time_base = t + step;
        nsteps = 7;
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&ts, data);
        return buffer[0];
    }
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    double buffer[2] __attribute__((aligned(16)));
    static void init()
    {
    }
    double compute(int64_t t, const void *data, int64_t step, const StreamBuffer &streambuffer)
    {
        if (!is_vector)
            return compute_scalar(t, data, streambuffer);
        if (nsteps > 0) {
            // We compute up to two numbers for ASIMD and have at most one extra...
            assert(nsteps == 1);
            if (time_base == t) {
                nsteps = 0;
                return buffer[1];
            }
        }
        auto ts = rel_time(t, streambuffer) + vld1q_f64(streambuffer.time) * (double)step;
        time_base = t + step;
        nsteps = 1;
        auto func = (void (*)(double*, const double*, const void*))ramp_func;
        func(buffer, (const double*)&ts, data);
        return buffer[0];
    }
#else
    static void init()
    {
    }
    double compute(int64_t t, const void *data, int64_t, const StreamBuffer &streambuffer)
    {
        return compute_scalar(t, data, streambuffer);
    }
#endif
};

#if NACS_CPU_X86 || NACS_CPU_X86_64
double (BCGen::DataStream::*BCGen::DataStream::compute_vector)(
    int64_t, const void*, int64_t, const StreamBuffer&);
#endif

struct BCGen::TTLPulse {
    uint8_t chn;
    uint32_t id;
    uint32_t val;
    int64_t time; // in sequence unit before sorting and 10ns unit after sorting
};

struct BCGen::TTLManager {
    // All time in sequence time units

    // The time it takes to react to channel turning off (`~off_val` -> `off_val`)
    // This is the time the channel needs to be turned off before when it needed to be off.
    int64_t off_delay = 0;
    // The time it takes to react to channel turning on (0 -> 1)
    // This is the time the channel needs to be turned on before when it needed to be on.
    int64_t on_delay = 0;
    // Minimum off time. Off interval shorter than this will be skipped.
    int64_t skip_time = 0;
    // Minimum on time. On time shorter than this will be extended
    // (the off pulse will be delayed so that it's at least this time after the on one).
    int64_t min_time = 0;

    bool off_val = false;

    operator bool() const
    {
        return off_delay || on_delay || skip_time || min_time;
    }

    void clear() const
    {
        pulses.clear();
    }

    void process(const BCGen &seq, uint8_t chn) const
    {
        // Note that start time is in sequence unit and not FPGA unit at this point.
        std::sort(pulses.begin(), pulses.end(), [&] (auto &p1, auto &p2) {
            if (p1.time < p2.time)
                return true;
            if (p1.time > p2.time)
                return false;
            return p1.id < p2.id;
        });
        auto start_it = seq.m_real_start_vals.find({ChnType::TTL, chn});

        // The current pulse.
        // If `cur.val == output.val`, this pulse has already been outputted,
        // which may happen if an off pulse was skipped,
        // or before any actual output pulse was added.
        // If `cur.val != output.val`, the pulse is not outputted yet
        // (can be merged with the next one).
        // In this case, we have `output.time < cur.time`
        // and the code that modifies the current pulse should make sure
        // that the previous value is outputted when needed.

        // If `cur.time < 0`, this is a startvalue that we don't need to output
        // and will also be assumed to have happend sufficiently long before
        // that we don't need to worry about skip time/min time
        TTLPulse cur;
        if (start_it == seq.m_real_start_vals.end()) {
            cur.val = false;
        }
        else {
            cur.val = start_it->second != 0;
        }
        cur.time = -1;
        cur.id = 0;
        cur.chn = chn;
        TTLPulse output = cur;

        for (auto &pulse: pulses) {
            if (cur.val == pulse.val)
                continue;
            auto &t = pulse.time;
            if (pulse.val == off_val) {
                // If a turn off needs to be started at t <= 0,
                // generate a pulse at t = 0 instead of updating the start value.
                // This makes sure the channel is as programmed by the user,
                // at the beginning of the sequence
                t = max(t - off_delay, 0);
                if (cur.time >= 0) {
                    // If we are turning off, we might need to extend the turn off time.
                    t = max(t, cur.time + min_time);
                }
            }
            else {
                assert(cur.val == off_val);
                t = max(t - on_delay, 0);
                // If we are turning on, we might want to skip the previous one
                if (cur.time >= 0 && t <= cur.time + skip_time) {
                    // We only need to skip if the current pulse is an actual output
                    // (i.e. `cur.time >= 0`) in which case
                    // `cur.val == output.val` can only happen after an off pulse is skipped
                    // and there should be`cur.val != off_val` so it cannot happen.
                    assert(cur.val != output.val);
                    cur = output;
                    continue;
                }
                if (t <= 0 && seq.first_bseq) {
                    assert(cur.time < 0);
                    assert(output.time < 0);
                    seq.m_real_start_vals.insert_or_assign({ChnType::TTL, chn}, !off_val);
                    output = pulse;
                    output.time = -1;
                    cur = output;
                    continue;
                }
            }
            if (cur.val != output.val) {
                // Flush the output
                seq.m_ttlpulses.push_back(cur);
                output = cur;
            }
            cur = pulse;
        }
        if (cur.time >= 0 && cur.val != output.val) {
            seq.m_ttlpulses.push_back(cur);
        }
    }

    mutable std::vector<TTLPulse> pulses{};
};

bool BCGen::preprocess_ttl_pulse(const TTLPulse &ttl_pulse) const
{
    // Note that start time is in sequence unit and not FPGA unit at this point.
    if (ttl_pulse.chn >= m_ttl_managers.size())
        return false;
    auto &ttl_mgr = m_ttl_managers[ttl_pulse.chn];
    if (!ttl_mgr)
        return false;
    ttl_mgr.pulses.push_back(ttl_pulse);
    return true;
}

void BCGen::preprocess_ttl_managers() const
{
    // Note that start time is in sequence unit and not FPGA unit at this point.
    auto nmgrs = m_ttl_managers.size();
    for (uint8_t chn = 0; chn < nmgrs; chn++) {
        auto &ttl_mgr = m_ttl_managers[chn];
        if (!ttl_mgr || ttl_mgr.pulses.empty())
            continue;
        ttl_mgr.process(*this, chn);
    }
}

NACS_EXPORT_ BCGen::BCGen()
{
    static Call init([&] { DataStream::init(); });
}

NACS_EXPORT() BCGen::~BCGen()
{
}

NACS_EXPORT() uint32_t BCGen::convert_value(ChnType type, double value)
{
    switch (type) {
    case ChnType::TTL:
        return value != 0;
    case ChnType::Freq:
        return dds_freq_to_mu(value);
    case ChnType::Amp:
        return dds_amp_to_mu(value);
    case ChnType::Phase:
        return dds_phase_to_mu(value);
    case ChnType::DAC:
        return dac_to_mu(value);
    default:
        return 0;
    }
}

NACS_EXPORT() double BCGen::compute_start_val(PulseType pulse_type,
                                              void (*ramp_func)(void),
                                              const void *data)
{
    return DataStream::compute_start_val(pulse_type, ramp_func, data);
}

NACS_EXPORT() void BCGen::add_ttl_manager(uint8_t chn, int64_t off_delay, int64_t on_delay,
                                          int64_t skip_time, int64_t min_time, bool off_val)
{
    if (chn >= ttl_mask.size() * 32)
        throw std::runtime_error("TTL Manager on invalid channel.");
    if (chn >= m_ttl_managers.size())
        m_ttl_managers.resize(chn + 1);
    m_ttl_managers[chn] = TTLManager{off_delay, on_delay, skip_time, min_time, off_val};
}

NACS_EXPORT() uint32_t BCGen::version() const
{
    // DMA instructions version
    if (dma)
        return 0;
    // Bytecode version
    return ttl_mask.size() > 1 ? 3 : 2;
}

void BCGen::populate_pulses(const HostSeq &host_seq) const
{
    m_pulses.clear();
    m_ttlpulses.clear();
    for (auto &ttl_mgr: m_ttl_managers)
        ttl_mgr.clear();
    for (auto &seq_pulse: seq_pulses) {
        auto chn_type = seq_pulse.chn_type;
        auto ramp_func = seq_pulse.ramp_func;
        // We cannot translate the time unit yet since that may affect sorting
        // as well as the `time == 0` check.
        auto time = host_seq.get_time(seq_pulse.time_id);
        // Pulses at t=0 will be merged into the initial values.
        // For the first subsequence, this will be set before the sequence starts
        // so we can ignore those pulses.
        if (seq_pulse.cond_id != uint32_t(-1) &&
            host_seq.get_value(seq_pulse.cond_id, seq_idx) == 0)
            continue;
        auto endvalue = convert_value(chn_type,
                                      host_seq.get_value(seq_pulse.endvalue_id, seq_idx));
        if (chn_type == ChnType::TTL) {
            // Merge t=0 pulse into start value
            if (time == 0 && first_bseq) {
                m_real_start_vals.insert_or_assign({ChnType::TTL, seq_pulse.chn}, endvalue);
                continue;
            }
            // Keep TTL pulses in a separate array since we know we don't have any ramps
            // and we'd like to offer slightly better timing guarantee for these.
            TTLPulse ttl_pulse{ .chn = seq_pulse.chn, .id = seq_pulse.id,
                .val = endvalue != 0, .time = time };
            if (!preprocess_ttl_pulse(ttl_pulse))
                m_ttlpulses.push_back(ttl_pulse);
            continue;
        }
        // The bytecode generation code will automatically omit it if the value isn't changed.
        auto len = seq_pulse.len_id == uint32_t(-1) ? 0 :
            convert_time(host_seq.get_value(seq_pulse.len_id, seq_idx));
        if (len < 0)
            len = 0;
        if (len == 0)
            ramp_func = nullptr;
        // Merge t=0 pulse into start value
        if (time == 0 && first_bseq) {
            if (!ramp_func) {
                m_real_start_vals.insert_or_assign({chn_type, seq_pulse.chn}, endvalue);
                continue;
            }
            auto start_val = DataStream::compute_start_val(seq_pulse.pulse_type, ramp_func,
                                                           host_seq.values.data());
            m_real_start_vals.insert_or_assign({chn_type, seq_pulse.chn},
                                               convert_value(chn_type, start_val));
        }
        m_pulses.push_back({
                .pulse_type = seq_pulse.pulse_type,
                .chn_type = chn_type,
                .chn = seq_pulse.chn,
                .need_endvalue = true,
                .id = seq_pulse.id,
                .endvalue = endvalue,
                .time = time,
                .len = len,
                .ramp_func = ramp_func
            });
    }
    preprocess_ttl_managers();
    std::sort(m_pulses.begin(), m_pulses.end(), [&] (auto &p1, auto &p2) {
        if (p1.time < p2.time)
            return true;
        if (p1.time > p2.time)
            return false;
        return p1.id < p2.id;
    });
    for (auto &pulse: m_pulses) {
        assume(pulse.time >= 0);
        pulse.time = convert_time(pulse.time);
    }
    std::sort(m_ttlpulses.begin(), m_ttlpulses.end(), [&] (auto &p1, auto &p2) {
        if (p1.time < p2.time)
            return true;
        if (p1.time > p2.time)
            return false;
        return p1.id < p2.id;
    });
    for (auto &pulse: m_ttlpulses) {
        assume(pulse.time >= 0);
        pulse.time = convert_time(pulse.time);
    }
}

void BCGen::merge_pulses() const
{
    ChnMap<int64_t> next_times;
    next_times.fill(INT64_MAX);
    int32_t npulses = (int32_t)m_pulses.size();
    int32_t to = npulses - 1;
    int32_t from = npulses - 1;
    for (; from >= 0; from--) {
        auto &pulse = m_pulses[from];
        auto &time = next_times[{pulse.chn_type, pulse.chn}];
        assert(pulse.time <= time);
        if (pulse.time == time)
            continue;
        if (pulse.len >= time - pulse.time) {
            // Truncated ramp, shorten length and no need to output endvalue.
            pulse.len = time - pulse.time;
            pulse.need_endvalue = false;
        }
        time = pulse.time;
        if (from != to)
            m_pulses[to] = std::move(pulse);
        to--;
    }
    if (to >= 0) {
        m_pulses.erase(m_pulses.begin(), m_pulses.begin() + (to + 1));
    }
}

namespace {

class TimeKeeper {
    static constexpr int32_t avg_window = 1024; // Use a power of 2
    // This is roughly the density of TTL pulses we can support
    // when the DDS ramps are on:
    //
    //   We can push over roughly one command every 40 cycles
    //   and each DDS pulse takes 50 cycles. The command buffer size is 4096.
    //   Ignoring the length of the TTL pulse,
    //   we can fit roughly 4096 DDS pulses and 1024 TTL pulses
    //   within a timeframe of 4096 * 50 cycles.
    static constexpr int32_t min_time_span = avg_window * 200;

public:
    int64_t next_available(int64_t time) const
    {
        if (m_npulses < avg_window)
            return time;
        auto idx = m_npulses % avg_window;
        auto prev_idx = (m_npulses - 1) % avg_window;
        // The time span between the last `avg_window` TTL pulses.
        // Note that `idx` is the one we'll override and is the oldest time.
        auto prev_time = m_pulse_times[prev_idx];
        auto cur_span = prev_time - m_pulse_times[idx];
        assert(cur_span > 0);
        auto update_time = [&] (bool cond, int64_t limit) {
            if (cond && time < limit) {
                time = limit;
            }
        };
        update_time(cur_span < min_time_span / 2, prev_time + 250);
        update_time(cur_span < min_time_span, prev_time + 150);
        update_time(cur_span < min_time_span * 2, prev_time + 70);
        update_time(cur_span < min_time_span * 4, prev_time + 20);
        update_time(cur_span < min_time_span * 8, prev_time + 5);
        update_time(true, prev_time + 1);
        return time;
    }
    void add_time(int64_t time)
    {
        // Ignore
        if (m_npulses != 0) {
            auto prev_time = m_pulse_times[(m_npulses - 1) % avg_window];
            // The user might've ignored the time limit.
            if (time < prev_time)
                return;
            assert(time != prev_time);
        }
        m_pulse_times[m_npulses % avg_window] = time;
        m_npulses++;
    }

private:
    size_t m_npulses = 0;
    int64_t m_pulse_times[avg_window];
};

}

void BCGen::merge_ttl_pulses() const
{
    uint32_t ttl_vals[NUM_TTL_BANKS] = {};
    for (auto it = m_real_start_vals.begin(), end = m_real_start_vals.end(); it != end;) {
        auto [chn, val] = *it;
        if (chn.first != ChnType::TTL) {
            ++it;
            continue;
        }
        auto bank = chn.second / 32;
        uint8_t ttlchn = chn.second % 32;
        if (bank < NUM_TTL_BANKS)
            ttl_vals[bank] = setBit(ttl_vals[bank], ttlchn, val != 0);
        it = m_real_start_vals.erase(it);
    }
    for (int bank = 0; bank < min(int(ttl_mask.size()), NUM_TTL_BANKS); bank++) {
        // Skip ttl banks that we are not using
        if (ttl_mask[bank] == 0)
            continue;
        m_real_start_vals.emplace(std::make_pair(ChnType::TTL, bank * 32),
                                  ttl_vals[bank]);
    }

    TimeKeeper time_keeper;
    uint32_t nttlpulses = (uint32_t)m_ttlpulses.size();

    struct BankStatus {
        int64_t prev_t{-1};
        uint32_t last_idx{0};
    } bank_statuses[NUM_TTL_BANKS] = {};

    auto check_time_conflict = [&] (int64_t time) {
        for (auto status: bank_statuses) {
            if (status.prev_t == time) {
                return false;
            }
        }
        return true;
    };

    uint32_t to = 0;
    uint32_t from = 0;
    for (;from < nttlpulses; from++) {
        auto &pulse = m_ttlpulses[from];
        auto chn = pulse.chn;
        auto bank = chn / 32;
        auto &ttl_val = ttl_vals[bank];
        auto new_ttl_val = setBit(ttl_val, chn % 32, pulse.val != 0);
        if (new_ttl_val == ttl_val)
            continue;
        ttl_val = new_ttl_val;
        auto &bank_status = bank_statuses[bank];
        auto time = pulse.time;

        if (bank_status.prev_t >= time) {
            // Shouldn't trigger before we actually added a TTL pulse.
            assert(to != 0);
            m_ttlpulses[bank_status.last_idx].val = new_ttl_val;
            continue;
        }

        // Bump `to` since we'll now add a new pulse one way or another
        // The DMA can output all TTL banks in every cycle without any throughput limit.
        auto available_time = dma ? time : time_keeper.next_available(time);

        if (available_time > time) {
            if (bank_status.prev_t >= 0 && time < bank_status.prev_t + 1000) {
                // We'll respect the time limit if we have output on this bank
                // in the past 10us (arbitrary limit).
                time = available_time;
            }
            else {
                // We are ignoring the throughput limit since it's caused by other banks
                // but we still need to make sure we don't schedule any output
                // at exactly the same time as another bank.
                // We do so by maintaining the invariance that for each bank,
                // there's no output on this bank in the range
                // [current_time, bank_status.prev_t)
                while (check_time_conflict(time)) {
                    time++;
                }
            }
        }
        if (!dma)
            time_keeper.add_time(time);

        auto &out_pulse = m_ttlpulses[to];
        out_pulse.chn = uint8_t(bank * 32);
        if (from != to)
            out_pulse.id = pulse.id;
        out_pulse.val = new_ttl_val;
        out_pulse.time = time;
        bank_status.prev_t = time;
        bank_status.last_idx = to;
        to++;
    }
    m_ttlpulses.resize(to);
    std::sort(m_ttlpulses.begin(), m_ttlpulses.end(), [&] (auto &p1, auto &p2) {
        if (p1.time < p2.time)
            return true;
        if (p1.time > p2.time)
            return false;
        return p1.id < p2.id;
    });
}

NACS_EXPORT() void BCGen::generate(const HostSeq &host_seq) const
{
    m_real_start_vals = start_vals;
    populate_pulses(host_seq);
    merge_pulses();
    merge_ttl_pulses();
    if (dma) {
        emit_dma(host_seq.values.data());
    }
    else {
        emit_bytecode(host_seq.values.data());
    }
}

static uint8_t pulse_mintime(ChnType chn_typ)
{
    if (chn_typ == ChnType::DAC)
        return PulseTime::DAC;
    if (chn_typ == ChnType::Freq)
        return PulseTime::DDSFreq;
    if (chn_typ == ChnType::Amp)
        return PulseTime::DDSAmp;
    if (chn_typ == ChnType::Phase)
        return PulseTime::DDSPhase;
    if (chn_typ == ChnType::Clock)
        return PulseTime::Clock;
    return PulseTime::Min2;
}

template<typename Inst>
struct BCGen::Writer {
    struct DDS {
        uint32_t freq = 0;
        uint16_t amp = 0;
        uint16_t phase = 0;
        bool freq_set = false;
        bool amp_set = false;
        bool phase_set = false;
    };
    struct DAC {
        uint16_t V = 0;
        bool set = false;
    };

    Writer(std::vector<uint8_t> &bytecode, uint8_t start_chn)
        : start_ttl_chn(start_chn),
          bytecode(bytecode)
    {
        if (Inst::version < 3) {
            assert(start_ttl_chn < 32);
        }
    }

    template<typename T>
    uint32_t write(T obj)
    {
        static_assert(std::is_trivial_v<T>);
        auto len = (uint32_t)bytecode.size();
        bytecode.resize(len + sizeof(obj));
        memcpy(&bytecode[len], &obj, sizeof(obj));
        return len;
    }

    template<typename T>
    uint32_t add_inst(T inst)
    {
        max_time_left = 0;
        return write<T>(inst);
    }

    void write_header(int64_t len_ns, const uint32_t *ttl_masks, uint32_t nmasks)
    {
        // Version number is in a different ZMQ message and won't be included here.
        // See `BCGen::version()`.
        write(len_ns);
        if (Inst::version < 3) {
            assert(nmasks == 1);
            write(setBit(*ttl_masks, start_ttl_chn, true));
        }
        else {
            assert(nmasks <= NUM_TTL_BANKS);
            assert(start_ttl_chn < nmasks * 32);
            write(nmasks);
            uint8_t start_ttl_bank = start_ttl_chn / 32;
            for (uint32_t bank = 0; bank < nmasks; bank++) {
                auto ttl_mask = ttl_masks[bank];
                if (bank == start_ttl_bank)
                    ttl_mask = setBit(ttl_mask, start_ttl_chn % 32, true);
                write(ttl_mask);
            }
        }
    }

    // States
    uint32_t cur_ttl[NUM_TTL_BANKS] = {0};
    bool ttl_set[NUM_TTL_BANKS] = {false};
    DDS dds[22] = {};
    DAC dac[4] = {};

    uint8_t start_ttl_chn;
    std::vector<uint8_t> &bytecode;

    // This is the time of the end of the last pulse.
    // We use the difference between this time and the time of the pulse
    // to decide how much wait time we need to add.
    int64_t cur_t = 0;
    uint8_t max_time_left = 0;
    ssize_t last_timed_inst = 0;

    // Increase the wait time encoded in the last instruction by `t`.
    // The caller should have checked that the time fits in the maximum time possible to be
    // encoded.
    void inc_last_time(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t b = bytecode[last_timed_inst];
        uint8_t op = b & 0xf;
        uint8_t tmask = 0xf0;
        if constexpr (Inst::version >= 3) {
            switch (op) {
            case ByteCode::TTLAll:
                tmask = 0x10;
                break;
            case ByteCode::TTL1_v3:
                break;
            case ByteCode::TTL3_v3:
                tmask = 0x30;
                break;
            default:
                abort();
            }
        }
        else {
            switch (op) {
            case ByteCode::TTLAll:
                break;
            case ByteCode::TTL2_v1:
                tmask = 0x30;
                break;
            case ByteCode::TTL5:
                tmask = 0x70;
                break;
            default:
                abort();
            }
        }
        t = uint8_t(t + ((b & tmask) >> 4));
        bytecode[last_timed_inst] = uint8_t((b & ~tmask) | ((t << 4) & tmask));
    }

    void add_wait_trigger(uint8_t chn, bool trig_raise, uint32_t timeout)
    {
        using namespace ByteCode;
        assert(timeout <= 0xffffff);
        add_inst(typename Inst::WaitTrigger{OpCode::WaitTrigger, 0, trig_raise,
                chn, timeout & 0xffffff});
    }

    // The `add***` functions below provides an abstraction to the bytecode and hides
    // the detail about the bytecode encoding. The code tries to use the most compact
    // encoding when available and can fallback to generic encoding if not.
    // This is the same level as the API of `ExeState`.
    // Note that the user of the writer should generally not call `add_wait` directly
    // but simply call the functions below to add pulses at specific times.
    void add_wait(int64_t dt)
    {
        using namespace ByteCode;
        assert(dt >= 0);
        if (dt == 0)
            return;
        cur_t += dt;
        if (dt <= max_time_left) {
            inc_last_time(uint8_t(dt));
            return;
        }
        uint16_t times[17];
        memset(times, 0, sizeof(times));
        for (int i = 15; i >= 0; i--) {
            if (dt <= uint16_t(2047 + max_time_left)) {
                assert(dt > max_time_left);
                if (max_time_left) {
                    dt -= max_time_left;
                    inc_last_time(max_time_left);
                }
                times[16] = uint16_t(dt);
                dt = 0;
                break;
            }
            if (i > 0) {
                // The first number not representable by the next one is
                // 0x10000 * 2^(3 * i - 3)
                // or
                // 0x2000 << (3 * i)
                int64_t tnext_max = int64_t(0x2000) << (3 * i);
                if (dt < tnext_max) {
                    continue;
                }
            }
            times[i] = uint16_t(dt >> (3 * i));
            dt -= int64_t(times[i]) << (3 * i);
            if (dt == 0) {
                break;
            }
            else if (dt <= max_time_left) {
                inc_last_time(uint8_t(dt));
                break;
            }
        }
        for (int i = 15; i >= 0; i--) {
            if (!times[i])
                continue;
            add_inst(typename Inst::Wait{OpCode::Wait, uint8_t(i & 0xf), times[i]});
        }
        if (times[16]) {
            add_inst(typename Inst::Wait2{OpCode::Wait2, 1, uint16_t(times[16] & 2047)});
        }
    }

    int add_ttl(int64_t t, uint32_t ttl, int bank)
    {
        using namespace ByteCode;
        if (ttl == cur_ttl[bank] && ttl_set[bank])
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::Min2;
        auto changes = ttl ^ cur_ttl[bank];
        if (!ttl_set[bank]) {
            ttl_set[bank] = true;
            changes = ttl;
        }
        else {
            assert(changes);
        }
        cur_ttl[bank] = ttl;
        auto nchgs = __builtin_popcountll(changes);
        // The case where nchgs == 0, i.e. initial setting to 0, will be handled by
        // the generic case below.
        if (nchgs == 1) {
            auto bit = __builtin_ffs(changes) - 1;
            if constexpr (Inst::version >= 3) {
                last_timed_inst = add_inst(typename Inst::TTL1{OpCode::TTL1_v3, 0,
                        uint8_t(bank & 7), uint8_t(bit & 0x1f)});
                max_time_left = 0xf;
            }
            else {
                last_timed_inst = add_inst(typename Inst::TTL2{OpCode::TTL2_v1, 0,
                        uint8_t(bit & 0x1f), uint8_t(bit & 0x1f)});
                max_time_left = 3;
            }
        }
        else if (nchgs == 2) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            if constexpr (Inst::version >= 3) {
                last_timed_inst = add_inst(typename Inst::TTL3{OpCode::TTL3_v3, 0,
                        uint8_t(bank & 7), uint8_t(bit1 & 0x1f), uint8_t(bit2 & 0x1f),
                        uint8_t(bit2 & 0x1f)});
                max_time_left = 3;
            }
            else {
                last_timed_inst = add_inst(typename Inst::TTL2{OpCode::TTL2_v1, 0,
                        uint8_t(bit1 & 0x1f), uint8_t(bit2 & 0x1f)});
                max_time_left = 3;
            }
        }
        else if (nchgs == 3) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            if constexpr (Inst::version >= 3) {
                last_timed_inst = add_inst(typename Inst::TTL3{OpCode::TTL3_v3, 0,
                        uint8_t(bank & 7), uint8_t(bit1 & 0x1f), uint8_t(bit2 & 0x1f),
                        uint8_t(bit3 & 0x1f)});
                max_time_left = 3;
            }
            else {
                add_inst(typename Inst::TTL4{OpCode::TTL4_v1, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit3 & 0x1f)});
            }
        }
        else if (nchgs == 4) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit3);
            auto bit4 = __builtin_ffs(changes) - 1;
            if constexpr (Inst::version >= 3) {
                add_inst(typename Inst::TTL5{OpCode::TTL5, uint8_t(bank & 7),
                        uint8_t(bit1 & 0x1f), uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit4 & 0x1f)});
            }
            else {
                add_inst(typename Inst::TTL4{OpCode::TTL4_v1, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f), uint8_t(bit4 & 0x1f)});
            }
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
            if constexpr (Inst::version >= 3) {
                add_inst(typename Inst::TTL5{OpCode::TTL5, uint8_t(bank & 7),
                        uint8_t(bit1 & 0x1f), uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit4 & 0x1f)});
            }
            else {
                last_timed_inst = add_inst(typename Inst::TTL5{OpCode::TTL5, 0, uint8_t(bit1 & 0x1f),
                        uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                        uint8_t(bit4 & 0x1f), uint8_t(bit5 & 0x1f)});
                max_time_left = 7;
            }
        }
        else {
            if constexpr (Inst::version >= 3) {
                last_timed_inst = add_inst(typename Inst::TTLAll{OpCode::TTLAll, 0,
                        uint8_t(bank & 7), ttl});
                max_time_left = 1;
            }
            else {
                last_timed_inst = add_inst(typename Inst::TTLAll{OpCode::TTLAll, 0, ttl});
                max_time_left = 15;
            }
        }
        return PulseTime::Min2;
    }

    int add_clock(int64_t t, uint8_t period)
    {
        using namespace ByteCode;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::Clock;
        add_inst(typename Inst::Clock{OpCode::Clock, 0, uint8_t(period - 1)});
        return PulseTime::Clock;
    }

    int add_freq(int64_t t, uint8_t chn, uint32_t freq)
    {
        using namespace ByteCode;
        if (dds[chn].freq_set && freq == dds[chn].freq)
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::DDSFreq;
        uint32_t dfreq = freq - dds[chn].freq;
        dds[chn].freq = freq;
        dds[chn].freq_set = true;
        if (dfreq <= 0x3f || dfreq >= 0xffffffc0) {
            add_inst(typename Inst::DDSDetFreq2{OpCode::DDSDetFreq2, uint8_t(chn & 0x1f),
                    uint8_t(dfreq & 0x7f)});
        }
        else if (dfreq <= 0x3fff || dfreq >= 0xffffc000) {
            add_inst(typename Inst::DDSDetFreq3{OpCode::DDSDetFreq3, uint8_t(chn & 0x1f),
                    uint16_t(dfreq & 0x7fff)});
        }
        else if (dfreq <= 0x003fffff || dfreq >= 0xffc00000) {
            add_inst(typename Inst::DDSDetFreq4{OpCode::DDSDetFreq4, uint8_t(chn & 0x1f),
                    uint32_t(dfreq & 0x7fffff)});
        }
        else {
            add_inst(typename Inst::DDSFreq{OpCode::DDSFreq, uint8_t(chn & 0x1f),
                    uint32_t(freq & 0x7fffffff)});
        }
        return PulseTime::DDSFreq;
    }

    int add_phase(int64_t t, uint8_t chn, uint16_t phase)
    {
        using namespace ByteCode;
        if (dds[chn].phase_set && phase == dds[chn].phase)
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::DDSPhase;
        dds[chn].phase = phase;
        dds[chn].phase_set = true;
        add_inst(typename Inst::DDSPhase{OpCode::DDSPhase, uint8_t(chn & 0x1f), 0, phase});
        return PulseTime::DDSPhase;
    }

    int add_amp(int64_t t, uint8_t chn, uint16_t amp)
    {
        using namespace ByteCode;
        if (dds[chn].amp_set && amp == dds[chn].amp)
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::DDSAmp;
        uint16_t damp = uint16_t(amp - dds[chn].amp);
        dds[chn].amp = amp;
        dds[chn].amp_set = true;
        if (damp <= 0x3f || damp >= 0xffc0) {
            add_inst(typename Inst::DDSDetAmp{OpCode::DDSDetAmp, uint8_t(chn & 0x1f),
                    uint8_t(damp & 0x7f)});
        }
        else {
            add_inst(typename Inst::DDSAmp{OpCode::DDSAmp, 0, uint8_t(chn & 0x1f),
                    uint16_t(amp & 0xfff)});
        }
        return PulseTime::DDSAmp;
    }

    int add_dac(int64_t t, uint8_t chn, uint16_t V)
    {
        using namespace ByteCode;
        if (dac[chn].set && V == dac[chn].V)
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::DAC;
        uint16_t dV = uint16_t(V - dac[chn].V);
        dac[chn].V = V;
        dac[chn].set = true;
        if (dV <= 0x1ff || dV >= 0xfe00) {
            add_inst(typename Inst::DACDet{OpCode::DACDet, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            add_inst(typename Inst::DAC{OpCode::DAC, 0, uint8_t(chn & 0x3), V});
        }
        return PulseTime::DAC;
    }

    int add_pulse(ChnType chn_typ, uint8_t chn, uint32_t val, int64_t t)
    {
        switch (chn_typ) {
        case ChnType::TTL:
            return add_ttl(t, val, chn);
        case ChnType::Clock:
            return add_clock(t, uint8_t(val));
        case ChnType::Freq:
            return add_freq(t, chn, val);
        case ChnType::Amp:
            return add_amp(t, chn, uint16_t(val));
        case ChnType::Phase:
            return add_phase(t, chn, uint16_t(val));
        case ChnType::DAC:
            return add_dac(t, chn, uint16_t(val));
        default:
            assert(chn_typ != ChnType::TTL);
            return 0;
        }
    }
};

namespace {

struct SinglePulse {
    ChnType typ;
    uint8_t chn;
    uint32_t val;
    uint8_t mintime() const
    {
        return pulse_mintime(typ);
    }
    static SinglePulse clock(uint8_t period)
    {
        return { ChnType::Clock, 0, period };
    }
    static SinglePulse ttl(uint32_t val, uint8_t bank)
    {
        return { ChnType::TTL, bank, val };
    }
};

struct PulseSequence {
    // Map from start time to pulse object
    std::map<int64_t,SinglePulse> single_pulses;
    auto add(int64_t t, SinglePulse pulse, int64_t ub=INT64_MAX)
    {
        auto it = single_pulses.upper_bound(t); // it->first > t
        if (it != single_pulses.begin()) {
            auto prev_it = it;
            --prev_it; // prev_it->first <= t
            auto prev_end = prev_it->first + prev_it->second.mintime();
            if (prev_end > ub)
                return single_pulses.end();
            t = std::max(t, prev_end);
        }
        auto mintime = pulse.mintime();
        while (it != single_pulses.end() && it->first < t + mintime) {
            t = it->first + it->second.mintime();
            if (t > ub)
                return single_pulses.end();
            ++it;
        }
        return single_pulses.emplace(t, pulse).first;
    }
    // Find the first empty slot or the first pulse that is on the same channel
    // as `pulse` no earlier than `t`.
    auto add_endval(int64_t t, SinglePulse pulse)
    {
        auto it = single_pulses.lower_bound(t); // it->first >= t
        if (it != single_pulses.begin()) {
            auto prev_it = it;
            --prev_it; // prev_it->first < t
            auto prev_end = prev_it->first + prev_it->second.mintime();
            t = std::max(t, prev_end);
        }
        auto mintime = pulse.mintime();
        while (true) {
            if (it == single_pulses.end() || it->first >= t + mintime)
                return single_pulses.emplace(t, pulse).first;
            if (it->second.typ == pulse.typ && it->second.chn == pulse.chn) {
                // Found a new pulse on the same channel,
                // no need to output the end value anymore.
                return it;
            }
            t = it->first + it->second.mintime();
            ++it;
        }
    }
};

}

void BCGen::emit_bytecode(const void *data) const
{
    bytecode.clear();
    if (version() >= 3) {
        _emit_bytecode<ByteCode::Inst_v3>(data);
    }
    else {
        _emit_bytecode<ByteCode::Inst_v1>(data);
    }
}

template<typename Inst>
void BCGen::_emit_bytecode(const void *data) const
{
    Writer<Inst> writer(bytecode, start_ttl_chn);
    writer.write_header(len_ns, &ttl_mask[0], uint32_t(ttl_mask.size()));

    for (int bank = 0; bank < NUM_TTL_BANKS; bank++) {
        auto ttl_start_it = m_real_start_vals.find({ChnType::TTL, bank * 32});
        if (ttl_start_it != m_real_start_vals.end()) {
            writer.cur_ttl[bank] = ttl_start_it->second;
        }
    }
    // Initialize channels
    bool has_init = false;
    ChnMap<uint32_t> chn_values;
    chn_values.fill(0);
    for (auto [chn, val]: m_real_start_vals) {
        if (chn.first == ChnType::TTL) {
            uint8_t bank = chn.second / 32;
            if (first_bseq) {
                has_init = true;
                writer.add_ttl(writer.cur_t, val, bank);
            }
            continue;
        }
        chn_values[chn] = val;
        if (first_bseq) {
            has_init = true;
            auto min_dt = writer.add_pulse(chn.first, chn.second, val, writer.cur_t);
            assert(min_dt >= 0);
            (void)min_dt;
        }
    }
    // Add some time margin after channel initialization
    // (which should only be non-empty for first bseq)
    if (has_init) {
        assert(first_bseq);
        writer.add_wait(2500);
    }
    if (trig_type != TrigType::None)
        writer.add_wait_trigger(trig_channel, trig_type == TrigType::Raise,
                                trig_timeout_cycle);

    PulseSequence pseq;
    // Note that this will also make sure the initial ttl value is set.
    auto set_start_ttl = [&] (int64_t t, bool val) {
        auto start_ttl_bank = uint8_t(start_ttl_chn / 32);
        auto base = writer.cur_ttl[start_ttl_bank];
        return pseq.add(t, SinglePulse::ttl(setBit(base, start_ttl_chn % 32, val),
                                            start_ttl_bank));
    };
    set_start_ttl(-int32_t(seq_delay), false);

    // 1. Clock out: these will always be outputted at the right time with the highest priority
    for (auto clock: clocks)
        pseq.add(convert_time(clock.time), SinglePulse::clock(clock.period));
    // Add the start_ttl on after the clock since we don't care about it's timing
    // that much.
    set_start_ttl(-int32_t(seq_delay) - 100, true);
    // 2. TTL out: after we make sure we can do the clock correct,
    //    we will try to respect ttl output time as much as possible.
    for (auto &pulse: m_ttlpulses)
        pseq.add(pulse.time, SinglePulse::ttl(pulse.val, pulse.chn / 32));

    // 3. DDS/DAC: output at the correct time if we can
    ChnMap<int64_t> last_times;
    // This is a value that is more negative than the largest mintime
    last_times.fill(-1000);
    std::vector<Pulse*> ramps;
    const StreamBuffer streambuffer(fpga_clock_div);
    ChnMap<DataStream> streams;
    for (auto &pulse: m_pulses) {
        auto &last_time = last_times[{pulse.chn_type, pulse.chn}];
        auto &chn_value = chn_values[{pulse.chn_type, pulse.chn}];
        auto &stream = streams[{pulse.chn_type, pulse.chn}];

        auto mintime = pulse_mintime(pulse.chn_type);
        auto last_finish_time = last_time + mintime;
        uint32_t new_value = pulse.endvalue;
        // Initialize stream to be used later
        if (pulse.ramp_func)
            stream.add_pulse(pulse);

        if (pulse.time <= last_finish_time) {
            // The last output on this channel is already in the future.
            // Simply update it to the value specified by the new pulse.
            auto pulse_it = pseq.single_pulses.find(last_time);
            assert(pulse_it != pseq.single_pulses.end());
            assert(pulse_it->second.typ == pulse.chn_type);
            assert(pulse_it->second.chn == pulse.chn);
            // The last output time is ahead of the pulse start time but may or may not
            // be ahead of the pulse end time. We might still need to schedule the ramp
            // for this pulse if there is one.
            if (pulse.ramp_func && pulse.time + pulse.len > last_finish_time) {
                auto t = std::max(last_time, pulse.time);
                auto start_val = stream.compute(t, data, 0, streambuffer);
                new_value = convert_value(pulse.chn_type, start_val);
                ramps.push_back(&pulse);
            }
            if (new_value != chn_value) {
                pulse_it->second.val = new_value;
                chn_value = new_value;
            }
            continue;
        }
        // First figure out the time we can do an output before we decide which value
        // to use.
        auto new_pulse_it = pseq.add(pulse.time, { pulse.chn_type, pulse.chn, new_value });
        auto new_time = new_pulse_it->first;
        if (pulse.ramp_func && new_time + mintime < pulse.time + pulse.len) {
            // If the output time isn't exactly at the pulse start time,
            // compute the actual value for the output time.
            auto start_val = stream.compute(new_time, data, 0, streambuffer);
            new_value = convert_value(pulse.chn_type, start_val);
            new_pulse_it->second.val = new_value;
            ramps.push_back(&pulse);
        }
        last_time = new_time;
        chn_value = new_value;
    }

    struct RampState {
        Pulse *pulse;
        decltype(pseq.single_pulses.begin()) last_it;
    };

    ChnMap<RampState> ramp_states;
    ramp_states.fill(RampState{ nullptr, pseq.single_pulses.end() });
    auto nramps = (int)ramps.size();
    auto ramp_idx = 0;
    int64_t cur_time = 0;
    int cur_chn = 0;
    auto last_progress_chn = ramp_states.chn_num - 1;
    auto next_chn = [&] { return cur_chn == ramp_states.chn_num - 1 ? 0 : cur_chn + 1; };
    int ramp_step = 0;
    while (true) {
        if (ramp_step == 0) {
            // No on-going ramp, jump ahead or exit
            if (ramp_idx >= nramps)
                break;
            cur_time = ramps[ramp_idx]->time;
        }
        for (; ramp_idx < nramps; ramp_idx++) {
            auto pulse = ramps[ramp_idx];
            if (pulse->time > cur_time)
                break;
            auto &ramp_state = ramp_states[{pulse->chn_type, pulse->chn}];
            if (!ramp_state.pulse)
                ramp_step += pulse_mintime(pulse->chn_type);
            ramp_state = { pulse, pseq.single_pulses.end() };
            streams[{pulse->chn_type, pulse->chn}].add_pulse(*pulse);
        }
        assert(ramp_step);
#ifndef NDEBUG
        {
            auto has_ramp = false;
            for (auto &state: ramp_states)
                has_ramp |= state.pulse != nullptr;
            assert(has_ramp);
        }
#endif
        // Find the next channel
        while (!ramp_states[cur_chn].pulse) {
            if (cur_chn == last_progress_chn)
                cur_time += 50;
            cur_chn = next_chn();
        }
        auto &ramp_state = ramp_states[cur_chn];
        auto pulse = ramp_state.pulse;
        auto endtime = pulse->time + pulse->len;
        auto mintime = pulse_mintime(pulse->chn_type);
        if (endtime > cur_time) {
            auto new_pulse_it = pseq.add(cur_time,
                                         { pulse->chn_type, pulse->chn, 0 },
                                         endtime);
            if (new_pulse_it != pseq.single_pulses.end()) {
                // Found a new slot for the ramp output
                auto new_time = new_pulse_it->first;
                auto val = streams[cur_chn].compute(new_time, data, ramp_step,
                                                    streambuffer);
                new_pulse_it->second.val = convert_value(pulse->chn_type, val);
                if (ramp_state.last_it != pseq.single_pulses.end() &&
                    ramp_state.last_it->second.val == new_pulse_it->second.val) {
                    // Repeated value, skipping
                    if (cur_chn == last_progress_chn)
                        cur_time += 50;
                    pseq.single_pulses.erase(new_pulse_it);
                }
                else {
                    last_progress_chn = cur_chn;
                    ramp_state.last_it = new_pulse_it;
                    cur_time = new_time + mintime;
                }
                cur_chn = next_chn();
                continue;
            }
        }
        last_progress_chn = cur_chn;
        // No more slot available
        if (pulse->need_endvalue) {
            if (ramp_state.last_it != pseq.single_pulses.end()) {
                ramp_state.last_it->second.val = pulse->endvalue;
            }
            else {
                pseq.add_endval(std::min(endtime, cur_time),
                                { pulse->chn_type, pulse->chn, pulse->endvalue });
            }
        }
        ramp_state.pulse = nullptr;
        ramp_step -= mintime;
        cur_chn = next_chn();
    }
    assert(!pseq.single_pulses.empty());
    // Shift cur_t to match where we schedule the first pulse to be.
    writer.cur_t = pseq.single_pulses.begin()->first;
    for (auto [t, pulse]: pseq.single_pulses)
        writer.add_pulse(pulse.typ, pulse.chn, pulse.val, t);

    // Wait for 1us at the end of the sequence.
    // This just adds a little spacing before the termination sequence.
    writer.add_wait(std::max(int64_t(0), -writer.cur_t) + 100);
}

/**
 * DMA instructions generation (version 0, see `DMA::Inst_v0`).
 *
 * The hardware executes all the output instructions before a wait in the same cycle
 * (a "group") with at most one output on each of the DDS buses, the DAC and the clock
 * in each group. A DDS or DAC output keeps its bus busy for a while and a new output
 * on a busy bus would corrupt the previous one.
 *
 * Unlike the bytecode, which executes one output at a time, the outputs on different
 * devices are scheduled independently and run concurrently.
 * * The TTL and clock outputs are done at exactly the requested time.
 * * Each of the DDS buses and the DAC are scheduled independently (`DMABus`).
 *   The DDS writes are started early so that the frequency update (FUD),
 *   which makes the new values take effect, happens at the requested time.
 *   Multiple DDS parameters (of any DDS's on the same bus) updated at the same time
 *   are written together with only the last write triggering the FUD.
 *
 * All the time in this part is in DMA cycles (5ns) unless otherwise noted.
 */

namespace {

// Each 10ns (the time unit of the bytecode) is 2 DMA cycles.
constexpr int64_t dma_cycle_per_10ns = 2;
// Extra cycles between two DAC outputs (safety margin).
// No margin is needed for the DDS since a new write can start right when the previous
// write finishes (the time reported by the server), the same as the two halves
// of a 32-bit write.
constexpr int64_t dma_dac_margin = 1;
// SPI clock divider for the DAC (DMA cycles per SPI clock edge minus 1).
// This is the same SPI clock as the DAC output of the bytecode.
constexpr uint16_t dma_dac_div = 1;
// Minimum time between the starts of two DAC outputs
// (18 bits with two clock edges per bit, plus the chip select).
constexpr int64_t dma_dac_time = 37 * (dma_dac_div + 1) + 2;
// The hardware runs the commands from the CPU (e.g. DDS writes) in waits longer than this.
// Writes without FUD must not be separated by a wait this long
// so that no FUD from the CPU can happen before the FUD of the DMA.
constexpr int64_t dma_long_wait = 128;
// Wait at the start of the sequence so that the hardware can buffer some instructions.
constexpr int64_t dma_start_wait = 200; // 1us
// Wait after initializing the channels (first basic sequence only), same as the bytecode.
constexpr int64_t dma_init_wait = 2500 * dma_cycle_per_10ns; // 25us
// Delay of the FUD output of the DDS writes relative to the TTL outputs
// done in the same group. Assumed to be 0.
constexpr int64_t dma_dds_fud_delay = 0;

// Convert the clock period (in 10ns, 0 for off, same as the input of the bytecode)
// to the period of the DMA clock output instruction.
// The DMA clock output toggles every `period + 1` cycles and is off for 0x1ff.
static uint16_t dma_clock_period(uint8_t period)
{
    auto p = uint8_t(period - 1);
    return uint16_t((p << 1) | 1);
}

// Which bus the DDS/DAC channel is on (0 and 1 for the DDS buses and 2 for the DAC),
// -1 for other channels.
static int dma_bus_id(ChnType type, uint8_t chn)
{
    switch (type) {
    case ChnType::Freq:
    case ChnType::Amp:
    case ChnType::Phase:
        assert(chn < 22);
        // DDS 0-10 are on bus 0 and 11-21 are on bus 1.
        return chn >= 11 ? 1 : 0;
    case ChnType::DAC:
        assert(chn < 4);
        return 2;
    default:
        return -1;
    }
}

struct DMADDSWrite {
    uint8_t dds_id; // DDS number on the bus
    uint8_t addr; // 16 bit register index (byte address / 2)
    bool is32;
    bool fud;
    uint32_t data;
};

// All the outputs at the same time.
struct DMAGroup {
    // TTL values of the banks in `ttl_set`
    uint32_t ttl[NUM_TTL_BANKS];
    uint8_t ttl_set = 0;
    static_assert(NUM_TTL_BANKS <= 8);
    bool has_clock = false;
    uint16_t clock_period = 0;
    bool has_dds[2] = {false, false};
    DMADDSWrite dds[2];
    bool has_dac = false;
    uint32_t dac_data = 0;
};

using DMAGroups = std::map<int64_t,DMAGroup>;

}

// Schedule the outputs on a DDS bus or the DAC.
//
// The outputs are done in batches. A batch is a series of back-to-back writes on the bus
// for any number of DDS's (DAC channels) on the bus. For the DDS's, only the last write
// of the batch triggers the frequency update (FUD), which is shared by all the DDS's
// on the bus, so that all the values written in the batch take effect at the same time,
// i.e. at the (rising) edge of the FUD. We refer to this time as the time of the batch.
// For the DAC (no FUD), the time of the batch is the start of the first write.
//
// * The values that must change at a specific time (pulses, start of ramps
//   and end values of ramps) are written in a batch at exactly that time if possible
//   (i.e. if the bus is available early enough) or as soon as possible otherwise.
// * The ramps are sampled (at the time of the batch) as fast as possible in the time
//   between these batches without delaying them, with all the ramps on the bus
//   in the same batch if possible.
// * Only the 16-bit registers that are changed are written. Two registers in the same
//   32-bit pair (frequency, or phase and amplitude) are written with a single 32-bit write.
struct BCGen::DMABus {
    struct Slot {
        Pulse *ramp = nullptr;
        DataStream stream{};
        // Value to be written in the next batch
        bool has_pending = false;
        uint32_t pending = 0;
    };
    // A DDS or a DAC channel
    struct Unit {
        // Last value written to each 16-bit register
        // (frequency low, frequency high, phase, amplitude for the DDS).
        bool has_cur[4] = {false, false, false, false};
        uint16_t cur[4] = {0, 0, 0, 0};
    };
    // A write to the hardware.
    struct Write {
        uint8_t unit;
        uint8_t reg; // Index of the first 16-bit register (see `Unit::cur`)
        bool is32;
        bool fud;
        uint32_t data;
    };
    // A change at a specific time.
    struct Event {
        int64_t time;
        unsigned slot;
        // Start of `pulse` or the end of the ramp `pulse` (setting the end value).
        bool ramp_end;
        Pulse *pulse;
    };

    DMABus(int bus_id, const DDSTiming &timing)
        : bus_id(bus_id),
          timing(timing)
    {
        if (is_dac()) {
            unit_nslots = 1;
            units.resize(4);
            retry_time = dma_dac_time + dma_dac_margin;
        }
        else {
            // A unit is a DDS with 3 slots (freq, amp, phase).
            unit_nslots = 3;
            units.resize(11);
            retry_time = timing.set16_fud;
        }
        slots.resize(units.size() * unit_nslots);
    }

    bool is_dac() const
    {
        return bus_id == 2;
    }
    unsigned slot_index(const Pulse &pulse) const
    {
        if (is_dac())
            return pulse.chn;
        unsigned param = (pulse.chn_type == ChnType::Freq ? 0 :
                          pulse.chn_type == ChnType::Amp ? 1 : 2);
        return (pulse.chn - bus_id * 11) * 3 + param;
    }
    // The registers of the unit set by the slot.
    uint8_t slot_regs(unsigned slot) const
    {
        if (is_dac())
            return 1;
        switch (slot % 3) {
        case 0: return 0x3; // frequency
        case 1: return 0x8; // amplitude
        default: return 0x4; // phase
        }
    }
    // The time from the start of the write until the next write on the bus can start.
    int64_t write_time(const Write &write) const
    {
        if (is_dac())
            return dma_dac_time + dma_dac_margin;
        if (write.is32)
            return write.fud ? timing.set32_fud : timing.set32;
        return write.fud ? timing.set16_fud : timing.set16;
    }
    // Total time of the writes.
    int64_t total_time(const std::vector<Write> &writes) const
    {
        int64_t res = 0;
        for (auto &write: writes)
            res += write_time(write);
        return res;
    }
    // The time from the start of the batch to the time the outputs take effect.
    int64_t lead_time(const std::vector<Write> &writes) const
    {
        if (writes.empty() || is_dac())
            return 0;
        return total_time(writes) - timing.fud_to_end + dma_dds_fud_delay;
    }

    // Add the writes for `unit` to `writes`, writing the registers in `regs` (bit mask).
    void add_writes(std::vector<Write> &writes, unsigned unit, uint8_t regs,
                    const uint16_t *vals) const
    {
        if (is_dac()) {
            if (regs & 1) {
                writes.push_back({ .unit = uint8_t(unit), .reg = 0, .is32 = false,
                        .fud = false, .data = (uint32_t(unit) << 16) | vals[0] });
            }
            return;
        }
        for (uint8_t reg = 0; reg < 4; reg += 2) {
            auto pair = (regs >> reg) & 3;
            if (pair == 3) {
                writes.push_back({ .unit = uint8_t(unit), .reg = reg, .is32 = true,
                        .fud = false,
                        .data = uint32_t(vals[reg]) | (uint32_t(vals[reg + 1]) << 16) });
            }
            else if (pair) {
                auto r = uint8_t(reg + (pair == 2));
                writes.push_back({ .unit = uint8_t(unit), .reg = r, .is32 = false,
                        .fud = false, .data = vals[r] });
            }
        }
    }
    void finish_writes(std::vector<Write> &writes) const
    {
        // Only the last write of the batch triggers the FUD (shared by all the DDS's
        // on the bus) so that all the values are updated at the same time.
        if (!is_dac() && !writes.empty()) {
            writes.back().fud = true;
        }
    }

    // The upper bound of the lead time of the batch writing the slots in `sel`.
    int64_t max_lead_time(const std::vector<unsigned> &sel)
    {
        if (is_dac())
            return 0;
        auto &writes = scratch_writes;
        writes.clear();
        const uint16_t vals[4] = {};
        for (size_t i = 0; i < sel.size();) {
            auto unit = sel[i] / unit_nslots;
            uint8_t regs = 0;
            for (; i < sel.size() && sel[i] / unit_nslots == unit; i++)
                regs = uint8_t(regs | slot_regs(sel[i]));
            add_writes(writes, unit, regs, vals);
        }
        finish_writes(writes);
        return lead_time(writes);
    }

    // Compute the writes for the slots in `sel` (sorted) for the batch at time `t`.
    // The ramps are sampled at `t` (the time of the DAC write for the DAC).
    // The ramps that are not running at the sampling time are ignored.
    void plan(std::vector<Write> &writes, const std::vector<unsigned> &sel,
              int64_t t, const void *data, const StreamBuffer &streambuffer)
    {
        writes.clear();
        for (size_t i = 0; i < sel.size();) {
            auto unit = sel[i] / unit_nslots;
            auto &u = units[unit];
            uint16_t vals[4] = {};
            uint8_t regs = 0;
            auto sample_t = t;
            if (is_dac())
                sample_t += int64_t(writes.size()) * (dma_dac_time + dma_dac_margin);
            for (; i < sel.size() && sel[i] / unit_nslots == unit; i++) {
                auto &slot = slots[sel[i]];
                uint32_t val;
                if (slot.has_pending) {
                    val = slot.pending;
                }
                else if (slot.ramp && slot.ramp->time <= sample_t &&
                         sample_t < slot.ramp->time + slot.ramp->len) {
                    auto type = slot.ramp->chn_type;
                    val = convert_value(type, slot.stream.compute(sample_t, data, retry_time,
                                                                  streambuffer));
                }
                else {
                    continue;
                }
                uint8_t slot_reg = slot_regs(sel[i]);
                if (is_dac()) {
                    vals[0] = uint16_t(val);
                }
                else if (slot_reg == 0x3) {
                    vals[0] = uint16_t(val);
                    vals[1] = uint16_t(val >> 16);
                }
                else if (slot_reg == 0x8) {
                    // Only the lower 12 bits of the amplitude are used, same as the bytecode.
                    vals[3] = uint16_t(val & 0xfff);
                }
                else {
                    vals[2] = uint16_t(val);
                }
                for (unsigned reg = 0; reg < 4; reg++) {
                    if (((slot_reg >> reg) & 1) &&
                        (!u.has_cur[reg] || u.cur[reg] != vals[reg])) {
                        regs = uint8_t(regs | (1 << reg));
                    }
                }
            }
            add_writes(writes, unit, regs, vals);
        }
        finish_writes(writes);
    }

    // Plan the batch for the slots in `sel` to start as soon as possible at or after `t`.
    // Returns the start time.
    int64_t plan_asap(std::vector<Write> &writes, const std::vector<unsigned> &sel,
                      int64_t t, int64_t &lead, const void *data,
                      const StreamBuffer &streambuffer)
    {
        // The registers that need to be written (and therefore the lead time)
        // depends on the values, which depends on the time of the batch for the ramps.
        // Start with the guess in `lead` and only increase the guess
        // so that this converges.
        for (int i = 0; ; i++) {
            plan(writes, sel, t + lead, data, streambuffer);
            auto new_lead = lead_time(writes);
            if (new_lead <= lead) {
                // Delay the start (if needed) to make the outputs happen
                // at the time the ramps were sampled.
                auto start = t + lead - new_lead;
                lead = new_lead;
                return start;
            }
            lead = i < 3 ? new_lead : std::max(new_lead, max_lead_time(sel));
        }
    }

    // Add the writes to the output starting at `t` and return the time
    // when the bus is available again.
    int64_t emit(DMAGroups &groups, const std::vector<Write> &writes, int64_t t)
    {
        for (auto &write: writes) {
            auto &group = groups[t];
            auto &u = units[write.unit];
            if (is_dac()) {
                assert(!group.has_dac);
                group.has_dac = true;
                group.dac_data = write.data;
                u.has_cur[0] = true;
                u.cur[0] = uint16_t(write.data);
            }
            else {
                assert(!group.has_dds[bus_id]);
                group.has_dds[bus_id] = true;
                // The first register is at 0x2c.
                group.dds[bus_id] = { .dds_id = write.unit, .addr = uint8_t(0x16 + write.reg),
                    .is32 = write.is32, .fud = write.fud, .data = write.data };
                u.has_cur[write.reg] = true;
                u.cur[write.reg] = uint16_t(write.data);
                if (write.is32) {
                    u.has_cur[write.reg + 1] = true;
                    u.cur[write.reg + 1] = uint16_t(write.data >> 16);
                }
            }
            t += write_time(write);
        }
        return t;
    }

    void apply_event(const Event &event)
    {
        auto &slot = slots[event.slot];
        auto pulse = event.pulse;
        if (event.ramp_end) {
            assert(slot.ramp == pulse);
            slot.ramp = nullptr;
            nramps--;
            slot.has_pending = true;
            slot.pending = pulse->endvalue;
        }
        else if (pulse->ramp_func) {
            if (!slot.ramp)
                nramps++;
            slot.ramp = pulse;
            slot.stream.add_pulse(*pulse);
            slot.has_pending = false;
        }
        else {
            if (slot.ramp) {
                slot.ramp = nullptr;
                nramps--;
            }
            slot.has_pending = true;
            slot.pending = pulse->endvalue;
        }
    }

    // Schedule all the `pulses` (sorted by time), with the bus available starting at `t`.
    // Returns the time when the bus is available after the last output.
    int64_t run(DMAGroups &groups, int64_t t, const void *data,
                const StreamBuffer &streambuffer)
    {
        std::vector<Event> events;
        events.reserve(pulses.size() * 2);
        for (auto pulse: pulses) {
            auto slot = slot_index(*pulse);
            events.push_back({ pulse->time, slot, false, pulse });
            // The truncated ramps (`!need_endvalue`) are ended by the next pulse
            // on the same channel.
            if (pulse->ramp_func && pulse->need_endvalue) {
                events.push_back({ pulse->time + pulse->len, slot, true, pulse });
            }
        }
        std::stable_sort(events.begin(), events.end(), [] (auto &e1, auto &e2) {
            return e1.time < e2.time;
        });
        pulses.clear();

        const auto nunits = unsigned(units.size());
        std::vector<Write> writes;
        std::vector<unsigned> sel;
        size_t idx = 0;
        while (true) {
            // Stop the (truncated) ramps that have ended.
            // The other ramps are ended by the event setting the end value.
            for (auto &slot: slots) {
                if (slot.ramp && !slot.ramp->need_endvalue &&
                    slot.ramp->time + slot.ramp->len <= t) {
                    slot.ramp = nullptr;
                    nramps--;
                }
            }
            bool has_event = idx < events.size();
            if (!has_event && !nramps)
                break;
            // The batch for the next event must start before `deadline`.
            int64_t deadline = INT64_MAX;
            int64_t event_t = INT64_MAX;
            if (has_event) {
                event_t = events[idx].time;
                sel.clear();
                for (auto i = idx; i < events.size() && events[i].time == event_t; i++)
                    sel.push_back(events[i].slot);
                std::sort(sel.begin(), sel.end());
                deadline = event_t - max_lead_time(sel);
            }
            if (nramps && t < deadline) {
                // Sample the ramps, starting from `next_ramp_unit`
                // and dropping the units from the end if the batch does not finish
                // before the deadline.
                auto &ramp_units = scratch_units;
                ramp_units.clear();
                for (unsigned i = 0; i < nunits; i++) {
                    auto unit = (next_ramp_unit + i) % nunits;
                    for (unsigned j = 0; j < unit_nslots; j++) {
                        if (slots[unit * unit_nslots + j].ramp) {
                            ramp_units.push_back(unit);
                            break;
                        }
                    }
                }
                bool done = false;
                bool changed = true;
                for (auto n = ramp_units.size(); n > 0; n--) {
                    sel.clear();
                    for (size_t i = 0; i < n; i++) {
                        auto unit = ramp_units[i];
                        for (unsigned j = 0; j < unit_nslots; j++) {
                            auto s = unit * unit_nslots + j;
                            if (slots[s].ramp) {
                                sel.push_back(s);
                            }
                        }
                    }
                    std::sort(sel.begin(), sel.end());
                    // The lead time of the ramps is usually the same as the last one.
                    auto start = plan_asap(writes, sel, t, ramp_lead, data, streambuffer);
                    if (writes.empty()) {
                        changed = n != ramp_units.size();
                        break;
                    }
                    if (start + total_time(writes) <= deadline) {
                        t = emit(groups, writes, start);
                        next_ramp_unit = (ramp_units[n - 1] + 1) % nunits;
                        done = true;
                        break;
                    }
                }
                if (done)
                    continue;
                if (!changed) {
                    // No change in any of the ramps, try again later.
                    t = std::min(t + retry_time, deadline);
                    continue;
                }
                // Not enough time for the ramps before the next event.
            }
            if (!has_event)
                continue;
            // Output the values for the event at the requested time if possible.
            sel.clear();
            for (; idx < events.size() && events[idx].time == event_t; idx++) {
                apply_event(events[idx]);
                sel.push_back(events[idx].slot);
            }
            std::sort(sel.begin(), sel.end());
            plan(writes, sel, event_t, data, streambuffer);
            auto lead = lead_time(writes);
            auto start = event_t - lead;
            if (start < t) {
                // Too late, output as soon as possible.
                // Starting from the lead time for the requested time makes sure
                // that the outputs never happen before the requested time
                // (and that all the ramps starting at this time are sampled).
                start = plan_asap(writes, sel, t, lead, data, streambuffer);
            }
            for (auto s: sel)
                slots[s].has_pending = false;
            if (!writes.empty()) {
                t = emit(groups, writes, start);
            }
            else {
                // No change, e.g. a ramp starting from the current value.
                // The ramps can only be sampled after this time.
                t = std::max(t, event_t);
            }
        }
        return t;
    }

    const int bus_id;
    const DDSTiming &timing;
    unsigned unit_nslots;
    int64_t retry_time;
    std::vector<Unit> units;
    std::vector<Slot> slots;
    // Input
    std::vector<Pulse*> pulses;
    int nramps = 0;
    unsigned next_ramp_unit = 0;
    // Lead time of the last batch of ramps.
    int64_t ramp_lead = 0;
    std::vector<Write> scratch_writes;
    std::vector<unsigned> scratch_units;
};

void BCGen::emit_dma(const void *data) const
{
    using namespace DMA;
    bytecode.clear();
    auto &code = bytecode;

    if (!dds_timing.set16 || !dds_timing.set16_fud || !dds_timing.set32 ||
        !dds_timing.set32_fud || !dds_timing.fud_to_end ||
        dds_timing.fud_to_end >= dds_timing.set16_fud)
        throw std::runtime_error("Invalid DDS timing for DMA.");
    if (dds_timing.set32 >= dma_long_wait)
        throw std::runtime_error("DDS write time too long for DMA.");
    auto nmasks = uint32_t(ttl_mask.size());
    if (nmasks < 1 || nmasks > NUM_TTL_BANKS)
        throw std::runtime_error("Invalid number of TTL banks.");
    if (start_ttl_chn < 0 || uint32_t(start_ttl_chn) >= nmasks * 32)
        throw std::runtime_error("Invalid start TTL channel.");

    auto write = [&] (auto obj) {
        static_assert(std::is_trivial_v<decltype(obj)>);
        auto len = code.size();
        code.resize(len + sizeof(obj));
        memcpy(&code[len], &obj, sizeof(obj));
    };

    // Header: [nbanks: 4B][ttl_mask: 4B x nbanks]
    // Version number is in a different ZMQ message and won't be included here.
    // See `BCGen::version()`.
    auto start_bank = start_ttl_chn / 32;
    auto start_bit = uint8_t(start_ttl_chn % 32);
    uint32_t masks[NUM_TTL_BANKS] = {};
    write(nmasks);
    for (uint32_t bank = 0; bank < nmasks; bank++) {
        auto mask = ttl_mask[bank];
        if (int(bank) == start_bank)
            mask = setBit(mask, start_bit, true);
        masks[bank] = mask;
        write(mask);
    }

    // Emitting instructions
    auto emit_wait = [&] (int64_t cycles) {
        assert(cycles >= 0);
        // A wait of `cycle` takes `cycle + 1` cycles.
        constexpr int64_t max_wait2 = int64_t(1) << 28;
        constexpr int64_t max_wait1 = int64_t(1) << 12;
        while (cycles > max_wait2) {
            write(Inst_v0::Wait2(uint32_t(max_wait2 - 1)));
            cycles -= max_wait2;
        }
        if (cycles > max_wait1) {
            write(Inst_v0::Wait2(uint32_t(cycles - 1)));
        }
        else if (cycles > 0) {
            write(Inst_v0::Wait1(uint16_t(cycles - 1)));
        }
    };
    // The TTL values in the hardware for the banks in `ttl_known`.
    uint32_t ttl_hw[NUM_TTL_BANKS] = {};
    bool ttl_known[NUM_TTL_BANKS] = {};
    // Emit the TTL changes of a group.
    // Each instruction sets all the channels of one 4-bit sub-bank (`TTLSet4`, 2 bytes),
    // or two 8-bit (`TTLSet16`, 4 bytes) or 16-bit (`TTLSet32`, 6 bytes) sub-banks,
    // which can be in any of the banks. Each sub-bank must be set at most once in a group.
    // The instructions are chosen to minimize the number (and then the total size)
    // of the instructions.
    auto emit_ttl = [&] (const DMAGroup &group) {
        uint32_t changes[NUM_TTL_BANKS] = {};
        for (int bank = 0; bank < NUM_TTL_BANKS; bank++) {
            if (!((group.ttl_set >> bank) & 1))
                continue;
            auto mask = masks[bank];
            // Skipping the banks without any channel controlled by the DMA is required
            // since the hardware ignores the high bits of the bank number
            // for the banks that don't exist.
            if (!mask)
                continue;
            // The value must not have any bits set outside of the mask
            // since the hardware does not mask the value.
            auto val = group.ttl[bank] & mask;
            changes[bank] = ttl_known[bank] ? (val ^ ttl_hw[bank]) : mask;
            ttl_hw[bank] = val;
            ttl_known[bank] = true;
        }
        // Dynamic programming over the 16-bit sub-banks with changes.
        // For each 16-bit sub-bank, the options are setting it as a 16-bit sub-bank
        // or setting each of its changed 8-bit sub-banks as an 8-bit sub-bank
        // or as (one or two) 4-bit sub-banks.
        // Since two 8-bit (16-bit) sub-banks shares one instruction,
        // the state is the parity of the number of the 8-bit and 16-bit sub-banks.
        // Cost is the number of instructions times 1024 plus the size in bytes.
        constexpr int nhalfs = NUM_TTL_BANKS * 2;
        struct Choice {
            int prev_state = -1;
            // Bit 0: 16-bit sub-bank, bit 1/2: the lower/higher 8-bit sub-bank
            // as an 8-bit sub-bank (otherwise 4-bit sub-banks).
            uint8_t opt = 0;
        };
        int64_t costs[4] = {0, INT64_MAX, INT64_MAX, INT64_MAX}; // state: hpar | bpar << 1
        Choice choices[nhalfs][4];
        int halfs[nhalfs];
        int nchanged = 0;
        for (int half = 0; half < nhalfs; half++) {
            auto change = uint16_t(changes[half / 2] >> ((half % 2) * 16));
            if (!change)
                continue;
            auto idx = nchanged++;
            halfs[idx] = half;
            uint8_t byte_changed = uint8_t(((change & 0xff) != 0) | ((change >> 8 != 0) << 1));
            int64_t new_costs[4] = {INT64_MAX, INT64_MAX, INT64_MAX, INT64_MAX};
            auto update = [&] (int state, int new_state, int64_t cost, uint8_t opt) {
                if (cost < new_costs[new_state]) {
                    new_costs[new_state] = cost;
                    choices[idx][new_state] = {state, opt};
                }
            };
            for (int state = 0; state < 4; state++) {
                auto cost = costs[state];
                if (cost == INT64_MAX)
                    continue;
                auto hpar = state & 1;
                auto bpar = state >> 1;
                // As a 16-bit sub-bank
                update(state, (hpar ^ 1) | (bpar << 1),
                       cost + (hpar ? 0 : 1024 + 6), 1);
                // As 8-bit or 4-bit sub-banks
                for (uint8_t bytes = 0; bytes < 4; bytes++) {
                    if (bytes & ~byte_changed)
                        continue;
                    auto c = cost;
                    auto bp = bpar;
                    for (int i = 0; i < 2; i++) {
                        if (!((byte_changed >> i) & 1))
                            continue;
                        if ((bytes >> i) & 1) {
                            c += bp ? 0 : 1024 + 4;
                            bp ^= 1;
                        }
                        else {
                            auto byte = (change >> (i * 8)) & 0xff;
                            c += (((byte & 0xf) != 0) + ((byte >> 4) != 0)) * (1024 + 2);
                        }
                    }
                    update(state, hpar | (bp << 1), c, uint8_t(bytes << 1));
                }
            }
            memcpy(costs, new_costs, sizeof(costs));
        }
        if (!nchanged)
            return;
        int state = 0;
        for (int s = 1; s < 4; s++) {
            if (costs[s] < costs[state]) {
                state = s;
            }
        }
        // Back track
        std::vector<uint8_t> sub16, sub8, sub4;
        for (int idx = nchanged - 1; idx >= 0; idx--) {
            auto &choice = choices[idx][state];
            auto half = halfs[idx];
            state = choice.prev_state;
            if (choice.opt & 1) {
                sub16.push_back(uint8_t(half));
                continue;
            }
            auto change = uint16_t(changes[half / 2] >> ((half % 2) * 16));
            for (int i = 0; i < 2; i++) {
                auto byte = (change >> (i * 8)) & 0xff;
                if (!byte)
                    continue;
                if ((choice.opt >> (i + 1)) & 1) {
                    sub8.push_back(uint8_t(half * 2 + i));
                    continue;
                }
                for (int j = 0; j < 2; j++) {
                    if ((byte >> (j * 4)) & 0xf) {
                        sub4.push_back(uint8_t(half * 4 + i * 2 + j));
                    }
                }
            }
        }
        assert(state == 0);
        auto sub_val = [&] (int idx, int width) {
            auto per_bank = 32 / width;
            auto val = ttl_hw[idx / per_bank] >> ((idx % per_bank) * width);
            return width == 32 ? val : (val & ((uint32_t(1) << width) - 1));
        };
        // Emit in the order of the sub-banks.
        std::reverse(sub16.begin(), sub16.end());
        std::reverse(sub8.begin(), sub8.end());
        std::reverse(sub4.begin(), sub4.end());
        for (auto sub: sub4)
            write(Inst_v0::TTLSet4(sub, uint8_t(sub_val(sub, 4))));
        // For an odd number of sub-banks, the last one is set twice
        // with the same value in the same instruction.
        for (size_t i = 0; i < sub8.size(); i += 2) {
            auto sub2 = i + 1 < sub8.size() ? sub8[i + 1] : sub8[i];
            write(Inst_v0::TTLSet16(sub8[i], uint8_t(sub_val(sub8[i], 8)),
                                    sub2, uint8_t(sub_val(sub2, 8))));
        }
        for (size_t i = 0; i < sub16.size(); i += 2) {
            auto sub2 = i + 1 < sub16.size() ? sub16[i + 1] : sub16[i];
            write(Inst_v0::TTLSet32(sub16[i], uint16_t(sub_val(sub16[i], 16)),
                                    sub2, uint16_t(sub_val(sub2, 16))));
        }
    };
    auto emit_outputs = [&] (const DMAGroup &group) {
        emit_ttl(group);
        if (group.has_clock)
            write(Inst_v0::ClockOut(group.clock_period));
        for (uint8_t bus = 0; bus < 2; bus++) {
            if (!group.has_dds[bus])
                continue;
            auto &dds = group.dds[bus];
            if (dds.is32) {
                write(Inst_v0::DDSSet32(bus, dds.dds_id, dds.fud, dds.addr, dds.data));
            }
            else {
                write(Inst_v0::DDSSet16(bus, dds.dds_id, dds.fud, dds.addr,
                                        uint16_t(dds.data)));
            }
        }
        if (group.has_dac) {
            write(Inst_v0::DAC(0, dma_dac_div, 0, 0, group.dac_data));
        }
    };
    // Emit the groups between `start` and `end` (exclusive).
    // The instructions emitted before must end with a wait.
    // The instructions emitted will end with a wait.
    auto emit_groups = [&] (const DMAGroups &groups, int64_t start, int64_t end) {
        auto cur = start;
        for (auto &[t, group]: groups) {
            assert(t >= cur);
            emit_wait(t - cur);
            emit_outputs(group);
            cur = t;
        }
        assert(end > cur || groups.empty());
        emit_wait(end - cur);
    };

    const StreamBuffer streambuffer(fpga_clock_div);
    DMABus buses[3] = {DMABus(0, dds_timing), DMABus(1, dds_timing),
        DMABus(2, dds_timing)};

    // Let the hardware buffer some instructions before running any output.
    emit_wait(dma_start_wait);

    // TTL values at the beginning of the sequence.
    uint32_t ttl_vals[NUM_TTL_BANKS] = {};
    for (int bank = 0; bank < NUM_TTL_BANKS; bank++) {
        auto it = m_real_start_vals.find({ChnType::TTL, bank * 32});
        if (it != m_real_start_vals.end()) {
            ttl_vals[bank] = it->second;
        }
    }

    // 1. Initialize all the channels for the first basic sequence.
    if (first_bseq) {
        DMAGroups groups;
        std::vector<Pulse> init_pulses;
        init_pulses.reserve(m_real_start_vals.size());
        for (auto [chn, val]: m_real_start_vals) {
            if (chn.first == ChnType::TTL) {
                // Only the bank values are left after `merge_ttl_pulses`.
                auto bank = chn.second / 32;
                auto &group = groups[0];
                group.ttl[bank] = val;
                group.ttl_set = uint8_t(group.ttl_set | (1 << bank));
                continue;
            }
            init_pulses.push_back({
                    .pulse_type = PulseType::Value,
                    .chn_type = chn.first,
                    .chn = chn.second,
                    .need_endvalue = true,
                    .id = 0,
                    .endvalue = val,
                    .time = 0,
                    .len = 0,
                    .ramp_func = nullptr
                });
        }
        for (auto &pulse: init_pulses) {
            auto bus_id = dma_bus_id(pulse.chn_type, pulse.chn);
            if (bus_id >= 0) {
                buses[bus_id].pulses.push_back(&pulse);
            }
        }
        int64_t init_end = 1;
        for (auto &bus: buses)
            init_end = std::max(init_end, bus.run(groups, 0, data, streambuffer));
        emit_groups(groups, 0, init_end + dma_init_wait);
    }

    // 2. Trigger
    // The sequence time starts after the trigger
    if (trig_type != TrigType::None) {
        write(Inst_v0::WaitTrig(trig_channel, trig_type == TrigType::Raise,
                                uint64_t(trig_timeout_cycle) * dma_cycle_per_10ns));
    }

    // 3. The sequence
    // The time here is the sequence time (in DMA cycles) and may start
    // before 0 for the start trigger and the clock outputs.
    DMAGroups groups;

    // TTL, including the start trigger at `-seq_delay`
    // (high for 1us before that, same as the bytecode).
    struct TTLEvent {
        int64_t time;
        int bank;
        uint32_t val;
        // 0 or 1 to set the start trigger, -1 for a normal TTL pulse.
        int start_trigger;
    };
    std::vector<TTLEvent> ttl_events;
    ttl_events.reserve(m_ttlpulses.size() + 2);
    auto start_trigger_lo = -int64_t(seq_delay) * dma_cycle_per_10ns;
    auto start_trigger_hi = start_trigger_lo - 100 * dma_cycle_per_10ns;
    ttl_events.push_back({start_trigger_hi, start_bank, 0, 1});
    ttl_events.push_back({start_trigger_lo, start_bank, 0, 0});
    for (auto &pulse: m_ttlpulses)
        ttl_events.push_back({pulse.time, pulse.chn / 32, pulse.val, -1});
    std::stable_sort(ttl_events.begin(), ttl_events.end(),
                     [] (auto &e1, auto &e2) { return e1.time < e2.time; });
    for (auto &event: ttl_events) {
        auto &val = ttl_vals[event.bank];
        if (event.start_trigger >= 0) {
            val = setBit(val, start_bit, event.start_trigger != 0);
        }
        else if (event.bank == start_bank) {
            // The TTL pulses don't include the start trigger channel.
            val = setBit(event.val, start_bit, (val >> start_bit) & 1);
        }
        else {
            val = event.val;
        }
        auto &group = groups[event.time];
        group.ttl[event.bank] = val;
        group.ttl_set = uint8_t(group.ttl_set | (1 << event.bank));
    }

    for (auto &clock: clocks) {
        auto &group = groups[convert_time(clock.time)];
        group.has_clock = true;
        group.clock_period = dma_clock_period(clock.period);
    }

    // DDS and DAC
    for (auto &pulse: m_pulses) {
        auto bus_id = dma_bus_id(pulse.chn_type, pulse.chn);
        if (bus_id >= 0) {
            buses[bus_id].pulses.push_back(&pulse);
        }
    }
    // Wait for all the outputs to finish (including the DDS and DAC outputs).
    // The DDS writes for the beginning of the sequence may start early,
    // but not before the start trigger so that the time between the trigger
    // and the start of the sequence is fixed.
    int64_t seq_end = (len_ns + 4) / 5;
    for (auto &bus: buses)
        seq_end = std::max(seq_end, bus.run(groups, start_trigger_hi, data, streambuffer));
    assert(!groups.empty());
    seq_end = std::max(seq_end, groups.rbegin()->first + 1);
    emit_groups(groups, std::min(int64_t(0), groups.begin()->first), seq_end);

    // 4. Finish
    if (!clocks.empty() && clock_tail_ns > 0) {
        // Same as the end of the bytecode sequence on molecube2,
        // output the clock for a while (10ms by default) after the sequence.
        // This is a hack that is believed to make the NI card happy.
        write(Inst_v0::ClockOut(dma_clock_period(10)));
        emit_wait(int64_t((clock_tail_ns + 4) / 5));
        write(Inst_v0::ClockOut(dma_clock_period(0)));
        // The sequence must end with a wait so that the last outputs are executed.
        emit_wait(10);
    }
}

}
