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
#include "utils.h"
#include "../../nacs-utils/processor.h"

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
        auto available_time = time_keeper.next_available(time);

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
    emit_bytecode(host_seq.values.data());
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
    ChnMap<uint32_t> chn_values;
    chn_values.fill(0);
    for (auto [chn, val]: m_real_start_vals) {
        if (chn.first == ChnType::TTL) {
            // TTL on the start trigger bank will be set during `start()`
            uint8_t bank = chn.second / 32;
            if (first_bseq && bank != start_ttl_chn / 32)
                writer.add_ttl(writer.cur_t, val, bank);
            continue;
        }
        chn_values[chn] = val;
        if (first_bseq) {
            auto min_dt = writer.add_pulse(chn.first, chn.second, val, writer.cur_t);
            assert(min_dt >= 0);
            (void)min_dt;
        }
    }

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
    // Shift cur_t so that the first pulse happens 25us after the initialization step.
    writer.cur_t = pseq.single_pulses.begin()->first - 2500;
    for (auto [t, pulse]: pseq.single_pulses)
        writer.add_pulse(pulse.typ, pulse.chn, pulse.val, t);

    // Wait for 1us at the end of the sequence.
    // This just adds a little spacing before the termination sequence.
    writer.add_wait(std::max(int64_t(0), -writer.cur_t) + 100);
}

}
