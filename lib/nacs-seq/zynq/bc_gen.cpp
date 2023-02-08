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
#include "../../nacs-utils/processor.h"

#include <string.h>

namespace NaCs::Seq::Zynq {

struct BCGen::Pulse {
    PulseType pulse_type;
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
    // 0 for non-ramps or pulses that are too short
    int64_t endtime; // in 10ns unit
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
        endtime = time + pulse.len;
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
        assert(t <= endtime);
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
        assert(t <= endtime);
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
        assert(t <= endtime);
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
    for (uint8_t chn = 0; chn < nmgrs && chn < 32; chn++) {
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
    case ChnType::Freq: {
        constexpr double factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;
        if (value <= 0)
            return 0;
        // The instruction only have 31 bits so it actually won't encode exactly 1.75 GHz.
        if (value > 1.7499999987776392e9) // last float that would round to 2^31 - 2
            return 0x7fffffff;
        return round<int32_t>(value * factor);
    }
    case ChnType::Amp: {
        constexpr double factor = 4095;
        if (value <= 0)
            return 0;
        if (value >= 1)
            return 4095;
        return round<int32_t>(value * factor);
    }
    case ChnType::Phase: {
        // Value is phase in degree
        value = fmod(value, 360);
        constexpr double phase_factor = (1 << 14) / 90.0;
        return round<uint16_t>(value * phase_factor);
    }
    case ChnType::DAC: {
        constexpr double factor = 65535 / 20.0;
        constexpr double offset = 10.0;
        if (value <= -10)
            return 0;
        if (value >= 10)
            return 0xffff;
        value += offset;
        return round<int32_t>(value * factor);
    }
    default:
        return 0;
    }
}

NACS_EXPORT() void BCGen::add_ttl_manager(uint8_t chn, int64_t off_delay, int64_t on_delay,
                                          int64_t skip_time, int64_t min_time, bool off_val)
{
    if (chn >= 32)
        throw std::runtime_error("TTL Manager on invalid channel.");
    if (chn >= m_ttl_managers.size())
        m_ttl_managers.resize(chn + 1);
    m_ttl_managers[chn] = TTLManager{off_delay, on_delay, skip_time, min_time, off_val};
}

NACS_EXPORT() uint32_t BCGen::version()
{
    // Bytecode version
    return 2;
}

void BCGen::populate_pulses(const HostSeq &host_seq) const
{
    for (auto &pulses: m_pulses)
        pulses.clear();
    m_ttlpulses.clear();
    for (auto &ttl_mgr: m_ttl_managers)
        ttl_mgr.clear();
    for (auto &seq_pulse: seq_pulses) {
        auto chn_type = seq_pulse.chn_type;
        auto ramp_func = seq_pulse.ramp_func;
        // We cannot translate the time unit yet since that may affect sorting
        // as well as the `time == 0` check.
        auto time = host_seq.get_time(seq_pulse.time_id);
        // Pulses at t=0 are already merged into the initial value.
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
        m_pulses[{chn_type, seq_pulse.chn}].push_back({
                .pulse_type = seq_pulse.pulse_type,
                .id = seq_pulse.id,
                .endvalue = endvalue,
                .time = time,
                .len = len,
                .ramp_func = ramp_func
            });
    }
    preprocess_ttl_managers();
    for (auto &pulses: m_pulses) {
        std::sort(pulses.begin(), pulses.end(), [&] (auto &p1, auto &p2) {
            if (p1.time < p2.time)
                return true;
            if (p1.time > p2.time)
                return false;
            return p1.id < p2.id;
        });
        for (auto &pulse: pulses) {
            assume(pulse.time >= 0);
            pulse.time = convert_time(pulse.time);
        }
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
    for (auto &pulses: m_pulses) {
        int64_t prev_time = -1;
        uint32_t to = 0;
        uint32_t from = 0;
        uint32_t npulses = (uint32_t)pulses.size();
        for (;from < npulses; from++, to++) {
            auto &pulse = pulses[from];
            if (prev_time == pulse.time) {
                to--;
                pulses[to] = pulse;
            }
            else if (from != to) {
                prev_time = pulse.time;
                pulses[to] = std::move(pulse);
            }
        }
        pulses.resize(to);
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
    bool add_ttl_time(int64_t time)
    {
        if (m_npulses < avg_window) {
            m_pulse_times[m_npulses++] = time;
            return true;
        }
        auto idx = m_npulses % avg_window;
        auto prev_idx = (m_npulses - 1) % avg_window;
        // The time span between the last `avg_window` TTL pulses.
        // Note that `idx` is the one we'll override and is the oldest time.
        auto prev_time = m_pulse_times[prev_idx];
        auto cur_span = prev_time - m_pulse_times[idx];
        assert(cur_span > 0);
        assert(time > prev_time);
        if (cur_span < min_time_span / 2 && time < prev_time + 250)
            return false;
        if (cur_span < min_time_span && time < prev_time + 150)
            return false;
        if (cur_span < min_time_span * 2 && time < prev_time + 70)
            return false;
        if (cur_span < min_time_span * 4 && time < prev_time + 20)
            return false;
        if (cur_span < min_time_span * 8 && time < prev_time + 5)
            return false;
        m_npulses++;
        m_pulse_times[idx] = time;
        return true;
    }

private:
    size_t m_npulses = 0;
    int64_t m_pulse_times[avg_window];
};

}

void BCGen::merge_ttl_pulses() const
{
    uint32_t ttl_val = 0;
    for (auto it = m_real_start_vals.begin(), end = m_real_start_vals.end(); it != end;) {
        auto [chn, val] = *it;
        if (chn.first != ChnType::TTL) {
            ++it;
            continue;
        }
        if (chn.second < 32)
            ttl_val = setBit(ttl_val, chn.second, val != 0);
        it = m_real_start_vals.erase(it);
    }
    // Note that unlike for DDS channels,
    // the initial TTL value computed here is always going to be used
    // when we set the startup trigger signal.
    // This also means that we can ignore ttl pulses that doesn't change the value
    // even if it's the first pulse in the sequence.
    m_real_start_vals.emplace(std::make_pair(ChnType::TTL, 0), ttl_val);

    TimeKeeper time_keeper;
    int64_t prev_ttl_t = -1;
    uint32_t to = 0;
    uint32_t from = 0;
    uint32_t nttlpulses = (uint32_t)m_ttlpulses.size();
    for (;from < nttlpulses; from++, to++) {
        auto &pulse = m_ttlpulses[from];
        auto new_ttl_val = setBit(ttl_val, pulse.chn, pulse.val != 0);
        if (new_ttl_val == ttl_val) {
            to--;
            continue;
        }
        ttl_val = new_ttl_val;
        if (prev_ttl_t == pulse.time || !time_keeper.add_ttl_time(pulse.time)) {
            // Neither condition should trigger before we actually added a TTL pulse.
            assert(to != 0);
            to--;
            m_ttlpulses[to].val = ttl_val;
        } else {
            prev_ttl_t = pulse.time;
            pulse.chn = 0;
            pulse.val = ttl_val;
            if (from != to) {
                m_ttlpulses[to] = std::move(pulse);
            }
        }
    }
    m_ttlpulses.resize(to);
}

NACS_EXPORT() void BCGen::generate(const HostSeq &host_seq) const
{
    m_real_start_vals = start_vals;
    populate_pulses(host_seq);
    merge_pulses();
    merge_ttl_pulses();
    emit_bytecode(host_seq.values.data());
}

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

    Writer(std::vector<uint8_t> &bytecode, int8_t start_chn)
        : start_ttl_chn(start_chn),
          bytecode(bytecode)
    {
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

    void write_header(int64_t len_ns, uint32_t ttl_mask)
    {
        // Version number is in a different ZMQ message and won't be included here.
        // See `BCGen::version()`.
        write(len_ns);
        write(setBit(ttl_mask, start_ttl_chn, true));
    }

    // States
    uint32_t cur_ttl = 0;
    bool ttl_set = false;
    DDS dds[22] = {};
    DAC dac[4] = {};

    int8_t start_ttl_chn;
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
            add_inst(Inst::Wait{OpCode::Wait, uint8_t(i & 0xf), times[i]});
        }
        if (times[16]) {
            add_inst(Inst::Wait2{OpCode::Wait2, 1, uint16_t(times[16] & 2047)});
        }
    }

    int add_ttl(int64_t t, uint32_t ttl)
    {
        using namespace ByteCode;
        if (ttl == cur_ttl && ttl_set)
            return 0;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::Min2;
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
            last_timed_inst = add_inst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit & 0x1f),
                    uint8_t(bit & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 2) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            last_timed_inst = add_inst(Inst::TTL2{OpCode::TTL2, 0, uint8_t(bit1 & 0x1f),
                    uint8_t(bit2 & 0x1f)});
            max_time_left = 3;
        }
        else if (nchgs == 3) {
            auto bit1 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit1);
            auto bit2 = __builtin_ffs(changes) - 1;
            changes = changes ^ (1 << bit2);
            auto bit3 = __builtin_ffs(changes) - 1;
            add_inst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
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
            add_inst(Inst::TTL4{OpCode::TTL4, uint8_t(bit1 & 0x1f),
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
            last_timed_inst = add_inst(Inst::TTL5{OpCode::TTL5, 0, uint8_t(bit1 & 0x1f),
                    uint8_t(bit2 & 0x1f), uint8_t(bit3 & 0x1f),
                    uint8_t(bit4 & 0x1f), uint8_t(bit5 & 0x1f)});
            max_time_left = 7;
        }
        else {
            last_timed_inst = add_inst(Inst::TTLAll{OpCode::TTLAll, 0, ttl});
            max_time_left = 15;
        }
        return PulseTime::Min2;
    }

    int add_ttl_single(int64_t t, uint8_t chn, bool val)
    {
        return add_ttl(t, setBit(cur_ttl, chn, val));
    }

    int add_clock(int64_t t, uint8_t period)
    {
        using namespace ByteCode;
        assert(t >= cur_t);
        add_wait(t - cur_t);
        cur_t += PulseTime::Clock;
        add_inst(Inst::Clock{OpCode::Clock, 0, uint8_t(period - 1)});
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
            add_inst(Inst::DDSDetFreq2{OpCode::DDSDetFreq2, uint8_t(chn & 0x1f),
                    uint8_t(dfreq & 0x7f)});
        }
        else if (dfreq <= 0x3fff || dfreq >= 0xffffc000) {
            add_inst(Inst::DDSDetFreq3{OpCode::DDSDetFreq3, uint8_t(chn & 0x1f),
                    uint16_t(dfreq & 0x7fff)});
        }
        else if (dfreq <= 0x003fffff || dfreq >= 0xffc00000) {
            add_inst(Inst::DDSDetFreq4{OpCode::DDSDetFreq4, uint8_t(chn & 0x1f),
                    uint32_t(dfreq & 0x7fffff)});
        }
        else {
            add_inst(Inst::DDSFreq{OpCode::DDSFreq, uint8_t(chn & 0x1f),
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
        add_inst(Inst::DDSPhase{OpCode::DDSPhase, uint8_t(chn & 0x1f), 0, phase});
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
            add_inst(Inst::DDSDetAmp{OpCode::DDSDetAmp, uint8_t(chn & 0x1f),
                    uint8_t(damp & 0x7f)});
        }
        else {
            add_inst(Inst::DDSAmp{OpCode::DDSAmp, 0, uint8_t(chn & 0x1f),
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
            add_inst(Inst::DACDet{OpCode::DACDet, uint8_t(chn & 0x3), uint16_t(dV & 0x3ff)});
        }
        else {
            add_inst(Inst::DAC{OpCode::DAC, 0, uint8_t(chn & 0x3), V});
        }
        return PulseTime::DAC;
    }

    void start()
    {
        // wait 25us
        auto t = cur_t + 2500;
        // Note that this will also make sure the initial ttl value is set.
        add_ttl_single(t, start_ttl_chn, true);
        // 1us
        t += 100;
        add_ttl_single(t, start_ttl_chn, false);
        cur_t = 0;
    }

    void end(int64_t t)
    {
        // 1us
        t += 100;
        assert(t >= cur_t);
        add_wait(t - cur_t);
    }

    uint8_t pulse_mintime(ChnType chn_typ)
    {
        if (chn_typ == ChnType::DAC)
            return PulseTime::DAC;
        if (chn_typ == ChnType::Freq)
            return PulseTime::DDSFreq;
        if (chn_typ == ChnType::Amp)
            return PulseTime::DDSAmp;
        if (chn_typ == ChnType::Phase)
            return PulseTime::DDSPhase;
        return PulseTime::Min2;
    }

    int add_pulse(ChnType chn_typ, uint8_t chn, uint32_t val, int64_t t)
    {
        switch (chn_typ) {
        case ChnType::TTL:
            return add_ttl(t, val);
        case ChnType::Freq:
            return add_freq(t, chn, val);
        case ChnType::Amp:
            return add_amp(t, chn, uint16_t(val));
        case ChnType::Phase:
            return add_phase(t, chn, uint16_t(val));
        case ChnType::DAC:
            return add_dac(t, chn, uint16_t(val));
        default:
            return 0;
        }
    }
};

void BCGen::emit_bytecode(const void *data) const
{
    bytecode.clear();
    Writer writer(bytecode, start_ttl_chn);
    writer.write_header(len_ns, ttl_mask);

    auto ttl_start_it = m_real_start_vals.find({ChnType::TTL, 0});
    if (ttl_start_it != m_real_start_vals.end())
        writer.cur_ttl = ttl_start_it->second;
    // Initialize channels
    if (first_bseq) {
        for (auto [chn, val]: m_real_start_vals) {
            // TTL will be set during `start()`
            if (chn.first == ChnType::TTL)
                continue;
            auto min_dt = writer.add_pulse(chn.first, chn.second, val, writer.cur_t);
            assert(min_dt >= 0);
            (void)min_dt;
        }
    }
    uint32_t clock_idx = 0;
    auto nclocks = (uint32_t)clocks.size();
    int64_t next_clock_time;
    if (clocks.empty()) {
        next_clock_time = INT64_MAX;
    }
    else {
        next_clock_time = convert_time(clocks.front().time);
    }
    int64_t time_offset = 0;
    auto real_seq_delay = seq_delay;
    // writer.cur_t == time_offset -> sequence time 0.
    if (next_clock_time < 0) {
        auto clock_offset = uint32_t(-next_clock_time);

        // If we want to start the clock right away (which we most likely will)
        // the time for the first clock command will be negative since
        // we want the first lowering edge to happen at the correct time.
        if (real_seq_delay <= clock_offset) {
            real_seq_delay = 0;
        }
        else {
            real_seq_delay -= clock_offset;
        }
        writer.start(); // cur_t == 0
        auto period = clocks.front().period;
        writer.add_clock(real_seq_delay, period);
        // cur_t == real_seq_delay + PulseTime::Clock
        // the clock output was added at `real_seq_delay` or `cur_t - PulseTime::Clock`
        clock_idx = 1;
        if (clock_idx < nclocks) {
            assume(clocks[clock_idx].time >= 0);
            next_clock_time = convert_time(clocks[clock_idx].time);
        }
        else {
            next_clock_time = INT64_MAX;
        }
        // Now align the start of the sequence time to get the correct t=0.
        // we should have `<input time for 1st clock> + time_offset == <cur_t for 1st clock>`
        // or `-clock_offset + time_offset == cur_t - PulseTime::Clock`
        // or `time_offset == cur_t - PulseTime::Clock + clock_offset`
        // Moreover, since we might have a pulse at t=0,
        // which would have a `cur_t` of `0 + time_offset`
        // we want to make sure `time_offset >= cur_t` so that we don't go back in time.
        if (PulseTime::Clock < clock_offset) {
            time_offset = writer.cur_t + clock_offset - PulseTime::Clock;
        }
        else {
            time_offset = writer.cur_t;
        }
    }
    else {
        // Start the sequence and restart timer.
        writer.start();
        time_offset = real_seq_delay;
    }

    // Now we have three kinds of pulses to deal with
    // 1. Clock out: these will always be outputted at the right time with the highest priority
    // 2. TTL out: after we make sure we can do the clock correct,
    //    we will try to respect ttl output time as much as possible.
    // 3. DDS/DAC: output at the correct time if we can

    const StreamBuffer streambuffer(fpga_clock_div);
    ChnMap<DataStream> streams;
    ChnMap<uint32_t> pulse_idxs;
    pulse_idxs.fill(0);

    int64_t next_ttl_time = m_ttlpulses.empty() ? INT64_MAX : m_ttlpulses.front().time;

    bool has_pulses = false;
    for (auto &pulses: m_pulses)
        has_pulses |= !pulses.empty();
    auto nttlpulses = (uint32_t)m_ttlpulses.size();
    // Index points to the item to be processed next.
    uint32_t ttlpulse_idx = 0;
    // Time between each ramp update. 50 for each DDS, 45 for each DAC
    uint32_t step_time = 0;

    // Record the last channel that could do an output before the deadline.
    // This way, when we are cut-off by a deadline, we can restart from where we left off.
    int chn_offset = 0;

    int64_t cur_time = 0;

    constexpr auto max_ramp_step = max(PulseTime::_DDS, PulseTime::DAC);
    auto should_continue = [&] {
        return clock_idx < nclocks || ttlpulse_idx < nttlpulses || has_pulses;
    };

    while (should_continue()) {
        // Record the time when we started.
        // If the channels are ramping too slowly we could end up not forwarding
        // the time at all in the loop otherwise.
        auto start_time = cur_time;
        auto next_critical_time = min(next_clock_time, next_ttl_time);

        // Whether there's a on-going ramp.
        // When this is true, we should not jump forward in time to the next pulse until
        // the current ramps are all finished.
        bool has_ramp = false;
        // This is used to determine how much we can jump forward
        // (together with next_critical_time) if there's no on-going ramps.
        // We need to update this whenever we set has_pulse without has_ramp.
        int64_t next_pulse_time = INT64_MAX;

        // First check if we can do any normal pulse output before the deadline.
        // We'll output at most one pulse per channel.
        // Also recompute if we have any pulses to handle in the next cycle.
        if (has_pulses) {
            has_pulses = false;

            int new_chn_offset = chn_offset;

            // If there might be space to do a output.
            for (int _chn_idx = 0; _chn_idx < m_pulses.chn_num; _chn_idx++) {
                int chn_idx = (_chn_idx + chn_offset) % m_pulses.chn_num;
                auto [chn_type, chn] = m_pulses.to_channel(chn_idx);
                auto &pulses = m_pulses[chn_idx];
                auto &stream = streams[chn_idx];
                auto &pulse_idx = pulse_idxs[chn_idx];
                auto mindt = writer.pulse_mintime(chn_type);
                if (cur_time + mindt > next_critical_time) {
                    if (pulse_idx < pulses.size()) {
                        // We don't have time to output anything
                        // but if a new ramp is starting we can setup our datastream.
                        auto &pulse = pulses[pulse_idx];
                        if (pulse.ramp_func && pulse.time <= cur_time + mindt) {
                            if (!stream.ramp_func)
                                step_time += mindt;
                            stream.add_pulse(pulse);
                            pulse_idx++;
                            has_ramp = true;
                        }
                        else if (stream.ramp_func) {
                            has_ramp = true;
                        }
                        else {
                            next_pulse_time = min(next_pulse_time, pulse.time);
                        }
                        has_pulses = true;
                    }
                    else if (stream.ramp_func) {
                        has_pulses = true;
                        has_ramp = true;
                    }
                    continue;
                }
                new_chn_offset = chn_idx + 1;
                // Before we compute any values, make sure the pulse is still the current one.
                // First check if it is replaced.
                if (pulse_idx < pulses.size()) {
                    if (pulses[pulse_idx].time < cur_time + mindt) {
                        // If we have many DDS pulses in this time range,
                        // find the last one...
                        while (pulse_idx + 1 < pulses.size() &&
                               pulses[pulse_idx + 1].time < cur_time + mindt)
                            pulse_idx++;
                        auto &pulse = pulses[pulse_idx];
                        assert(pulse.time < cur_time + mindt);
                        pulse_idx++;
                        if (pulse.ramp_func) {
                            if (!stream.ramp_func)
                                step_time += mindt;
                            stream.add_pulse(pulse);
                        }
                        else {
                            // For non-ramp pulses, simply output its value.
                            if (pulse_idx < pulses.size()) {
                                next_pulse_time = min(next_pulse_time, pulses[pulse_idx].time);
                                has_pulses = true;
                            }
                            // If the time is slightly later,
                            // see if we have enough time to move it to the correct time.
                            auto t = cur_time;
                            if (pulse.time > t)
                                t = min(pulse.time, next_critical_time - mindt);
                            // The `add_pulse` might be a no-op
                            // in which case we do not need to forward the `cur_time`.
                            // Most importantly, if the first pulse in the sequence is skipped,
                            // `writer.cur_t` might not have be increased
                            // to be larger than `time_offset` yet.
                            if (writer.add_pulse(chn_type, chn, pulse.endvalue,
                                                 t + time_offset)) {
                                assert(writer.cur_t >= time_offset);
                                cur_time = max(cur_time, writer.cur_t - time_offset);
                            }
                            if (stream.ramp_func)
                                step_time -= mindt;
                            assert(step_time <= streams.chn_num * max_ramp_step);
                            stream.clear();
                            continue;
                        }
                    }
                }
                // No current ramp, move to the next channel.
                if (!stream.ramp_func) {
                    if (pulse_idx < pulses.size()) {
                        next_pulse_time = min(next_pulse_time, pulses[pulse_idx].time);
                        has_pulses = true;
                    }
                    continue;
                }
                // Next, check if the current one is finished.
                if (stream.endtime <= cur_time + mindt) {
                    if (writer.add_pulse(chn_type, chn, stream.endvalue,
                                         cur_time + time_offset)) {
                        assert(writer.cur_t >= time_offset);
                        cur_time = max(cur_time, writer.cur_t - time_offset);
                    }
                    step_time -= mindt;
                    assert(step_time <= streams.chn_num * max_ramp_step);
                    stream.clear();
                    if (pulse_idx < pulses.size()) {
                        next_pulse_time = min(next_pulse_time, pulses[pulse_idx].time);
                        has_pulses = true;
                    }
                    continue;
                }
                // Now we have a on-going ramp, output one sample on it.
                has_ramp = true;
                has_pulses = true;
                auto v = convert_value(chn_type, stream.compute(max(cur_time, stream.time),
                                                                data, step_time, streambuffer));
                if (writer.add_pulse(chn_type, chn, v, cur_time + time_offset)) {
                    assert(writer.cur_t >= time_offset);
                    cur_time = max(cur_time, writer.cur_t - time_offset);
                }
            }
            chn_offset = new_chn_offset;
        }

        // Code below assume we still need to do something,
        // break out of the loop if we have nothing left to do.
        if (!should_continue())
            break;
        // If there's no pulse for a while, we can forward the time.
        if (!has_ramp) {
            assert(min(next_critical_time, next_pulse_time) != INT64_MAX);
            cur_time = max(cur_time, min(next_critical_time, next_pulse_time));
        }
        else if (start_time == cur_time) {
            // We have ramps but didn't output anything.
            // This is either because we hit a deadline,
            // or because none of the ramps changed anything.
            // Make sure we are alway forwarding the time.
            cur_time = min(next_critical_time, cur_time + max_ramp_step);
        }
        // We can run another loop without having to check TTL or clocks.
        if (cur_time + max_ramp_step <= next_critical_time)
            continue;
        // Now we can deal with TTL and Clock pulses.
        if (next_ttl_time < next_clock_time) {
            static_assert(PulseTime::Min2 == 1);
            // Try outputting a TTL pulse. Since its length is 1,
            // we'll always have enough time to output the clock at the correct time
            // as long as it's after the TTL in time.

            // Since `next_ttl_time < next_clock_time <= INT64_MAX`
            // we know that we have a current ttl pulse.
            assert(next_ttl_time == m_ttlpulses[ttlpulse_idx].time);
            if (writer.add_ttl(next_ttl_time + time_offset, m_ttlpulses[ttlpulse_idx].val)) {
                assert(writer.cur_t >= time_offset);
                cur_time = max(cur_time, writer.cur_t - time_offset);
            }
            ttlpulse_idx++;
            if (ttlpulse_idx < nttlpulses) {
                next_ttl_time = m_ttlpulses[ttlpulse_idx].time;
            }
            else {
                next_ttl_time = INT64_MAX;
            }
            continue;
        }
        if (clock_idx < nclocks) {
            assert(next_clock_time == convert_time(clocks[clock_idx].time));
            writer.add_clock(next_clock_time + time_offset, clocks[clock_idx].period);
            assert(writer.cur_t >= time_offset);
            cur_time = max(cur_time, writer.cur_t - time_offset);
            clock_idx++;
            if (clock_idx < nclocks) {
                assume(clocks[clock_idx].time >= 0);
                next_clock_time = convert_time(clocks[clock_idx].time);
            }
            else {
                next_clock_time = INT64_MAX;
            }
        }
    }
    // Wait for 1us at the end of the sequence.
    // This just adds a little spacing before the termination sequence.
    writer.end(cur_time + time_offset);
}

}
