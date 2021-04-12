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

#include "data_gen.h"

#include "../../nacs-utils/number.h"
#include "../../nacs-utils/processor.h"

#include <atomic>

#include <assert.h>

namespace NaCs::Seq::NiDAQ {

namespace {

// This needs to be consistent with the compiler of the ramp functions.
struct StreamBuffer {
    StreamBuffer(uint32_t step_size)
        : step_size(step_size)
    {
#if NACS_CPU_X86 || NACS_CPU_X86_64
        for (int i = 0; i < 8; i++) {
            static_assert(sizeof(time) == 8 * sizeof(double));
            time[i] = i * step_size;
        }
#elif NACS_CPU_AARCH64
        for (int i = 0; i < 2; i++) {
            static_assert(sizeof(time) == 2 * sizeof(double));
            time[i] = i * step_size;
        }
#endif
    }
#if NACS_CPU_X86 || NACS_CPU_X86_64
    double time[8] __attribute__((aligned(64)));
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    double time[2] __attribute__((aligned(16)));
#endif
    uint32_t step_size;
};

// Clamping values to [-10, 10]
// Note that we only do this for scalar values.
// We assume the ramp functions have clamping implemented in them since it's more efficient
static double bound_output(double value)
{
    return std::isnan(value) ? 0 : bound(-10, value, 10);
}

}

struct DataGen::DataStream {
    static void ramp_scalar(double *__restrict__ data, size_t begin, size_t end,
                            int64_t toffset, const Pulse &pulse, const StreamBuffer &streambuffer)
    {
        auto func = (double (*)(double))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        for (auto i = begin; i < end; i++) {
            data[i] = func(double(i * step_size - start_time + toffset));
        }
    }
    static void fill_scalar(double *__restrict__ data, size_t begin, size_t end, double value)
    {
        value = bound_output(value);
        for (auto i = begin; i < end; i++) {
            data[i] = value;
        }
    }

    // This needs to be consistent with the compiler of the ramp functions.
#if NACS_CPU_X86 || NACS_CPU_X86_64
    using ramp_func_t = void (*)(double*, size_t, size_t, int64_t,
                                 const Pulse&, const StreamBuffer&);
    static std::atomic<ramp_func_t> _ramp_vector;
    static void ramp(double *__restrict__ data, size_t begin, size_t end,
                     int64_t toffset, const Pulse &pulse, const StreamBuffer &streambuffer)
    {
        if (pulse.pulse_type == PulseType::Scalar) {
            ramp_scalar(data, begin, end, toffset, pulse, streambuffer);
            return;
        }
        // Relaxed load is OK here since we only care about the pointer itself
        // but not anything ordered with it.
        _ramp_vector.load(std::memory_order_relaxed)(data, begin, end, toffset,
                                                     pulse, streambuffer);
    }
    static void ramp_vector_sse(double *__restrict__ data, size_t begin, size_t end,
                                int64_t toffset, const Pulse &pulse,
                                const StreamBuffer &streambuffer)
    {
        auto func = (void (*)(double*, const double*))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        auto i = begin;
        for (; i + 1 < end; i += 2) {
            auto ts __attribute__((aligned(16))) =
                double(i * step_size - start_time + toffset) + _mm_load_pd(streambuffer.time);
            func(&data[i], (const double*)&ts);
        }
        if (i < end) {
            assert(i == end - 1);
            auto ts __attribute__((aligned(16))) =
                double(i * step_size - start_time + toffset) + _mm_load_pd(streambuffer.time);
            double res[2] __attribute__((aligned(16)));
            func(res, (const double*)&ts);
            data[i] = res[0];
        }
    }
    __attribute__((target("avx")))
    static void ramp_vector_avx(double *__restrict__ data, size_t begin, size_t end,
                                int64_t toffset, const Pulse &pulse,
                                const StreamBuffer &streambuffer)
    {
        auto func = (void (*)(double*, const double*))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        auto i = begin;
        for (; i + 3 < end; i += 4) {
            auto ts __attribute__((aligned(32))) =
                double(i * step_size - start_time + toffset) + _mm256_load_pd(streambuffer.time);
            func(&data[i], (const double*)&ts);
        }
        if (i < end) {
            assert(i > end - 4);
            auto ts __attribute__((aligned(32))) =
                double(i * step_size - start_time + toffset) + _mm256_load_pd(streambuffer.time);
            double res[4] __attribute__((aligned(32)));
            func(res, (const double*)&ts);
            assume(end - i < 4);
            memcpy(&data[i], &res[0], (end - i) * sizeof(double));
        }
    }
    __attribute__((target("avx2")))
    static void ramp_vector_avx2(double *__restrict__ data, size_t begin, size_t end,
                                 int64_t toffset, const Pulse &pulse,
                                 const StreamBuffer &streambuffer)
    {
        // Exactly the same as the AVX version but can use `vbroadcast` instructions
        auto func = (void (*)(double*, const double*))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        auto i = begin;
        for (; i + 3 < end; i += 4) {
            auto ts __attribute__((aligned(32))) =
                double(i * step_size - start_time + toffset) + _mm256_load_pd(streambuffer.time);
            func(&data[i], (const double*)&ts);
        }
        if (i < end) {
            assert(i > end - 4);
            auto ts __attribute__((aligned(32))) =
                double(i * step_size - start_time + toffset) + _mm256_load_pd(streambuffer.time);
            double res[4] __attribute__((aligned(32)));
            func(res, (const double*)&ts);
            assume(end - i < 4);
            memcpy(&data[i], &res[0], (end - i) * sizeof(double));
        }
    }
    __attribute__((target("avx512f,avx512dq")))
    static void ramp_vector_avx512(double *__restrict__ data, size_t begin, size_t end,
                                   int64_t toffset, const Pulse &pulse,
                                   const StreamBuffer &streambuffer)
    {
        auto func = (void (*)(double*, const double*))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        auto i = begin;
        for (; i + 7 < end; i += 8) {
            auto ts __attribute__((aligned(64))) =
                double(i * step_size - start_time + toffset) + _mm512_load_pd(streambuffer.time);
            func(&data[i], (const double*)&ts);
        }
        if (i < end) {
            assert(i > end - 8);
            auto ts __attribute__((aligned(64))) =
                double(i * step_size - start_time + toffset) + _mm512_load_pd(streambuffer.time);
            double res[8] __attribute__((aligned(64)));
            func(res, (const double*)&ts);
            assume(end - i < 8);
            memcpy(&data[i], &res[0], (end - i) * sizeof(double));
        }
    }
    static void ramp_vector_dispatch(double *__restrict__ data, size_t begin, size_t end,
                                     int64_t toffset, const Pulse &pulse,
                                     const StreamBuffer &streambuffer)
    {
        // Lazy dispatch on first call
        // In principle, this function can be called multiple times
        // if we have more than one thread calling this but that's OK
        // since different execution should have the same result.
        auto &host = CPUInfo::get_host();
        if (host.test_feature(X86::Feature::avx512f) &&
            host.test_feature(X86::Feature::avx512dq)) {
            _ramp_vector.store(&ramp_vector_avx512, std::memory_order_relaxed);
        }
        else if (host.test_feature(X86::Feature::avx2)) {
            _ramp_vector.store(&ramp_vector_avx2, std::memory_order_relaxed);
        }
        else if (host.test_feature(X86::Feature::avx)) {
            _ramp_vector.store(&ramp_vector_avx, std::memory_order_relaxed);
        }
        else {
            _ramp_vector.store(&ramp_vector_sse, std::memory_order_relaxed);
        }
        _ramp_vector.load(std::memory_order_relaxed)(data, begin, end, toffset,
                                                     pulse, streambuffer);
    }
    using fill_func_t = void (*)(double*, size_t, size_t, double);
    static std::atomic<fill_func_t> _fill_vector;
    static void fill(double *__restrict__ data, size_t begin, size_t end, double value)
    {
        _fill_vector.load(std::memory_order_relaxed)(data, begin, end, value);
    }
    // GCC somehow generates slightly suboptimal code for AVX and AVX2 version.
    // Ref https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100088
    __attribute__((target("avx")))
    static void fill_vector_avx(double *__restrict__ data,
                                size_t begin, size_t end, double value)
    {
        value = bound_output(value);
        for (auto i = begin; i < end; i++) {
            data[i] = value;
        }
    }
    __attribute__((target("avx2")))
    static void fill_vector_avx2(double *__restrict__ data,
                                 size_t begin, size_t end, double value)
    {
        value = bound_output(value);
        for (auto i = begin; i < end; i++) {
            data[i] = value;
        }
    }
    __attribute__((target("avx512f,avx512dq")))
    static void fill_vector_avx512(double *__restrict__ data,
                                   size_t begin, size_t end, double value)
    {
        value = bound_output(value);
        for (auto i = begin; i < end; i++) {
            data[i] = value;
        }
    }
    static void fill_vector_dispatch(double *__restrict__ data,
                                     size_t begin, size_t end, double value)
    {
        // Lazy dispatch on first call
        // In principle, this function can be called multiple times
        // if we have more than one thread calling this but that's OK
        // since different execution should have the same result.
        auto &host = CPUInfo::get_host();
        if (host.test_feature(X86::Feature::avx512f) &&
            host.test_feature(X86::Feature::avx512dq)) {
            _fill_vector.store(&fill_vector_avx512, std::memory_order_relaxed);
        }
        else if (host.test_feature(X86::Feature::avx2)) {
            _fill_vector.store(&fill_vector_avx2, std::memory_order_relaxed);
        }
        else if (host.test_feature(X86::Feature::avx)) {
            _fill_vector.store(&fill_vector_avx, std::memory_order_relaxed);
        }
        else {
            _fill_vector.store(&fill_scalar, std::memory_order_relaxed);
        }
        _fill_vector.load(std::memory_order_relaxed)(data, begin, end, value);
    }
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    static void ramp(double *__restrict__ data, size_t begin, size_t end,
                     int64_t toffset, const Pulse &pulse, const StreamBuffer &streambuffer)
    {
        if (pulse.pulse_type == PulseType::Scalar) {
            ramp_scalar(data, begin, end, pulse, streambuffer);
            return;
        }
        auto func = (void (*)(double*, const double*))pulse.ramp_func;
        auto step_size = streambuffer.step_size;
        auto start_time = pulse.time;
        auto i = begin;
        for (; i + 1 < end; i += 2) {
            auto ts __attribute__((aligned(16))) =
                double(i * step_size - start_time + toffset) + vld1q_f64(streambuffer.time);
            func(&data[i], (const double*)&ts);
        }
        if (i < end) {
            assert(i == end - 1);
            auto ts __attribute__((aligned(16))) =
                double(i * step_size - start_time + toffset) + vld1q_f64(streambuffer.time);
            double res[2] __attribute__((aligned(16)));
            func(res, (const double*)&ts);
            data[i] = res[0];
        }
    }
    static void fill(double *__restrict__ data, size_t begin, size_t end, double value)
    {
        fill_scalar(data, begin, end, value);
    }
#else
    static void ramp(double *__restrict__ data, size_t begin, size_t end,
                     int64_t toffset, const Pulse &pulse, const StreamBuffer &streambuffer)
    {
        assert(pulse.pulse_type == PulseType::Scalar);
        ramp_scalar(data, begin, end, toffset, pulse, streambuffer);
    }
    static void fill(double *__restrict__ data, size_t begin, size_t end, double value)
    {
        fill_scalar(data, begin, end, value);
    }
#endif
};

#if NACS_CPU_X86 || NACS_CPU_X86_64
std::atomic<DataGen::DataStream::ramp_func_t>
DataGen::DataStream::_ramp_vector = ramp_vector_dispatch;
std::atomic<DataGen::DataStream::fill_func_t>
DataGen::DataStream::_fill_vector = fill_vector_dispatch;
#endif

NACS_EXPORT_ DataGen::DataGen()
{
}

NACS_EXPORT() DataGen::~DataGen()
{
}

NACS_EXPORT() int64_t DataGen::get_time(uint32_t idx) const
{
    assert(types[idx] == HostSeq::Type::Int64);
    return values[idx].i64;
}

NACS_EXPORT() double DataGen::get_value(uint32_t idx) const
{
    auto type = types[idx];
    auto value = values[idx];
    switch (type) {
    case HostSeq::Type::Bool:
        return value.b;
    case HostSeq::Type::Int32:
        return value.i32;
    case HostSeq::Type::Float64:
        return value.f64;
    case HostSeq::Type::Int64:
        return (double)value.i64;
    default:
        return 0;
    }
}

NACS_EXPORT() std::vector<DataGen::Pulse> &DataGen::get_pulses(uint32_t chn)
{
    if (chn >= m_pulses.size())
        m_pulses.resize(chn + 1);
    return m_pulses[chn];
}

NACS_EXPORT() void DataGen::compute_times()
{
    active_times.clear();
    int64_t max_step = 0;
    for (auto &pulses: m_pulses) {
        for (auto &pulse: pulses) {
            auto time = get_time(pulse.time_id);
            auto len = pulse.len_id == uint32_t(-1) ? 0 :
                round<int64_t>(get_value(pulse.len_id));
            if (len < 0)
                len = 0;
            // This is the first time point after the start of the pulse.
            // The pulse has no effect on anything before this point
            // The `(x - 1) / size + 1` trick only works for positive numbers
            // so we explicitly check for `time == 0` and return 0 instead.
            auto start_step = time == 0 ? 0 : ((time - 1) / step_size + 1);
            // This is the first time point after the end of the pulse.
            // The value at this point (assuming not being overwriten by another pulse)
            // is `endvalue`
            auto end_step = len == 0 ? start_step : ((time + len - 1) / step_size + 1);
            if (pulse.cond_id == uint32_t(-1) || get_value(pulse.cond_id) != 0) {
                // Note that in principle we could have a short pulse
                // interrupting a longer one and the time range for the longer one
                // would not actually be needed.
                // However, including them is still correct
                // and I don't expect that to be a case that is worth optimizing for...
                if (end_step > max_step)
                    max_step = end_step;
                active_times.emplace_back(start_step, end_step);
                pulse.time = time;
                pulse.len = len;
                pulse.start_step = start_step;
                pulse.end_step = end_step;
            }
            else {
                // Sort disabled pulse to the very end so that they won't accidentally
                // interrupt any normal ramps
                pulse.time = INT64_MAX;
                pulse.len = 0;
                pulse.start_step = INT64_MAX;
                pulse.end_step = INT64_MAX;
            }
        }
        std::sort(pulses.begin(), pulses.end(), [&] (auto &p1, auto &p2) {
            if (p1.time < p2.time)
                return true;
            if (p1.time > p2.time)
                return false;
            return p1.id < p2.id;
        });
    }
    // Output up to 1000 clock cycles (2 ms for 500 kHz rate) up front.
    // When the NI DAQ card does not recieved any clock or trigger it will time out.
    // If it recieved some clocks but not enough continuously it may prematurely
    // stop the sequence without reporting any error.
    // From experiment, outputing 1000 cycles up front seems to fix that problem.
    if (max_step <= 1000) {
        active_times.clear();
        active_times.push_back({0, max_step});
        nsamples = size_t(max_step + 1);
    }
    else {
        assert(!active_times.empty());
        active_times.push_back({0, 1000});
        std::sort(active_times.begin(), active_times.end(), [&] (auto &p1, auto &p2) {
            return p1 < p2;
        });
        nsamples = 1001;
        uint32_t from = 1;
        uint32_t ntimes = (uint32_t)active_times.size();
        auto *prev_time = &active_times[0];
        for (;from < ntimes; from++) {
            auto &time = active_times[from];
            assert(prev_time->first <= time.first);
            // Merge nearby clock pulses to avoid too frequent on-off.
            // the 20 cycle here is a pretty arbitrary cutoff
            // (note that it should at least be 1 since
            //  we record active time as closed range and a range ending on `n`
            //  with the next one starting on `n + 1` is continuous).
            if (prev_time->second + 20 >= time.first) {
                if (time.second > prev_time->second) {
                    nsamples += size_t(time.second - prev_time->second);
                    prev_time->second = time.second;
                }
            }
            else {
                prev_time++;
                if (&time != prev_time)
                    *prev_time = time;
                nsamples += size_t(time.second - time.first + 1);
            }
        }
        active_times.resize((prev_time - &active_times[0]) + 1);
    }
}

NACS_EXPORT() void DataGen::generate_data()
{
    data.resize(nchns * nsamples);
    const StreamBuffer streambuffer(step_size);
    auto nactive_times = active_times.size();
    for (uint32_t i = 0; i < nchns; i++) {
        auto chn_data = &data[i * nsamples];
        auto &pulses = m_pulses[i];
        auto value = start_values[i];
        size_t data_idx = 0; // Next data to be filled.
        size_t active_time_idx = 0; // Current active_time to look at
        size_t total_active_time = 0; // Total steps that we've looked at
        int64_t toffset = 0;
        auto npulses = pulses.size();
        auto find_time_idx = [&] (int64_t step) {
            assert(active_times[active_time_idx].first <= step);
            while (active_times[active_time_idx].second < step) {
                auto active_time = active_times[active_time_idx];
                total_active_time += active_time.second - active_time.first + 1;
                assert(active_time_idx + 1 < nactive_times);
                toffset += (active_times[active_time_idx + 1].first -
                            active_time.second - 1) * step_size;
                (void)nactive_times;
                active_time_idx++;
                assert(active_times[active_time_idx].first <= step);
            }
            return total_active_time + (step - active_times[active_time_idx].first);
        };
        for (uint32_t pulse_idx = 0; pulse_idx < npulses; pulse_idx++) {
            auto &pulse = pulses[pulse_idx];
            auto cond = pulse.cond_id == uint32_t(-1) || get_value(pulse.cond_id) != 0;
            if (!cond)
                continue;
            // For each pulse, fill all the previous values with the current value.
            auto start_step = pulse.start_step;
            auto start_idx = find_time_idx(start_step);
            assert(data_idx <= start_idx);
            if (data_idx < start_idx)
                DataStream::fill(chn_data, data_idx, start_idx, value);
            value = get_value(pulse.endvalue_id);
            data_idx = start_idx;

            // Then look at the next pulse to figure out the range of data
            // we need to fill in using the ramp.
            auto end_step = pulse.end_step;
            if (pulse_idx + 1 < npulses)
                end_step = min(end_step, pulses[pulse_idx + 1].start_step);
            if (start_step >= end_step) // no actual output necessary
                continue;
            assert(active_times[active_time_idx].second >= end_step);
            // We do not include the `end_step` here
            // since that will be filled with `endvalue` instead.
            auto nsteps = end_step - start_step;
            DataStream::ramp(chn_data, data_idx, data_idx + nsteps, toffset,
                             pulse, streambuffer);
            data_idx = data_idx + nsteps;
        }
        if (data_idx < nsamples) {
            DataStream::fill(chn_data, data_idx, nsamples, value);
        }
    }
}

}
