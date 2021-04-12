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

#define CATCH_CONFIG_MAIN

#include "../../lib/nacs-seq/nidaq/data_gen.h"

#include "../../lib/nacs-utils/number.h"
#include "../../lib/nacs-utils/processor.h"

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

using namespace NaCs;
using namespace NaCs::Seq;

using active_time_t = std::vector<std::pair<int64_t,int64_t>>;

static unsigned vector_size = [] {
#if NACS_CPU_X86 || NACS_CPU_X86_64
    auto &host = CPUInfo::get_host();
    if (host.test_feature(X86::Feature::avx512f) &&
        host.test_feature(X86::Feature::avx512dq)) {
        return 8;
    }
    else if (host.test_feature(X86::Feature::avx)) {
        return 4;
    }
    else {
        return 2;
    }
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    return 2;
#else
    return 1;
#endif
}();

static void check_data_range(NiDAQ::DataGen &data_gen, uint32_t chn,
                             size_t begin, size_t end, double val)
{
    REQUIRE(begin <= end);
    REQUIRE(end <= data_gen.nsamples);
    REQUIRE(chn < data_gen.nchns);
    REQUIRE(data_gen.nsamples * data_gen.nchns == data_gen.data.size());
    auto offset = chn * data_gen.nsamples;
    for (size_t i = begin; i < end; i++) {
        INFO(i);
        REQUIRE(data_gen.data[i + offset] == val);
    }
}

template<typename Func>
static void check_data_range(NiDAQ::DataGen &data_gen, uint32_t chn,
                             size_t begin, size_t end, Func &&func)
{
    REQUIRE(begin <= end);
    REQUIRE(end <= data_gen.nsamples);
    REQUIRE(chn < data_gen.nchns);
    REQUIRE(data_gen.nsamples * data_gen.nchns == data_gen.data.size());
    auto offset = chn * data_gen.nsamples;
    for (size_t i = begin; i < end; i++) {
        INFO(i);
        REQUIRE(data_gen.data[i + offset] == func(i));
    }
}

TEST_CASE("Empty sequence") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.0);
    data_gen.start_values.push_back(2.0);
    data_gen.nchns = 2;
    data_gen.get_pulses(0);
    data_gen.get_pulses(1);
    data_gen.step_size = 1000;
    data_gen.values = nullptr;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 0}}));
    REQUIRE(data_gen.nsamples == 1);
    data_gen.generate_data();
    REQUIRE(data_gen.data == (std::vector<double>{1.0, 2.0}));
}

TEST_CASE("Short sequence") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 300001 }, { .f64 = -2.3 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 2, 2,
                         uint32_t(-1), 3, uint32_t(-1), nullptr);
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 301}}));
    REQUIRE(data_gen.nsamples == 302);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 302, 1.9);
    check_data_range(data_gen, 1, 0, 301, 2.6);
    check_data_range(data_gen, 1, 301, 302, -2.3);
}

TEST_CASE("Short sequence with ramp") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 300400 }, { .f64 = -2.3 }, { .f64 = 100800 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Scalar, 2, 2,
                         4, 3, uint32_t(-1), (void(*)(void))(double(*)(double))[] (double t) {
                             REQUIRE(t <= 100800);
                             return t / 1000;
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 402}}));
    REQUIRE(data_gen.nsamples == 403);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 403, 1.9);
    check_data_range(data_gen, 1, 0, 301, 2.6);
    check_data_range(data_gen, 1, 301, 402, [&] (auto t) {
        return (double(t) * 1000 - 300400) / 1000;
    });
    check_data_range(data_gen, 1, 402, 403, -2.3);
}

TEST_CASE("Short sequence with ramp (vector)") {
    if (vector_size < 2)
        return;
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 300400 }, { .f64 = -2.3 }, { .f64 = 100800 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Vector, 2, 2,
                         4, 3, uint32_t(-1),
                         (void(*)(void))(void(*)(double*, const double*))
                         [] (double *out, const double *in) {
                             REQUIRE(in[0] <= 100800);
                             for (unsigned i = 0; i < vector_size; i++) {
                                 out[i] = in[i] / 1000;
                             }
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 402}}));
    REQUIRE(data_gen.nsamples == 403);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 403, 1.9);
    check_data_range(data_gen, 1, 0, 301, 2.6);
    check_data_range(data_gen, 1, 301, 402, [&] (auto t) {
        return (double(t) * 1000 - 300400) / 1000;
    });
    check_data_range(data_gen, 1, 402, 403, -2.3);
}

TEST_CASE("Start time") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000001 }, { .f64 = -2.3 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 2, 2,
                         uint32_t(-1), 3, uint32_t(-1), nullptr);
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3001}}));
    REQUIRE(data_gen.nsamples == 1002);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1002, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1002, -2.3);
}

TEST_CASE("Extend start time") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 1010001 }, { .f64 = -2.3 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 2, 2,
                         uint32_t(-1), 3, uint32_t(-1), nullptr);
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1011}}));
    REQUIRE(data_gen.nsamples == 1012);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1012, 1.9);
    check_data_range(data_gen, 1, 0, 1011, 2.6);
    check_data_range(data_gen, 1, 1011, 1012, -2.3);
}

TEST_CASE("Extend start time with ramp") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 990500 }, { .f64 = -2.3 }, { .f64 = 20500 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Scalar, 2, 2,
                         4, 3, uint32_t(-1), (void(*)(void))(double(*)(double))[] (double t) {
                             REQUIRE(t <= 20500);
                             return t / 1000;
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1011}}));
    REQUIRE(data_gen.nsamples == 1012);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1012, 1.9);
    check_data_range(data_gen, 1, 0, 991, 2.6);
    check_data_range(data_gen, 1, 991, 1011, [&] (auto t) {
        return (double(t) * 1000 - 990500) / 1000;
    });
    check_data_range(data_gen, 1, 1011, 1012, -2.3);
}

TEST_CASE("Extend start time with ramp (vector)") {
    if (vector_size < 2)
        return;
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 990500 }, { .f64 = -2.3 }, { .f64 = 20500 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Vector, 2, 2,
                         4, 3, uint32_t(-1),
                         (void(*)(void))(void(*)(double*, const double*))
                         [] (double *out, const double *in) {
                             REQUIRE(in[0] <= 20500);
                             for (unsigned i = 0; i < vector_size; i++) {
                                 out[i] = in[i] / 1000;
                             }
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1011}}));
    REQUIRE(data_gen.nsamples == 1012);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1012, 1.9);
    check_data_range(data_gen, 1, 0, 991, 2.6);
    check_data_range(data_gen, 1, 991, 1011, [&] (auto t) {
        return (double(t) * 1000 - 990500) / 1000;
    });
    check_data_range(data_gen, 1, 1011, 1012, -2.3);
}

TEST_CASE("Ramp") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, uint32_t(-1), nullptr);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Scalar, 2, 2,
                         4, 3, uint32_t(-1), (void(*)(void))(double(*)(double))[] (double t) {
                             REQUIRE(t <= 100950);
                             return t / 1000;
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1051, 1.9);
    check_data_range(data_gen, 0, 1051, 1103, 2.0);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1102, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1102, 1103, -3.9);
}

TEST_CASE("Ramp (vector)") {
    if (vector_size < 2)
        return;
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, uint32_t(-1), nullptr);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Vector, 2, 2,
                         4, 3, uint32_t(-1),
                         (void(*)(void))(void(*)(double*, const double*))
                         [] (double *out, const double *in) {
                             REQUIRE(in[0] <= 100950);
                             for (unsigned i = 0; i < vector_size; i++) {
                                 out[i] = in[i] / 1000;
                             }
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1051, 1.9);
    check_data_range(data_gen, 0, 1051, 1103, 2.0);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1102, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1102, 1103, -3.9);
}

TEST_CASE("Condition interrupt ramp") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }, { .b = true }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Bool};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, 7, nullptr);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Scalar, 2, 2,
                         4, 3, uint32_t(-1), (void(*)(void))(double(*)(double))[] (double t) {
                             REQUIRE(t <= 100950);
                             return t / 1000;
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;

    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1051, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1051, 1103, 2.0);

    values[7].b = false;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1102, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1102, 1103, -3.9);
}

TEST_CASE("Condition interrupt ramp (vector)") {
    if (vector_size < 2)
        return;
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }, { .b = true }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Bool};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, 7, nullptr);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Vector, 2, 2,
                         4, 3, uint32_t(-1),
                         (void(*)(void))(void(*)(double*, const double*))
                         [] (double *out, const double *in) {
                             REQUIRE(in[0] <= 100950);
                             for (unsigned i = 0; i < vector_size; i++) {
                                 out[i] = in[i] / 1000;
                             }
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;

    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1051, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1051, 1103, 2.0);

    values[7].b = false;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1102, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1102, 1103, -3.9);
}

TEST_CASE("Conditional ramp") {
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }, { .b = true }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Bool};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, uint32_t(-1), nullptr);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Scalar, 2, 2,
                         4, 3, 7, (void(*)(void))(double(*)(double))[] (double t) {
                             REQUIRE(t <= 100950);
                             return t / 1000;
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;

    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1051, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1051, 1103, 2.0);

    values[7].b = false;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3051, 3051}}));
    REQUIRE(data_gen.nsamples == 1002);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1002, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1002, 2.0);
}

TEST_CASE("Conditional ramp (vector)") {
    if (vector_size < 2)
        return;
    NiDAQ::DataGen data_gen;
    data_gen.start_values.push_back(1.2);
    data_gen.start_values.push_back(2.6);
    data_gen.nchns = 2;
    HostSeq::Value values[] = {{ .i64 = 200000 }, { .f64 = 1.9 },
                               { .i64 = 3000100 }, { .f64 = -3.9 }, { .f64 = 100950 },
                               { .i64 = 3050600 }, { .f64 = 2.0 }, { .b = true }};
    data_gen.types = {HostSeq::Type::Int64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Float64,
        HostSeq::Type::Int64, HostSeq::Type::Float64, HostSeq::Type::Bool};
    auto &pulses1 = data_gen.get_pulses(0);
    pulses1.emplace_back(NiDAQ::DataGen::PulseType::Value, 3, 0,
                         uint32_t(-1), 1, uint32_t(-1), nullptr);
    auto &pulses2 = data_gen.get_pulses(1);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Value, 1, 5,
                         uint32_t(-1), 6, uint32_t(-1), nullptr);
    pulses2.emplace_back(NiDAQ::DataGen::PulseType::Vector, 2, 2,
                         4, 3, 7,
                         (void(*)(void))(void(*)(double*, const double*))
                         [] (double *out, const double *in) {
                             REQUIRE(in[0] <= 100950);
                             for (unsigned i = 0; i < vector_size; i++) {
                                 out[i] = in[i] / 1000;
                             }
                         });
    data_gen.step_size = 1000;
    data_gen.values = values;

    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3001, 3102}}));
    REQUIRE(data_gen.nsamples == 1103);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1103, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1051, [&] (auto t) {
        return (double(t - 1000) * 1000 - 100) / 1000;
    });
    check_data_range(data_gen, 1, 1051, 1103, 2.0);

    values[7].b = false;
    data_gen.compute_times();
    REQUIRE(data_gen.active_times == (active_time_t{{0, 1000}, {3051, 3051}}));
    REQUIRE(data_gen.nsamples == 1002);
    data_gen.generate_data();
    check_data_range(data_gen, 0, 0, 200, 1.2);
    check_data_range(data_gen, 0, 201, 1002, 1.9);
    check_data_range(data_gen, 1, 0, 1001, 2.6);
    check_data_range(data_gen, 1, 1001, 1002, 2.0);
}
