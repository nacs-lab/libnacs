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

#ifndef __NACS_SEQ_NIDAQ_DATA_GEN_H__
#define __NACS_SEQ_NIDAQ_DATA_GEN_H__

#include "../host_seq.h"

#include <vector>

#include <stdlib.h>

namespace NaCs::Seq::NiDAQ {

// The channels are identified by an 0-based index
// which may correspond to a physical device and channel.
// It may not correspond to any real channel number.
// The caller is responsible for mapping this to the real channel numbers.
class DataGen {
    struct DataStream;

public:
    enum class PulseType : uint8_t {
        Value,
        Scalar,
        Vector,
    };
    struct Pulse {
        Pulse() = default;
        Pulse(PulseType pulse_type, uint32_t id, uint32_t time_id,
              uint32_t len_id, uint32_t endvalue_id, uint32_t cond_id,
              void (*ramp_func)(void))
            : pulse_type(pulse_type), id(id), time_id(time_id),
              len_id(len_id), endvalue_id(endvalue_id),
              cond_id(cond_id), ramp_func(ramp_func)
        {}
        PulseType pulse_type;

        uint32_t id;
        uint32_t time_id;
        uint32_t len_id; // may be -1, `len == -1` <=> no-ramp pulse
        uint32_t endvalue_id;
        uint32_t cond_id;
        // double (*)(double in) for scalar function
        // void (*)(double *out, const double *in) for vector function
        void (*ramp_func)(void);

        int64_t time = 0; // in sequence time unit
        int64_t len = 0; // in sequence time unit
        int64_t start_step = 0;
        int64_t end_step = 0;
        friend class DataGen;
    };

    DataGen();
    ~DataGen();

    std::vector<Pulse> &get_pulses(uint32_t chn);
    // Sort pulses and generate active_times
    void compute_times();
    // Compute times must be called before generate data
    void generate_data();

    int64_t get_time(uint32_t idx) const;
    double get_value(uint32_t idx) const;

    // Inputs
    std::vector<double> start_values;
    uint32_t nchns;
    uint32_t step_size; // How many sequence clock cycle per refresh step
    std::vector<HostSeq::Type> types;
    HostSeq::Value *values;

    // Outputs
    std::vector<std::pair<int64_t,int64_t>> active_times;
    size_t nsamples;
    std::vector<double> data;
private:

    std::vector<std::vector<Pulse>> m_pulses;
};

}

#endif
