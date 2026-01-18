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

#ifndef __NACS_SEQ_ZYNQ_BC_GEN_H__
#define __NACS_SEQ_ZYNQ_BC_GEN_H__

#include "../host_seq.h"

#include <nacs-utils/number.h>

#include <map>
#include <set>
#include <vector>

#include <assert.h>

namespace NaCs::Seq::Zynq {

class BCGen {
    struct Pulse;
    struct TTLPulse;
    struct TTLManager;
    template<typename Inst> struct Writer;
    struct DataStream;

public:
    enum class ChnType : uint8_t {
        TTL = 1,
        Freq = 2,
        Amp = 3,
        DAC = 4,
        Phase = 5,
        Clock = 6,
    };
    enum class PulseType : uint8_t {
        Value,
        Scalar,
        Vector,
    };
    struct SeqPulse {
        ChnType chn_type;
        PulseType pulse_type;
        uint8_t chn;

        uint32_t id;
        uint32_t time_id;
        uint32_t len_id; // may be -1, `len == -1` <=> no-ramp pulse
        uint32_t endvalue_id;
        uint32_t cond_id; // may be -1
        // double (*)(double in, const void *data) for scalar function
        // void (*)(double *out, const double *in, const void *data) for vector function
        // Use reference passing for vector function since GCC doesn't support
        // vector calling convention on windows
        // and clang doesn't support per-function vector calling convention.
        void (*ramp_func)(void);
    };
    struct Clock {
        int64_t time; // in sequence time unit
        uint8_t period;
    };

    BCGen();
    ~BCGen();
    void generate(const HostSeq &host_seq) const;
    static uint32_t convert_value(ChnType type, double value);
    static double compute_start_val(PulseType pulse_type, void (*ramp_func)(void),
                                    const void *data);
    void add_ttl_manager(uint8_t chn, int64_t off_delay, int64_t on_delay,
                         int64_t skip_time, int64_t min_time, bool off_val);
    uint32_t version() const;

    // Inputs
    uint32_t seq_idx;
    uint32_t seq_delay; // in unit of 10ns.
    std::vector<SeqPulse> seq_pulses;
    // Assume this includes all the channels needed.
    // For `first_bseq`, this is needed for all channels.
    // For non-`first_bseq`, this is needed only for TTL pulses.
    // In both cases, this is expected to be the value without correcting for t=0 pulses.
    std::map<std::pair<ChnType,uint8_t>,uint32_t> start_vals;
    std::vector<Clock> clocks;
    bool first_bseq;

    int8_t start_ttl_chn;
    std::vector<uint32_t> ttl_mask{0};
    // How many sequence cycle corresponds to one FPGA clock cycle.
    uint32_t fpga_clock_div;
    int64_t len_ns;

    // Output
    mutable std::vector<uint8_t> bytecode;

private:
    template<typename T>
    struct ChnMap {
        static constexpr int chn_num = 22 + 22 + 4 + 22;
        T channels[chn_num];

        template<typename T2>
        void fill(T2 v)
        {
            for (auto &chn: channels) {
                chn = v;
            }
        }
        static std::pair<BCGen::ChnType,uint8_t> to_channel(int idx)
        {
            if (idx < 22)
                return {BCGen::ChnType::Freq, idx};
            idx -= 22;
            if (idx < 22)
                return {BCGen::ChnType::Amp, idx};
            idx -= 22;
            if (idx < 4)
                return {BCGen::ChnType::DAC, idx};
            idx -= 4;
            return {BCGen::ChnType::Phase, idx};
        }
        static int to_index(std::pair<BCGen::ChnType,uint8_t> chn)
        {
            if (chn.first == BCGen::ChnType::Freq) {
                assert(chn.second < 22);
                return chn.second;
            }
            else if (chn.first == BCGen::ChnType::Amp) {
                assert(chn.second < 22);
                return chn.second + 22;
            }
            else if (chn.first == BCGen::ChnType::DAC) {
                assert(chn.second < 4);
                return chn.second + 22 * 2;
            }
            else if (chn.first == BCGen::ChnType::Phase) {
                assert(chn.second < 22);
                return chn.second + 22 * 2 + 4;
            }
            assert(false && "Invalid channel type.");
            abort();
        }
        const T &operator[](int idx) const
        {
            assert(idx < chn_num);
            return channels[idx];
        }
        T &operator[](int idx)
        {
            assert(idx < chn_num);
            return channels[idx];
        }
        const T &operator[](std::pair<BCGen::ChnType,uint8_t> chn) const
        {
            return this->operator[](to_index(chn));
        }
        T &operator[](std::pair<BCGen::ChnType,uint8_t> chn)
        {
            return this->operator[](to_index(chn));
        }
        T *begin()
        {
            return channels;
        }
        T *end()
        {
            return &channels[chn_num];
        }
        const T *begin() const
        {
            return channels;
        }
        const T *end() const
        {
            return &channels[chn_num];
        }
    };

    int64_t convert_time(int64_t seq_time) const
    {
        seq_time += seq_time >= 0 ? int32_t(fpga_clock_div) / 2 : -int32_t(fpga_clock_div) / 2;
        return seq_time / fpga_clock_div;
    }
    int64_t convert_time(double seq_time) const
    {
        return round<int64_t>(seq_time / fpga_clock_div);
    }
    void populate_pulses(const HostSeq &host_seq) const;
    bool preprocess_ttl_pulse(const TTLPulse &ttl_pulse) const;
    void preprocess_ttl_managers() const;
    void merge_pulses() const;
    void merge_ttl_pulses() const;

    void emit_bytecode(const void *data) const;

    template<typename Inst>
    void _emit_bytecode(const void *data) const;

    std::vector<TTLManager> m_ttl_managers;
    mutable std::vector<Pulse> m_pulses;
    mutable std::vector<TTLPulse> m_ttlpulses;
    // This is basically a copy of `start_vals` which we will mutate during generation.
    // We need to keep all the inputs untouched in order to be able to do multiple generations.
    mutable std::map<std::pair<ChnType,uint8_t>,uint32_t> m_real_start_vals;
};

}

#endif
