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

#ifndef __NACS_SEQ_ZYNQ_EXEHELPER_P_H__
#define __NACS_SEQ_ZYNQ_EXEHELPER_P_H__

#include <ostream>

#include "pulse_time.h"
#include "bc_gen.h"

// This file contains a few classes/functions that are useful for both bytecode
// and command list.

namespace NaCs::Seq::Zynq {

namespace {

struct Printer {
    void ttl(uint32_t ttl, uint64_t t)
    {
        stm << "ttl=0x" << std::hex << ttl << " t=0x" << t << std::dec << std::endl;
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        stm << "ttl(" << int(chn) << ")=" << int(val)
            << " t=0x" << std::hex << t << std::dec << std::endl;
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        stm << "freq(" << int(chn) << ")=0x"
            << std::hex << freq << std::dec << std::endl;
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        stm << "amp(" << int(chn) << ")=0x"
            << std::hex << amp << std::dec << std::endl;
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        stm << "phase(" << int(chn) << ")=0x"
            << std::hex << phase << std::dec << std::endl;
    }
    void dds_detphase(uint8_t chn, uint16_t detphase)
    {
        stm << "phase(" << int(chn) << ")+=0x"
            << std::hex << detphase << std::dec << std::endl;
    }
    void dds_reset(uint8_t chn)
    {
        stm << "reset(" << int(chn) << ")" << std::endl;
    }
    void dac(uint8_t chn, uint16_t V)
    {
        stm << "dac(" << int(chn) << ")=0x" << std::hex << V << std::dec << std::endl;
    }
    void wait(uint64_t t)
    {
        stm << "wait(0x" << std::hex << t << ")" << std::dec << std::endl;
    }
    void clock(uint8_t period)
    {
        stm << "clock(" << int(period) << ")" << std::endl;
    }
    std::ostream &stm;
};

struct ArgPrinter {
    void ttl(uint32_t ttl, uint64_t t)
    {
        printf("ttl: %u, %lu\n", ttl, t);
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        printf("ttl1: %u, %d, %lu\n", chn, val, t);
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        printf("dds_freq: %u, %u\n", chn, freq);
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        printf("dds_amp: %u, %u\n", chn, amp);
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        printf("dds_phase: %u, %u\n", chn, phase);
    }
    void dds_detphase(uint8_t chn, uint16_t detphase)
    {
        printf("dds_detphase: %u, %u\n", chn, detphase);
    }
    void dds_reset(uint8_t chn)
    {
        printf("dds_reset: %u\n", chn);
    }
    void dac(uint8_t chn, uint16_t V)
    {
        printf("dac: %u, %u\n", chn, V);
    }
    void wait(uint64_t t)
    {
        printf("wait: %lu\n", t);
    }
    void clock(uint8_t period)
    {
        printf("clock: %u\n", period);
    }
};

struct TimeKeeper {
    void ttl(uint32_t, uint64_t t)
    {
        total_t += t;
    }
    void ttl1(uint8_t, bool, uint64_t t)
    {
        total_t += t;
    }
    void dds_freq(uint8_t, uint32_t)
    {
        total_t += PulseTime::DDSFreq;
    }
    void dds_amp(uint8_t, uint16_t)
    {
        total_t += PulseTime::DDSAmp;
    }
    void dds_phase(uint8_t, uint16_t)
    {
        total_t += PulseTime::DDSPhase;
    }
    void dds_detphase(uint8_t, uint16_t)
    {
        total_t += PulseTime::DDSPhase;
    }
    void dds_reset(uint8_t)
    {
        total_t += PulseTime::DDSReset;
    }
    void dac(uint8_t, uint16_t)
    {
        total_t += PulseTime::DAC;
    }
    void wait(uint64_t t)
    {
        total_t += t;
    }
    void clock(uint8_t)
    {
        total_t += PulseTime::Clock;
    }
    uint64_t total_t = 0;
};
struct PulseCollector : TimeKeeper {

    PulseCollector(BCGen::ChnType type, uint8_t num) :
        m_chntype(type),
        chn_num(num)
    {
    }

    void ttl(uint32_t ttl, uint64_t t)
    {
        if (first_ttl) {
            cur_ttl_val = ttl;
            first_ttl = false;
            return;
        }
        if (!seq_started) {
            // Look for a lowering edge on ttl 0
            if (get_bit(cur_ttl_val, 0) != 0 && get_bit(ttl, 0) == 0) {
                seq_started = true;
            }
            else {
                cur_ttl_val = ttl;
                return;
            }
        }
        if (m_chntype == BCGen::ChnType::TTL){
            ts.push_back(total_t);
            vals.push_back(uint32_t(get_bit(ttl, chn_num) > 0));
        }
        cur_ttl_val = ttl;
        TimeKeeper::ttl(ttl, t); // Advances the time
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        if (first_ttl) {
            if (val)
                cur_ttl_val = 1 << chn;
            first_ttl = false;
            return;
        }
        if (!seq_started) {
            // check for lowering edge on chn 0
            if (chn != 0 || val || get_bit(cur_ttl_val, 0) == 0) {
                // cannot be start condition, just update ttls
                if (val)
                    cur_ttl_val = (1 << chn) | cur_ttl_val;
                else
                    cur_ttl_val = (~(1 << chn)) & cur_ttl_val;
                return;
            }
            else {
                // must be start condition
                seq_started = true;
            }
        }
        if (val)
            cur_ttl_val = (1 << chn) | cur_ttl_val;
        else
            cur_ttl_val = (~(1 << chn)) & cur_ttl_val;
        if (m_chntype == BCGen::ChnType::TTL){
            ts.push_back(total_t);
            vals.push_back(uint32_t(get_bit(cur_ttl_val, chn_num) > 0));
        }
        TimeKeeper::ttl1(chn, val, t);
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        if (m_chntype == BCGen::ChnType::Freq && chn_num == chn) {
            ts.push_back(total_t);
            vals.push_back(double(freq) * 3.5e9/((uint64_t(1) << 16) * (uint64_t(1) << 16))); // 3.5e9 / 2^32 is conversion factor
        }
        if (seq_started)
            TimeKeeper::dds_freq(chn, freq);
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        if (m_chntype == BCGen::ChnType::Amp && chn_num == chn) {
            ts.push_back(total_t);
            vals.push_back(double(amp) / 4096.0);
        }
        if (seq_started)
            TimeKeeper::dds_amp(chn, amp);
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        //if (m_chntype == BCGen::ChnType::Phase && chn_num == chn) {
        //    ts.push_back(total_t);
        //    vals.push_back(double(phase) * 90.0 / (1 << 14)); // in degrees
        //}
        if (seq_started)
            TimeKeeper::dds_phase(chn, phase);
    }
    void dds_detphase(uint8_t chn, uint16_t phase)
    {
        //if (m_chntype == BCGen::ChnType::DDSPhase && chn_num == chn) {
        //    ts.push_back(total_t);
        //    if (vals.size() == 0)
        //        vals.push_back(double(phase) * 90.0 / (1 << 14)); // in degrees
        //    else
        //        vals.push_back(vals.back() + double(phase) * 90.0 / (1 << 14));
        //}
        if (seq_started)
            TimeKeeper::dds_detphase(chn, phase);
    }
    void dds_reset(uint8_t chn)
    {
        if ((m_chntype == BCGen::ChnType::Freq || m_chntype == BCGen::ChnType::Amp) && chn_num == chn) {
            ts.push_back(total_t);
            vals.push_back(0);
        }
        if (seq_started)
            TimeKeeper::dds_reset(chn);
    }
    void dac(uint8_t chn, uint16_t V)
    {
        if (seq_started) {
            TimeKeeper::dac(chn, V);
        }
    }
    void wait(uint64_t t)
    {
        if (seq_started) {
            TimeKeeper::wait(t);
        }
    }
    void clock(uint8_t period)
    {
        if (seq_started) {
            TimeKeeper::clock(period);
        }
    }

    uint32_t get_bit(uint32_t val, uint32_t idx)
    {
        // Returns 0 if 0, otherwise returns non zero value, but not necessarily 1
        return (val & (1 << idx));
    }

    std::vector<uint64_t> ts;
    std::vector<double> vals;
    BCGen::ChnType m_chntype;
    uint8_t chn_num;
    bool seq_started = false;
    bool first_ttl = true;
    uint32_t cur_ttl_val = 0;
};

}

} // NaCs::Seq::Zynq

#endif // __NACS_SEQ_ZYNQ_EXEHELPER_P_H__
