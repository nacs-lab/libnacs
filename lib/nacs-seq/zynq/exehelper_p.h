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
#include "utils.h"

// This file contains a few classes/functions that are useful for both bytecode
// and command list.

namespace NaCs::Seq::Zynq {

namespace {

struct Printer {
    void ttl(uint32_t ttl, uint64_t t, int bank=0)
    {
        if (bank == 0) {
            stm << "ttl";
        }
        else {
            stm << "ttl[" << bank << "]=";
        }
        stm << "0x" << std::hex << ttl << " t=";
        if (dec) {
            stm << std::dec << t << std::endl;
        }
        else {
            stm << "0x" << std::hex << t << std::dec << std::endl;
        }
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        stm << "ttl(" << int(chn) << ")=" << int(val) << " t=";
        if (dec) {
            stm << std::dec << t << std::endl;
        }
        else {
            stm << "0x" << std::hex << t << std::dec << std::endl;
        }
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        stm << "freq(" << int(chn) << ")=";
        if (dec) {
            stm << dds_freq_from_mu(freq) << std::endl;
        }
        else {
            stm << "0x" << std::hex << freq << std::dec << std::endl;
        }
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        stm << "amp(" << int(chn) << ")=";
        if (dec) {
            stm << dds_amp_from_mu(amp) << std::endl;
        }
        else {
            stm << "0x" << std::hex << amp << std::dec << std::endl;
        }
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        stm << "phase(" << int(chn) << ")=";
        if (dec) {
            stm << dds_phase_from_mu(phase) << std::endl;
        }
        else {
            stm << "0x" << std::hex << phase << std::dec << std::endl;
        }
    }
    void dds_detphase(uint8_t chn, uint16_t detphase)
    {
        stm << "phase(" << int(chn) << ")";
        if (dec) {
            auto sdet = int16_t(detphase);
            if (sdet < 0) {
                stm << "-=" << -double(sdet) / double(1 << 16) << std::endl;
            }
            else {
                stm << "+=" << double(sdet) / double(1 << 16) << std::endl;
            }
        }
        else {
            stm << "+=0x" << std::hex << detphase << std::dec << std::endl;
        }
    }
    void dds_reset(uint8_t chn)
    {
        stm << "reset(" << int(chn) << ")" << std::endl;
    }
    void dac(uint8_t chn, uint16_t V)
    {
        stm << "dac(" << int(chn) << ")=";
        if (dec) {
            stm << V * (20.0 / 65535) - 10.0 << std::endl;
        }
        else {
            stm << "0x" << std::hex << V << std::dec << std::endl;
        }
    }
    void wait(uint64_t t)
    {
        if (dec) {
            stm << "wait(" << t << ")" << std::endl;
        }
        else {
            stm << "wait(0x" << std::hex << t << ")" << std::dec << std::endl;
        }
    }
    void clock(uint8_t period)
    {
        stm << "clock(" << int(period) << ")" << std::endl;
    }
    std::ostream &stm;
    bool dec;
};

struct TimeKeeper {
    void ttl(uint32_t, uint64_t t, int=0)
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

}

} // NaCs::Seq::Zynq

#endif // __NACS_SEQ_ZYNQ_EXEHELPER_P_H__
