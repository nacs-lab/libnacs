/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "cmdlist.h"

namespace NaCs {
namespace Seq {
namespace CmdList {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len)
{
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b;
        assert(op < 10);
        auto len = cmd_size[op];
        i += len;
        count++;
    }
    return count;
}

namespace {

struct Printer {
    void ttl(uint32_t ttl, uint64_t t)
    {
        stm << "TTL: val=" << std::hex << ttl << std::dec << " t=" << t << std::endl;
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        stm << "TTL1: chn=" << int(chn) << " val=" << int(val) << " t=" << t << std::endl;
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        stm << "DDS Freq: chn=" << int(chn)
            << " freq=" << std::hex << freq << std::dec << std::endl;
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        stm << "DDS Amp: chn=" << int(chn)
            << " amp=" << std::hex << amp << std::dec << std::endl;
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        stm << "DDS Phase: chn=" << int(chn)
            << " phase=" << std::hex << phase << std::dec << std::endl;
    }
    void dds_detphase(uint8_t chn, uint16_t detphase)
    {
        stm << "DDS Det Phase: chn=" << int(chn)
            << " detphase=" << std::hex << detphase << std::dec << std::endl;
    }
    void dds_reset(uint8_t chn)
    {
        stm << "DDS Reset: chn=" << int(chn) << std::endl;
    }
    void dac(uint8_t chn, uint16_t V)
    {
        stm << "DAC: chn=" << int(chn) << " V=" << std::hex << V << std::dec << std::endl;
    }
    void wait(uint64_t t)
    {
        stm << "Wait: t=" << t << std::endl;
    }
    void clock(uint8_t period)
    {
        stm << "Clock: period=" << int(period) << std::endl;
    }
    std::ostream &stm;
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
        total_t += 50;
    }
    void dds_amp(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_phase(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_detphase(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_reset(uint8_t)
    {
        total_t += 50;
    }
    void dac(uint8_t, uint16_t)
    {
        total_t += 45;
    }
    void wait(uint64_t t)
    {
        total_t += t;
    }
    void clock(uint8_t)
    {
        total_t += 5;
    }
    uint64_t total_t = 0;
};

}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len)
{
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

}
}
}
