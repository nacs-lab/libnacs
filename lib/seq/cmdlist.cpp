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

#include "exehelper_p.h"

#include "cmdlist.h"

#include "../utils/streams.h"

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

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len)
{
    print(stm, code, code_len, 0);
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask)
{
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
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

namespace {

/**
 * The writer class maintains all the states for cmdlist generation.
 */
class Writer {
    // States
    uint32_t all_ttl_mask;

    uint8_t max_time_left = 0;
    ssize_t last_timed_inst = 0;

    buff_ostream &stm;

    template<typename Inst>
    ssize_t addInst(Inst inst)
    {
        auto len = stm.tellg();
        stm.write((char*)&inst, sizeof(inst));
        max_time_left = 0;
        return (ssize_t)len;
    }

    // Increase the wait time encoded in the last instruction by `t`.
    // The caller should have checked that the time fits in the maximum time possible to be
    // encoded.
    void incLastTime(uint8_t t)
    {
        assert(t <= max_time_left);
        max_time_left = uint8_t(max_time_left - t);
        uint8_t op = stm[last_timed_inst];
        uint8_t *tp;
        uint8_t tmask;
        int tshift;
        if (op == 1) {
            // TTL1
            tmask = 0xc0;
            tp = (uint8_t*)&stm[last_timed_inst + 1];
            tshift = 6;
        }
        else {
            assert(0 && "Invalid command to increase time.");
            abort();
        }
        uint8_t b = *tp;
        t = uint8_t(t + ((b & tmask) >> tshift));
        *tp = uint8_t((b & ~tmask) | ((t << tshift) & tmask));
    }

public:
    // The `add***` functions below provides an abstraction to the cmdlist and hides
    // the detail about the cmdlist encoding.
    // This is the same level as the API of `ExeState`.
    void addTTL(uint32_t ttl)
    {
        all_ttl_mask = uint32_t(-1);
        addInst(Inst::TTLAll{OpCode::TTLAll, ttl});
    }

    void addTTL1(uint8_t chn, bool val)
    {
        chn = chn & 0x1f;
        all_ttl_mask = all_ttl_mask | (1 << chn);
        last_timed_inst = addInst(Inst::TTL1{OpCode::TTL1, 0, val, chn});
        max_time_left = 3;
    }

    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        if (dt <= max_time_left) {
            incLastTime(uint8_t(dt));
            return;
        }
        if (max_time_left) {
            dt -= max_time_left;
            incLastTime(max_time_left);
        }
        addInst(Inst::Wait{OpCode::Wait, dt});
    }

    void addClock(uint8_t period)
    {
        addInst(Inst::Clock{OpCode::Clock, period});
    }

    void addDDSFreq(uint8_t chn, uint32_t freq)
    {
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        addInst(Inst::DDSFreq{OpCode::DDSFreq, chn, freq});
    }

    void addDDSAmp(uint8_t chn, uint16_t amp)
    {
        if (amp > 4095)
            amp = 4095;
        addInst(Inst::DDSAmp{OpCode::DDSAmp, chn, amp});
    }

    void addDDSPhase(uint8_t chn, uint16_t phase)
    {
        addInst(Inst::DDSPhase{OpCode::DDSPhase, chn, phase});
    }

    void addDDSDetPhase(uint8_t chn, uint16_t det_phase)
    {
        addInst(Inst::DDSDetPhase{OpCode::DDSDetPhase, chn, det_phase});
    }

    void addDDSReset(uint8_t chn)
    {
        addInst(Inst::DDSReset{OpCode::DDSReset, chn});
    }

    void addDAC(uint8_t chn, uint16_t amp)
    {
        addInst(Inst::DAC{OpCode::DAC, chn, amp});
    }

    Writer(buff_ostream &stm, uint32_t ttl_mask)
        : all_ttl_mask(ttl_mask),
          stm(stm)
    {}

    uint32_t get_ttl_mask() const
    {
        return all_ttl_mask;
    }
};

}

}
}
}
