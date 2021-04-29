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

#include "exehelper_p.h"
#include "parser_p.h"

#include "../utils/streams.h"

#include <cctype>
#include <iomanip>
#include <type_traits>

#include <assert.h>

namespace NaCs::Seq::CmdList {

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

    uint64_t max_time_left = 0;
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
    void incLastTime(uint64_t t)
    {
        assert(t <= max_time_left);
        max_time_left = max_time_left - t;
        uint8_t op = stm[last_timed_inst];
        if (op == 1) {
            assert(t <= 3);
            // TTL1
            uint8_t b = stm[last_timed_inst + 1];
            uint8_t tb = uint8_t(t + (b & 0x3));
            stm[last_timed_inst + 1] = uint8_t((b & ~0x3) | (tb & 0x3));
        }
        else if (op == 2) {
            uint64_t ti;
            memcpy(&ti, &stm[last_timed_inst + 1], 8);
            ti += t;
            memcpy(&stm[last_timed_inst + 1], &ti, 8);
        }
        else {
            assert(0 && "Invalid command to increase time.");
            abort();
        }
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
        all_ttl_mask = all_ttl_mask | (1 << chn);
        last_timed_inst = addInst(Inst::TTL1{OpCode::TTL1, 0, val, uint8_t(chn & 0x1f)});
        max_time_left = 3;
    }

    void addWait(uint64_t dt)
    {
        if (dt == 0)
            return;
        if (dt <= max_time_left) {
            incLastTime(dt);
            return;
        }
        if (max_time_left) {
            dt -= max_time_left;
            incLastTime(max_time_left);
        }
        last_timed_inst = addInst(Inst::Wait{OpCode::Wait, dt});
        max_time_left = UINT64_MAX - dt;
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

struct Parser : ParserBase {
    using ParserBase::ParserBase;

    void parse_ttl(Writer &writer)
    {
        skip_whitespace();
        auto c0 = peek();
        if (c0 == '=') {
            writer.addTTL(read_ttlall());
        }
        else if (c0 == '(') {
            auto res = read_ttl1();
            writer.addTTL1(res.first, res.second);
        }
        else {
            syntax_error("Invalid ttl command: expecting `(` or `=`", colno + 1);
        }
        writer.addWait(read_ttlwait());
    }

    bool parse_cmd(Writer &writer)
    {
        skip_whitespace();
        auto nres = read_name();
        if (nres.second == -1)
            syntax_error("Expecting command name", colno + 1);
        if (nres.first == "ttl") {
            parse_ttl(writer);
        }
        else if (nres.first == "wait") {
            writer.addWait(read_waitcmd());
        }
        else if (nres.first == "freq") {
            auto res = read_freqcmd();
            writer.addDDSFreq(res.first, res.second);
        }
        else if (nres.first == "amp") {
            auto res = read_ampcmd();
            writer.addDDSAmp(res.first, res.second);
        }
        else if (nres.first == "phase") {
            auto res = read_phasecmd();
            if (res.second.first) {
                writer.addDDSDetPhase(res.first, res.second.second);
            }
            else {
                writer.addDDSPhase(res.first, res.second.second);
            }
        }
        else if (nres.first == "reset") {
            writer.addDDSReset(read_ddschn("reset"));
        }
        else if (nres.first == "dac") {
            auto res = read_daccmd();
            writer.addDAC(res.first, res.second);
        }
        else if (nres.first == "clock") {
            writer.addClock(read_clockcmd());
        }
        else {
            syntax_error("Unknown command name", -1, nres.second + 1, colno);
        }
        if (!checked_next_line())
            return false;
        return skip_comments();
    }
};

}

NACS_EXPORT() uint32_t parse(buff_ostream &ostm, std::istream &istm)
{
    Parser parser(istm);
    bool cont;
    uint32_t ttl_mask;
    std::tie(cont, ttl_mask) = parser.read_ttlmask();
    if (!cont)
        return ttl_mask;
    Writer writer(ostm, ttl_mask);
    while (parser.parse_cmd(writer)) {
    }
    return writer.get_ttl_mask();
}

}
