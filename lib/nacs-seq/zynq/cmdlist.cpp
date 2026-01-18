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

#include "cmdlist.h"

#include "exehelper_p.h"
#include "parser_p.h"

#include "../../nacs-utils/streams.h"

#include <array>
#include <cctype>
#include <iomanip>
#include <type_traits>

#include <assert.h>

namespace NaCs::Seq::Zynq::CmdList {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len, uint32_t version)
{
    if (version == 0 || version > 3)
        throw std::runtime_error("Invalid CmdList version number.");
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b;
        assert(op < 10);
        auto len = version == 3 ? Inst_v3::cmd_size[op] : Inst_v1::cmd_size[op];
        i += len;
        count++;
    }
    return count;
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask, uint32_t version, bool dec)
{
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
    Printer printer{stm, dec};
    ExeState state;
    state.run(printer, code, code_len, version);
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         const std::vector<uint32_t> &ttl_masks, uint32_t version,
                         bool dec)
{
    stm << "ttl_mask=" << std::hex;
    bool isfirst = true;
    for (auto ttl_mask: ttl_masks) {
        stm << (isfirst ? "0x" : " 0x") << ttl_mask;
        isfirst = false;
    }
    stm << std::dec << std::endl;
    Printer printer{stm, dec};
    ExeState state;
    state.run(printer, code, code_len, version);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len, uint32_t version)
{
    TimeKeeper keeper;
    ExeState state;
    state.run(keeper, code, code_len, version);
    return keeper.total_t;
}

namespace {

/**
 * The writer class maintains all the states for cmdlist generation.
 */
template<typename Inst>
class Writer {
    // States
    std::array<uint32_t,NUM_TTL_BANKS> all_ttl_mask;

    uint64_t max_time_left = 0;
    ssize_t last_timed_inst = 0;

    buff_ostream &stm;

    template<typename _Inst>
    ssize_t addInst(_Inst inst)
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
        if (op == OpCode::Wait) {
            uint64_t ti;
            memcpy(&ti, &stm[last_timed_inst + 1], 8);
            ti += t;
            memcpy(&stm[last_timed_inst + 1], &ti, 8);
            return;
        }
        if constexpr (Inst::version == 1) {
            if (op == OpCode::TTL1) {
                assert(t <= 0x3);
                uint8_t b = stm[last_timed_inst + 1];
                uint8_t tb = uint8_t(t + (b & 0x3));
                stm[last_timed_inst + 1] = uint8_t((b & ~0x3) | (tb & 0x3));
                return;
            }
        }
        else {
            static_assert(Inst::version == 3);
            if (op == OpCode::TTLAll) {
                assert(t <= 0x1f);
                uint8_t b = stm[last_timed_inst + 1];
                uint8_t tb = uint8_t(t + (b & 0x1f));
                stm[last_timed_inst + 1] = uint8_t((b & ~0x1f) | (tb & 0x1f));
                return;
            }
            if (op == OpCode::TTL1) {
                assert(t <= 0x7f);
                uint8_t b = stm[last_timed_inst + 1];
                uint8_t tb = uint8_t(t + (b & 0x7f));
                stm[last_timed_inst + 1] = uint8_t((b & ~0x7f) | (tb & 0x7f));
                return;
            }
        }
        assert(0 && "Invalid command to increase time.");
        abort();
    }

public:
    static constexpr auto inst_version = Inst::version;

    // The `add***` functions below provides an abstraction to the cmdlist and hides
    // the detail about the cmdlist encoding.
    // This is the same level as the API of `ExeState`.
    void addTTL(uint32_t ttl, int bank)
    {
        assert(bank < NUM_TTL_BANKS);
        all_ttl_mask[bank] = uint32_t(-1);
        if constexpr (Inst::version == 1) {
            assert(bank == 0);
            addInst(typename Inst::TTLAll{OpCode::TTLAll, ttl});
        }
        else {
            static_assert(Inst::version == 3);
            last_timed_inst = addInst(typename Inst::TTLAll{OpCode::TTLAll, 0,
                    uint8_t(bank & (NUM_TTL_BANKS - 1)), ttl});
            max_time_left = 0x1f;
        }
    }

    void addTTL1(uint8_t chn, bool val)
    {
        int bank = chn >> 5;
        assert(bank < NUM_TTL_BANKS);
        all_ttl_mask[bank] = all_ttl_mask[bank] | (1 << (chn & 0x1f));
        if constexpr (Inst::version == 1) {
            assert(bank == 0);
            last_timed_inst = addInst(typename Inst::TTL1{OpCode::TTL1, 0, val,
                    uint8_t(chn & 0x1f)});
            max_time_left = 3;
        }
        else {
            static_assert(Inst::version == 3);
            last_timed_inst = addInst(typename Inst::TTL1{OpCode::TTL1, 0, val, uint8_t(chn)});
            max_time_left = 0x7f;
        }
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
        last_timed_inst = addInst(typename Inst::Wait{OpCode::Wait, dt});
        max_time_left = UINT64_MAX - dt;
    }

    void addClock(uint8_t period)
    {
        addInst(typename Inst::Clock{OpCode::Clock, period});
    }

    void addDDSFreq(uint8_t chn, uint32_t freq)
    {
        if (freq > 0x7fffffff)
            freq = 0x7fffffff;
        addInst(typename Inst::DDSFreq{OpCode::DDSFreq, chn, freq});
    }

    void addDDSAmp(uint8_t chn, uint16_t amp)
    {
        if (amp > 4095)
            amp = 4095;
        addInst(typename Inst::DDSAmp{OpCode::DDSAmp, chn, amp});
    }

    void addDDSPhase(uint8_t chn, uint16_t phase)
    {
        addInst(typename Inst::DDSPhase{OpCode::DDSPhase, chn, phase});
    }

    void addDDSDetPhase(uint8_t chn, uint16_t det_phase)
    {
        addInst(typename Inst::DDSDetPhase{OpCode::DDSDetPhase, chn, det_phase});
    }

    void addDDSReset(uint8_t chn)
    {
        addInst(typename Inst::DDSReset{OpCode::DDSReset, chn});
    }

    void addDAC(uint8_t chn, uint16_t amp)
    {
        addInst(typename Inst::DAC{OpCode::DAC, chn, amp});
    }

    Writer(buff_ostream &stm, const std::array<uint32_t,NUM_TTL_BANKS> &ttl_mask)
        : all_ttl_mask(ttl_mask),
          stm(stm)
    {}

    const std::array<uint32_t,NUM_TTL_BANKS> &get_ttl_mask() const
    {
        return all_ttl_mask;
    }
};

struct Parser : ParserBase {
    using ParserBase::ParserBase;

    template<typename Writer>
    void parse_ttl(Writer &writer)
    {
        skip_whitespace();
        auto c0 = peek();
        if (c0 == '=') {
            writer.addTTL(read_ttlall0(), 0);
        }
        else if (c0 == '(') {
            auto res = read_ttl1(Writer::inst_version >= 3 ? 8 : 1);
            writer.addTTL1(res.first, res.second);
        }
        else if (c0 == '[') {
            auto [ttl, bank] = read_ttlall(Writer::inst_version >= 3 ? 8 : 1);
            writer.addTTL(ttl, bank);
        }
        else {
            syntax_error("Invalid ttl command: expecting `(`, `[` or `=`", colno + 1);
        }
        writer.addWait(read_ttlwait());
    }

    template<typename Writer>
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

NACS_EXPORT() SeqMetadata parse(buff_ostream &ostm, std::istream &istm, uint32_t version)
{
    Parser parser(istm);
    if (version == 0 || version > 3)
        throw std::runtime_error("Invalid CmdList version number.");
    if (version >= 2)
        parser.min_time = PulseTime::Min2;
    SeqMetadata metadata{version, true};
    bool cont;
    std::tie(cont, metadata.ttl_masks) = parser.read_ttlmask(version >= 3 ? 8 : 1);
    auto nmasks = metadata.ttl_masks.size();
    assert(nmasks <= NUM_TTL_BANKS);
    if (!cont)
        return metadata;
    std::array<uint32_t,NUM_TTL_BANKS> ttl_masks{};
    for (size_t bank = 0; bank < nmasks; bank++)
        ttl_masks[bank] = metadata.ttl_masks[bank];
    if (version <= 2) {
        Writer<Inst_v1> writer(ostm, ttl_masks);
        while (parser.parse_cmd(writer)) {
        }
        // Write the updated mask from writer to metadata return
        metadata.ttl_masks.resize(1);
        metadata.ttl_masks[0] = writer.get_ttl_mask()[0];
    }
    else {
        Writer<Inst_v3> writer(ostm, ttl_masks);
        while (parser.parse_cmd(writer)) {
        }
        // Write the updated mask from writer to metadata return
        metadata.ttl_masks.resize(NUM_TTL_BANKS);
        const auto &new_ttl_masks = writer.get_ttl_mask();
        for (size_t bank = 0; bank < NUM_TTL_BANKS; bank++) {
            metadata.ttl_masks[bank] = new_ttl_masks[bank];
        }
    }
    return metadata;
}

}
