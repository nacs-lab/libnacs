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

#include "exehelper_p.h"

#include "bytecode.h"

#include <assert.h>
#include <math.h>

namespace NaCs::Seq::Zynq::ByteCode {

NACS_EXPORT() size_t count(const uint8_t *code, size_t code_len, uint32_t version)
{
    if (version == 0 || version > 2)
        throw std::runtime_error("Invalid CmdList version number.");
    // Do not use the `ExeState` helper to avoid fusing of instructions.
    size_t count = 0;
    for (size_t i = 0; i < code_len;) {
        uint8_t b = code[i];
        uint8_t op = b & 0xf;
        assert(op < 14);
        auto len = inst_size[op];
        i += len;
        count++;
    }
    return count;
}

NACS_EXPORT() void print(std::ostream &stm, const uint8_t *code, size_t code_len,
                         uint32_t ttl_mask, uint32_t version)
{
    if (version == 0 || version > 2)
        throw std::runtime_error("Invalid CmdList version number.");
    if (ttl_mask)
        stm << "ttl_mask=0x" << std::hex << ttl_mask << std::dec << std::endl;
    Printer printer{stm};
    ExeState state;
    if (version >= 2)
        state.min_time = PulseTime::Min2;
    state.run(printer, code, code_len);
}

NACS_EXPORT() uint64_t total_time(const uint8_t *code, size_t code_len, uint32_t version)
{
    if (version == 0 || version > 2)
        throw std::runtime_error("Invalid CmdList version number.");
    TimeKeeper keeper;
    ExeState state;
    if (version >= 2)
        state.min_time = PulseTime::Min2;
    state.run(keeper, code, code_len);
    return keeper.total_t;
}

}