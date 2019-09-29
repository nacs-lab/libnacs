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

#include "../../lib/nacs-utils/timer.h"
#include "../../lib/nacs-seq/zynq/bytecode.h"
#include "../../lib/nacs-seq/zynq/seq.h"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace NaCs;
using namespace NaCs::Seq::Zynq;

int main(int argc, char **argv)
{
    assert(argc >= 2);

    std::ifstream istm(argv[1], std::ios::binary);
    istm.seekg(0, std::ios::end);
    auto filesize = (size_t)istm.tellg();
    istm.seekg(0, std::ios::beg);
    std::vector<uint32_t> data(filesize / 4);
    istm.read((char*)data.data(), filesize);
    Timer timer;
    auto code = Sequence::fromBinary(data.data(), data.size())
        .toByteCode(nullptr);
    timer.print();
    size_t code_len;
    auto code2 = Sequence::fromBinary(data.data(), data.size())
        .toByteCode(&code_len, nullptr);

    if (argc >= 3) {
        std::ofstream ostm(argv[2], std::ios::binary);
        ostm.write((const char*)&code[0], code.size());
    }

    assert(code_len == code.size());
    assert(memcmp(code2, &code[0], code_len) == 0);
    free(code2);

    // ByteCode::print(std::cout, code);

    std::cout << ByteCode::count(code) << std::endl;
    std::cout << code.size() << std::endl;

    return 0;
}
