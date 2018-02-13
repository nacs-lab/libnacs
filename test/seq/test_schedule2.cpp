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

#ifdef NDEBUG
#  undef NDEBUG
#endif

#include <nacs-utils/timer.h>
#include <nacs-utils/base64.h>
#include <nacs-seq/pulser.h>

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace NaCs;

int main(int argc, char **argv)
{
    assert(argc >= 2);

    std::ifstream istm(argv[1]);
    std::string data(std::istreambuf_iterator<char>(istm), {});
    assert(Base64::validate((const uint8_t*)data.data(), data.size()));
    tic();
    auto code = Seq::Sequence::fromBase64((const uint8_t*)data.data(),
                                          data.size()).toByteCode(nullptr);
    printToc();
    size_t code_len;
    auto code2 = Seq::Sequence::fromBase64((const uint8_t*)data.data(),
                                           data.size()).toByteCode(&code_len, nullptr);

    if (argc >= 3) {
        std::ofstream ostm(argv[2]);
        ostm.write((const char*)&code[0], code.size());
    }

    assert(code_len == code.size());
    assert(memcmp(code2, &code[0], code_len) == 0);
    free(code2);

    // Seq::ByteCode::print(std::cout, code);

    std::cout << Seq::ByteCode::count(code) << std::endl;
    std::cout << code.size() << std::endl;

    return 0;
}
