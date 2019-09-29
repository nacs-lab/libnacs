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

#include "../lib/nacs-seq/zynq/seq.h"

#include <iostream>
#include <fstream>

using namespace NaCs;

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cout << "ERROR: wrong number of arguments." << std::endl;
        return 1;
    }

    std::ifstream istm(argv[1], std::ios::binary);
    istm.seekg(0, std::ios::end);
    auto filesize = (size_t)istm.tellg();
    istm.seekg(0, std::ios::beg);
    std::vector<uint32_t> data(filesize / 4);
    istm.read((char*)data.data(), filesize);
    Seq::Zynq::Sequence::dumpBinary(std::cout, data.data(), data.size());

    return 0;
}
