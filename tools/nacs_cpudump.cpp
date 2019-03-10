/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../lib/utils/processor.h"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace NaCs;

void dump(const CPUInfo &info)
{
    info.dump_llvm(std::cout);
    std::cout << "ID: ";
    info.dump(std::cout);
    std::cout << std::endl;
}

int main(int argc, char **argv)
{
    if (argc > 3)
        std::cerr << "ERROR: wrong number of arguments." << std::endl;
    if (argc == 3) {
        dump(*CPUInfo::create(argv[1], argv[2]));
        return 0;
    }
    if (argc == 2) {
        if (strcmp(argv[1], "host") == 0 || strcmp(argv[1], "native") == 0) {
            dump(CPUInfo::get_host());
            return 0;
        }
        dump(*CPUInfo::create(argv[1], "generic"));
        return 0;
    }
    dump(CPUInfo::get_host());

    return 0;
}
