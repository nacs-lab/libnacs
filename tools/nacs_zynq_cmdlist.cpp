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

#include "../lib/nacs-seq/zynq/cmdlist.h"
#include "../lib/nacs-utils/log.h"
#include "../lib/nacs-utils/errors.h"
#include "../lib/nacs-utils/streams.h"

#include <iostream>
#include <fstream>

using namespace NaCs;

int parse(int argc, char **argv)
{
    if (argc < 1) {
        Log::error("No input file specified.\n");
        return 1;
    }
    else if (argc < 2) {
        Log::error("No output file specified.\n");
        return 1;
    }
    else if (argc != 2) {
        Log::error("Wrong number of arguments.\n");
        return 1;
    }
    std::ifstream istm(argv[0]);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    // Hard code for now
    uint32_t ver = 2;

    uvector_ostream vstm;
    Seq::Zynq::SeqMetadata meta;
    try {
        meta = Seq::Zynq::CmdList::parse(vstm, istm, ver);
    }
    catch (const SyntaxError &err) {
        std::cerr << err;
        return 1;
    }
    std::ofstream ostm(argv[1], std::ios::binary);
    if (!ostm) {
        Log::error("Cannot open output file.\n");
        return 1;
    }
    ostm.write((char*)&ver, 4);
    auto v = vstm.get_buf();
    uint64_t len_ns = Seq::Zynq::CmdList::total_time(v.data(), v.size(), ver) * 10;
    ostm.write((char*)&len_ns, 8);
    ostm.write((char*)&meta.ttl_masks[0], 4);
    ostm.write((char*)v.data(), v.size());
    return 0;
}

int print(int argc, char **argv)
{
    std::ostream *stm;
    std::unique_ptr<std::ostream> ustm;
    if (argc < 1) {
        Log::error("No input file specified.\n");
        return 1;
    }
    else if (argc < 2) {
        stm = &std::cout;
    }
    else if (argc == 2) {
        ustm = std::make_unique<std::ofstream>(argv[1]);
        if (!ustm->good()) {
            Log::error("Cannot open output file.\n");
            return 1;
        }
        stm = ustm.get();
    }
    else {
        Log::error("Wrong number of arguments.\n");
        return 1;
    }
    std::ifstream istm(argv[0], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    std::string str(std::istreambuf_iterator<char>(istm), {});
    if (str.size() < 16) {
        Log::error("Cmd list too short.\n");
        return 1;
    }
    auto str_data = (const uint8_t*)str.data();
    auto str_sz = str.size();

    uint32_t ver;
    memcpy(&ver, str_data, 4);
    str_data += 4;
    str_sz -= 4;
    if (ver == 0 || ver > 2) {
        Log::error("Wrong cmd list file version.\n");
        return 1;
    }

    uint64_t len_ns;
    memcpy(&len_ns, str_data, 8);
    str_data += 8;
    str_sz -= 8;

    uint32_t ttl_mask;
    memcpy(&ttl_mask, str_data, 4);
    str_data += 4;
    str_sz -= 4;

    *stm << "# " << len_ns << " ns" << std::endl;
    Seq::Zynq::CmdList::print(*stm, str_data, str_sz, ttl_mask, ver);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        Log::error("No action specified.\n");
        return 1;
    }
    if (strcmp(argv[1], "parse") == 0) {
        return parse(argc - 2, argv + 2);
    }
    else if (strcmp(argv[1], "print") == 0) {
        return print(argc - 2, argv + 2);
    }
    else {
        Log::error("Unknown action: %s.\n", argv[1]);
        return 1;
    }
}
