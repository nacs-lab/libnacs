/*************************************************************************
 *   Copyright (c) 2016 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#define CATCH_CONFIG_MAIN

#include <../../lib/nacs-utils/thread.h>
#include <../../lib/nacs-utils/timer.h>
#include <../../lib/nacs-utils/utils.h>

#include <thread>
#include <iostream>

#include <catch2/catch.hpp>
#include <stdint.h>

using namespace NaCs;

constexpr size_t pipe_sz = 4096 * 128;
struct Pipe : DataPipe<uint64_t> {
    Pipe()
        : DataPipe<uint64_t>(buff, pipe_sz)
    {}
private:
    uint64_t buff[pipe_sz];
};

void push_data(Pipe &pipe, size_t n)
{
    size_t i = 0;
    while (true) {
        size_t sz;
        uint64_t *ptr;
        while (true) {
            ptr = pipe.get_write_ptr(&sz);
            if (sz <= 7 && sz > 0)
                pipe.sync_writer();
            sz &= ~(size_t)7;
            if (sz != 0)
                break;
            CPU::pause();
        }
        for (size_t j = 0; j < sz; j++) {
            ptr[j] = i;
            i++;
        }
        // This can produce more data then n
        pipe.wrote_size(sz);
        CPU::wake();
        if (i >= n) {
            break;
        }
    }
}

void check_data(Pipe &pipe, size_t n)
{
    size_t i = 0;
    while (true) {
        size_t sz;
        const uint64_t *ptr;
        while (true) {
            ptr = pipe.get_read_ptr(&sz);
            if (sz <= 7 && sz > 0)
                pipe.sync_reader();
            sz &= ~(size_t)7;
            if (sz != 0)
                break;
            CPU::pause();
        }
        for (size_t j = 0; j < sz; j++) {
            REQUIRE(ptr[j] == i);
            i++;
        }
        // This can produce more data then n
        pipe.read_size(sz);
        CPU::wake();
        if (i >= n) {
            break;
        }
    }
}

TEST_CASE("DataPipe") {
#ifdef  __OPTIMIZE__
    static constexpr size_t n = 4096ul * 4096ul * 100ul;
#else
    static constexpr size_t n = 4096ul * 4096ul * 10ul;
#endif
    std::unique_ptr<Pipe> pipe(new Pipe);
    std::thread writer{[&] {
            push_data(*pipe, n);
        }};
    auto t1 = getTime();
    check_data(*pipe, n);
    auto t2 = getTime();
    writer.join();
    std::cout << "Time per data = " << double(t2 - t1) / n << " (ns)" << std::endl;
}
