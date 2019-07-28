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

#include <../../lib/utils/thread.h>
#include <../../lib/utils/timer.h>
#include <../../lib/utils/utils.h>

#include <thread>
#include <iostream>

#include <assert.h>
#include <stdint.h>

using namespace NaCs;

constexpr uint32_t pipe_sz = 4096 * 128;
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
        uint32_t sz;
        uint64_t *ptr;
        while (true) {
            ptr = pipe.get_write_ptr(&sz);
            if (sz <= 7 && sz > 0)
                pipe.sync_writer();
            sz &= ~(uint32_t)7;
            if (sz != 0)
                break;
            CPU::pause();
        }
        for (uint32_t j = 0; j < sz; j++) {
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
        uint32_t sz;
        const uint64_t *ptr;
        while (true) {
            ptr = pipe.get_read_ptr(&sz);
            if (sz <= 7 && sz > 0)
                pipe.sync_reader();
            sz &= ~(uint32_t)7;
            if (sz != 0)
                break;
            CPU::pause();
        }
        for (uint32_t j = 0; j < sz; j++) {
            assert(ptr[j] == i);
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

int main()
{
    static constexpr size_t n = 4096ul * 4096ul * 1000ul;
    Pipe pipe;
    std::thread writer{[&] {
            push_data(pipe, n);
        }};
    auto t1 = getTime();
    check_data(pipe, n);
    auto t2 = getTime();
    writer.join();
    std::cout << "Time per data = " << double(t2 - t1) / n << " (ns)" << std::endl;
    return 0;
}
