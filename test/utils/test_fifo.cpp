/*************************************************************************
 *   Copyright (c) 2015 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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

#include <nacs-utils/timer.h>
#include <nacs-utils/container.h>

#include <mutex>
#include <atomic>
#include <iostream>
#include <thread>
#include <chrono>
#include <vector>

#include <assert.h>

static constexpr size_t N = 1000000;

using namespace NaCs;
using namespace std::literals;

template<size_t sz, typename F>
static inline void
fifo_list_try_pusher(F &fifo)
{
    std::this_thread::yield();
    size_t buff[sz];
    for (size_t i = 0;i < N;) {
        auto end = min(sz, N - i);
        for (size_t j = 0;j < end;j++) {
            buff[j] = i + j;
        }
        size_t pushed;
        while ((pushed = fifo.tryPush(buff, 1)) == 0) {
            std::this_thread::yield();
        }
        i += pushed;
    }
}

template<size_t sz, typename F>
static inline void
fifo_list_pusher(F &fifo)
{
    std::this_thread::yield();
    size_t buff[sz];
    for (size_t i = 0;i < N;) {
        auto end = min(sz, N - i);
        for (size_t j = 0;j < end;j++) {
            buff[j] = i + j;
        }
        fifo.push(buff, end);
        i += end;
    }
}

template<size_t sz, typename F>
static inline void
fifo_poper(F &fifo)
{
    size_t buff[sz];
    for (size_t i = 0;i < N;) {
        size_t size;
        while (!(size = fifo.pop(buff, sz)))
            std::this_thread::yield();
        for (size_t j = 0;j < size;j++) {
            assert(buff[j] == i + j);
        }
        i += size;
    }
}

template<size_t sz, typename F>
static inline void
fifo_try_poper(F &fifo)
{
    size_t buff[sz];
    for (size_t i = 0;i < N;) {
        size_t size;
        while (!(size = fifo.tryPop(buff, sz)))
            std::this_thread::yield();
        for (size_t j = 0;j < size;j++) {
            assert(buff[j] == i + j);
        }
        i += size;
    }
}

template<typename Lock, size_t init_size=256>
static inline void
test_fifo()
{
    FIFO<size_t, Lock> fifo(init_size);
    assert(fifo.size() == 0);
    assert(fifo.capacity() == init_size);
    assert(fifo.spaceLeft() == init_size);

    tic();
    std::thread t1([&] {
            std::this_thread::yield();
            for (size_t i = 0;i < N;i++) {
                while (!fifo.tryPush(i)) {
                    std::this_thread::yield();
                }
            }
        });
    for (size_t i = 0;i < N;i++) {
        while (!fifo.size()) {
            std::this_thread::yield();
        }
        auto res = fifo.pop();
        assert(res == i);
    }
    t1.join();
    auto write_time = toc();

    std::cout << "Average time per write: "
              << double(write_time) / double(N) / 1e3 << " us"
              << std::endl;

    tic();
    std::thread t3([&] {
            fifo_list_try_pusher<1>(fifo);
        });
    fifo_poper<3>(fifo);
    t3.join();
    auto write_time3 = toc();

    std::cout << "Average time per write: "
              << double(write_time3) / double(N) / 1e3 << " us"
              << std::endl;

    tic();
    std::thread t4([&] {
            fifo_list_try_pusher<3>(fifo);
        });
    fifo_poper<1>(fifo);
    t4.join();
    auto write_time4 = toc();

    std::cout << "Average time per write: "
              << double(write_time4) / double(N) / 1e3 << " us"
              << std::endl;

    assert(fifo.capacity() == init_size);

    tic();
    std::thread t2([&] {
            std::this_thread::yield();
            for (size_t i = 0;i < N;i++) {
                fifo.push(i);
            }
        });
    size_t i = 0;
    while (i < N) {
        size_t size = fifo.size();
        if (!size) {
            std::this_thread::yield();
            continue;
        }
        for (size_t j = 0;j < size;j++) {
            assert(fifo.pop() == i++);
        }
    }
    t2.join();
    auto write_time2 = toc();

    std::cout << "Average time per write: "
              << double(write_time2) / double(N) / 1e3 << " us"
              << std::endl;
    std::cout << "Total size: " << fifo.capacity() << std::endl;

    tic();
    std::thread t5([&] {
            fifo_list_pusher<1>(fifo);
        });
    fifo_try_poper<3>(fifo);
    t5.join();
    auto write_time5 = toc();

    std::cout << "Average time per write: "
              << double(write_time5) / double(N) / 1e3 << " us"
              << std::endl;
    std::cout << "Total size: " << fifo.capacity() << std::endl;

    tic();
    std::thread t6([&] {
            fifo_list_pusher<3>(fifo);
        });
    fifo_try_poper<1>(fifo);
    t6.join();
    auto write_time6 = toc();
    std::cout << "Average time per write: "
              << double(write_time6) / double(N) / 1e3 << " us"
              << std::endl;
    std::cout << "Total size: " << fifo.capacity() << std::endl;
}

int
main()
{
    std::cout << "<SpinLock, 256>: " << std::endl;
    test_fifo<SpinLock, 256>();
    std::cout << "<SpinLock, 8192>: " << std::endl;
    test_fifo<SpinLock, 8192>();
    std::cout << "<std::mutex, 256>: " << std::endl;
    test_fifo<std::mutex, 256>();
    std::cout << "<std::mutex, 8192>: " << std::endl;
    test_fifo<std::mutex, 8192>();
    return 0;
}
