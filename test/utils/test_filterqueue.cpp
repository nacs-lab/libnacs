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

#include "../../lib/utils/container.h"

#include <assert.h>

#include <thread>
#include <random>

using namespace NaCs;

std::default_random_engine rng{std::random_device()()};
std::uniform_int_distribution<int> dist{0, 10};

void test(int n)
{
    FilterQueue<std::atomic<int>> queue;
    {
        auto peekres = queue.peek();
        assert(!peekres.first);
        assert(!peekres.second);
        assert(!queue.pop());
    }
    std::thread tfilter([&] {
            for (int i = 0; i < n; i++) {
                std::atomic<int> *p;
                while ((p = queue.get_filter()) == nullptr)
                    CPU::pause();
                assert(p == queue.get_filter());
                assert(p->load(std::memory_order_relaxed) == i);
                p->store(i * 2 + 1, std::memory_order_relaxed);
                queue.forward_filter();
                CPU::wake();
            }
        });
    int nwrite = 0;
    int nread = 0;
    auto check_it = [&] () {
        int i0 = nread;
        for (auto p: queue) {
            int i = i0++;
            int v = p->load(std::memory_order_relaxed);
            assert(v == (i * 2 + 1) || v == i);
        }
    };
    while (nread < n) {
        bool do_read = true;
        int det = nwrite - nread;
        if (det == 0) {
            do_read = false;
        }
        else if (nwrite < n && det <= 10) {
            do_read = dist(rng) < det;
        }
        check_it();
        if (do_read) {
            std::atomic<int> *p;
            auto exp = queue.peek();
            assert(exp.first);
            if (exp.second) {
                p = queue.pop();
                assert(p);
            }
            else {
                while ((p = queue.pop()) == nullptr) {
                    CPU::pause();
                }
            }
            assert(exp.first == p);
            assert(p->load(std::memory_order_relaxed) == 2 * nread + 1);
            delete p;
            nread++;
        }
        else {
            if (det)
                assert(queue.peek().first);
            queue.push(new std::atomic<int>(nwrite));
            CPU::wake();
            nwrite++;
        }
    }
    tfilter.join();
}

int main()
{
    for (int i = 0; i < 10; i++)
        test(10000);
    return 0;
}
