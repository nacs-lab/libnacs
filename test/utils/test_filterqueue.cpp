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

#include <nacs-utils/container.h>

#include <assert.h>

#include <thread>
#include <random>

using namespace NaCs;

std::default_random_engine rng{std::random_device()()};
std::uniform_int_distribution<int> dist{0, 10};

void test(int n)
{
    FilterQueue<int> queue;
    std::thread tfilter([&] {
            for (int i = 0; i < n; i++) {
                int *p;
                while ((p = queue.get_filter()) == nullptr)
                    std::this_thread::yield(); // pause
                assert(p == queue.get_filter());
                assert(*p == i);
                *p = i * 2 + 1;
                queue.forward_filter();
            }
        });
    int nwrite = 0;
    int nread = 0;
    while (nread < n) {
        bool do_read = true;
        int det = nwrite - nread;
        if (det == 0) {
            do_read = false;
        }
        else if (nwrite < n && det <= 10) {
            do_read = dist(rng) < det;
        }
        if (do_read) {
            int *p;
            while ((p = queue.pop()) == nullptr)
                std::this_thread::yield(); // pause
            assert(*p == 2 * nread + 1);
            delete p;
            nread++;
        }
        else {
            queue.push(new int(nwrite));
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
