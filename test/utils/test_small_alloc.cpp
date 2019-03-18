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

#include "test_helpers.h"

#include "../../lib/utils/mem.h"
#include "../../lib/utils/timer.h"

#include <assert.h>

#include <random>

using namespace NaCs;

template<size_t nstatic>
class AllocatorTester {
    static constexpr auto padding = 8192;
    static constexpr auto ncounters = nstatic ? nstatic : 5;
    using Eltype = UnmovableIncrement;
    using Allocator = SmallAllocator<Eltype,nstatic>;
    struct Allocation {
        Eltype *obj;
        int pool_idx;
        int counter_idx;
    };
public:
    AllocatorTester()
        : m_palloc((char*)malloc(padding * 2 + sizeof(Allocator))),
          m_pend(m_palloc + padding * 2 + sizeof(Allocator)),
          m_allocator(new ((char*)m_palloc + padding) Allocator),
          m_base((char*)m_allocator->alloc()),
          m_rng(std::random_device()())
    {
        memset(m_allocated, 0, sizeof(m_allocated));
        m_allocated[nstatic] = true; // dummy
        memset(m_counters, 0, sizeof(m_counters));
        memset(m_counters2, 0, sizeof(m_counters2));
        m_allocator->free((Eltype*)m_base);
        if (nstatic) {
            assert(m_palloc + padding <= m_base);
            assert(m_pend - padding > m_base);
        }
    }
    void test(size_t nrun)
    {
        for (size_t i = 0; i < nrun; i++)
            test_iteration();
        while (!m_allocations.empty()) {
            free();
        }
    }

private:
    int check_alloc(char *p)
    {
        if (m_palloc > p || p >= m_pend) {
            // Check that all slots are full
            for (auto b: m_allocated)
                assert(b);
            return -1;
        }
        assert(nstatic);
        assert(p >= m_base);
        assert(p < m_base + sizeof(Eltype) * nstatic);
        auto offset = p - m_base;
        assert(offset % sizeof(Eltype) == 0);
        auto idx = offset / sizeof(Eltype);
        assert(!m_allocated[idx]);
        m_allocated[idx] = true;
        return int(idx);
    }
    void alloc()
    {
        auto counter_idx = m_counter_dist(m_rng);
        auto obj = m_allocator->alloc(counter_idx < 0 ? nullptr : &m_counters[counter_idx]);
        auto pool_idx = check_alloc((char*)obj);
        m_allocations.push_back(Allocation{obj, pool_idx, counter_idx});
    }
    void check_counters()
    {
        for (size_t i = 0; i < nstatic; i++) {
            assert(m_counters[i] == m_counters2[i]);
        }
    }
    void free()
    {
        assert(!m_allocations.empty());
        check_counters();
        auto to_free =
            std::uniform_int_distribution<size_t>{0, m_allocations.size() - 1}(m_rng);
        auto alloc = m_allocations[to_free];
        if (alloc.pool_idx >= 0) {
            assert(m_base);
            assert((char*)alloc.obj == m_base + sizeof(Eltype) * alloc.pool_idx);
            assert(m_allocated[alloc.pool_idx]);
            m_allocated[alloc.pool_idx] = false;
        }
        m_allocator->free(alloc.obj);
        if (alloc.counter_idx >= 0) {
            assert(m_counters[alloc.counter_idx] == m_counters2[alloc.counter_idx] + 1);
            m_counters2[alloc.counter_idx]++;
        }
        check_counters();
        m_allocations.erase(m_allocations.begin() + to_free);
    }
    void test_iteration()
    {
        auto action = m_action_dist(m_rng);
        if (action >= m_allocations.size()) {
            alloc();
        }
        else {
            free();
        }
    }
    // Allocate the SmallAllocator with 8kB of padding in the front and end
    // in order to detect buffer overflow error.
    char *const m_palloc;
    char *const m_pend;
    Allocator *const m_allocator;
    char *const m_base;
    bool m_allocated[nstatic + 1]; // To suppress the warning when `nstatic == 0`
    int m_counters[ncounters];
    int m_counters2[ncounters];
    std::default_random_engine m_rng;
    std::uniform_int_distribution<int> m_counter_dist{-1, ncounters - 1};
    std::uniform_int_distribution<size_t> m_action_dist{0, (nstatic ? nstatic * 2 : 10) - 1};
    std::vector<Allocation> m_allocations;
};

template<size_t nstatic>
NACS_NOINLINE void benchmark(size_t nlive, size_t ncycle)
{
    std::vector<int*> pointers(nlive, nullptr);
    SmallAllocator<int,nstatic> allocator;
    Timer timer;
    for (size_t i = 0; i < ncycle; i++) {
        size_t idx = i % nlive;
        if (auto old = pointers[idx])
            allocator.free(old);
        pointers[idx] = allocator.alloc(1);
        asm volatile ("" :: "r"(pointers[idx]) : "memory");
        asm volatile ("" :: "r"(&pointers[0]) : "memory");
    }
    timer.print();
    for (auto ptr: pointers) {
        if (ptr) {
            allocator.free(ptr);
        }
    }
}

int main()
{
    AllocatorTester<0> tester0;
    tester0.test(100000);
    AllocatorTester<10> tester;
    tester.test(100000);
    AllocatorTester<32> tester2;
    tester2.test(100000);
    AllocatorTester<1023> tester3;
    tester3.test(10000000);

    benchmark<0>(10, 10000000);
    benchmark<4>(10, 10000000);
    benchmark<8>(10, 10000000);
    benchmark<16>(10, 10000000);

    benchmark<0>(100, 10000000);
    benchmark<64>(100, 10000000);
    benchmark<128>(100, 10000000);

    return 0;
}
