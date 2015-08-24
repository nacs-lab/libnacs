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
#include <nacs-utils/thread.h>

#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <iostream>
#include <condition_variable>
#include <thread>
#include <cmath>
#include <chrono>

using namespace std::literals;

static constexpr int N = 1000000;

static inline void
tocPerCycle(int n=N)
{
    std::cout << "Time: " << double(NaCs::toc()) / double(n) / 1e3
              << " us" << std::endl;
}

class DummyLock;

template<typename Lock>
struct test_lock {
    inline void
    operator()()
    {
        Lock lock;
        NaCs::tic();
        for (int i = 0;i < N;i++) {
            std::lock_guard<Lock> locker(lock);
        }
        tocPerCycle();
    }
};

template<typename Func, typename Lock, int n>
class RepeatRunner {
    Func &m_func;
    Lock &m_lock;
public:
    RepeatRunner(Func &func, Lock &lock)
        : m_func(func),
          m_lock(lock)
    {}
    inline void
    operator()()
    {
        for (int i = 0;i < n;i++) {
            std::lock_guard<Lock> locker(m_lock);
            m_func();
        }
    }
};

template<int n=N, int _nthread=2, typename Lock, typename Func>
static inline void
test_n_threads(Func &&func, Lock &lock)
{
    NaCs::tic();
    RepeatRunner<Func, Lock, n> runner(func, lock);
    constexpr static int nthread =
        std::is_same<DummyLock, Lock>() ? 1 : _nthread;
    std::thread threads[nthread];
    for (auto &t: threads) {
        t = std::thread(runner);
    }
    for (auto &t: threads) {
        t.join();
    }
    tocPerCycle(n * nthread);
}

template<typename Lock>
struct test_thread_lock {
    inline void
    operator()()
    {
        Lock lock;
        std::atomic<int> ai;
        volatile double f = 1.2;
        test_n_threads([&] {
                ++ai;
                f = std::cos(f);
            }, lock);
        test_n_threads([&] {
                for (int j = 0;j < 10;j++) {
                    f = std::cos(f);
                }
            }, lock);
        volatile int *volatile ptr = nullptr;
        test_n_threads<N / 100>([&] {
                if (ptr) {
                    delete ptr;
                    ptr = nullptr;
                } else {
                    ptr = new int[1024 * 1024];
                }
            }, lock);
        test_n_threads<N / 100>([&] {
                std::this_thread::sleep_for(1us);
            }, lock);
        test_n_threads<N / 100, 6>([&] {
                for (int j = 0;j < 100;j++) {
                    f = std::cos(f);
                    f = f + 1;
                    f = std::sin(f);
                }
            }, lock);
    }
};

static void __attribute__((noinline))
f()
{
    asm("");
}

static inline void
f_inline()
{
    asm("");
}

class DummyLock {
public:
    inline void
    lock() const
    {}
    inline void
    unlock() const
    {}
};

template<bool yield=false>
class SpinLock {
    std::atomic_bool m_spin;
public:
    SpinLock()
        : m_spin(false)
    {}
    inline void
    lock()
    {
        while (m_spin.exchange(true)) {
            if (yield) {
                std::this_thread::yield();
            }
        }
    }
    inline void
    unlock()
    {
        m_spin = false;
    }
};

class PthreadMutex {
    pthread_mutex_t m_lock;
public:
    inline
    PthreadMutex()
        : m_lock()
    {
        pthread_mutex_init(&m_lock, nullptr);
    }
    inline
    ~PthreadMutex()
    {
        pthread_mutex_destroy(&m_lock);
    }
    inline void
    lock()
    {
        pthread_mutex_lock(&m_lock);
    }
    inline void
    unlock()
    {
        pthread_mutex_unlock(&m_lock);
    }
};

template<typename Lock=std::mutex>
struct uniqueMutex: public std::unique_lock<Lock> {
    uniqueMutex() :
        std::unique_lock<Lock>(),
        m_lock()
    {
        std::unique_lock<Lock>::operator=(std::unique_lock<Lock>(m_lock));
    }
private:
    Lock m_lock;
};

template<typename CondVar, typename Locker>
static void
test_cond_var()
{
    CondVar cv;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_one();
    }
    tocPerCycle();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_all();
    }
    tocPerCycle();
    std::thread thread([&] {
            Locker locker;
            for (int i = 0;i < N * 2;i++) {
                cv.wait(locker);
                std::this_thread::yield();
            }
        });
    thread.detach();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_one();
    }
    tocPerCycle();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_all();
    }
    tocPerCycle();
}

template<template<typename...> class Tester>
static inline void
locks_tester()
{
    std::cout << "std::mutex" << std::endl;
    Tester<std::mutex>()();
    std::cout << "PthreadMutex" << std::endl;
    Tester<PthreadMutex>()();
    std::cout << "SpinLock<false>" << std::endl;
    Tester<SpinLock<false> >()();
    std::cout << "SpinLock<true>" << std::endl;
    Tester<NaCs::SpinLock>()();
    std::cout << "DummyLock" << std::endl;
    Tester<DummyLock>()();
}

int
main()
{
    std::cout << "Condition Variable, unique_lock<mutex>" << std::endl;
    test_cond_var<std::condition_variable, uniqueMutex<> >();

    locks_tester<test_lock>();
    locks_tester<test_thread_lock>();

    std::cout << "function call" << std::endl;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f();
    }
    tocPerCycle();

    std::cout << "inlined function call" << std::endl;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f_inline();
    }
    tocPerCycle();
    return 0;
}
