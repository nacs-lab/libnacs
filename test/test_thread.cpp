#include <nacs-utils/timer.h>
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
            lock.lock();
            lock.unlock();
        }
        tocPerCycle();
        NaCs::tic();
        for (int i = 0;i < N;i++) {
            std::lock_guard<Lock> locker(lock);
        }
        tocPerCycle();
        NaCs::tic();
        for (int i = 0;i < N;i++) {
            std::unique_lock<Lock> locker(lock);
        }
        tocPerCycle();

        std::unique_lock<Lock> unique_locker(lock, std::defer_lock);
        NaCs::tic();
        for (int i = 0;i < N;i++) {
            unique_locker.lock();
            unique_locker.unlock();
        }
        tocPerCycle();
    }
};

template<typename Func, typename Lock, int n>
struct n_runner {
    Func &m_func;
    Lock &m_lock;
    n_runner(Func &func, Lock &lock)
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

template<typename Lock, int n=N>
struct thread_tester {
    template<typename Func>
    inline void
    operator()(Func &&func, Lock &lock)
    {
        NaCs::tic();
        n_runner<Func, Lock, n> runner(func, lock);
        std::thread t1(runner);
        std::thread t2(runner);
        t1.join();
        t2.join();
        tocPerCycle(n * 2);
    }
};

template<int n>
struct thread_tester<DummyLock, n> {
    template<typename Func>
    inline void
    operator()(Func &&func, DummyLock &lock)
    {
        NaCs::tic();
        n_runner<Func, DummyLock, n> runner(func, lock);
        std::thread(runner).join();
        tocPerCycle(n);
    }
};

template<typename Lock>
struct test_thread_lock {
    inline void
    operator()()
    {
        Lock lock;
        std::atomic<int> ai;
        volatile double f = 1.2;
        thread_tester<Lock>()([&] {
                ++ai;
                f = std::cos(f);
            }, lock);
        thread_tester<Lock>()([&] {
                for (int j = 0;j < 10;j++) {
                    f = std::cos(f);
                }
            }, lock);
        volatile int *volatile ptr = nullptr;
        thread_tester<Lock, N / 100>()([&] {
                if (ptr) {
                    delete ptr;
                    ptr = nullptr;
                } else {
                    ptr = new int[1024 * 1024];
                }
            }, lock);
        thread_tester<Lock, N / 100>()([&] {
                std::this_thread::sleep_for(1us);
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
    PthreadMutex()
        : m_lock(PTHREAD_MUTEX_INITIALIZER)
    {}
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
    Tester<SpinLock<true> >()();
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
