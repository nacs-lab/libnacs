#include <nacs-utils/timer.h>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <iostream>
#include <condition_variable>
#include <thread>

static constexpr int N = 1000000;

static inline void
tocPerCycle()
{
    std::cout << "Time: " << double(NaCs::toc()) / double(N) / 1e3
              << " us" << std::endl;
}

template<typename Lock>
static inline void
test_lock()
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

int
main()
{
    std::cout << "Condition Variable, unique_lock<mutex>" << std::endl;
    test_cond_var<std::condition_variable, uniqueMutex<> >();
    std::cout << "Condition Variable Any, unique_lock<mutex>" << std::endl;
    test_cond_var<std::condition_variable_any, uniqueMutex<> >();
    std::cout << "Condition Variable Any, mutex" << std::endl;
    test_cond_var<std::condition_variable_any, std::mutex>();
    std::cout << "Condition Variable Any, unique_lock<DummyLock>" << std::endl;
    test_cond_var<std::condition_variable_any, uniqueMutex<DummyLock> >();
    std::cout << "Condition Variable Any, DummyLock" << std::endl;
    test_cond_var<std::condition_variable_any, DummyLock>();

    std::cout << "std::mutex" << std::endl;
    test_lock<std::mutex>();
    std::cout << "std::timed_mutex" << std::endl;
    test_lock<std::timed_mutex>();
    std::cout << "std::recursive_mutex" << std::endl;
    test_lock<std::recursive_mutex>();
    std::cout << "std::recursive_timed_mutex" << std::endl;
    test_lock<std::recursive_timed_mutex>();
    std::cout << "std::shared_timed_mutex" << std::endl;
    test_lock<std::shared_timed_mutex>();

    std::cout << "pthread mutex" << std::endl;
    pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        pthread_mutex_lock(&lock);
        pthread_mutex_unlock(&lock);
    }
    tocPerCycle();

    std::cout << "Simple spin lock" << std::endl;
    std::atomic_bool spin(false);
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        while (spin.exchange(true)) {
        }
        spin = false;
    }
    tocPerCycle();

    std::cout << "DummyLock" << std::endl;
    test_lock<DummyLock>();

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
