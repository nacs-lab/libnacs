#include <nacs-utils/timer.h>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <iostream>
#include <condition_variable>
#include <thread>

static constexpr int N = 1000000;

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
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        std::lock_guard<Lock> locker(lock);
    }
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        std::unique_lock<Lock> locker(lock);
    }
    NaCs::printToc();

    std::unique_lock<Lock> unique_locker(lock, std::defer_lock);
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        unique_locker.lock();
        unique_locker.unlock();
    }
    NaCs::printToc();
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

static void
test_cond_var()
{
    std::condition_variable cv;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_one();
    }
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_all();
    }
    NaCs::printToc();
    std::thread thread([&] {
            std::mutex lock;
            std::unique_lock<std::mutex> locker(lock);
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
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_all();
    }
    NaCs::printToc();
}

int
main()
{
    std::cout << "Condition Variable" << std::endl;
    test_cond_var();

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
    std::cout << "DummyLock" << std::endl;
    test_lock<DummyLock>();

    std::cout << "function call" << std::endl;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f();
    }
    NaCs::printToc();

    std::cout << "inlined function call" << std::endl;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f_inline();
    }
    NaCs::printToc();

    std::cout << "pthread mutex" << std::endl;
    pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        pthread_mutex_lock(&lock);
        pthread_mutex_unlock(&lock);
    }
    NaCs::printToc();

    std::cout << "Simple spin lock" << std::endl;
    std::atomic_bool spin(false);
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        while (spin.exchange(true)) {
        }
        spin = false;
    }
    NaCs::printToc();
    return 0;
}
