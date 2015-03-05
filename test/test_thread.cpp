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
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        cv.notify_all();
    }
    NaCs::printToc();
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
    return 0;
}
