#include <nacs-utils/timer.h>
#include <mutex>
#include <shared_mutex>
#include <atomic>

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

int
main()
{
    test_lock<std::mutex>();
    test_lock<std::timed_mutex>();
    test_lock<std::recursive_mutex>();
    test_lock<std::recursive_timed_mutex>();
    test_lock<std::shared_timed_mutex>();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f();
    }
    NaCs::printToc();
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        f_inline();
    }
    NaCs::printToc();

    pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    NaCs::tic();
    for (int i = 0;i < N;i++) {
        pthread_mutex_lock(&lock);
        pthread_mutex_unlock(&lock);
    }
    NaCs::printToc();

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
