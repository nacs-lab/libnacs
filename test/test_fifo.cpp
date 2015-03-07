#include <nacs-utils/timer.h>
#include <nacs-utils/container.h>

#include <mutex>
#include <atomic>
#include <iostream>
#include <thread>
#include <chrono>
#include <vector>

#include <assert.h>

static constexpr int N = 10000000;

using namespace NaCs;
using namespace std::literals;

template<typename Lock, int init_size=256>
static inline void
test_fifo()
{
    FIFO<int> fifo(init_size);
    assert(fifo.size() == 0);
    assert(fifo.capacity() == init_size);
    assert(fifo.spaceLeft() == init_size);

    tic();
    std::thread t1([&] {
            std::this_thread::yield();
            for (int i = 0;i < N;i++) {
                while (!fifo.tryPush(i)) {
                    std::this_thread::yield();
                }
            }
        });
    for (int i = 0;i < N;i++) {
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

    assert(fifo.capacity() == init_size);

    tic();
    std::thread t2([&] {
            std::this_thread::yield();
            for (int i = 0;i < N;i++) {
                fifo.push(i);
            }
        });
    int i = 0;
    while (i < N) {
        int size = fifo.size();
        if (!size) {
            std::this_thread::yield();
            continue;
        }
        for (int j = 0;j < size;j++) {
            assert(fifo.pop() == i++);
        }
    }
    t2.join();
    auto write_time2 = toc();

    std::cout << "Average time per write: "
              << double(write_time2) / double(N) / 1e3 << " us"
              << std::endl;
    std::cout << "Total size: " << fifo.capacity() << std::endl;
}

int
main()
{
    test_fifo<SpinLock, 256>();
    test_fifo<SpinLock, 8192>();
    test_fifo<std::mutex, 256>();
    test_fifo<std::mutex, 8192>();
    return 0;
}
