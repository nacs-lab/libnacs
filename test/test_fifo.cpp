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

int
main()
{
    FIFO<int> fifo(256);
    assert(fifo.size() == 0);
    assert(fifo.capacity() == 256);
    assert(fifo.spaceLeft() == 256);

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

    assert(fifo.capacity() == 256);

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

    return 0;
}
