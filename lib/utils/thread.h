/*************************************************************************
 *   Copyright (c) 2015 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef _NACS_UTILS_THREAD_H_
#define _NACS_UTILS_THREAD_H_

#include "utils.h"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <thread>
#include <mutex>

namespace NaCs {

template<bool real, typename Lock>
struct CondLock {
    CondLock(Lock&)
    {}
};

template<typename Lock>
struct CondLock<true, Lock> : std::lock_guard<Lock> {
    using std::lock_guard<Lock>::lock_guard;
};

/**
 * A straightforward implementation of a spin lock
 *
 * According to the benchmark (see test/test_thread.cpp) this lock is
 * faster than the stdlib lock (pthread_mutex_t or std::mutex) both without
 * and with certain level of lock contention and on both ARM and x86_64
 * (Haswell).
 *
 * As one would expect, the spin lock can slow down other (especially CPU
 * intensive) tasks running on the system. Yielding the CPU when failed to
 * aquire the lock improves the situation by a lot but this effect can still
 * be seen when the total busy threads is more than twice the number of
 * hardware threads (also see the benchmark). Therefore, it is important to
 * make sure to only use it only for performance critical lock and minimize
 * the critical region protected by the lock.
 *
 * This class should meet the requirements for a
 * [Mutex](http://en.cppreference.com/w/cpp/concept/Mutex) so it should
 * be easy to swap it with another implantation in the future.
 */
class SpinLock {
    std::atomic_bool m_spin;
public:
    SpinLock()
        : m_spin(false)
    {}
    inline bool
    try_lock()
    {
        return !m_spin.exchange(true, std::memory_order_acquire);
    }
    inline void
    lock()
    {
        while (!try_lock()) {
            CPU::pause();
        }
    }
    inline void
    unlock()
    {
        m_spin.store(false, std::memory_order_release);
        CPU::wake();
    }
};

/**
 * A high throughput pipe designed to pass a stream of data between two thread,
 * one producer and one consumer.
 * In the expected usecase, both read and write blocking must be very short so
 * no locking is implemented and both the producer and the consumer are expected to busy wait
 * before the blocking is cleared.
 *
 * Synchronization are done with the two atomic variables `m_read_ptr` and `m_write_ptr`.
 * Since the reader and writer are expected to read and write in small blocks
 * there will be a lot of updates to these variables.
 * Caches are added for both the reader and the writer to use their local copy if the
 * upper/lower bound provided by it is good enough.
 * Fields are aligned to cache boundary to make sure each atomic variables/caches
 * are in their own cacheline.
 */
template<typename T>
class DataPipe {
public:
    DataPipe(T *buff, size_t buff_sz, size_t max_block_sz=4096 / sizeof(T))
        : m_reader_cache(buff, buff_sz, max_block_sz),
          m_writer_cache(buff, buff_sz, max_block_sz)
    {
        if ((buff_sz & (buff_sz - 1)) != 0) {
            throw std::invalid_argument("Buffer size must be a power of 2");
        }
    }
    // Interface for the consumer of the data
    // This function does thte same thing as `get_write_buff`.
    // Prefered to be called by the reader.
    inline const T *get_read_buff(size_t *sz=nullptr) const
    {
        if (sz)
            *sz = m_reader_cache.buff_sz;
        return m_reader_cache.buff;
    }
    inline const T *get_read_ptr(size_t *sz)
    {
        auto buff_sz = m_reader_cache.buff_sz;
        auto read_p = m_reader_cache.read;
        if (read_p == m_reader_cache.write)
            sync_reader();
        auto write_p = m_reader_cache.write;
        *sz = std::min(get_avail_sz(read_p, write_p, buff_sz),
                       m_reader_cache.max_block_sz);
        return &m_reader_cache.buff[read_p];
    }
    inline void read_size(size_t sz)
    {
        auto buff_sz = m_reader_cache.buff_sz;
        m_reader_cache.read = (m_reader_cache.read + sz) & (buff_sz - 1);
        m_read_ptr.store(m_reader_cache.read, std::memory_order_release);
    }
    // Interface for the generator of the data
    // This function does thte same thing as `get_read_buff`.
    // Prefered to be called by the writer.
    inline const T *get_write_buff(size_t *sz=nullptr) const
    {
        if (sz)
            *sz = m_writer_cache.buff_sz;
        return m_writer_cache.buff;
    }
    inline T *get_write_ptr(size_t *sz)
    {
        auto buff_sz = m_writer_cache.buff_sz;
        auto write_p = m_writer_cache.write;
        if (write_p == ((m_writer_cache.read - 1) & (buff_sz - 1)))
            sync_writer();
        auto read_p = (m_writer_cache.read - 1) & (buff_sz - 1);
        *sz = std::min(get_avail_sz(write_p, read_p, buff_sz),
                       m_writer_cache.max_block_sz);
        return &m_writer_cache.buff[write_p];
    }
    inline void wrote_size(size_t sz)
    {
        auto buff_sz = m_writer_cache.buff_sz;
        m_writer_cache.write = (m_writer_cache.write + sz) & (buff_sz - 1);
        m_write_ptr.store(m_writer_cache.write, std::memory_order_release);
    }
    // Return `true` if the reader is no more than `sz` behind the writer
    // Useful to limit buffering to ensure that the writer can change setting with low
    // latency.
    inline bool check_reader(size_t sz)
    {
        auto buff_sz = m_writer_cache.buff_sz;
        auto write_p = m_writer_cache.write;
        auto read_p = m_writer_cache.read;
        if ((write_p - read_p) & (buff_sz - 1) <= sz)
            return true;
        sync_writer();
        read_p = m_writer_cache.read;
        return (write_p - read_p) & (buff_sz - 1) <= sz;
    }

    // The reader/writer should call these functions
    // if the available size is none zero but too small.
    // This will force the reader/writer cache to be flushed even though
    // there are non-zero available space.
    //
    // Calling this in other cases (on the right thread) will not affect correctness
    // but can have negative effect on performance.
    inline void sync_reader()
    {
        m_reader_cache.write = m_write_ptr.load(std::memory_order_acquire);
    }
    inline void sync_writer()
    {
        m_writer_cache.read = m_read_ptr.load(std::memory_order_acquire);
    }

private:
    inline size_t get_avail_sz(size_t start, size_t end, size_t buff_sz)
    {
        if (start < end)
            return end - start;
        if (start > end)
            return buff_sz - start;
        return 0;
    }

    struct Cache {
        T *const buff;
        const size_t buff_sz;
        const size_t max_block_sz;
        size_t write = 0;
        size_t read = 0;
        Cache(T *buff, size_t buff_sz, size_t max_block_sz)
            : buff(buff),
              buff_sz(buff_sz),
              max_block_sz(max_block_sz)
        {
        }
    };
    // Read pointer, write by reader, read by writer
    std::atomic<size_t> m_read_ptr __attribute__ ((aligned(64))) {0};
    // Cached value of the read/write pointers for the reader
    // This value should be in the cache of the reader and can be read/write with minimum
    // cache conflict.
    // The reader guarantees that `m_reader_cache.read == m_read_ptr`
    Cache m_reader_cache __attribute__ ((aligned(64)));
    // Write pointer, write by writer, read by reader
    // m_read_ptr == m_write_ptr means that the buffer is empty.
    std::atomic<size_t> m_write_ptr __attribute__ ((aligned(64))) {0};
    // Cached value of the read/write pointers for the writer
    // This value should be in the cache of the writer and can be read/write with minimum
    // cache conflict.
    // The writer guarantees that `m_writer_cache.write == m_write_ptr`
    Cache m_writer_cache __attribute__ ((aligned(64)));
} __attribute__((aligned(64)));


}

#endif
