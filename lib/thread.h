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

#ifndef _NACS_UTILS_THREAD_H_
#define _NACS_UTILS_THREAD_H_

#include <atomic>
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
        return m_spin.exchange(true);
    }
    inline void
    lock()
    {
        while (m_spin.exchange(true)) {
            // It is not a issue if we go to sleep and yielding execution
            // here seems to make the timing in the present of lock contention
            // shorter and more consistent.
            std::this_thread::yield();
        }
    }
    inline void
    unlock()
    {
        m_spin = false;
    }
};

}

#endif
