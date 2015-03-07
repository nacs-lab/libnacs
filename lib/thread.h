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

namespace NaCs {

/**
 * A straightforward implementation of a spin lock
 *
 * According to the benchmark (see test/test_thread.cpp) this lock is
 * faster than the stdlib lock (pthread_mutex_t or std::mutex) both with
 * and without lock contention and on both ARM and x86_64 (Haswell).
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
