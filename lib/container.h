/*************************************************************************
 *   Copyright (c) 2013 - 2015 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_UTILS_CONTAINER_H__
#define __NACS_UTILS_CONTAINER_H__

#include "thread.h"
#include "number.h"

namespace NaCs {

template<typename T, typename Lock=SpinLock>
class FIFO {
    FIFO(const FIFO&) = delete;
    void operator=(const FIFO&) = delete;
    // m_alloc is always a power of 2
    size_t m_alloc;
    T *m_buff;
    mutable Lock m_lock;
    size_t m_read_p;
    size_t m_write_p;
    inline void
    doPush(T &v)
    {
        m_buff[m_write_p & (m_alloc - 1)] = v;
        ++m_write_p;
    }
public:
    inline
    FIFO(size_t size=256)
        : m_alloc(1 << max(getBits(size) - 1, 2)),
          m_buff((T*)malloc(sizeof(T) * m_alloc)),
          m_lock(),
          m_read_p(0),
          m_write_p(0)
    {
    }
    inline size_t
    capacity() const
    {
        return m_alloc;
    }
    template<bool lock=true>
    inline size_t
    size() const
    {
        CondLock<lock, Lock> locker(m_lock);
        return (m_write_p - m_read_p) & (m_alloc - 1);
    }
    template<bool lock=true>
    inline size_t
    spaceLeft() const
    {
        CondLock<lock, Lock> locker(m_lock);
        return capacity() - size<false>();
    }
    template<bool lock=true>
    inline bool
    tryPush(T &v)
    {
        CondLock<lock, Lock> locker(m_lock);
        if (spaceLeft<false>() <= 1)
            return false;
        doPush(v);
        return true;
    }
    template<bool lock=true>
    inline void
    push(T &v)
    {
        CondLock<lock, Lock> locker(m_lock);
        if (spaceLeft<false>() <= 1) {
            auto mask = m_alloc - 1;
            m_read_p &= mask;
            m_write_p &= mask;
            m_buff = (T*)realloc(m_buff, sizeof(T) * 2 * m_alloc);
            if (m_read_p > m_write_p) {
                if (m_write_p)
                    memcpy(m_buff + m_alloc, m_buff, m_write_p * sizeof(T));
                m_write_p += m_alloc;
            }
            m_alloc *= 2;
        }
        doPush(v);
    }
    template<bool lock=true>
    inline T
    pop()
    {
        CondLock<lock, Lock> locker(m_lock);
        T v = m_buff[m_read_p & (m_alloc - 1)];
        ++m_read_p;
        return v;
    }
};

}

#endif
