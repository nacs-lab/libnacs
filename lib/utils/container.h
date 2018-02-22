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

#include <memory>
#include <type_traits>

namespace NaCs {

// Wrapping an arbitrary pointer/object's lifetime
class AnyPtr : std::unique_ptr<void,void(*)(void*)> {
    // C++20
    template<typename T>
    using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
    template<typename T>
    struct Destructor {
        static void deleter(void *p)
        {
            delete (T*)p;
        }
    };
public:
    using std::unique_ptr<void,void(*)(void*)>::unique_ptr;
    template<typename T, class=std::enable_if_t<!std::is_same<T,void>::value>>
    AnyPtr(T *v)
        : unique_ptr((void*)v, Destructor<T>::deleter)
    {
    }
    template<typename T,
             class=std::enable_if_t<!std::is_lvalue_reference<T>::value &&
                                    !std::is_same<remove_cvref_t<T>,AnyPtr>::value &&
                                    !std::is_pointer<T>::value>>
    AnyPtr(T &&v)
        : AnyPtr(new T(std::move(v)))
    {
    }
};

template<typename T, typename Lock=SpinLock>
class FIFO {
    // Not supported by stdc++ yet.
    // static_assert(std::is_trivially_copyable<T>(), "");
    FIFO(const FIFO&) = delete;
    void operator=(const FIFO&) = delete;
    // m_alloc is always a power of 2
    size_t m_alloc;
    T *m_buff;
    mutable Lock m_lock;
    size_t m_read_p; // always less than m_alloc
    size_t m_write_p; // always less than m_alloc
    inline void
    doPush(const T *v)
    {
        // There has to be enough space in the buffer before calling this
        // function (i.e. `spaceLeft() > 1`).
        memcpy(m_buff + m_write_p, v, sizeof(T));
        m_write_p = (m_write_p + 1) & (m_alloc - 1);
    }
    inline void
    doPush(const T *v, size_t len)
    {
        // There has to be enough space in the buffer before calling this
        // function (i.e. `spaceLeft() > len`).
        const auto start_p = m_write_p;
        m_write_p = (m_write_p + len) & (m_alloc - 1);
        if (m_write_p > start_p || m_write_p == 0) {
            memcpy(m_buff + start_p, v, len * sizeof(T));
        } else {
            const auto start_len = len - m_write_p;
            memcpy(m_buff + start_p, v, start_len * sizeof(T));
            memcpy(m_buff, v + start_len, m_write_p * sizeof(T));
        }
    }
    inline void
    doubleSize()
    {
        m_buff = (T*)realloc(m_buff, sizeof(T) * 2 * m_alloc);
        if (m_read_p > m_write_p) {
            if (m_write_p)
                memcpy(m_buff + m_alloc, m_buff, m_write_p * sizeof(T));
            m_write_p += m_alloc;
        }
        m_alloc *= 2;
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
    inline size_t
    tryPush(const T *v, size_t len=1)
    {
        CondLock<lock, Lock> locker(m_lock);
        len = min(spaceLeft<false>() - 1, len);
        if (len > 0)
            doPush(v, len);
        return len;
    }
    template<bool lock=true>
    inline bool
    tryPush(const T &v)
    {
        CondLock<lock, Lock> locker(m_lock);
        if (spaceLeft<false>() <= 1)
            return false;
        doPush(&v);
        return true;
    }
    template<bool lock=true>
    inline void
    push(const T *v, size_t len=1)
    {
        CondLock<lock, Lock> locker(m_lock);
        while (spaceLeft<false>() <= len) {
            doubleSize();
        }
        doPush(v, len);
    }
    template<bool lock=true>
    inline void
    push(const T &v)
    {
        CondLock<lock, Lock> locker(m_lock);
        if (spaceLeft<false>() <= 1) {
            doubleSize();
        }
        doPush(&v);
    }
    template<bool lock=true>
    inline T
    pop()
    {
        // This function has no protection for overflowing
        CondLock<lock, Lock> locker(m_lock);
        T v = m_buff[m_read_p];
        m_read_p = (m_read_p + 1) & (m_alloc - 1);
        return v;
    }
    template<bool lock=true>
    inline size_t
    pop(T *v, size_t len=1)
    {
        CondLock<lock, Lock> locker(m_lock);
        len = min(size<false>(), len);
        if (!len)
            return 0;
        const auto start_p = m_read_p;
        m_read_p = (m_read_p + len) & (m_alloc - 1);

        if (m_read_p > start_p || m_read_p == 0) {
            memcpy(v, m_buff + start_p, len * sizeof(T));
        } else {
            const auto start_len = len - m_read_p;
            memcpy(v, m_buff + start_p, start_len * sizeof(T));
            memcpy(v + start_len, m_buff, m_read_p * sizeof(T));
        }
        return len;
    }
    inline size_t
    tryPop(T *v, size_t len=1)
    {
        std::unique_lock<Lock> locker(m_lock, std::defer_lock);
        if (!locker.try_lock()) {
            return 0;
        }
        return pop<false>(v, len);
    }
};

}

#endif
