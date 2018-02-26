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
#include "mem.h"

#include <memory>
#include <type_traits>

#include <assert.h>

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

/**
 * This is a queue from one thread (producer) to another thread (filter) and back.
 * Based on Herb Sutter's version in
 * http://www.drdobbs.com/parallel/writing-lock-free-code-a-corrected-queue/210604448
 * Wrapped in an API that's suitable for passing the data to another thread an back.
 * Assuming atomics on pointers are lockless, this should be wait-free,
 * and have all allocations done only on the producer thread.
 *
 * The items are chained in a singly linked list. The list is never empty.
 * The head pointer points to the head of the list where items are removed
 * from the list by the producer and the tail pointer points to the tail
 * of the list where new items are added by the producer.
 * A peak pointer always points to somewhere in between and represents the points
 *
 * The head and the tail pointers belong to the producer thread and may not be
 * read by the filter thread.
 * The peak pointer belongs to the filter thread but can also be read by the
 * producer thread.
 * Read(write) on the peak pointer by the producer(filter)
 * should have acquire(release) ordering. Since the pointer is only writen by
 * the filter thread, reading from the filter thread does not need to be ordered.
 *
 * All three pointers can only be moved forward in the list.
 *
 * The filter thread may never change the list link.
 * It may access any elements following the peak pointer,
 * except for the item pointed to directly by the peak pointer,
 * which only allow reading of the next pointer from the filter thread.
 *
 * The producer thread may add new element to the end of the list
 * and remove element from anywhere above the peak pointer.
 * The producer thread may access any elements in the list.
 *
 * In another word, the filter thread owns the lists up to
 * but not including the peak pointer (from the tail),
 * moving the peak pointer forward hands the list element to the producer thread.
 * Access of elements are synchronized by the user.
 *
 * With the above constraints, the defined operations are,
 *
 * Producer:
 * * Push element
 * * Remove element (assumed to be before peak point)
 *
 * Filter:
 * * Read peak point
 * * Forward peak point
 */
template<typename T>
class FilterQueue {
    struct Item;
    using atomic_ptr = std::atomic<Item*>;
    struct Item {
        T *obj;
        atomic_ptr next;
        Item(T *obj)
            : obj(obj),
              next(nullptr)
        {
        }
    };
public:
    // For producer thread
    T *pop()
    {
        auto peak = m_peak.load(std::memory_order_acquire);
        auto head = m_head;
        while (true) {
            auto res = head->obj;
            if (head == peak) {
                head->obj = nullptr;
                return res;
            }
            m_head = head->next.load(std::memory_order_relaxed);
            auto new_head = m_head;
            m_alloc.free(head);
            if (res)
                return res;
            head = new_head;
        }
    }
    void push(T *p)
    {
        assert(p);
        auto item = m_alloc.alloc(p);
        auto tail = m_tail;
        m_tail = item;
        tail->next.store(item, std::memory_order_release);
    }

    // For filter thread
    T *get_filter()
    {
        if (auto item = get_peak())
            return item->obj;
        return nullptr;
    }
    void forward_filter()
    {
        // `get_filter` must have returned a non `NULL`
        assert(m_peak_cache);
        m_peak.store(m_peak_cache, std::memory_order_release);
        m_peak_cache = m_peak_cache->next.load(std::memory_order_acquire);
    }

    FilterQueue()
        : m_head(m_alloc.alloc(nullptr)),
          m_tail(m_head),
          m_peak(m_head)
    {
    }

private:
    inline Item *get_peak()
    {
        if (!m_peak_cache) {
            m_peak_cache = m_peak.load(std::memory_order_relaxed)
                ->next.load(std::memory_order_acquire);
        }
        return m_peak_cache;
    }
    // Accessed by producer only
    SmallAllocator<Item,32> m_alloc;
    Item *m_head;
    Item *m_tail;
    // Accessed by both threads
    atomic_ptr m_peak __attribute__ ((aligned(64)));
    // Accessed by filter thread only
    Item *m_peak_cache __attribute__ ((aligned(64))) = nullptr; // for filter
};

}

#endif
