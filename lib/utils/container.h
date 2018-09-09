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
#include <utility>

#include <assert.h>

namespace NaCs {

// Wrapping an arbitrary pointer/object's lifetime
// This is basically a `unique_ptr` with the type info hidden.
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
 * A mid pointer always points to somewhere in between and represents the points
 *
 * The head and the tail pointers belong to the producer thread and may not be
 * read by the filter thread.
 * The mid pointer belongs to the filter thread but can also be read by the
 * producer thread.
 * Read(write) on the mid pointer by the producer(filter)
 * should have acquire(release) ordering. Since the pointer is only writen by
 * the filter thread, reading from the filter thread does not need to be ordered.
 *
 * All three pointers can only be moved forward in the list.
 *
 * The filter thread may never change the list link.
 * It may access any elements following the mid pointer,
 * except for the item pointed to directly by the mid pointer,
 * which only allow reading of the next pointer from the filter thread.
 *
 * The producer thread may add new element to the end of the list
 * and remove element from anywhere above the mid pointer.
 * The producer thread may access any elements in the list.
 *
 * In another word, the filter thread owns the lists up to
 * but not including the mid pointer (from the tail),
 * moving the mid pointer forward hands the list element to the producer thread.
 * Access of elements are synchronized by the user.
 *
 * With the above constraints, the defined operations are,
 *
 * Producer:
 * * Push element
 * * Remove element (assumed to be before mid point)
 * * Peak at the head
 *
 * Filter:
 * * Read filter(mid) point
 * * Forward filter(mid) point
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
    inline Item *get_mid()
    {
        if (!m_mid_cache) {
            m_mid_cache = m_mid.load(std::memory_order_relaxed)
                ->next.load(std::memory_order_acquire);
        }
        return m_mid_cache;
    }
    inline Item *forward_head(Item *old_head)
    {
        m_head = old_head->next.load(std::memory_order_relaxed);
        auto new_head = m_head;
        m_alloc.free(old_head);
        return new_head;
    }
public:
    // For producer thread
    /**
     * Remove a filtered element from the head of the list.
     * If there isn't any filtered element left, return NULL.
     */
    T *pop()
    {
        // Only one NULL obj on the head is allowed.
        auto mid = m_mid.load(std::memory_order_acquire);
        auto head = m_head;
        auto res = head->obj;
        if (head == mid) {
            // Empty
            head->obj = nullptr;
            return res;
        }
        head = forward_head(head);
        if (res)
            return res;
        res = head->obj;
        assert(res);
        if (head != mid) {
            forward_head(head);
        }
        else {
            head->obj = nullptr;
        }
        return res;
    }
    /**
     * Push one element into the tail of the list to be processed.
     */
    void push(T *p)
    {
        assert(p);
        auto item = m_alloc.alloc(p);
        auto tail = m_tail;
        m_tail = item;
        tail->next.store(item, std::memory_order_release);
    }
    /**
     * Peak at the head of the list. Return the element and whether the element is already
     * processed (filtered). If the list is empty, return `{NULL, false}`.
     */
    std::pair<T*,bool> peak()
    {
        auto head = m_head;
        if (auto res = head->obj)
            return {res, true};
        auto mid = m_mid.load(std::memory_order_acquire);
        auto next = head->next.load(std::memory_order_relaxed);
        if (head == mid)
            return {next->obj, false};
        if (!next)
            return {nullptr, false};
        m_head = next;
        m_alloc.free(head);
        assert(next->obj);
        return {next->obj, true};
    }

    // For filter thread
    /**
     * Get the current item to be processed by the filter.
     */
    T *get_filter()
    {
        if (auto item = get_mid())
            return item->obj;
        return nullptr;
    }
    /**
     * Finish processing the current item.
     */
    void forward_filter()
    {
        // `get_filter` must have returned a non `NULL`
        assert(m_mid_cache);
        m_mid.store(m_mid_cache, std::memory_order_release);
        m_mid_cache = m_mid_cache->next.load(std::memory_order_acquire);
    }

    FilterQueue()
        : m_head(m_alloc.alloc(nullptr)),
          m_tail(m_head),
          m_mid(m_head)
    {
    }

private:
    // Accessed by producer only
    SmallAllocator<Item,32> m_alloc;
    Item *m_head;
    Item *m_tail;
    // Accessed by both threads
    atomic_ptr m_mid __attribute__ ((aligned(64)));
    // Accessed by filter thread only
    Item *m_mid_cache __attribute__ ((aligned(64))) = nullptr; // for filter
};

// A fixed size queue of maximum size `n`
template<typename T, uint32_t n>
class FixedQueue {
    struct alignas(T) Ceil {
        char buff[sizeof(T)];
    };
public:
    template<typename T2>
    void push(T2 &&v)
    {
        auto idx = (m_ptr + m_sz) % n;
        m_sz++;
        assert(m_sz <= n);
        auto ptr = &m_buf[idx];
        new (&ptr->buff) T(std::forward<T2>(v));
    }
    T pop()
    {
        auto ptr = &front();
        T v(std::move(*ptr));
        ptr->~T();
        m_ptr = (m_ptr + 1) % n;
        m_sz -= 1;
        return v;
    }
    const T &front() const
    {
        return *reinterpret_cast<const T*>(&m_buf[m_ptr]);
    }
    T &front()
    {
        return *reinterpret_cast<T*>(&m_buf[m_ptr]);
    }
    uint32_t size() const
    {
        return m_sz;
    }
    bool empty() const
    {
        return m_sz == 0;
    }
    bool full() const
    {
        return m_sz == n;
    }
private:
    Ceil m_buf[n];
    uint32_t m_ptr{0};
    uint32_t m_sz{0};
};

}

#endif
