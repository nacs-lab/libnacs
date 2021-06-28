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

#ifndef __NACS_KERNEL_DMA_BUFFER_H__
#define __NACS_KERNEL_DMA_BUFFER_H__

#include <stddef.h>

#include <type_traits>

namespace NaCs::Kernel {

class DMABufferBase {
    void resize_slow(size_t);
protected:
    void *m_buff;
    size_t m_size;
public:
    DMABufferBase(size_t=4096);
    void resize(size_t);
    void release();
    void send(size_t);
    ~DMABufferBase();
};

inline void
DMABufferBase::resize(size_t new_size)
{
    if (new_size <= m_size)
        return;
    resize_slow(new_size);
}

template<typename T>
class DMABuffer : private DMABufferBase {
    size_t m_len;
public:
    DMABuffer(size_t len=0, size_t reserve_size=4096 / sizeof(T))
        : DMABufferBase(reserve_size * sizeof(T)),
          m_len(0)
    {
        static_assert(std::is_pod<T>()(), "");
        resize(len);
    }
    void
    size_hint(size_t len)
    {
        DMABufferBase::resize(len * sizeof(T));
    }
    void
    resize(size_t len)
    {
        m_len = len;
        size_hint(len);
    }
    void
    release()
    {
        m_len = 0;
        DMABufferBase::release();
    }
    void
    send()
    {
        DMABufferBase::send(m_len * sizeof(T));
        m_len = 0;
    }
    size_t
    length() const
    {
        return m_len;
    }
    const T&
    operator[](size_t i) const
    {
        return ((T*)m_buff)[i];
    }
    T&
    operator[](size_t i)
    {
        return ((T*)m_buff)[i];
    }
    void
    push_back(const T &v)
    {
        size_t idx = m_len;
        resize(m_len + 1);
        (*this)[idx] = v;
    }
};

}

#endif
