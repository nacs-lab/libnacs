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

#include "dma_buffer.h"
#include "devctl_p.h"

#include <nacs-utils/utils.h>
#include <nacs-utils/number.h>

namespace NaCs {

using namespace NaCs;

namespace Kernel {

static const auto page_size = sysconf(_SC_PAGESIZE);

NACS_EXPORT(kernel)
DMABufferBase::DMABufferBase(size_t size)
    : m_buff(nullptr),
      m_size(0)
{
    resize(size);
}

NACS_EXPORT(kernel) void
DMABufferBase::resize_slow(size_t size)
{
    size = alignTo((off_t)size, page_size);
    if (!m_buff) {
        m_buff = allocDmaBuffer(size);
    } else {
        m_buff = reallocDmaBuffer(m_buff, m_size, size);
    }
    m_size = size;
}

NACS_EXPORT(kernel)
DMABufferBase::~DMABufferBase()
{
    release();
}

NACS_EXPORT(kernel) void
DMABufferBase::release()
{
    if (m_buff) {
        freeDmaBuffer(m_buff, m_size);
        m_buff = nullptr;
        m_size = 0;
    }
}

NACS_EXPORT(kernel) void
DMABufferBase::send(size_t len)
{
    sendDmaBuffer(m_buff, len);
    m_buff = nullptr;
    m_size = 0;
}

}
}
