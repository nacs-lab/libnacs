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

#include "devctl_p.h"
#include "device_p.h"

#include <nacs-utils/utils.h>
#include <nacs-utils/fd_utils.h>

#include <stropts.h>
#include <sys/mman.h>

namespace NaCs {
namespace Kernel {

static const auto page_size = sysconf(_SC_PAGESIZE);

NACS_EXPORT knacs_version_t
getDriverVersion()
{
    knacs_version_t ver;
    checkErrno(ioctl(getFD(), KNACS_GET_VERSION, &ver));
    return ver;
}

NACS_EXPORT void*
mapPulseCtrl()
{
    return mapFile(getFD(), 0, 32 * 4);
}

void*
allocDmaBuffer(size_t len)
{
    return mapFile(getFD(), page_size, len);
}

void*
reallocDmaBuffer(void *buff, size_t old_size, size_t new_size)
{
    return mremap(buff, old_size, new_size, MREMAP_MAYMOVE);
}

void
freeDmaBuffer(void *buff, size_t old_size)
{
    munmap(buff, old_size);
}

void
sendDmaBuffer(void *buff, size_t len)
{
    knacs_dma_buffer_t kernel_buff = {
        (unsigned long)len,
        buff
    };
    checkErrno(ioctl(getFD(), KNACS_SEND_DMA_BUFFER, &kernel_buff));
}

}
}
