/*************************************************************************
 *   Copyright (c) 2015 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "devctl.h"

#include "device_p.h"

#include "../nacs-utils/utils.h"
#include "../nacs-utils/mem.h"
#include "../nacs-utils/fd_utils.h"

#include <sys/ioctl.h>
#include <sys/mman.h>

namespace NaCs::Kernel {

NACS_EXPORT() knacs_version_t
getDriverVersion()
{
    knacs_version_t ver;
    checkErrno(ioctl(getFD(), KNACS_GET_VERSION, &ver));
    return ver;
}

NACS_EXPORT() void*
mapPulseCtrl()
{
    return mapFile(getFD(), 0, 32 * 4);
}

NACS_EXPORT() void *allocOCMBuffer(size_t len)
{
    return mapFile(getFD(), page_size, len);
}

NACS_EXPORT() void freeOCMBuffer(void *buff, size_t old_size)
{
    munmap(buff, old_size);
}

NACS_EXPORT() void *allocDMABuffer(size_t len)
{
    return mapFile(getFD(), page_size * 2, len);
}

NACS_EXPORT() void freeDMABuffer(void *buff, size_t old_size)
{
    munmap(buff, old_size);
}

}
