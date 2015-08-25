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

#include "device.h"
#include "device_p.h"

#include <nacs-utils/utils.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

namespace NaCs {
namespace Kernel {

static int knacs_fd = -1;

NACS_EXPORT void
init()
{
    init("/dev/knacs");
}

NACS_EXPORT void
init(const char *name)
{
    init(open(name, O_RDWR | O_SYNC));
}

NACS_EXPORT void
init(int fd)
{
    knacs_fd = fd;
}

NACS_EXPORT bool
initialized()
{
    return knacs_fd >= 0;
}

int
getFD()
{
    static int fd = knacs_fd ? knacs_fd : ([] {
            init();
            return knacs_fd;
        })();
    return fd;
}

}
}
