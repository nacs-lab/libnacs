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

#include "../../lib/kernel/devctl.h"
#include "../../lib/kernel/device.h"

#include <stdio.h>

using namespace NaCs;

int
main()
{
    // Initialization should happen automatically, explicitly check for it
    // in this test.
    Kernel::init();
    if (!Kernel::initialized()) {
        fprintf(stderr, "Open KNaCs device failed\n");
        return -1;
    }
    auto ver = Kernel::getDriverVersion();
    printf("KNaCs Version: %d.%d\n", ver.major, ver.minor);
    return 0;
}
