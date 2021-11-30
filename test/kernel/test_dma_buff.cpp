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

#include "../../lib/nacs-kernel/devctl.h"

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

using namespace NaCs;

int
main()
{
    auto ptr = Kernel::allocOCMBuffer(4096);
    auto ptr2 = Kernel::allocDMABuffer(4096);

    printf("ocm: %p\n", ptr);
    printf("dma: %p\n", ptr2);

    Kernel::freeOCMBuffer(ptr, 4096);
    Kernel::freeDMABuffer(ptr2, 4096);

    return 0;
}
