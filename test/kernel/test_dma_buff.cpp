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

#include <nacs-kernel/devctl.h>

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

using namespace NaCs;

int
main()
{
    volatile uint32_t *buff = (uint32_t*)Kernel::allocDmaBuffer(4096);
    printf("Buffer allocated\n");
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);
    buff[0] = 1;
    buff[1023] = 2;
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);

    buff = (uint32_t*)Kernel::rellocDmaBuffer((void*)buff, 4096, 4096 * 2);
    printf("Buffer re-allocated\n");
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);
    printf("Buff[1024] = %d\n", buff[1024]);
    printf("Buff[2047] = %d\n", buff[2047]);
    buff[1024] = 3;
    buff[2047] = 4;
    printf("Buff[1024] = %d\n", buff[1024]);
    printf("Buff[2047] = %d\n", buff[2047]);

    volatile uint32_t *buff2 = (uint32_t*)Kernel::allocDmaBuffer(4096);
    buff2[1] = 2;

    Kernel::freeDmaBuffer((void*)buff, 4096 * 2);
    Kernel::freeDmaBuffer((void*)buff2, 4096);

    return 0;
}
