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

#include "../../lib/kernel/dma_buffer.h"

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

using namespace NaCs;

int
main()
{
    Kernel::DMABuffer<uint32_t> buff(1024);
    printf("Buffer created %p\n", &buff[0]);
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);
    buff[0] = 1;
    buff[1023] = 2;
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);

    buff.resize(2048);

    printf("Buffer resized, %p\n", &buff[0]);
    printf("Buff[0] = %d\n", buff[0]);
    printf("Buff[1023] = %d\n", buff[1023]);
    printf("Buff[1024] = %d\n", buff[1024]);
    printf("Buff[2047] = %d\n", buff[2047]);
    buff[1024] = 3;
    buff[2047] = 4;
    printf("Buff[1024] = %d\n", buff[1024]);
    printf("Buff[2047] = %d\n", buff[2047]);

    Kernel::DMABuffer<uint32_t> buff2;
    printf("Buffer 2 created %p\n", &buff2[0]);
    buff2.push_back(2);

    buff.send();
    printf("Buffer 1 sent\n");
    buff2.send();
    printf("Buffer 2 sent\n");

    buff.release();
    printf("Buffer 1 release\n");
    buff2.release();
    printf("Buffer 2 release\n");

    return 0;
}
