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

#include "../../lib/nacs-kernel/devctl.h"

#include <stdio.h>
#include <stdint.h>

using namespace NaCs;

int
main()
{
    volatile uint32_t *regs = (uint32_t*)Kernel::mapPulseCtrl();
    printf("Pulse controller mapped at %p\n", regs);
    for (int i = 0;i < 32;i++) {
        printf("Reg[%d] = %d\n", i, int(regs[i]));
        regs[i] = i * 2 - 1;
    }
    for (int i = 0;i < 32;i++) {
        printf("Reg[%d] = %d\n", i, int(regs[i]));
    }
    return 0;
}
