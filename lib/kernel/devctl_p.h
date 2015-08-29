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

#include "devctl.h"
#include <stddef.h>

#ifndef __NACS_KERNEL_DEVCTL_P_H__
#define __NACS_KERNEL_DEVCTL_P_H__

namespace NaCs {
namespace Kernel {

void *allocDmaBuffer(size_t);
void *reallocDmaBuffer(void*, size_t old_size, size_t new_size);
void freeDmaBuffer(void *, size_t old_size);

}
}

#endif
