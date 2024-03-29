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

#include <knacs/knacs.h>

#include <stdlib.h>

#ifndef __NACS_KERNEL_DEVCTL_H__
#define __NACS_KERNEL_DEVCTL_H__

namespace NaCs::Kernel {

knacs_version_t getDriverVersion();
void *mapPulseCtrl();

void *allocOCMBuffer(size_t);
void freeOCMBuffer(void*, size_t);

void *allocDMABuffer(size_t);
void freeDMABuffer(void*, size_t);

}

#endif
