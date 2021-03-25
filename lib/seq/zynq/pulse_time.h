/*************************************************************************
 *   Copyright (c) 2016 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include <nacs-utils/utils.h>

#ifndef __NACS_SEQ_ZYNQ_PULSE_TIME_H__
#define __NACS_SEQ_ZYNQ_PULSE_TIME_H__

namespace NaCs::Seq::Zynq {

namespace PulseTime {
static constexpr uint8_t Min = 3;
static constexpr uint8_t _DDS = 50;
static constexpr uint8_t DDSFreq = _DDS;
static constexpr uint8_t DDSAmp = _DDS;
static constexpr uint8_t DDSPhase = _DDS;
static constexpr uint8_t DDSReset = _DDS;
static constexpr uint8_t Clear = 5;
static constexpr uint8_t LoopBack = 5;
static constexpr uint8_t Clock = 5;
static constexpr uint8_t DAC = 45;
};

}

#endif
