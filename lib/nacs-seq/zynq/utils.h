/*************************************************************************
 *   Copyright (c) 2026 - 2026 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_SEQ_ZYNQ_UTILS_H__
#define __NACS_SEQ_ZYNQ_UTILS_H__

#include "../../nacs-utils/number.h"

namespace NaCs::Seq::Zynq {

enum class ChnType : uint8_t {
    TTL = 1,
    Freq = 2,
    Amp = 3,
    DAC = 4,
    Phase = 5,
    Clock = 6,
};

static inline uint32_t dds_freq_to_mu(double freq)
{
    constexpr double factor = 1.0 * (1 << 16) * (1 << 16) / 3.5e9;
    if (freq <= 0)
        return 0;
    // The instruction only have 31 bits so it actually won't encode exactly 1.75 GHz.
    if (freq > 1.7499999987776392e9) // last float that would round to 2^31 - 2
        return 0x7fffffff;
    return round<int32_t>(freq * factor);
}

static constexpr inline double dds_freq_from_mu(uint32_t freq)
{
    constexpr double factor = 3.5e9 / (1 << 16) / (1 << 16);
    return double(freq) * factor;
}

static inline uint16_t dds_amp_to_mu(double amp)
{
    constexpr double factor = 4095;
    if (amp <= 0)
        return 0;
    if (amp >= 1)
        return 4095;
    return (uint16_t)round<int32_t>(amp * factor);
}

static constexpr inline double dds_amp_from_mu(uint32_t amp)
{
    return double(amp) / 4095;
}

static inline uint32_t dds_phase_to_mu(double phase)
{
    // phase in [0, 1]
    return round<int32_t>(phase * (1 << 16)) & ((1 << 16) - 1);
}

static constexpr inline double dds_phase_from_mu(uint32_t phase)
{
    return double(phase) / double(1 << 16);
}

static inline uint16_t dac_to_mu(double V)
{
    constexpr double factor = 65535 / 20.0;
    constexpr double offset = 10.0;
    if (V <= -10)
        return 0;
    if (V >= 10)
        return 0xffff;
    V += offset;
    return (uint16_t)round<int32_t>(V * factor);
}

static constexpr inline double dac_from_mu(uint32_t V)
{
    return V * (20.0 / 65535) - 10.0;
}

}

#endif
