/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include "timing.h"

#include <nacs-utils/utils.h>

namespace NaCs {
namespace Seq {
namespace Time {

NACS_EXPORT(seq) Keeper::Keeper(const Constraints &cons)
    : m_cons(cons),
      m_plens(cons.avg_window - 1)
{
    reset();
}

NACS_EXPORT(seq) void Keeper::addPulse(uint64_t dt)
{
    size_t pidx = m_npulses % (m_cons.avg_window - 1);
    m_total_len = m_total_len + dt - m_plens[pidx];
    m_plens[pidx] = dt;
}

NACS_EXPORT(seq) uint64_t Keeper::minDt(uint64_t min_dt) const
{
    auto window_sz = m_cons.avg_dt * m_cons.avg_window;
    min_dt = std::max<uint64_t>(min_dt, 1);
    if (m_npulses < m_cons.avg_window - 1 || m_total_len >= window_sz)
        return min_dt;
    return std::max(window_sz - m_total_len, min_dt);
}

NACS_EXPORT(seq) void Keeper::reset(void)
{
    m_npulses = 0;
    m_total_len = 0;
    memset(m_plens.data(), 0, m_plens.size() * sizeof(uint64_t));
}

}
}
}
