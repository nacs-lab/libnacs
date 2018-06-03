/*************************************************************************
 *   Copyright (c) 2013 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "number.h"

namespace NaCs {

NACS_EXPORT() double linearInterpolate(double x, double x0, double dx,
                                       uint32_t npoints, const double *points)
{
    // Offset
    x -= x0;
    // Scale
    if (unlikely(x <= 0)) {
        return points[0];
    }
    else if (unlikely(x >= dx)) {
        return points[npoints - 1];
    }
    x = x * (npoints - 1) / dx;
    double lof = 0.0;
    x = modf(x, &lof);
    uint32_t lo = (uint32_t)lof;
    double vlo = points[lo];
    if (x == 0)
        return vlo;
    double vhi = points[lo + 1];
    return x * vhi + (1 - x) * vlo;
}

}
