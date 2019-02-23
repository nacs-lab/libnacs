/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/utils/number.h"
#include "../../lib/utils/timer.h"
#include <assert.h>
#include <math.h>

using namespace NaCs;

bool approx(double a, double b)
{
    double diff = abs(a - b);
    double avg = abs(a + b) / 2;
    return diff < 2e-8 || diff / avg < 2e-8;
}

int main()
{
    const double points[] = {0, 0.1, 0.2, 0.6};
    assert(approx(linearInterpolate(0, 2, 3, 4, points), 0));
    assert(approx(linearInterpolate(2, 2, 3, 4, points), 0));
    assert(approx(linearInterpolate(2.3, 2, 3, 4, points), 0.03));
    assert(approx(linearInterpolate(3.4, 2, 3, 4, points), 0.14));
    assert(approx(linearInterpolate(4.5, 2, 3, 4, points), 0.4));
    assert(approx(linearInterpolate(5, 2, 3, 4, points), 0.6));
    assert(approx(linearInterpolate(7, 2, 3, 4, points), 0.6));
    return 0;
}
