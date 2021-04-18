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

#include "winpath_helper.h"

#include "../../lib/utils/dlload.h"
#include "../../lib/utils/number.h"
#include <assert.h>
#include <math.h>

using namespace NaCs;

NACS_EXPORT() int global_symbol = 1;

int main()
{
    auto hdl = DL::open("libopenlibm", DL::GLOBAL | DL::NOW);
    assert(hdl && "Unable to load openlibm");
    auto psin = (double(*)(double))DL::sym(hdl, "sin");
    assert(psin && "Cannot find sin in openlibm");
    assert(approx(sin(1.5), psin(1.5)));
    DL::close(hdl);

    auto *pgv = DL::sym(nullptr, "global_symbol");
    assert(pgv == (void*)&global_symbol);

    return 0;
}
