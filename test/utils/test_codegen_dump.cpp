/*************************************************************************
 *   Copyright (c) 2016 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "codegen_helper.h"

#include "../../lib/utils/number.h"
#include "../../lib/utils/streams.h"
#include "../../lib/utils/timer.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <math.h>

int main(int, char **argv)
{
    TestCtx ctx;
    IR::Builder builder(IR::Type::Int32, {});
    builder.createRet(builder.getConstInt(0));
    ctx.get_llvm_test(builder.get()).get_ptr(argv[1]);
    return 0;
}
