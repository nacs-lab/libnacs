/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "ir.h"

#ifndef __NACS_UTILS_IR_P_H__
#define __NACS_UTILS_IR_P_H__

namespace NaCs {
namespace IR {

typedef double (*fptr_t)(double);

enum class BuiltinType : uint8_t {
    Invalid,
    F64_F64,
    F64_F64F64,
    F64_F64F64F64,
    F64_F64I32,
    F64_I32F64,
};

const char *getBuiltinSymbol(Builtins id);
BuiltinType getBuiltinType(Builtins id);

#ifndef NACS_HAS_EXP10
NACS_EXPORT() double nacs_exp10(double x);
#endif

}
}

#endif
