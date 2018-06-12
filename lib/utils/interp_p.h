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
#include "number.h"

#ifndef __NACS_UTILS_INTERP_P_H__
#define __NACS_UTILS_INTERP_P_H__

namespace NaCs {
namespace IR {

// Functions that are used both by the interpreter as well as for constant propagation.
static NACS_UNUSED TagVal evalAdd(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() + val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() + val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static NACS_UNUSED TagVal evalSub(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() - val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() - val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static NACS_UNUSED TagVal evalMul(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() * val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() * val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static NACS_UNUSED TagVal evalFDiv(TagVal val1, TagVal val2)
{
    return val1.get<double>() / val2.get<double>();
}

static NACS_UNUSED TagVal evalCmp(CmpType cmptyp, TagVal val1, TagVal val2)
{
    double v1 = val1.get<double>();
    double v2 = val2.get<double>();
    switch (cmptyp) {
    case CmpType::eq:
        return v1 == v2;
    case CmpType::gt:
        return v1 > v2;
    case CmpType::ge:
        return v1 >= v2;
    case CmpType::lt:
        return v1 < v2;
    case CmpType::le:
        return v1 <= v2;
    case CmpType::ne:
        return v1 != v2;
    default:
        return false;
    }
}

ExeContext *get_interp_context();

}
}

#endif
