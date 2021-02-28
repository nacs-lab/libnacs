/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "error.h"

#include "../nacs-utils/utils.h"

namespace NaCs::Seq {

Error::Error(Type type, uint16_t code, Type type1, uint64_t id1, const char *what)
    : std::runtime_error(what),
      type(type),
      type1(type1),
      code(code),
      id1(id1)
{
}

Error::Error(Type type, uint16_t code, Type type1, uint64_t id1, const std::string &what)
    : std::runtime_error(what),
      type(type),
      type1(type1),
      code(code),
      id1(id1)
{
}

Error::Error(Type type, uint16_t code, Type type1, uint64_t id1,
             Type type2, uint64_t id2, const char *what)
    : std::runtime_error(what),
      type(type),
      type1(type1),
      type2(type2),
      code(code),
      id1(id1),
      id2(id2)
{
}

Error::Error(Type type, uint16_t code, Type type1, uint64_t id1,
             Type type2, uint64_t id2, const std::string &what)
    : std::runtime_error(what),
      type(type),
      type1(type1),
      type2(type2),
      code(code),
      id1(id1),
      id2(id2)
{
}

Error::Error(const Error &other)
    : std::runtime_error(other),
      type(other.type),
      type1(other.type1),
      type2(other.type2),
      code(other.code),
      id1(other.id1),
      id2(other.id2)
{
}

Error::~Error()
{
}

}
