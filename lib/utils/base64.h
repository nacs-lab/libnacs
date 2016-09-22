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

#include "utils.h"

#include <vector>

#ifndef __NACS_UTILS_BASE64_H__
#define __NACS_UTILS_BASE64_H__

namespace NaCs {
namespace Base64 {

size_t encode_len(const uint8_t *data, size_t len);
void encode(uint8_t *dest, const uint8_t *data, size_t len);

size_t decode_len(const uint8_t *data, size_t len);
void decode(uint8_t *dest, const uint8_t *data, size_t len);
bool validate(const uint8_t *data, size_t len);

std::vector<uint8_t> encode(const std::vector<uint8_t> &data);
std::vector<uint8_t> decode(const std::vector<uint8_t> &data);
bool validate(const std::vector<uint8_t> &data);

}
}

#endif
