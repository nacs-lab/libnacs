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

#include "base64.h"
#include <assert.h>

namespace NaCs {
namespace Base64 {

NACS_EXPORT(utils) size_t encode_len(const uint8_t*, size_t len)
{
    size_t nbits2 = len * 4;
    size_t nbytes = nbits2 / 3;
    switch (nbits2 % 3) {
    case 1:
        return nbytes + 3;
    case 2:
        return nbytes + 2;
    default:
        return nbytes;
    }
}

static const uint8_t encode_table[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

NACS_EXPORT(utils) void encode(uint8_t *dest, const uint8_t *data, size_t len)
{
    size_t i = 0;
    for (;i + 3 <= len;i += 3) {
        uint8_t d0 = data[i];
        uint8_t d1 = data[i + 1];
        uint8_t d2 = data[i + 2];
        dest[0] = encode_table[d0 >> 2];
        dest[1] = encode_table[((d0 & 3) << 4) | (d1 >> 4)];
        dest[2] = encode_table[((d1 & 0xf) << 2) | (d2 >> 6)];
        dest[3] = encode_table[d2 & 0x3f];
        dest += 4;
    }
    switch (len - i) {
    case 1: {
        uint8_t d0 = data[i];
        dest[0] = encode_table[d0 >> 2];
        dest[1] = encode_table[(d0 & 3) << 4];
        dest[2] = '=';
        dest[3] = '=';
        break;
    }
    case 2: {
        uint8_t d0 = data[i];
        uint8_t d1 = data[i + 1];
        dest[0] = encode_table[d0 >> 2];
        dest[1] = encode_table[((d0 & 3) << 4) | (d1 >> 4)];
        dest[2] = encode_table[(d1 & 0xf) << 2];
        dest[3] = '=';
        break;
    }
    default:
        break;
    }
}

static constexpr uint8_t decode_table[256] = {
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 62, 64, 64, 64, 63,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 64, 64, 64, 64, 64, 64,
    64,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 64, 64, 64, 64, 64,
    64, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
};

NACS_EXPORT(utils) size_t decode_len(const uint8_t *data, size_t len)
{
    assert(len % 4 == 0);
    if (len < 4)
        return 0;
    const uint8_t *end = data + len;
    size_t nmax = len / 4 * 3;
    if (end[-2] == '=') {
        return nmax - 2;
    } else if (end[-1] == '=') {
        return nmax - 1;
    } else {
        return nmax;
    }
}

NACS_EXPORT(utils) void decode(uint8_t *dest, const uint8_t *data, size_t len)
{
    if (len < 4)
        return;
    assert(len % 4 == 0);
    size_t i = 0;
    for (;i < len - 4;i += 4) {
        uint8_t c0 = data[i];
        uint8_t c1 = data[i + 1];
        uint8_t c2 = data[i + 2];
        uint8_t c3 = data[i + 3];
        uint8_t d0 = decode_table[c0];
        uint8_t d1 = decode_table[c1];
        uint8_t d2 = decode_table[c2];
        uint8_t d3 = decode_table[c3];
        dest[0] = uint8_t(d0 << 2) | uint8_t(d1 >> 4);
        dest[1] = uint8_t((d1 & 0xf) << 4) | uint8_t(d2 >> 2);
        dest[2] = uint8_t((d2 & 3) << 6) | d3;
        dest += 3;
    }
    uint8_t c0 = data[i];
    uint8_t c1 = data[i + 1];
    uint8_t c2 = data[i + 2];
    uint8_t c3 = data[i + 3];
    uint8_t d0 = decode_table[c0];
    uint8_t d1 = decode_table[c1];
    uint8_t d2 = decode_table[c2];
    uint8_t d3 = decode_table[c3];
    dest[0] = uint8_t(d0 << 2) | uint8_t(d1 >> 4);
    if (c2 != '=') {
        dest[1] = uint8_t((d1 & 0xf) << 4) | uint8_t(d2 >> 2);
        if (c3 != '=') {
            dest[2] = uint8_t((d2 & 3) << 6) | d3;
        }
    }
}

NACS_EXPORT(utils) std::vector<uint8_t> encode(const std::vector<uint8_t> &data)
{
    std::vector<uint8_t> res(encode_len(data.data(), data.size()));
    encode(res.data(), data.data(), data.size());
    return res;
}

NACS_EXPORT(utils) std::vector<uint8_t> decode(const std::vector<uint8_t> &data)
{
    std::vector<uint8_t> res(decode_len(data.data(), data.size()));
    decode(res.data(), data.data(), data.size());
    return res;
}

NACS_EXPORT(utils) bool validate(const uint8_t *data, size_t len)
{
    if (len % 4 != 0)
        return false;
    if (len == 0)
        return true;
    for (size_t i = 0;i + 2 < len;i++) {
        uint8_t d = data[i];
        if (decode_table[d] == 64) {
            return false;
        }
    }
    // This doesn't catch if the final bits are not 0
    const uint8_t *end = data + len;
    if (end[-2] == '=') {
        return end[-1] == '=';
    } else if (decode_table[end[-2]] == 64) {
        return false;
    } else if (end[-1] == '=') {
        return true;
    } else {
        return decode_table[end[-1]] != 64;
    }
}

NACS_EXPORT(utils) bool validate(const std::vector<uint8_t> &data)
{
    return validate(data.data(), data.size());
}

}
}
