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

#include <nacs-utils/base64.h>
#include <assert.h>

#include <iostream>

using namespace NaCs;

static void show_vector_u8(const std::vector<uint8_t> &vect)
{
    for (char c: vect)
        std::cout << c;
    std::cout << std::endl;
}

int main()
{
    std::vector<uint8_t> data1 = {'a'};
    std::vector<uint8_t> data2 = {'a', 'b'};
    std::vector<uint8_t> data3 = {'a', 'b', 'c'};
    std::vector<uint8_t> data4 = {'a', 'b', 'c', 'd'};
    auto encode1 = Base64::encode(data1);
    auto encode2 = Base64::encode(data2);
    auto encode3 = Base64::encode(data3);
    auto encode4 = Base64::encode(data4);
    show_vector_u8(encode1);
    show_vector_u8(encode2);
    show_vector_u8(encode3);
    show_vector_u8(encode4);
    Base64::validate(encode1);
    Base64::validate(encode2);
    Base64::validate(encode3);
    Base64::validate(encode4);
    auto decode1 = Base64::decode(encode1);
    auto decode2 = Base64::decode(encode2);
    auto decode3 = Base64::decode(encode3);
    auto decode4 = Base64::decode(encode4);
    show_vector_u8(decode1);
    show_vector_u8(decode2);
    show_vector_u8(decode3);
    show_vector_u8(decode4);
    return 0;
}
