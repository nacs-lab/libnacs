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

#include "../error_helper.h"

#include <../../lib/utils/mem.h>
#include <../../lib/utils/streams.h>

#include <iostream>

#include <assert.h>
#include <stdint.h>

using namespace NaCs;

struct ostream : uvector_ostream {
    using uvector_ostream::write;
    template<typename T>
    void write(T v)
    {
        static_assert(std::is_trivial_v<T>);
        write((const char*)&v, sizeof(T));
    }
    void write_string(const char *s)
    {
        write(s, strlen(s) + 1);
    }
};

int main()
{
    ostream stm;
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(32);
    stm.write<uint32_t>(98);
    stm.write<uint8_t>(2);
    stm.write<uint64_t>(((uint64_t(356 * 854) * 2034) << 32) +
                        8745273984ll); // 0x24dbcdb209424a80
    stm.write_string("hello");
    stm.write<uint16_t>(56898);
    stm.write_string("test");
    stm.write<int16_t>(-98);
    stm.write<int16_t>(-97);
    stm.write<int16_t>(-86);
    stm.write<int16_t>(56);
    stm.write<int16_t>(12);
    stm.write<int32_t>(-9);
    stm.write<int8_t>(1);
    stm.write<int8_t>(2);
    stm.write<int8_t>(3);

    auto buf = stm.get_buf();
    Mem::Reader reader(buf.data(), buf.size());

    assert(reader.read<uint32_t>() == 1);
    assert(reader.read<uint32_t>() == 32);
    assert(reader.read<uint32_t>() == 98);
    assert(reader.read<uint8_t>() == 2);
    assert(reader.read<uint64_t>() == 0x24dbcdb209424a80ull);
    auto [str1, len1] = reader.read_string();
    assert(strcmp(str1, "hello") == 0);
    assert(len1 == strlen("hello"));
    assert(reader.read<uint16_t>() == 56898);
    auto [str2, len2] = reader.read_string();
    assert(strcmp(str2, "test") == 0);
    assert(len2 == strlen("test"));
    auto p16 = reader.read_array<int16_t>(5);
    assert(p16[0] == -98);
    assert(p16[1] == -97);
    assert(p16[2] == -86);
    assert(p16[3] == 56);
    assert(p16[4] == 12);
    assert(reader.read<int32_t>() == -9);

    auto err = expect_error<std::overflow_error>([&] {
        reader.read<int32_t>();
    });
    assert(strcmp(err.what(), "Data terminates unexpectedly") == 0);
    err = expect_error<std::overflow_error>([&] {
        reader.read_string();
    });
    assert(strcmp(err.what(), "Data terminates unexpectedly") == 0);
    err = expect_error<std::overflow_error>([&] {
        reader.read_array<int8_t>(4);
    });
    assert(strcmp(err.what(), "Data terminates unexpectedly") == 0);
    auto p8 = reader.read_array<int8_t>(3);
    assert(p8[0] == 1);
    assert(p8[1] == 2);
    assert(p8[2] == 3);

    return 0;
}
