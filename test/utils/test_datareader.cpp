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

#define CATCH_CONFIG_MAIN

#include "../error_helper.h"

#include <../../lib/nacs-utils/mem.h>
#include <../../lib/nacs-utils/streams.h>

#include <iostream>

#include <stdint.h>

using namespace NaCs;
using namespace std::literals::string_literals;

struct ostream : uvector_ostream {
    using uvector_ostream::write;
    template<typename T>
    void write(T v)
    {
        STATIC_REQUIRE(std::is_trivial_v<T>);
        write((const char*)&v, sizeof(T));
    }
    void write_string(const char *s)
    {
        write(s, strlen(s) + 1);
    }
};

TEST_CASE("DataReader") {
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

    REQUIRE(reader.read<uint32_t>() == 1);
    REQUIRE(reader.read<uint32_t>() == 32);
    REQUIRE(reader.read<uint32_t>() == 98);
    REQUIRE(reader.read<uint8_t>() == 2);
    REQUIRE(reader.read<uint64_t>() == 0x24dbcdb209424a80ull);
    auto [str1, len1] = reader.read_string();
    REQUIRE(str1 == "hello"s);
    REQUIRE(len1 == strlen("hello"));
    REQUIRE(reader.read<uint16_t>() == 56898);
    auto [str2, len2] = reader.read_string();
    REQUIRE(str2 == "test"s);
    REQUIRE(len2 == strlen("test"));
    auto p16 = reader.read_array<int16_t>(5);
    REQUIRE(p16[0] == -98);
    REQUIRE(p16[1] == -97);
    REQUIRE(p16[2] == -86);
    REQUIRE(p16[3] == 56);
    REQUIRE(p16[4] == 12);
    REQUIRE(reader.read<int32_t>() == -9);

    auto err = expect_error<std::overflow_error>([&] {
        reader.read<int32_t>();
    });
    REQUIRE(err.what() == "Data terminates unexpectedly"s);
    err = expect_error<std::overflow_error>([&] {
        reader.read_string();
    });
    REQUIRE(err.what() == "Data terminates unexpectedly"s);
    err = expect_error<std::overflow_error>([&] {
        reader.read_array<int8_t>(4);
    });
    REQUIRE(err.what() == "Data terminates unexpectedly"s);
    auto p8 = reader.read_array<int8_t>(3);
    REQUIRE(p8[0] == 1);
    REQUIRE(p8[1] == 2);
    REQUIRE(p8[2] == 3);
}
