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

#include "../../lib/nacs-utils/dlload.h"
#include "../../lib/nacs-utils/number.h"

#include "nacs_lib_names.h"

#include <zmq.h>

#include <catch2/catch.hpp>

using namespace NaCs;

NACS_EXPORT() int global_symbol = 1;

TEST_CASE("dl") {
    // Use libzmq since it's a required dependency.
    auto hdl = DL::open(NACS_ZMQ_NAME, DL::GLOBAL | DL::NOW);
    REQUIRE(hdl);
    auto pversion = (void(*)(int*, int*, int*))DL::sym(hdl, "zmq_version");
    REQUIRE(pversion);
    int major = -1;
    int minor = -1;
    int patch = -1;
    pversion(&major, &minor, &patch);
    REQUIRE(major == ZMQ_VERSION_MAJOR);
    REQUIRE(minor == ZMQ_VERSION_MINOR);
    REQUIRE(patch == ZMQ_VERSION_PATCH);
    DL::close(hdl);

    auto *pgv = DL::sym(nullptr, "global_symbol");
    REQUIRE(pgv == (void*)&global_symbol);
}
