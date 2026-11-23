/*************************************************************************
 *   Copyright (c) 2026 - 2026 Yichao Yu <yyc1992@gmail.com>             *
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

// The Catch2 implementation and `main` are expensive to compile
// (~15 s and ~1 GB of memory with optimization and debug info)
// so compile them once here instead of in every test.
#define CATCH_CONFIG_MAIN

// catch.hpp uses `std::nothrow` without including `<new>`.
// This is otherwise only included by other headers indirectly.
#include <new>

#include <catch2/catch.hpp>
