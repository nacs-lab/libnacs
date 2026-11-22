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

#ifndef _NACS_TEST_BLOCK_HELPER_H_
#define _NACS_TEST_BLOCK_HELPER_H_

#include "../lib/nacs-utils/macros.h"

namespace {

// Runs the callable passed in immediately.
struct TestBlockRunner {
    template<typename F>
    void operator+(F &&f) const
    {
        f();
    }
};

}

// Run the following block in a separate non-inlined function.
// This keeps each function small, which reduces the compile time
// (in particular for the variable tracking for debug info) significantly.
// Usage: `TEST_BLOCK { ... };`
#define TEST_BLOCK TestBlockRunner{} + [&] () NACS_NOINLINE

#endif
