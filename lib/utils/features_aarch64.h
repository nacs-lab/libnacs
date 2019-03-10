/*************************************************************************
 *   Copyright (c) 2017 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

// AArch64 features definition
// hwcap
NACS_FEATURE_DEF(crypto, 3, 0)
NACS_FEATURE_DEF(crc, 7, 0)
NACS_FEATURE_DEF(lse, 8, 0) // ARMv8.1-Atomics
NACS_FEATURE_DEF(fullfp16, 9, 0)
NACS_FEATURE_DEF(rdm, 12, 0) // ARMv8.1-SIMD
NACS_FEATURE_DEF(jscvt, 13, UINT32_MAX) // Linux Kernel HWCAP name
NACS_FEATURE_DEF(fcma, 14, UINT32_MAX) // Linux Kernel HWCAP name
NACS_FEATURE_DEF(rcpc, 15, 60000)
NACS_FEATURE_DEF(dcpop, 16, UINT32_MAX) // Linux Kernel HWCAP name
// NACS_FEATURE_DEF(dotprod, ???, 60000) // ARMv8.2-DotProd
// NACS_FEATURE_DEF(ras, ???, 0)
// NACS_FEATURE_DEF(sve, ???, UINT32_MAX)

// hwcap2
// NACS_FEATURE_DEF(?, 32 + ?, 0)

// custom bits to match llvm model
NACS_FEATURE_DEF(v8_1a, 32 * 2 + 0, 0)
NACS_FEATURE_DEF(v8_2a, 32 * 2 + 1, 0)
NACS_FEATURE_DEF(v8_3a, 32 * 2 + 2, 60000)
// NACS_FEATURE_DEF(v8_4a, 32 * 2 + 3, ???)
