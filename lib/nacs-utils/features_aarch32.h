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

// Copy values `from arch/arm/include/uapi/asm/hwcap.h` from linux kernel source tree
// and match LLVM names.

// LLVM features in `llvm/lib/Target/ARM/ARM.td`

// AArch32 features definition
// hwcap
NACS_FEATURE_DEF(neon, 12, 0)
NACS_FEATURE_DEF(vfp3, 13, 0)
// NACS_FEATURE_DEF(vfpv3d16, 14, 0) // d16
NACS_FEATURE_DEF(vfp4, 16, 0)
NACS_FEATURE_DEF_NAME(hwdiv_arm, 17, 0, "hwdiv-arm")
NACS_FEATURE_DEF(hwdiv, 18, 0)
NACS_FEATURE_DEF(d32, 19, 0) // -d16

// hwcap2
NACS_FEATURE_DEF(crypto, 32 + 0, 0)
NACS_FEATURE_DEF(crc, 32 + 4, 0)
// NACS_FEATURE_DEF(ras, 32 + ???, 0)
// NACS_FEATURE_DEF(fullfp16, 32 + ???, 0)

// custom bits to match llvm model
NACS_FEATURE_DEF(aclass, 32 * 2 + 0, 0)
NACS_FEATURE_DEF(rclass, 32 * 2 + 1, 0)
NACS_FEATURE_DEF(mclass, 32 * 2 + 2, 0)
NACS_FEATURE_DEF(v7, 32 * 2 + 3, 0)
NACS_FEATURE_DEF(v8, 32 * 2 + 4, 0)
NACS_FEATURE_DEF(v8_1a, 32 * 2 + 5, 0)
NACS_FEATURE_DEF(v8_2a, 32 * 2 + 6, 0)
NACS_FEATURE_DEF(v8_3a, 32 * 2 + 7, 0)
NACS_FEATURE_DEF(v8_m_main, 32 * 2 + 8, 0)
NACS_FEATURE_DEF(v8_4a, 32 * 2 + 9, 0)
NACS_FEATURE_DEF(v8_5a, 32 * 2 + 10, 0)
NACS_FEATURE_DEF(v8_6a, 32 * 2 + 11, 110000)
