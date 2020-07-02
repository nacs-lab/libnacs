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

// X86 features definition
// EAX=1: ECX
NACS_FEATURE_DEF(sse3, 0, 0)
NACS_FEATURE_DEF(pclmul, 1, 0)
NACS_FEATURE_DEF(ssse3, 9, 0)
NACS_FEATURE_DEF(fma, 12, 0)
NACS_FEATURE_DEF(cx16, 13, 0)
NACS_FEATURE_DEF_NAME(sse41, 19, 0, "sse4.1")
NACS_FEATURE_DEF_NAME(sse42, 20, 0, "sse4.2")
NACS_FEATURE_DEF(movbe, 22, 0)
NACS_FEATURE_DEF(popcnt, 23, 0)
NACS_FEATURE_DEF(aes, 25, 0)
NACS_FEATURE_DEF(xsave, 26, 0)
NACS_FEATURE_DEF(avx, 28, 0)
NACS_FEATURE_DEF(f16c, 29, 0)
NACS_FEATURE_DEF(rdrnd, 30, 0)

// EAX=1: EDX
// NACS_FEATURE_DEF(, 32 + ?, ????)

// EAX=7,ECX=0: EBX
NACS_FEATURE_DEF(fsgsbase, 32 * 2 + 0, 0)
// NACS_FEATURE_DEF(sgx, 32 * 2 + 2, 0) // Disable for now since it's very hard to detect
NACS_FEATURE_DEF(bmi, 32 * 2 + 3, 0)
// NACS_FEATURE_DEF(hle, 32 * 2 + 4, 0) // Not used and gone in LLVM 5.0
NACS_FEATURE_DEF(avx2, 32 * 2 + 5, 0)
NACS_FEATURE_DEF(bmi2, 32 * 2 + 8, 0)
// NACS_FEATURE_DEF(invpcid, 32 * 2 + 10, 0) // Priviledged instruction
NACS_FEATURE_DEF(rtm, 32 * 2 + 11, 0)
// NACS_FEATURE_DEF(mpx, 32 * 2 + 14, 0) // Deprecated in LLVM 10.0
NACS_FEATURE_DEF(avx512f, 32 * 2 + 16, 0)
NACS_FEATURE_DEF(avx512dq, 32 * 2 + 17, 0)
NACS_FEATURE_DEF(rdseed, 32 * 2 + 18, 0)
NACS_FEATURE_DEF(adx, 32 * 2 + 19, 0)
// NACS_FEATURE_DEF(smap, 32 * 2 + 20, 0) // Not used and gone in LLVM 5.0
NACS_FEATURE_DEF(avx512ifma, 32 * 2 + 21, 0)
// NACS_FEATURE_DEF(pcommit, 32 * 2 + 22, 0) // Deprecated
NACS_FEATURE_DEF(clflushopt, 32 * 2 + 23, 0)
NACS_FEATURE_DEF(clwb, 32 * 2 + 24, 0)
NACS_FEATURE_DEF(avx512pf, 32 * 2 + 26, 0)
NACS_FEATURE_DEF(avx512er, 32 * 2 + 27, 0)
NACS_FEATURE_DEF(avx512cd, 32 * 2 + 28, 0)
NACS_FEATURE_DEF(sha, 32 * 2 + 29, 0)
NACS_FEATURE_DEF(avx512bw, 32 * 2 + 30, 0)
NACS_FEATURE_DEF(avx512vl, 32 * 2 + 31, 0)

// EAX=7,ECX=0: ECX
NACS_FEATURE_DEF(prefetchwt1, 32 * 3 + 0, 0)
NACS_FEATURE_DEF(avx512vbmi, 32 * 3 + 1, 0)
NACS_FEATURE_DEF(pku, 32 * 3 + 4, 0) // ospke
NACS_FEATURE_DEF(waitpkg, 32 * 3 + 5, 0)
NACS_FEATURE_DEF(avx512vbmi2, 32 * 3 + 6, 0)
NACS_FEATURE_DEF(shstk, 32 * 3 + 7, 0)
NACS_FEATURE_DEF(gfni, 32 * 3 + 8, 0)
NACS_FEATURE_DEF(vaes, 32 * 3 + 9, 0)
NACS_FEATURE_DEF(vpclmulqdq, 32 * 3 + 10, 0)
NACS_FEATURE_DEF(avx512vnni, 32 * 3 + 11, 0)
NACS_FEATURE_DEF(avx512bitalg, 32 * 3 + 12, 0)
NACS_FEATURE_DEF(avx512vpopcntdq, 32 * 3 + 14, 0)
NACS_FEATURE_DEF(rdpid, 32 * 3 + 22, 0)
NACS_FEATURE_DEF(cldemote, 32 * 3 + 25, 0)
NACS_FEATURE_DEF(movdiri, 32 * 3 + 27, 0)
NACS_FEATURE_DEF(movdir64b, 32 * 3 + 28, 0)
NACS_FEATURE_DEF(enqcmd, 32 * 3 + 29, 90000)

// EAX=7,ECX=0: EDX
// NACS_FEATURE_DEF(avx5124vnniw, 32 * 4 + 2, ?????)
// NACS_FEATURE_DEF(avx5124fmaps, 32 * 4 + 3, ?????)
NACS_FEATURE_DEF(avx512vp2intersect, 32 * 4 + 8, 90000)
NACS_FEATURE_DEF(serialize, 32 * 4 + 14, 110000)
NACS_FEATURE_DEF(tsxldtrk, 32 * 4 + 16, 110000)
NACS_FEATURE_DEF(pconfig, 32 * 4 + 18, 0)

// EAX=0x80000001: ECX
// ignore sahf on 32bit x86 since it is required
NACS_FEATURE_DEF(sahf, 32 * 5 + 0, 0)
NACS_FEATURE_DEF(lzcnt, 32 * 5 + 5, 0)
NACS_FEATURE_DEF(sse4a, 32 * 5 + 6, 0)
NACS_FEATURE_DEF(prfchw, 32 * 5 + 8, 0)
NACS_FEATURE_DEF(xop, 32 * 5 + 11, 0)
NACS_FEATURE_DEF(lwp, 32 * 5 + 15, 0)
NACS_FEATURE_DEF(fma4, 32 * 5 + 16, 0)
NACS_FEATURE_DEF(tbm, 32 * 5 + 21, 0)
NACS_FEATURE_DEF(mwaitx, 32 * 5 + 29, 0)

// EAX=0x80000001: EDX
// 3dnow is here but we don't care...
// NACS_FEATURE_DEF(, 32 * 6 + ?, ?????)

// EAX=0xd: EAX
NACS_FEATURE_DEF(xsaveopt, 32 * 7 + 0, 0)
NACS_FEATURE_DEF(xsavec, 32 * 7 + 1, 0)
NACS_FEATURE_DEF(xsaves, 32 * 7 + 3, 0)

// EAX=0x80000008: EBX
NACS_FEATURE_DEF(clzero, 32 * 8 + 0, 0)
NACS_FEATURE_DEF(wbnoinvd, 32 * 8 + 9, 0)

// EAX=7,ECX=1: EAX
NACS_FEATURE_DEF(avx512bf16, 32 * 9 + 5, 90000)

// EAX=0x14,ECX=0: EBX
NACS_FEATURE_DEF(ptwrite, 32 * 10 + 4, 0)
