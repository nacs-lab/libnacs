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

#include "macros.h"

#if (defined(__x86_64__) || defined(__x86_64)) && \
        defined(__GNUC__) && !defined(__clang__) && NACS_OS_WINDOWS
        .text
        .def    _ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32
        .scl    2
        .type   32
        .endef
        .section    .rdata,"dr"
        .p2align    5
.Ldouble1_v4:
        .quad   4607182418800017408     # double 1
        .quad   4607182418800017408     # double 1
        .quad   4607182418800017408     # double 1
        .quad   4607182418800017408     # double 1
        .text
        .globl  _ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32
        .p2align    4, 0x90
        // NaCs::linearInterpolate2_sse2(double __vector(2), unsigned int, double const*)
        .seh_proc   _ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32
_ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32:
        // The wrapper function below `_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd`,
        // assumes that this function does not touch `%rcx` and does not use the stack.
        // `_ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64`
        // assumes that this function does not touch `%rcx`.
        .seh_endprologue
        pxor        %xmm5, %xmm5
        comisd      %xmm0, %xmm5
        jnb         .L2sse2_v1_undflw
        movsd       .Ldouble1_v4(%rip), %xmm2
        leal        -1(%rdx), %eax
        comisd      %xmm2, %xmm0
        jnb         .L2sse2_v1_ovrflw
        pxor        %xmm3, %xmm3
        pxor        %xmm1, %xmm1
        cvtsi2sdq   %rax, %xmm3
        mulsd       %xmm0, %xmm3
        cvttsd2si   %xmm3, %eax
        cvtsi2sd    %eax, %xmm1
        cltq
        movsd       8(%r8,%rax,8), %xmm4
        subsd       %xmm1, %xmm3
        mulsd       %xmm3, %xmm4
        subsd       %xmm3, %xmm2
        mulsd       (%r8,%rax,8), %xmm2
        addsd       %xmm2, %xmm4
.L2sse2_v2:
        unpckhpd    %xmm0, %xmm0
        comisd      %xmm0, %xmm5
        jnb         .L2sse2_v2_undflw
.L2sse2_v2_not_undflw:
        movsd       .Ldouble1_v4(%rip), %xmm2
        subl        $1, %edx
        comisd      %xmm2, %xmm0
        jnb         .L2sse2_v2_ovrflw
        pxor        %xmm1, %xmm1
        cvtsi2sdq   %rdx, %xmm1
        mulsd       %xmm0, %xmm1
        pxor        %xmm0, %xmm0
        cvttsd2si   %xmm1, %eax
        cvtsi2sd    %eax, %xmm0
        cltq
        subsd       %xmm0, %xmm1
        movsd       8(%r8,%rax,8), %xmm0
        subsd       %xmm1, %xmm2
        mulsd       %xmm1, %xmm0
        movsd       (%r8,%rax,8), %xmm1
        mulsd       %xmm2, %xmm1
        addsd       %xmm0, %xmm1
        movapd      %xmm4, %xmm0
        unpcklpd    %xmm1, %xmm0
        ret
        .p2align 4,,10
        .p2align 3
.L2sse2_v1_ovrflw:
        unpckhpd    %xmm0, %xmm0
        movsd       (%r8,%rax,8), %xmm4
        comisd      %xmm0, %xmm5
        jb          .L2sse2_v2_not_undflw
.L2sse2_v2_undflw:
        movsd       (%r8), %xmm1
        movapd      %xmm4, %xmm0
        unpcklpd    %xmm1, %xmm0
        ret
        .p2align 4,,10
        .p2align 3
.L2sse2_v2_ovrflw:
        movsd       (%r8,%rdx,8), %xmm1
        movapd      %xmm4, %xmm0
        unpcklpd    %xmm1, %xmm0
        ret
        .p2align 4,,10
        .p2align 3
.L2sse2_v1_undflw:
        movsd       (%r8), %xmm4
        jmp         .L2sse2_v2
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32
        .p2align    4, 0x90
        // NaCs::linearInterpolate2_avx2(double __vector(2), unsigned int, double const*)
        .seh_proc   _ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32
_ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32:
        // The wrapper function below `_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd`,
        // and `_ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64`
        // assumes that this function does not touch `%rcx`.
        // Using the scratch area
        vmovaps     %xmm6, 8(%rsp)
        .seh_savexmm %xmm6, 8
        .seh_endprologue
        vmovapd     .Ldouble1_v4(%rip), %xmm2
        vxorpd      %xmm3, %xmm3, %xmm3
        addl        $-1, %edx
        vmovddup    (%r8,%rdx,8), %xmm1
        vmaxpd      %xmm3, %xmm0, %xmm3
        vxorpd      %xmm4, %xmm4, %xmm4
        vcvtsi2sdq  %rdx, %xmm4, %xmm4
        vmovddup    %xmm4, %xmm4
        vmulpd      %xmm3, %xmm4, %xmm3
        vcvttpd2dq  %xmm3, %xmm4
        vcvtdq2pd   %xmm4, %xmm5
        vcmpltpd    %xmm2, %xmm0, %xmm0
        vsubpd      %xmm5, %xmm3, %xmm3
        vmovapd     %xmm1, %xmm5
        vmovapd     %xmm0, %xmm6
        vgatherdpd  %xmm6, (%r8,%xmm4,8), %xmm5
        vgatherdpd  %xmm0, 8(%r8,%xmm4,8), %xmm1
        vsubpd      %xmm3, %xmm2, %xmm0
        vmulpd      %xmm0, %xmm5, %xmm0
        vfmadd231pd %xmm1, %xmm3, %xmm0
        vmovaps     8(%rsp), %xmm6
        ret
        .seh_endproc

        .def    _ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48
        .p2align    4, 0x90
        // NaCs::linearInterpolate4_avx(double __vector(4), unsigned int, double const*)
        .seh_proc _ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48
_ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48:
        // The wrapper function below `_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd`,
        // and `_ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112`
        // assumes that this function does not touch `%r11` and does not use the stack.
        .seh_endprologue
        lea           -0x1(%rdx), %r9d
        vxorpd        %xmm1, %xmm1, %xmm1
        vmovapd       .Ldouble1_v4(%rip), %ymm5
        vcmpltpd      %ymm5, %ymm0, %ymm4
        vmaxpd        %ymm1, %ymm0, %ymm0
        vxorpd        %xmm1, %xmm1, %xmm1
        vcvtsi2sd     %r9, %xmm1, %xmm1
        vmovddup      %xmm1, %xmm1
        vinsertf128   $1, %xmm1, %ymm1, %ymm1
        vmulpd        %ymm0, %ymm1, %ymm1
        vandps        %ymm4, %ymm1, %ymm1
        vcvttpd2dq    %ymm1, %xmm2
        vpextrd       $2, %xmm2, %edx
        vcvtdq2pd     %xmm2, %ymm0
        vmovd         %xmm2, %ecx
        movslq        %edx, %rdx
        vpextrd       $1, %xmm2, %r10d
        vsubpd        %ymm0, %ymm1, %ymm1
        movslq        %ecx, %rcx
        vpextrd       $3, %xmm2, %eax
        vmovsd        8(%r8,%rdx,8), %xmm2
        movslq        %r10d, %r10
        vmovsd        (%r8,%rdx,8), %xmm3
        cltq
        vmovhpd       8(%r8,%rax,8), %xmm2, %xmm0
        vmovsd        8(%r8,%rcx,8), %xmm2
        vmovhpd       8(%r8,%r10,8), %xmm2, %xmm2
        vinsertf128   $1, %xmm0, %ymm2, %ymm2
        vmovhpd       (%r8,%rax,8), %xmm3, %xmm0
        vmovsd        (%r8,%rcx,8), %xmm3
        vmulpd        %ymm1, %ymm2, %ymm2
        vsubpd        %ymm1, %ymm5, %ymm1
        vmovhpd       (%r8,%r10,8), %xmm3, %xmm3
        vinsertf128   $1, %xmm0, %ymm3, %ymm3
        vbroadcastsd  (%r8,%r9,8), %ymm0
        vmulpd        %ymm1, %ymm3, %ymm1
        vaddpd        %ymm1, %ymm2, %ymm1
        vblendvpd     %ymm4, %ymm1, %ymm0, %ymm0
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48
        .p2align    4, 0x90
        // NaCs::linearInterpolate4_avx2(double __vector(4), unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48
_ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48:
        // The wrapper function below `_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd`
        // and `_ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112`
        // assumes that this function does not touch `%rcx`.
        // Using the scratch area
        vmovaps      %xmm6, 8(%rsp)
        .seh_savexmm %xmm6, 8
        .seh_endprologue
        vmovapd      .Ldouble1_v4(%rip), %ymm2
        vxorpd       %xmm3, %xmm3, %xmm3
        addl         $-1, %edx
        vbroadcastsd (%r8,%rdx,8), %ymm1
        vmaxpd       %ymm3, %ymm0, %ymm3
        vxorpd       %xmm4, %xmm4, %xmm4
        vcvtsi2sdq   %rdx, %xmm4, %xmm4
        vbroadcastsd %xmm4, %ymm4
        vmulpd       %ymm3, %ymm4, %ymm3
        vcvttpd2dqy  %ymm3, %xmm4
        vcvtdq2pd    %xmm4, %ymm5
        vcmpltpd     %ymm2, %ymm0, %ymm0
        vsubpd       %ymm5, %ymm3, %ymm3
        vmovapd      %ymm1, %ymm5
        vmovapd      %ymm0, %ymm6
        vgatherdpd   %ymm6, (%r8,%xmm4,8), %ymm5
        vgatherdpd   %ymm0, 8(%r8,%xmm4,8), %ymm1
        vsubpd       %ymm3, %ymm2, %ymm0
        vmulpd       %ymm0, %ymm5, %ymm0
        vfmadd231pd  %ymm1, %ymm3, %ymm0
        vmovaps      8(%rsp), %xmm6
        ret
        .seh_endproc

// Do not enable for avx512 for now since the calling convention is currently the same as
// non-vectorcall as of Feb. 2019. The computer on which we need the asm implementation
// also doesn't support avx512 anyway so it doesn't really matter for now...
//         .def    _ZN4NaCs26linearInterpolate8_avx512fEDv8_djPKd@@80
//         .scl    2
//         .type   32
//         .endef
//         .globl  _ZN4NaCs26linearInterpolate8_avx512fEDv8_djPKd@@80
//         .p2align    4, 0x90
//         // NaCs::linearInterpolate8_avx512f(double __vector(8), unsigned int, double const*)
//         .seh_proc _ZN4NaCs26linearInterpolate8_avx512fEDv8_djPKd@@80
// _ZN4NaCs26linearInterpolate8_avx512fEDv8_djPKd@@80:
//         .seh_endprologue
//         vxorpd        %xmm1, %xmm1, %xmm1
//         leal          -1(%rdx), %eax
//         vbroadcastsd  .Ldouble1_v4(%rip), %zmm2
//         vbroadcastsd  (%r8,%rax,8), %zmm5
//         vcmpltpd      %zmm2, %zmm0, %k1
//         vmaxpd        %zmm1, %zmm0, %zmm0
//         vxorpd        %xmm1, %xmm1, %xmm1
//         vcvtsi2sdq    %rax, %xmm1, %xmm1
//         kmovb         %k1, %k2
//         vbroadcastsd  %xmm1, %zmm1
//         vmulpd        %zmm0, %zmm1, %zmm1
//         vcvttpd2dq    %zmm1, %ymm4
//         vcvtdq2pd     %ymm4, %zmm3
//         vsubpd        %zmm3, %zmm1, %zmm3
//         vmovapd       %zmm5, %zmm1
//         vgatherdpd    (%r8,%ymm4,8), %zmm5{%k2}
//         vgatherdpd    8(%r8,%ymm4,8), %zmm1{%k1}
//         vsubpd        %zmm3, %zmm2, %zmm0
//         vmulpd        %zmm5, %zmm0, %zmm0
//         vfmadd231pd   %zmm1, %zmm3, %zmm0
//         ret
//         .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_sse2(double __vector(2), double __vector(2), double __vector(2), unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64
_ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64:
        // The wrapper function below `_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd`
        // assumes that this function does not touch `%rcx`.
        .seh_endprologue
        movq    40(%rsp), %r8
        subpd   %xmm1, %xmm0
        divpd   %xmm2, %xmm0
        movl    %r9d, %edx
        jmp     _ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_avx2(double __vector(2), double __vector(2), double __vector(2), unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64
_ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64:
        // The wrapper function below `_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd`
        // assumes that this function does not touch `%rcx`.
        .seh_endprologue
        movq    40(%rsp), %r8
        vsubpd  %xmm1, %xmm0, %xmm0
        vdivpd  %xmm2, %xmm0, %xmm0
        movl    %r9d, %edx
        jmp     _ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32
        .seh_endproc

        .def    _ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx(double __vector(4), double __vector(4), double __vector(4), unsigned int, double const*)
        .seh_proc _ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112
_ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112:
        // The wrapper function below `_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd`
        // assumes that this function does not touch `%r11`.
        .seh_endprologue
        movq    40(%rsp), %r8
        vsubpd  %ymm1, %ymm0, %ymm0
        vdivpd  %ymm2, %ymm0, %ymm0
        movl    %r9d, %edx
        jmp     _ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx2(double __vector(4), double __vector(4), double __vector(4), unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112
_ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112:
        // The wrapper function below `_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd`
        // assumes that this function does not touch `%rcx`.
        .seh_endprologue
        movq    40(%rsp), %r8
        vsubpd  %ymm1, %ymm0, %ymm0
        vdivpd  %ymm2, %ymm0, %ymm0
        movl    %r9d, %edx
        jmp     _ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_sse2(double __vector(2)&, double __vector(2) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd
_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd:
        // Assumes that `_ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32`
        // does not touch `%rcx` and doesn't need the stack
        // subq    $8, %rsp
        // .seh_stackalloc 8
        .seh_endprologue
        movaps  (%rdx), %xmm0
        movl    %r8d, %edx
        movq    %r9, %r8
        callq   _ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32
        movaps  %xmm0, (%rcx)
        // addq    $8, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_sse2(double __vector(2)&, double __vector(2) const&, double __vector(2) const&, double __vector(2) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd
_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd:
        // Assumes that `_ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64`
        // does not touch `%rcx`
        subq    $56, %rsp
        .seh_stackalloc 56
        .seh_endprologue
        movl    96(%rsp), %eax
        movq    104(%rsp), %r10
        movaps  (%rdx), %xmm0
        movaps  (%r8), %xmm1
        movaps  (%r9), %xmm2
        movq    %r10, 32(%rsp)
        movl    %eax, %r9d
        callq   _ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64
        movaps  %xmm0, (%rcx)
        addq    $56, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_avx2(double __vector(2)&, double __vector(2) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd
_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd:
        // Assumes that `_ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32`
        // does not touch `%rcx`.
        subq    $40, %rsp
        .seh_stackalloc 40
        .seh_endprologue
        vmovaps (%rdx), %xmm0
        movl    %r8d, %edx
        movq    %r9, %r8
        callq   _ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32
        vmovaps %xmm0, (%rcx)
        addq    $40, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate2_avx2(double __vector(2)&, double __vector(2) const&, double __vector(2) const&, double __vector(2) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd
_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd:
        // Assumes that `_ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64`
        // does not touch `%r10`
        subq    $56, %rsp
        .seh_stackalloc 56
        .seh_endprologue
        movl    96(%rsp), %eax
        movq    104(%rsp), %r10
        vmovaps  (%rdx), %xmm0
        vmovaps  (%r8), %xmm1
        vmovaps  (%r9), %xmm2
        movq    %r10, 32(%rsp)
        movl    %eax, %r9d
        callq   _ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64
        vmovaps  %xmm0, (%rcx)
        addq    $56, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx(double __vector(4)&, double __vector(4) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd
_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd:
        // Assumes that `_ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48`
        // does not touch `%r11` and does not use the stack.
        // Also assumes that the caller knows about avx so we don't need `vzeroupper`
        // subq    $40, %rsp
        // .seh_stackalloc 40
        .seh_endprologue
        movq    %rcx, %r11
        vmovaps (%rdx), %ymm0
        movl    %r8d, %edx
        movq    %r9, %r8
        callq   _ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48
        vmovaps %ymm0, (%r11)
        // addq    $40, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx(double __vector(4)&, double __vector(4) const&, double __vector(4) const&, double __vector(4) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd
_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd:
        // Assumes that `_ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112`
        // does not touch `%r11`.
        // Also assumes that the caller knows about avx so we don't need `vzeroupper`
        subq    $56, %rsp
        .seh_stackalloc 56
        .seh_endprologue
        movq    %rcx, %r11
        movl    96(%rsp), %eax
        movq    104(%rsp), %rcx
        vmovaps (%rdx), %ymm0
        vmovaps (%r8), %ymm1
        vmovaps (%r9), %ymm2
        movq    %rcx, 32(%rsp)
        movl    %eax, %r9d
        callq   _ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112
        vmovaps %ymm0, (%r11)
        addq    $56, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx2(double __vector(4)&, double __vector(4) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd
_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd:
        // Assumes that `_ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48`
        // does not touch `%rcx`.
        // Also assumes that the caller knows about avx so we don't need `vzeroupper`
        subq    $40, %rsp
        .seh_stackalloc 40
        .seh_endprologue
        vmovaps (%rdx), %ymm0
        movl    %r8d, %edx
        movq    %r9, %r8
        callq   _ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48
        vmovaps %ymm0, (%rcx)
        addq    $40, %rsp
        ret
        .seh_endproc

        .def    _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd
        .scl    2
        .type   32
        .endef
        .globl  _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd
        .p2align        4, 0x90
        // NaCs::linearInterpolate4_avx2(double __vector(4)&, double __vector(4) const&, double __vector(4) const&, double __vector(4) const&, unsigned int, double const*)
        .seh_proc _ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd
_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd:
        subq    $56, %rsp
        .seh_stackalloc 56
        .seh_endprologue
        movl    96(%rsp), %eax
        movq    104(%rsp), %r11
        vmovaps (%rdx), %ymm0
        vmovaps (%r8), %ymm1
        vmovaps (%r9), %ymm2
        movq    %r11, 32(%rsp)
        movl    %eax, %r9d
        callq   _ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112
        vmovaps %ymm0, (%rcx)
        addq    $56, %rsp
        ret
        .seh_endproc

        .section    .drectve,"yn"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_sse2EDv2_djPKd@@32"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_avx2EDv2_djPKd@@32"
        .ascii  " -export:_ZN4NaCs22linearInterpolate4_avxEDv4_djPKd@@48"
        .ascii  " -export:_ZN4NaCs23linearInterpolate4_avx2EDv4_djPKd@@48"
        // .ascii  " -export:_ZN4NaCs26linearInterpolate8_avx512fEDv8_djPKd@@80"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_sse2EDv2_dS0_S0_jPKd@@64"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_avx2EDv2_dS0_S0_jPKd@@64"
        .ascii  " -export:_ZN4NaCs22linearInterpolate4_avxEDv4_dS0_S0_jPKd@@112"
        .ascii  " -export:_ZN4NaCs23linearInterpolate4_avx2EDv4_dS0_S0_jPKd@@112"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_jPKd"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_sse2ERDv2_dRKS0_S3_S3_jPKd"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_jPKd"
        .ascii  " -export:_ZN4NaCs23linearInterpolate2_avx2ERDv2_dRKS0_S3_S3_jPKd"
        .ascii  " -export:_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_jPKd"
        .ascii  " -export:_ZN4NaCs22linearInterpolate4_avxERDv4_dRKS0_S3_S3_jPKd"
        .ascii  " -export:_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_jPKd"
        .ascii  " -export:_ZN4NaCs23linearInterpolate4_avx2ERDv4_dRKS0_S3_S3_jPKd"
#endif
