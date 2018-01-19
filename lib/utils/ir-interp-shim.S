/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

// Linux only for now
#if defined(__x86_64__) || defined(__x86_64)
        .text
        .hidden nacs_exefunc_cb
        .globl  nacs_exefunc_cb
        .type   nacs_exefunc_cb, @function
nacs_exefunc_cb:
        .cfi_startproc
        pushq   %rbp
        .cfi_def_cfa_offset 16
        .cfi_offset 6, -16
        movq    %rsp, %rbp
        .cfi_def_cfa_register 6
        // slot sizes (aligned to 16 bytes)
        movl    (%rdi), %r11d
        subq    %r11, %rsp
        // Load nargs
        movl    4(%rdi), %r11d
        testq   %r11, %r11
        jle     .Lskip
.Lloop:
        movsbq  4(%rdi,%r11,4), %rax // Load arguments specifiers
        cmpq    $0, %rax
        cmoveq  %rsi, %rax
        je      .Lint_end
        cmpq    $1, %rax
        cmoveq  %rcx, %rax
        je      .Lint_end
        cmpq    $2, %rax
        cmoveq  %rdx, %rax
        je      .Lint_end
        cmpq    $3, %rax
        cmoveq  %r8, %rax
        je      .Lint_end
        cmpq    $4, %rax
        cmoveq  %r9, %rax
        je      .Lint_end

        cmpq    $5, %rax
        movsd   %xmm0, %xmm8
        je      .Lfloat_end
        cmpq    $6, %rax
        movsd   %xmm1, %xmm8
        je      .Lfloat_end
        cmpq    $7, %rax
        movsd   %xmm2, %xmm8
        je      .Lfloat_end
        cmpq    $8, %rax
        movsd   %xmm3, %xmm8
        je      .Lfloat_end
        cmpq    $9, %rax
        movsd   %xmm4, %xmm8
        je      .Lfloat_end
        cmpq    $10, %rax
        movsd   %xmm5, %xmm8
        je      .Lfloat_end
        cmpq    $11, %rax
        movsd   %xmm6, %xmm8
        je      .Lfloat_end
        subq    $12, %rax
        movsd   %xmm7, %xmm8
        je      .Lfloat_end

        movq    (%rbp,%rax), %rax
.Lint_end:
        movq    %rax, -8(%rsp,%r11,8)
        jmp     .Lloop_end
.Lfloat_end:
        movsd   %xmm8, -8(%rsp,%r11,8)
.Lloop_end:
        subq    $1, %r11
        jne     .Lloop
.Lskip:
        movq    %rsp, %rsi
        call    nacs_exefunc_real
        // No need to be conditional here...
        movq    %rdx, %xmm0
        movq    %rdx, %rax
        movq    %rbp, %rsp
        popq    %rbp
        .cfi_def_cfa 7, 8
        retq
        .cfi_endproc

        .size   nacs_exefunc_cb, . - nacs_exefunc_cb

#elif defined(__i386) || defined(__i386__)
#define _CPU_X86_
#elif defined(__aarch64__)
#define _CPU_AARCH64_
#elif defined(__arm__)
#define _CPU_ARM_
#endif
