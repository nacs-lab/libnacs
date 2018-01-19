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
        push    %rbp
        mov     %rsp, %rbp
        // slot sizes (aligned to 16 bytes)
        mov     (%rdi), %r11d
        sub     %r11, %rsp
        // Load nargs
        mov     4(%rdi), %r11d
        test    %r11, %r11
        jle     .Lskip
.Lloop:
        mov     4(%rdi,%r11,4), %eax // Load arguments specifiers
        cmp     $0, %rax
        cmove   %rsi, %rax
        je      .Lint_end
        cmp     $1, %rax
        cmove   %rcx, %rax
        je      .Lint_end
        cmp     $2, %rax
        cmove   %rdx, %rax
        je      .Lint_end
        cmp     $3, %rax
        cmove   %r8, %rax
        je      .Lint_end
        cmp     $4, %rax
        cmove   %r9, %rax
        je      .Lint_end

        cmp     $5, %rax
        movsd   %xmm0, %xmm8
        je      .Lfloat_end
        cmp     $6, %rax
        movsd   %xmm1, %xmm8
        je      .Lfloat_end
        cmp     $7, %rax
        movsd   %xmm2, %xmm8
        je      .Lfloat_end
        cmp     $8, %rax
        movsd   %xmm3, %xmm8
        je      .Lfloat_end
        cmp     $9, %rax
        movsd   %xmm4, %xmm8
        je      .Lfloat_end
        cmp     $10, %rax
        movsd   %xmm5, %xmm8
        je      .Lfloat_end
        cmp     $11, %rax
        movsd   %xmm6, %xmm8
        je      .Lfloat_end
        sub     $12, %rax
        movsd   %xmm7, %xmm8
        je      .Lfloat_end

        mov     (%rbp,%rax), %rax
.Lint_end:
        mov     %rax, -8(%rsp,%r11,8)
        jmp     .Lloop_end
.Lfloat_end:
        movsd   %xmm8, -8(%rsp,%r11,8)
.Lloop_end:
        sub     $1, %r11
        jne     .Lloop
.Lskip:
        mov     %rsp, %rsi
        call    nacs_exefunc_real
        // No need to be conditional here...
        movq    %rdx, %xmm0
        mov     %rdx, %rax
        mov     %rbp, %rsp
        pop     %rbp
        ret
        .cfi_endproc

        .size   nacs_exefunc_cb, . - nacs_exefunc_cb

#elif defined(__i386) || defined(__i386__)
#define _CPU_X86_
#elif defined(__aarch64__)
#define _CPU_AARCH64_
#elif defined(__arm__)
#define _CPU_ARM_
#endif
