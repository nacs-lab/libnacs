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
        cmpq    $12, %rax
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
        .text
        .hidden nacs_exefunc_cb
        .globl  nacs_exefunc_cb
        .type   nacs_exefunc_cb, @function
nacs_exefunc_cb: // sp: 12/16
        .cfi_startproc
        pushl   %ebp // sp: 8/16
        .cfi_def_cfa_offset 8
        .cfi_offset 6, -8
        movl    %esp, %ebp
        .cfi_def_cfa_register 6
        movl    8(%esp), %eax  // arg1 (data)
        // slot sizes (aligned to 16 bytes)
        subl    (%eax), %esp // sp: 8/16
        // Load nargs
        movl    4(%eax), %ecx
        testl   %ecx, %ecx
        jle     .Lskip
.Lloop:
        movl    4(%eax,%ecx,4), %edx // Load arguments specifiers
        cmpl    $0, %edx
        js      .Lfloat
        movl    (%ebp,%edx), %edx
        movl    %edx, -8(%esp,%ecx,8)
        jmp     .Lloop_end

.Lfloat:
        andl    $0x7fffffff, %edx
        movsd   (%ebp,%edx), %xmm0
        movsd   %xmm0, -8(%esp,%ecx,8)

.Lloop_end:
        subl    $1, %ecx
        jne     .Lloop

.Lskip:
        movl    %esp, %ecx
        subl    $12, %esp // sp: 12/16
        movl    %esp, %edx
        pushl   %ecx // sp: 8/16
        pushl   %eax // sp: 4/16
        pushl   %edx
        call    nacs_exefunc_real
        // No need to be conditional here...
        movl    12(%esp), %eax
        fldl    12(%esp)
        movl    %ebp, %esp
        popl    %ebp
        .cfi_def_cfa 7, 8
        retl
        .cfi_endproc

        .size   nacs_exefunc_cb, . - nacs_exefunc_cb
#elif defined(__aarch64__)
        .text
        .hidden nacs_exefunc_cb
        .globl  nacs_exefunc_cb
        .type   nacs_exefunc_cb, #function
nacs_exefunc_cb:
        .cfi_startproc
        stp     x29, x30, [sp, #-16]!
        .cfi_def_cfa_offset 16
        .cfi_offset 29, -16
        .cfi_offset 30, -8
        mov     x29, sp
        .cfi_def_cfa_register 29
        ldp     w9, w10, [x0]
        // slot sizes (aligned to 16 bytes)
        sub     sp, sp, x9
        cmp     x10, #1
        b.lt    .Lskip
.Lloop:
        add     x9, x0, x10, lsl #2
        ldur    w9, [x9, #-8] // Load arguments specifiers
        cmp     x9, #0
        mov     x1, x11
        b.eq    .Lint_end
        cmp     x9, #1
        mov     x2, x11
        b.eq    .Lint_end
        cmp     x9, #2
        mov     x3, x11
        b.eq    .Lint_end
        cmp     x9, #3
        mov     x4, x11
        b.eq    .Lint_end
        cmp     x9, #4
        mov     x5, x11
        b.eq    .Lint_end
        cmp     x9, #5
        mov     x6, x11
        b.eq    .Lint_end
        cmp     x9, #6
        mov     x7, x11
        b.eq    .Lint_end

        cmp     x9, #7
        mov     v0.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #8
        mov     v1.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #9
        mov     v2.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #10
        mov     v3.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #11
        mov     v4.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #12
        mov     v5.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #13
        mov     v6.16b, v16.16b
        b.eq    .Lfloat_end
        cmp     x9, #14
        mov     v7.16b, v16.16b
        b.eq    .Lfloat_end

        add     x9, x29, x9
        ldur    x11, [x9]
.Lint_end:
        add     x9, sp, x10, lsl #3
        stur    x11, [x9, #-8]
        b       .Lloop_end
.Lfloat_end:
        add     x9, sp, x10, lsl #3
        str     d16, [x9, #-8]
.Lloop_end:
        subs    x10, x10, #1
        b.ne    .Lloop
.Lskip:
        mov     x1, sp
        sub     sp, sp, #16
        mov     x8, sp
        bl      nacs_exefunc_real
        // No need to be conditional here...
        ldur    x0, [sp, #8]
        ldr     d0, [sp, #8]
        mov     sp, x29
        ldp     x29, x30, [sp], #16
        .cfi_def_cfa 31, 0
        ret
        .cfi_endproc

        .size   nacs_exefunc_cb, . - nacs_exefunc_cb

#elif defined(__arm__)

#else
#  error "Unsupported architecture"
#endif
