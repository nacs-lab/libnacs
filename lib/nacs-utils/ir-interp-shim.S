/*************************************************************************
 *   Copyright (c) 2016 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#if NACS_OS_DARWIN
# define CFUNC(name) _##name
#else
# define CFUNC(name) name
#endif

#if defined(__x86_64__) || defined(__x86_64)
#  if !NACS_OS_WINDOWS
        .text
#    if NACS_OS_DARWIN
        .private_extern CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
#    else
        .hidden CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
        .type   CFUNC(nacs_exefunc_cb), @function
#    endif
CFUNC(nacs_exefunc_cb):
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
        leaq    .labels_cases(%rip), %r10
.Lloop:
        movslq  4(%rdi,%r11,4), %rax // Load arguments specifiers
        cmpq    $12, %rax
        ja      .Lspill_start
        jmp     *(%r10,%rax,8)

.Lcase0:
        movq    %rsi, %rax
        jmp     .Lint_end
.Lcase1:
        movq    %rdx, %rax
        jmp     .Lint_end
.Lcase2:
        movq    %rcx, %rax
        jmp     .Lint_end
.Lcase3:
        movq    %r8, %rax
        jmp     .Lint_end
.Lcase4:
        movq    %r9, %rax
        jmp     .Lint_end
.Lcase5:
        movsd   %xmm0, %xmm8
        jmp     .Lfloat_end
.Lcase6:
        movsd   %xmm1, %xmm8
        jmp     .Lfloat_end
.Lcase7:
        movsd   %xmm2, %xmm8
        jmp     .Lfloat_end
.Lcase8:
        movsd   %xmm3, %xmm8
        jmp     .Lfloat_end
.Lcase9:
        movsd   %xmm4, %xmm8
        jmp     .Lfloat_end
.Lcase10:
        movsd   %xmm5, %xmm8
        jmp     .Lfloat_end
.Lcase11:
        movsd   %xmm6, %xmm8
        jmp     .Lfloat_end
.Lcase12:
        movsd   %xmm7, %xmm8
        jmp     .Lfloat_end

.Lspill_start:
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
        call    CFUNC(nacs_exefunc_real)
        // No need to be conditional here...
        movq    %rdx, %xmm0
        movq    %rdx, %rax
        movq    %rbp, %rsp
        popq    %rbp
        .cfi_def_cfa 7, 8
        retq
        .cfi_endproc

#    if !NACS_OS_DARWIN
        .size   CFUNC(nacs_exefunc_cb), . - CFUNC(nacs_exefunc_cb)

        .section        .data.rel.ro.local,"aw",@progbits
        .align  32
        .type   labels_cases, @object
        .size   labels_cases, 104
#    else
        .section __DATA,__const
        .p2align 3
#    endif
.labels_cases:
        .quad   .Lcase0
        .quad   .Lcase1
        .quad   .Lcase2
        .quad   .Lcase3
        .quad   .Lcase4
        .quad   .Lcase5
        .quad   .Lcase6
        .quad   .Lcase7
        .quad   .Lcase8
        .quad   .Lcase9
        .quad   .Lcase10
        .quad   .Lcase11
        .quad   .Lcase12
#  else
        .text
        .p2align        4,0x90
        .globl          CFUNC(nacs_exefunc_cb)
        .section        .drectve
        .ascii  " -export:"
        .ascii  "nacs_exefunc_cb"
        .section        .text
        .def    CFUNC(nacs_exefunc_cb)
        .scl    2
        .type   32
        .endef
CFUNC(nacs_exefunc_cb):
        .cfi_startproc
        pushq   %rbp
        .cfi_def_cfa_offset 16
        .cfi_offset 6, -16
        movq    %rsp, %rbp
        .cfi_def_cfa_register 6
        // slot sizes (aligned to 16 bytes)
        movl    (%rcx), %r11d
        subq    %r11, %rsp
        // Load nargs
        movl    4(%rcx), %r11d
        testq   %r11, %r11
        jle     .Lskip
.Lloop:
        movslq  4(%rcx,%r11,4), %rax // Load arguments specifiers
        cmpq    $0, %rax
        cmoveq  %rdx, %rax
        je      .Lint_end
        cmpq    $1, %rax
        cmoveq  %r8, %rax
        je      .Lint_end
        cmpq    $2, %rax
        cmoveq  %r9, %rax
        je      .Lint_end

        cmpq    $3, %rax
        movsd   %xmm1, %xmm4
        je      .Lfloat_end
        cmpq    $4, %rax
        movsd   %xmm2, %xmm4
        je      .Lfloat_end
        cmpq    $5, %rax
        movsd   %xmm3, %xmm4
        je      .Lfloat_end

        movq    (%rbp,%rax), %rax
.Lint_end:
        movq    %rax, -8(%rsp,%r11,8)
        jmp     .Lloop_end
.Lfloat_end:
        movsd   %xmm4, -8(%rsp,%r11,8)
.Lloop_end:
        subq    $1, %r11
        jne     .Lloop
.Lskip:
        movq    %rsp, %r8
        movq    %rcx, %rdx
        subq    $16, %rsp
        movq    %rsp, %rcx
        subq    $32, %rsp
        call    CFUNC(nacs_exefunc_real)
        // No need to be conditional here...
        movq    40(%rsp), %xmm0
        movq    40(%rsp), %rax
        movq    %rbp, %rsp
        popq    %rbp
        .cfi_def_cfa 7, 8
        retq
        .cfi_endproc

        .end
#  endif
#elif defined(__i386) || defined(__i386__)
        // Linux only for now
        .text
        .hidden CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
        .type   CFUNC(nacs_exefunc_cb), @function
CFUNC(nacs_exefunc_cb): // sp: 12/16
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
        call    CFUNC(nacs_exefunc_real)
        // No need to be conditional here...
        movl    12(%esp), %eax
        fldl    12(%esp)
        movl    %ebp, %esp
        popl    %ebp
        .cfi_def_cfa 7, 8
        retl
        .cfi_endproc

        .size   CFUNC(nacs_exefunc_cb), . - CFUNC(nacs_exefunc_cb)
#elif defined(__aarch64__)
        .text
#  if NACS_OS_DARWIN
        .private_extern CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
        .p2align 3
#  else
        .hidden CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
        .type   CFUNC(nacs_exefunc_cb), #function
#  endif
CFUNC(nacs_exefunc_cb):
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
        ldur    w9, [x9, #4] // Load arguments specifiers
        cmp     x9, #0
        mov     x11, x1
        b.eq    .Lint_end
        cmp     x9, #1
        mov     x11, x2
        b.eq    .Lint_end
        cmp     x9, #2
        mov     x11, x3
        b.eq    .Lint_end
        cmp     x9, #3
        mov     x11, x4
        b.eq    .Lint_end
        cmp     x9, #4
        mov     x11, x5
        b.eq    .Lint_end
        cmp     x9, #5
        mov     x11, x6
        b.eq    .Lint_end
        cmp     x9, #6
        mov     x11, x7
        b.eq    .Lint_end

        cmp     x9, #7
        mov     v16.16b, v0.16b
        b.eq    .Lfloat_end
        cmp     x9, #8
        mov     v16.16b, v1.16b
        b.eq    .Lfloat_end
        cmp     x9, #9
        mov     v16.16b, v2.16b
        b.eq    .Lfloat_end
        cmp     x9, #10
        mov     v16.16b, v3.16b
        b.eq    .Lfloat_end
        cmp     x9, #11
        mov     v16.16b, v4.16b
        b.eq    .Lfloat_end
        cmp     x9, #12
        mov     v16.16b, v5.16b
        b.eq    .Lfloat_end
        cmp     x9, #13
        mov     v16.16b, v6.16b
        b.eq    .Lfloat_end
        cmp     x9, #14
        mov     v16.16b, v7.16b
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
        bl      CFUNC(nacs_exefunc_real)
        // No need to be conditional here...
        mov     x0, x1
        fmov    d0, x1
        mov     sp, x29
        ldp     x29, x30, [sp], #16
        .cfi_def_cfa 31, 0
        ret
        .cfi_endproc

#  if !NACS_OS_DARWIN
        .size   CFUNC(nacs_exefunc_cb), . - CFUNC(nacs_exefunc_cb)
#  endif
#elif defined(__arm__)
        .text
        .hidden CFUNC(nacs_exefunc_cb)
        .globl  CFUNC(nacs_exefunc_cb)
        .type   CFUNC(nacs_exefunc_cb), #function
        .arm
        .fpu vfp
CFUNC(nacs_exefunc_cb):
        .fnstart
        push    {r8, r10, fp, lr}
        .save   {r8, r10, fp, lr}
        vpush.64 {d8}
        .vsave  {d8}
        mov     fp, sp
        .setfp  fp, sp
        ldm     r0, {r10, r12}
        // slot sizes (aligned to 16 bytes)
        sub     sp, sp, r10
        cmp     r12, #1
        blt     .Lskip
.Lloop:
        add     r10, r0, r12, lsl #2
        add     r8, sp, r12, lsl #3
        ldr     r10, [r10, #4] // Load arguments specifiers

        cmp     r10, #11
        // Each case below uses 2 instructions (8 bytes)
        // Add that offset to the PC direction.
        // The `pc` value during this instruction is 8 bytes ahead
        // so a `r10` of 0 will skip the next instruction.
        addlo   pc, pc, r10, lsl #3
        b       .Lspill_start

        mov     r10, r1
        b       .Lint_end
        mov     r10, r2
        b       .Lint_end
        mov     r10, r3
        b       .Lint_end

        vmov.f64 d8, d0
        b       .Lfloat_end
        vmov.f64 d8, d1
        b       .Lfloat_end
        vmov.f64 d8, d2
        b       .Lfloat_end
        vmov.f64 d8, d3
        b       .Lfloat_end
        vmov.f64 d8, d4
        b       .Lfloat_end
        vmov.f64 d8, d5
        b       .Lfloat_end
        vmov.f64 d8, d6
        b       .Lfloat_end
        vmov.f64 d8, d7
        b       .Lfloat_end

.Lspill_start:
        cmp     r10, #0
        blt     .Lfloat_spill
        add     r10, fp, r10
        ldr     r10, [r10]
.Lint_end:
        str     r10, [r8, #-8]
        b       .Lloop_end
.Lfloat_spill:
        bic     r10, r10, #0x80000000
        add     r10, fp, r10
        vldr    d8, [r10]
.Lfloat_end:
        vstr    d8, [r8, #-8]
.Lloop_end:
        subs    r12, r12, #1
        bne    .Lloop
.Lskip:
        mov     r2, sp
        mov     r1, r0
        sub     sp, sp, #24 // Make sure the stack is aligned
        mov     r0, sp
        bl      CFUNC(nacs_exefunc_real)
        // No need to be conditional here...
        ldr     r0, [sp, #8]
        vldr.64 d0, [sp, #8]
        mov     sp, fp
        vpop.64 {d8}
        pop     {r8, r10, fp, pc}
        .fnend

        .size   CFUNC(nacs_exefunc_cb), . - CFUNC(nacs_exefunc_cb)
#else
#  error "Unsupported architecture"
#endif
