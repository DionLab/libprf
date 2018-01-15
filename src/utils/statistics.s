# mark_description "Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 14.0 Build 20140422";
# mark_description "-I../../build -I../include -restrict -std=c99 -msse4.1 -O3 -S";
	.file "statistics.c"
	.text
..TXTST0:
# -- Begin  _mm_malloc
# mark_begin;
       .align    16,0x90
_mm_malloc:
# parameter 1: %rdi
# parameter 2: %rsi
..B1.1:                         # Preds ..B1.0
..___tag_value__mm_malloc.1:                                    #39.1
        pushq     %rsi                                          #39.1
..___tag_value__mm_malloc.3:                                    #
        movq      %rdi, %rdx                                    #39.1
        movl      $16, %esi                                     #45.7
        lea       (%rsp), %rdi                                  #45.7
        call      posix_memalign                                #45.7
                                # LOE rbx rbp r12 r13 r14 r15 eax
..B1.2:                         # Preds ..B1.1
        testl     %eax, %eax                                    #45.49
        jne       ..B1.4        # Prob 50%                      #45.49
                                # LOE rbx rbp r12 r13 r14 r15
..B1.3:                         # Preds ..B1.2
        movq      (%rsp), %rax                                  #46.12
        popq      %rcx                                          #46.12
..___tag_value__mm_malloc.4:                                    #
        ret                                                     #46.12
..___tag_value__mm_malloc.5:                                    #
                                # LOE
..B1.4:                         # Preds ..B1.2
        xorl      %eax, %eax                                    #48.12
        popq      %rcx                                          #48.12
..___tag_value__mm_malloc.6:                                    #
        ret                                                     #48.12
        .align    16,0x90
..___tag_value__mm_malloc.7:                                    #
                                # LOE
# mark_end;
	.type	_mm_malloc,@function
	.size	_mm_malloc,.-_mm_malloc
	.data
# -- End  _mm_malloc
	.text
# -- Begin  computeContengency
# mark_begin;
       .align    16,0x90
	.globl computeContengency
computeContengency:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %ecx
# parameter 5: %r8
# parameter 6: %r9
..B2.1:                         # Preds ..B2.0
..___tag_value_computeContengency.8:                            #14.1
        pushq     %rbp                                          #14.1
..___tag_value_computeContengency.10:                           #
        movq      %rsp, %rbp                                    #14.1
..___tag_value_computeContengency.11:                           #
        subq      $64, %rsp                                     #14.1
        movq      %r13, -24(%rbp)                               #14.1
..___tag_value_computeContengency.13:                           #
        movq      %rdi, %r13                                    #14.1
        movq      %r12, -32(%rbp)                               #14.1
..___tag_value_computeContengency.14:                           #
        movq      %rsi, %r12                                    #14.1
        movq      %r13, %rax                                    #16.26
        orq       %r12, %rax                                    #16.26
        movq      %rbx, -40(%rbp)                               #14.1
..___tag_value_computeContengency.15:                           #
        movq      %r8, %rbx                                     #14.1
        movq      %r15, -48(%rbp)                               #14.1
..___tag_value_computeContengency.16:                           #
        movl      %ecx, %r15d                                   #14.1
        movq      %r14, -56(%rbp)                               #14.1
..___tag_value_computeContengency.17:                           #
        movq      %r9, %r14                                     #14.1
        movq      %rdx, -16(%rbp)                               #14.1
        testq     $15, %rax                                     #16.26
        jne       ..B2.23       # Prob 28%                      #16.26
                                # LOE rbx r12 r13 r14 r15d
..B2.2:                         # Preds ..B2.1
        movq      %r14, %rdx                                    #21.71
        lea       -64(%rbp), %rdi                               #21.71
        shlq      $4, %rdx                                      #21.71
        movl      $16, %esi                                     #21.71
        call      posix_memalign                                #21.71
                                # LOE rbx r12 r13 r14 eax r15d
..B2.3:                         # Preds ..B2.2
        testl     %eax, %eax                                    #21.71
        jne       ..B2.21       # Prob 50%                      #21.71
                                # LOE rbx r12 r13 r14 r15d
..B2.4:                         # Preds ..B2.3
        movq      -64(%rbp), %rdx                               #21.71
        testq     %rdx, %rdx                                    #23.17
        je        ..B2.22       # Prob 2%                       #23.17
                                # LOE rdx rbx r12 r13 r14 r15d
..B2.5:                         # Preds ..B2.4
        xorl      $8, %r15d                                     #26.38
        movq      %rbx, %r8                                     #24.46
        movdqa    .L_2il0floatpacket.12(%rip), %xmm0            #25.38
        lea       63(,%r14,8), %rax                             #30.26
        andq      $-16, %r8                                     #24.46
        movd      %r15d, %xmm1                                  #26.38
        punpcklbw %xmm1, %xmm1                                  #26.38
        punpcklbw %xmm1, %xmm1                                  #26.38
        pshufd    $0, %xmm1, %xmm3                              #26.38
        addq      $15, %rax                                     #30.26
        andq      $-16, %rax                                    #30.26
        subq      %rax, %rsp                                    #30.26
        movq      %rsp, %rax                                    #30.26
                                # LOE rax rdx rbx r8 r12 r13 r14 xmm0 xmm3
..B2.28:                        # Preds ..B2.5
        andq      $-64, %rax                                    #31.85
        xorl      %ecx, %ecx                                    #33.13
        pxor      %xmm1, %xmm1                                  #34.36
                                # LOE rax rdx rcx rbx r8 r12 r13 r14 xmm0 xmm1 xmm3
..B2.6:                         # Preds ..B2.6 ..B2.28
        movdqa    %xmm1, (%rax,%rcx,8)                          #36.32
        addq      $2, %rcx                                      #37.5
        cmpq      %r14, %rcx                                    #38.17
        jb        ..B2.6        # Prob 82%                      #38.17
                                # LOE rax rdx rcx rbx r8 r12 r13 r14 xmm0 xmm1 xmm3
..B2.7:                         # Preds ..B2.6
        movq      $0, -8(%rax,%r14,8)                           #39.42
        movl      $16, %edi                                     #41.25
        cmpq      $16, %r8                                      #42.18
        jbe       ..B2.13       # Prob 10%                      #42.18
                                # LOE rax rdx rbx rdi r8 r12 r13 r14 xmm0 xmm3
..B2.8:                         # Preds ..B2.7
        movq      -16(%rbp), %rcx                               #
                                # LOE rax rdx rcx rbx rdi r8 r12 r13 r14 xmm0 xmm3
..B2.9:                         # Preds ..B2.11 ..B2.8
        movdqa    -16(%rdi,%r13), %xmm2                         #43.36
        xorl      %esi, %esi                                    #47.22
        movdqa    -16(%rdi,%r12), %xmm1                         #44.36
        pxor      %xmm0, %xmm2                                  #43.36
        pcmpeqb   %xmm3, %xmm1                                  #44.36
        addq      $16, %rdi                                     #45.4
                                # LOE rax rdx rcx rbx rsi rdi r8 r12 r13 r14 xmm0 xmm1 xmm2 xmm3
..B2.10:                        # Preds ..B2.10 ..B2.9
        movzbl    (%rsi,%rcx), %r9d                             #49.53
        movdqa    %xmm2, %xmm7                                  #50.39
        movdqa    %xmm1, %xmm6                                  #52.49
        movdqa    %xmm1, %xmm8                                  #53.49
        movd      %r9d, %xmm4                                   #49.53
        punpcklbw %xmm4, %xmm4                                  #49.53
        punpcklbw %xmm4, %xmm4                                  #49.53
        pshufd    $0, %xmm4, %xmm5                              #49.53
        pxor      %xmm0, %xmm5                                  #49.39
        pcmpgtb   %xmm5, %xmm7                                  #50.39
        pand      %xmm7, %xmm6                                  #52.49
        pandn     %xmm7, %xmm8                                  #53.49
        pmovmskb  %xmm6, %r10d                                  #52.31
        pmovmskb  %xmm8, %r11d                                  #53.31
        popcnt    %r15d, %r10d                                  #55.0
        popcnt    %r9d, %r11d                                   #56.0
        addl      %r15d, (%rax,%rsi,8)                          #57.5
        addl      %r9d, 4(%rax,%rsi,8)                          #58.5
        incq      %rsi                                          #60.15
        cmpq      %r14, %rsi                                    #60.28
        jb        ..B2.10       # Prob 82%                      #60.28
                                # LOE rax rdx rcx rbx rsi rdi r8 r12 r13 r14 xmm0 xmm1 xmm2 xmm3
..B2.11:                        # Preds ..B2.10
        cmpq      %r8, %rdi                                     #42.18
        jb        ..B2.9        # Prob 82%                      #42.18
                                # LOE rax rdx rcx rbx rdi r8 r12 r13 r14 xmm0 xmm3
..B2.13:                        # Preds ..B2.11 ..B2.7
        subq      %r8, %rdi                                     #63.23
        je        ..B2.17       # Prob 10%                      #64.7
                                # LOE rax rdx rbx rdi r8 r12 r13 r14 xmm0 xmm3
..B2.14:                        # Preds ..B2.13
        movdqu    -16(%r8,%r13), %xmm2                          #65.36
        movl      %edi, %ecx                                    #67.57
        movl      $1, %esi                                      #67.57
        movdqa    -16(%r8,%r12), %xmm1                          #66.36
        pxor      %xmm0, %xmm2                                  #65.36
        shll      %cl, %esi                                     #67.57
        xorl      %ecx, %ecx                                    #69.22
        movq      -16(%rbp), %r12                               #69.22
        decl      %esi                                          #67.65
        pcmpeqb   %xmm3, %xmm1                                  #66.36
                                # LOE rax rdx rcx rbx r12 r14 esi xmm0 xmm1 xmm2
..B2.15:                        # Preds ..B2.15 ..B2.14
        movzbl    (%rcx,%r12), %edi                             #71.53
        movdqa    %xmm2, %xmm6                                  #72.39
        movdqa    %xmm1, %xmm5                                  #74.49
        movdqa    %xmm1, %xmm7                                  #75.49
        movd      %edi, %xmm3                                   #71.53
        punpcklbw %xmm3, %xmm3                                  #71.53
        punpcklbw %xmm3, %xmm3                                  #71.53
        pshufd    $0, %xmm3, %xmm4                              #71.53
        pxor      %xmm0, %xmm4                                  #71.39
        pcmpgtb   %xmm4, %xmm6                                  #72.39
        pand      %xmm6, %xmm5                                  #74.49
        pandn     %xmm6, %xmm7                                  #75.49
        pmovmskb  %xmm5, %r8d                                   #74.31
        pmovmskb  %xmm7, %r9d                                   #75.31
        andl      %esi, %r8d                                    #76.5
        andl      %esi, %r9d                                    #77.5
        popcnt    %r10d, %r8d                                   #79.0
        popcnt    %r11d, %r9d                                   #80.0
        addl      %r10d, (%rax,%rcx,8)                          #81.5
        addl      %r11d, 4(%rax,%rcx,8)                         #82.5
        incq      %rcx                                          #84.15
        cmpq      %r14, %rcx                                    #84.28
        jb        ..B2.15       # Prob 82%                      #84.28
                                # LOE rax rdx rcx rbx r12 r14 esi xmm0 xmm1 xmm2
..B2.17:                        # Preds ..B2.15 ..B2.13
        xorl      %esi, %esi                                    #87.15
        movq      %rdx, %rcx                                    #
        testq     %r14, %r14                                    #87.20
        jbe       ..B2.22       # Prob 10%                      #87.20
                                # LOE rax rdx rcx rbx rsi r14
..B2.19:                        # Preds ..B2.17 ..B2.19
        movl      4(%rax,%rsi,8), %r9d                          #90.20
        movq      %rbx, %rdi                                    #88.20
        movl      %r9d, 8(%rcx)                                 #90.4
        negq      %r9                                           #91.20
        movl      (%rax,%rsi,8), %r8d                           #88.28
        incq      %rsi                                          #87.32
        subq      %r8, %rdi                                     #88.20
        addq      %rbx, %r9                                     #91.20
        movl      %edi, (%rcx)                                  #88.4
        movl      %r8d, 4(%rcx)                                 #89.4
        movl      %r9d, 12(%rcx)                                #91.4
        addq      $16, %rcx                                     #87.32
        cmpq      %r14, %rsi                                    #87.20
        jb        ..B2.19       # Prob 82%                      #87.20
        jmp       ..B2.22       # Prob 100%                     #87.20
                                # LOE rax rdx rcx rbx rsi r14
..B2.21:                        # Preds ..B2.3
        xorl      %edx, %edx                                    #21.71
                                # LOE rdx
..B2.22:                        # Preds ..B2.19 ..B2.21 ..B2.17 ..B2.4
        movq      %rdx, %rax                                    #96.9
        movq      -40(%rbp), %rbx                               #96.9
..___tag_value_computeContengency.18:                           #
        movq      -32(%rbp), %r12                               #96.9
..___tag_value_computeContengency.19:                           #
        movq      -24(%rbp), %r13                               #96.9
..___tag_value_computeContengency.20:                           #
        movq      -56(%rbp), %r14                               #96.9
..___tag_value_computeContengency.21:                           #
        movq      -48(%rbp), %r15                               #96.9
..___tag_value_computeContengency.22:                           #
        movq      %rbp, %rsp                                    #96.9
        popq      %rbp                                          #96.9
..___tag_value_computeContengency.23:                           #
        ret                                                     #96.9
..___tag_value_computeContengency.24:                           #
                                # LOE
..B2.23:                        # Preds ..B2.1
        movl      $.L_2__STRING.0, %edi                         #17.3
        movq      stderr(%rip), %rsi                            #17.3
        call      fputs                                         #17.3
                                # LOE
..B2.24:                        # Preds ..B2.23
        xorl      %eax, %eax                                    #18.10
        movq      -40(%rbp), %rbx                               #18.10
..___tag_value_computeContengency.30:                           #
        movq      -32(%rbp), %r12                               #18.10
..___tag_value_computeContengency.31:                           #
        movq      -24(%rbp), %r13                               #18.10
..___tag_value_computeContengency.32:                           #
        movq      -56(%rbp), %r14                               #18.10
..___tag_value_computeContengency.33:                           #
        movq      -48(%rbp), %r15                               #18.10
..___tag_value_computeContengency.34:                           #
        movq      %rbp, %rsp                                    #18.10
        popq      %rbp                                          #18.10
..___tag_value_computeContengency.35:                           #
        ret                                                     #18.10
        .align    16,0x90
..___tag_value_computeContengency.36:                           #
                                # LOE
# mark_end;
	.type	computeContengency,@function
	.size	computeContengency,.-computeContengency
	.data
# -- End  computeContengency
	.text
# -- Begin  computeHistogram
# mark_begin;
       .align    16,0x90
	.globl computeHistogram
computeHistogram:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B3.1:                         # Preds ..B3.0
..___tag_value_computeHistogram.37:                             #101.1
        lea       15(%rsi), %rcx                                #102.31
        pxor      %xmm0, %xmm0                                  #102.31
        xorl      %eax, %eax                                    #102.31
        movq      %xmm0, (%rsi)                                 #102.31
        andq      $-16, %rcx                                    #102.31
        movl      %eax, 8(%rsi)                                 #102.31
        movw      %ax, 12(%rsi)                                 #102.31
        movb      %al, 14(%rsi)                                 #102.31
        movl      $1008, %eax                                   #102.31
        movq      %xmm0, 1008(%rsi)                             #102.31
        movq      %xmm0, 1016(%rsi)                             #102.31
                                # LOE rax rdx rcx rbx rbp rsi rdi r12 r13 r14 r15 xmm0
..B3.12:                        # Preds ..B3.12 ..B3.1
        movaps    %xmm0, -16(%rcx,%rax)                         #102.31
        movaps    %xmm0, -32(%rcx,%rax)                         #102.31
        movaps    %xmm0, -48(%rcx,%rax)                         #102.31
        movaps    %xmm0, -64(%rcx,%rax)                         #102.31
        movaps    %xmm0, -80(%rcx,%rax)                         #102.31
        movaps    %xmm0, -96(%rcx,%rax)                         #102.31
        movaps    %xmm0, -112(%rcx,%rax)                        #102.31
        subq      $112, %rax                                    #102.31
        jne       ..B3.12       # Prob 88%                      #102.31
                                # LOE rax rdx rcx rbx rbp rsi rdi r12 r13 r14 r15 xmm0
..B3.2:                         # Preds ..B3.12
        testq     %rdx, %rdx                                    #104.21
        jbe       ..B3.9        # Prob 50%                      #104.21
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15
..B3.3:                         # Preds ..B3.2
        movq      %rdx, %rax                                    #104.2
        movl      $1, %ecx                                      #104.2
        shrq      $1, %rax                                      #104.2
        xorl      %r9d, %r9d                                    #104.2
        testq     %rax, %rax                                    #104.2
        jbe       ..B3.7        # Prob 10%                      #104.2
                                # LOE rax rdx rcx rbx rbp rsi rdi r9 r12 r13 r14 r15
..B3.4:                         # Preds ..B3.3
        movzbl    (%rdx,%rdi), %r8d                             #105.13
        movl      (%rsi,%r8,4), %ecx                            #105.3
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r12 r13 r14 r15 ecx
..B3.5:                         # Preds ..B3.5 ..B3.4
        incq      %r9                                           #104.2
        addl      $2, %ecx                                      #105.3
        cmpq      %rax, %r9                                     #104.2
        jb        ..B3.5        # Prob 63%                      #104.2
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r12 r13 r14 r15 ecx
..B3.6:                         # Preds ..B3.5
        movl      %ecx, (%rsi,%r8,4)                            #105.3
        lea       1(,%r9,2), %rcx                               #104.2
                                # LOE rdx rcx rbx rbp rsi rdi r12 r13 r14 r15
..B3.7:                         # Preds ..B3.6 ..B3.3
        decq      %rcx                                          #104.2
        cmpq      %rcx, %rdx                                    #104.2
        jbe       ..B3.9        # Prob 10%                      #104.2
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15
..B3.8:                         # Preds ..B3.7
        movzbl    (%rdx,%rdi), %eax                             #105.13
        incl      (%rsi,%rax,4)                                 #105.3
                                # LOE rbx rbp r12 r13 r14 r15
..B3.9:                         # Preds ..B3.2 ..B3.7 ..B3.8
        ret                                                     #107.1
        .align    16,0x90
..___tag_value_computeHistogram.39:                             #
                                # LOE
# mark_end;
	.type	computeHistogram,@function
	.size	computeHistogram,.-computeHistogram
	.data
# -- End  computeHistogram
	.section .rodata, "a"
	.align 16
	.align 16
.L_2il0floatpacket.12:
	.long	0x08080808,0x08080808,0x08080808,0x08080808
	.type	.L_2il0floatpacket.12,@object
	.size	.L_2il0floatpacket.12,16
	.section .rodata.str1.32, "aMS",@progbits,1
	.align 32
	.align 32
.L_2__STRING.0:
	.long	1667590211
	.long	1752440939
	.long	1818304613
	.long	1701734249
	.long	1953391981
	.long	543584032
	.long	1634890337
	.long	1918967929
	.long	1701672295
	.long	544437358
	.long	1646292852
	.long	909189221
	.long	1954112032
	.long	1818304613
	.long	1701734249
	.long	663908
	.type	.L_2__STRING.0,@object
	.size	.L_2__STRING.0,64
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
	.4byte 0x00000014
	.8byte 0x7801000100000000
	.8byte 0x0000019008070c10
	.4byte 0x00000000
	.4byte 0x00000034
	.4byte 0x0000001c
	.8byte ..___tag_value__mm_malloc.1
	.8byte ..___tag_value__mm_malloc.7-..___tag_value__mm_malloc.1
	.byte 0x04
	.4byte ..___tag_value__mm_malloc.3-..___tag_value__mm_malloc.1
	.2byte 0x100e
	.byte 0x04
	.4byte ..___tag_value__mm_malloc.4-..___tag_value__mm_malloc.3
	.2byte 0x080e
	.byte 0x04
	.4byte ..___tag_value__mm_malloc.5-..___tag_value__mm_malloc.4
	.2byte 0x100e
	.byte 0x04
	.4byte ..___tag_value__mm_malloc.6-..___tag_value__mm_malloc.5
	.4byte 0x0000080e
	.2byte 0x0000
	.4byte 0x000000a4
	.4byte 0x00000054
	.8byte ..___tag_value_computeContengency.8
	.8byte ..___tag_value_computeContengency.36-..___tag_value_computeContengency.8
	.byte 0x04
	.4byte ..___tag_value_computeContengency.10-..___tag_value_computeContengency.8
	.2byte 0x100e
	.byte 0x04
	.4byte ..___tag_value_computeContengency.11-..___tag_value_computeContengency.10
	.4byte 0x8610060c
	.2byte 0x0402
	.4byte ..___tag_value_computeContengency.13-..___tag_value_computeContengency.11
	.2byte 0x058d
	.byte 0x04
	.4byte ..___tag_value_computeContengency.14-..___tag_value_computeContengency.13
	.2byte 0x068c
	.byte 0x04
	.4byte ..___tag_value_computeContengency.15-..___tag_value_computeContengency.14
	.2byte 0x0783
	.byte 0x04
	.4byte ..___tag_value_computeContengency.16-..___tag_value_computeContengency.15
	.2byte 0x088f
	.byte 0x04
	.4byte ..___tag_value_computeContengency.17-..___tag_value_computeContengency.16
	.2byte 0x098e
	.byte 0x04
	.4byte ..___tag_value_computeContengency.18-..___tag_value_computeContengency.17
	.2byte 0x04c3
	.4byte ..___tag_value_computeContengency.19-..___tag_value_computeContengency.18
	.2byte 0x04cc
	.4byte ..___tag_value_computeContengency.20-..___tag_value_computeContengency.19
	.2byte 0x04cd
	.4byte ..___tag_value_computeContengency.21-..___tag_value_computeContengency.20
	.2byte 0x04ce
	.4byte ..___tag_value_computeContengency.22-..___tag_value_computeContengency.21
	.2byte 0x04cf
	.4byte ..___tag_value_computeContengency.23-..___tag_value_computeContengency.22
	.2byte 0x04c6
	.4byte ..___tag_value_computeContengency.24-..___tag_value_computeContengency.23
	.8byte 0x058d068c02860783
	.4byte 0x088f098e
	.byte 0x04
	.4byte ..___tag_value_computeContengency.30-..___tag_value_computeContengency.24
	.2byte 0x04c3
	.4byte ..___tag_value_computeContengency.31-..___tag_value_computeContengency.30
	.2byte 0x04cc
	.4byte ..___tag_value_computeContengency.32-..___tag_value_computeContengency.31
	.2byte 0x04cd
	.4byte ..___tag_value_computeContengency.33-..___tag_value_computeContengency.32
	.2byte 0x04ce
	.4byte ..___tag_value_computeContengency.34-..___tag_value_computeContengency.33
	.2byte 0x04cf
	.4byte ..___tag_value_computeContengency.35-..___tag_value_computeContengency.34
	.4byte 0x000000c6
	.4byte 0x00000014
	.4byte 0x000000fc
	.8byte ..___tag_value_computeHistogram.37
	.8byte ..___tag_value_computeHistogram.39-..___tag_value_computeHistogram.37
# End
