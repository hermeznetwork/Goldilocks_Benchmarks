        global  gl_mul
        global  gl_mmul
        global  gl_fromm
        global  gl_tom
        global  gl_add
        global  gl_sub
        section .text
        DEFAULT REL
gl_mul:
        mov     rax, rdi                ; result (rax) initially holds x
        mul     rsi                     ; is x less than y?
        div     qword [q]                   ; if so, set result to y
        mov     rax, rdx                ; if so, set result to z
        ret   

gl_mmul:
        xor     r10, r10
        mov     rax, rdi                
        mul     rsi     
        mov     r8, rdx
        mov     r9, rax                
        mul     qword [mm]       
        mul     qword [q]          
        add     rax, r9
        adc     rdx, r8                     ; rmul
        cmovc   r10, qword [cq]     
        add     rdx, r10
        mov     rax, rdx
        ret

gl_fromm:
        xor     r10, r10
        mov     rax, rdi
        mov     r9, rax                
        mul     qword [mm]       
        mul     qword [q]          
        add     rax, r9
        adc     rdx, r10                    ; rmul
        cmovc   r10, qword [cq]     
        add     rdx, r10
        mov     rax, rdx
        ret

gl_tom:
        xor     r10, r10
        mov     rax, rdi                
        mul     qword [r2]     
        mov     r8, rdx
        mov     r9, rax                
        mul     qword [mm]       
        mul     qword [q]          
        add     rax, r9
        adc     rdx, r8                     ; adc
        cmovc   r10, qword [cq]     
        add     rdx, r10
        mov     rax, rdx
        ret

gl_add:
        xor     r10, r10
        mov     rax, rdi
        add     rax, rsi
        cmovc   r10, qword [cq]     
        add     rax, r10
        ret

gl_sub:
        xor     r10, r10
        mov     rax, rdi
        sub     rax, rsi
        cmovc   r10, qword [q]     
        add     rax, r10
        ret


        section .data
q       dq      0xFFFFFFFF00000001
mm      dq      0xFFFFFFFeFFFFFFFF
cq      dq      0x00000000FFFFFFFF
r2      dq      0xFFFFFFFe00000001