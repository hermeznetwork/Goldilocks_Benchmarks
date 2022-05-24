#ifndef POSEIDON_GOLDILOCKS_OPT
#define POSEIDON_GOLDILOCKS_OPT

#include <inttypes.h>
#include <stdint.h> // for uint64_t
#include "poseidon_goldilocks_opt_constants.hpp"

#define SPONGE_WIDTH 12
#define MAX_WIDTH 12
#define HALF_N_FULL_ROUNDS 4
#define N_FULL_ROUNDS_TOTAL (2 * HALF_N_FULL_ROUNDS)
#define N_PARTIAL_ROUNDS 22
#define N_ROUNDS (N_FULL_ROUNDS_TOTAL + N_PARTIAL_ROUNDS)

#define GOLDILOCKS_PRIME 0xFFFFFFFF00000001ULL

#define ASM 1

#define uint128_t __uint128_t

typedef unsigned __int128 uint128_t;

extern "C" uint64_t gl_fromm(uint64_t a);
extern "C" uint64_t gl_tom(uint64_t a);

class Poseidon_goldilocks_opt
{
private:
    const static uint64_t Q;
    const static uint64_t MM;
    const static uint64_t CQ;
    const static uint64_t R2;
    uint64_t sum_counter = 0;
    uint64_t mul_counter = 0;

    inline void static pow7(uint64_t &x)
    {
#if ASM == 1
        uint64_t x2 = gl_mmul(x, x);
        uint64_t x3 = gl_mmul(x, x2);
        uint64_t x4 = gl_mmul(x2, x2);

        x = gl_mmul(x3, x4);
#else
        uint128_t x2 = ((uint128_t)x * (uint128_t)x) % GOLDILOCKS_PRIME;
        uint128_t x4 = (x2 * x2) % GOLDILOCKS_PRIME;
        uint128_t x3 = ((uint128_t)x * x2) % GOLDILOCKS_PRIME;
        x = (x3 * x4) % GOLDILOCKS_PRIME;
#endif
    };

    inline static uint64_t gl_add(uint64_t in1, uint64_t in2)
    {
        /*
                xor     r10, r10
                mov     rax, rdi
                add     rax, rsi
                cmovc   r10, qword [cq]
                add     rax, r10
                ret
        */
        uint64_t res;
        __asm__("xor   %%r10, %%r10\n\t"
                "mov   %1, %0\n\t"
                "add   %2, %0\n\t"
                "cmovc %3, %%r10\n\t"
                "add   %%r10, %0\n\t"
                : "=&a"(res)
                : "r"(in1), "r"(in2), "m"(CQ)
                : "%r10");
        return res;
    };

    inline static uint64_t gl_mul(uint64_t a, uint64_t b)
    {

        uint64_t r;
        uint64_t q;
        uint64_t m = GOLDILOCKS_PRIME;
        __asm__(
            "mulq   %3\n\t"
            "divq   %4\n\t"
            : "=a"(r), "=&d"(q)
            : "a"(a), "rm"(b), "rm"(m));
        return q;
    };

    inline static uint64_t gl_mmul(uint64_t in1, uint64_t in2)
    {
        uint64_t res;
        __asm__("xor   %%r10, %%r10\n\t"
                "mov   %1, %%rax\n\t"
                "mul   %2\n\t"
                "mov   %%rdx, %%r8\n\t"
                "mov   %%rax, %%r9\n\t"
                "mulq   %3\n\t"
                "mulq   %4\n\t"
                "add    %%r9, %%rax\n\t"
                "adc    %%r8, %%rdx\n\t"
                "cmovc %5, %%r10\n\t"
                "add   %%r10, %%rdx\n\t"
                : "=&d"(res)
                : "r"(in1), "r"(in2), "m"(MM), "m"(Q), "m"(CQ)
                : "%rax", "%r8", "%r9", "%r10");
        return res;
    }

public:
    void static hash(uint64_t (&input)[SPONGE_WIDTH]);
    void static print();
};

#endif // POSEIDON_GOLDILOCKS_OPT