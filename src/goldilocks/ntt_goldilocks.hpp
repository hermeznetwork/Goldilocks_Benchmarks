#ifndef GOLDILOCKS
#define GOLDILOCKS

#include <inttypes.h>
#include <stdint.h> // for uint64_t
#include <omp.h>
#include <iostream>
#include <cassert>
#include <gmp.h>

#define GOLDILOCKS_PRIME 0xFFFFFFFF00000001ULL
#define uint128_t __uint128_t

#define CACHESIZE 1 << 18

extern "C" uint64_t gl_fromm(uint64_t a);
extern "C" uint64_t gl_tom(uint64_t a);
// extern "C" uint64_t gl_mmul2(uint64_t a, uint64_t b);
// extern "C" uint64_t gl_sub(uint64_t a, uint64_t b);
//  extern "C" uint64_t gl_add(uint64_t a, uint64_t b);

typedef unsigned __int128 uint128_t;

class Goldilocks
{

private:
    const static uint64_t Q;
    const static uint64_t MM;
    const static uint64_t CQ;
    const static uint64_t R2;

    u_int32_t s;
    uint64_t nqr;
    uint64_t *roots;
    uint64_t *powTwoInv;
    u_int32_t nThreads;

public:
    Goldilocks(u_int64_t maxDomainSize, u_int32_t _nThreads = 0)
    {
        nThreads = _nThreads == 0 ? omp_get_max_threads() : _nThreads;

        u_int32_t domainPow = log2(maxDomainSize);

        mpz_t m_qm1d2;
        mpz_t m_q;
        mpz_t m_nqr;
        mpz_t m_aux;
        mpz_init(m_qm1d2);
        mpz_init(m_q);
        mpz_init(m_nqr);
        mpz_init(m_aux);

        u_int64_t negone = GOLDILOCKS_PRIME - 1;

        mpz_import(m_aux, 1, 1, sizeof(u_int64_t), 0, 0, &negone);
        mpz_add_ui(m_q, m_aux, 1);
        mpz_fdiv_q_2exp(m_qm1d2, m_aux, 1);

        mpz_set_ui(m_nqr, 2);
        mpz_powm(m_aux, m_nqr, m_qm1d2, m_q);
        while (mpz_cmp_ui(m_aux, 1) == 0)
        {
            mpz_add_ui(m_nqr, m_nqr, 1);
            mpz_powm(m_aux, m_nqr, m_qm1d2, m_q);
        }

        s = 1;
        mpz_set(m_aux, m_qm1d2);
        while ((!mpz_tstbit(m_aux, 0)) && (s < domainPow))
        {
            mpz_fdiv_q_2exp(m_aux, m_aux, 1);
            s++;
        }

        nqr = mpz_get_ui(m_nqr);

        if (s < domainPow)
        {
            throw std::range_error("Domain size too big for the curve");
        }

        uint64_t nRoots = 1LL << s;

        roots = new uint64_t[nRoots];
        powTwoInv = new uint64_t[s + 1];

        roots[0] = gl_tom(1);
        powTwoInv[0] = gl_tom(1);

        if (nRoots > 1)
        {
            mpz_powm(m_aux, m_nqr, m_aux, m_q);
            roots[1] = gl_tom(mpz_get_ui(m_aux));

            mpz_set_ui(m_aux, 2);
            mpz_invert(m_aux, m_aux, m_q);
            powTwoInv[1] = gl_tom(mpz_get_ui(m_aux));
        }

        for (uint64_t i = 2; i < nRoots; i++)
        {
            // roots[i] = gl_mmul(roots[i - 1], roots[1]);

            roots[i] = gl_mmul2(roots[i - 1], roots[1]);
            // assert(roots[i] == aux);
        }

        uint64_t aux = gl_mmul2(roots[nRoots - 1], roots[1]);

        assert(gl_fromm(aux) == 1);

        for (uint64_t i = 2; i <= s; i++)
        {
            powTwoInv[i] = gl_mmul2(powTwoInv[i - 1], powTwoInv[1]);
        }

        mpz_clear(m_qm1d2);
        mpz_clear(m_q);
        mpz_clear(m_nqr);
        mpz_clear(m_aux);
    };
    void ntt(u_int64_t *a, u_int64_t n);
    void reversePermutation(u_int64_t *dst, u_int64_t *a, u_int64_t n);
    u_int32_t log2(u_int64_t n);

    void printVector(u_int64_t *a, u_int64_t n);

    inline static uint64_t gl_add(uint64_t in1, uint64_t in2)
    {
        uint64_t res;
        __asm__("mov   %1, %0\n\t"
                "add   %2, %0\n\t"
                "jnc  1f\n\t"
                "add   %3, %0\n\t"
                "1: \n\t"
                : "=&a"(res)
                : "r"(in1), "r"(in2), "m"(CQ)
                :);
        return res;
    };
    inline static uint64_t gl_sub(uint64_t in1, uint64_t in2)
    {
        uint64_t res;
        __asm__("mov   %1, %0\n\t"
                "sub   %2, %0\n\t"
                "jnc  1f\n\t"
                "add   %3, %0\n\t"
                "1: \n\t"
                : "=&a"(res)
                : "r"(in1), "r"(in2), "m"(Q)
                :);
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
        // printf("r: %lu,q: %lu\n", r, q);
        return q;
    };

    inline static uint64_t gl_mmul2(uint64_t in1, uint64_t in2)
    {
        uint64_t res;
        __asm__("mov   %1, %%rax\n\t"
                "mul   %2\n\t"
                "mov   %%rdx, %%r8\n\t"
                "mov   %%rax, %%r9\n\t"
                "mulq   %3\n\t"
                "mulq   %4\n\t"
                "add    %%r9, %%rax\n\t"
                "adc    %%r8, %%rdx\n\t"
                "jnc  1f\n\t"
                "add   %5, %%rdx\n\t"
                "1:"
                : "=&d"(res)
                : "r"(in1), "r"(in2), "m"(MM), "m"(Q), "m"(CQ)
                : "%rax", "%r8", "%r9");
        return res;
    }

    inline uint64_t &root(u_int32_t domainPow, u_int64_t idx)
    {
        return roots[idx << (s - domainPow)];
    }
    void shuffle(u_int64_t *dst, u_int64_t *src, uint64_t n, uint64_t s);
    void intt(u_int64_t *a, u_int64_t n);

    void traspose(
        u_int64_t *dst,
        u_int64_t *src,
        uint64_t srcRowSize,
        uint64_t srcX,
        uint64_t srcWidth,
        uint64_t srcY,
        uint64_t srcHeight,
        uint64_t dstRowSize,
        uint64_t dstX,
        uint64_t dstY);
};

#endif // GOLDILOCKS