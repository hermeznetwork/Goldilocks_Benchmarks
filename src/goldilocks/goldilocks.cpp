#include "goldilocks.hpp"
#include <cstring> //memset
#include "poseidon_goldilocks_constants.hpp"
#include "poseidon_goldilocks_opt_constants.hpp"

const uint64_t Goldilocks::Q = 0xFFFFFFFF00000001LL;
const uint64_t Goldilocks::MM = 0xFFFFFFFeFFFFFFFFLL;
const uint64_t Goldilocks::CQ = 0x00000000FFFFFFFFLL;
const uint64_t Goldilocks::R2 = 0xFFFFFFFe00000001LL;

static inline u_int64_t BR(u_int64_t x, u_int64_t domainPow)
{
    x = (x >> 16) | (x << 16);
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    return (((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1)) >> (32 - domainPow);
}

void Goldilocks::ntt(u_int64_t *_a, u_int64_t n)
{
    u_int64_t *aux_a = new u_int64_t[n];
    u_int64_t *a = _a;
    u_int64_t *a2 = aux_a;
    u_int64_t *tmp;
    reversePermutation(a2, a, n);
    tmp = a2;
    a2 = a;
    a = tmp;

    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    u_int64_t maxBatchPow = s / 4;

    u_int64_t batchSize = 1 << maxBatchPow;
    u_int64_t nBatches = n / batchSize;

    for (u_int64_t s = 1; s <= domainPow; s += maxBatchPow)
    {
        u_int64_t sInc = s + maxBatchPow <= domainPow ? maxBatchPow : domainPow - s + 1;
        omp_set_dynamic(0);
        omp_set_num_threads(nThreads);
#pragma omp parallel for
        for (u_int64_t b = 0; b < nBatches; b++)
        {
            u_int64_t rs = s - 1;
            uint64_t re = domainPow - 1;
            uint64_t rb = 1 << rs;
            uint64_t rm = (1 << (re - rs)) - 1;

            for (u_int64_t si = 0; si < sInc; si++)
            {

                u_int64_t m = 1 << (s + si);
                u_int64_t mdiv2 = m >> 1;
                u_int64_t mdiv2i = 1 << si;
                u_int64_t mi = mdiv2i * 2;
                for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                {
                    u_int64_t t;
                    u_int64_t u;
                    u_int64_t ki = b * batchSize + (i / mdiv2i) * mi;
                    u_int64_t ji = i % mdiv2i;

                    u_int64_t j = (b * batchSize / 2 + i);
                    j = (j & rm) * rb + (j >> (re - rs));
                    j = j % mdiv2;

                    t = gl_mmul(root(s + si, j), a[ki + ji + mdiv2i]);
                    u = a[ki + ji];
                    a[ki + ji] = gl_add(t, u);
                    a[ki + ji + mdiv2i] = gl_sub(u, t);
                }
                /* DOEN'T WORK, it tries to compute the roots dinamicaly
                 u_int64_t m = 1 << (s + si);
                 u_int64_t mdiv2 = m >> 1;
                 u_int64_t mdiv2i = 1 << si;
                 u_int64_t mi = mdiv2i * 2;

                 u_int64_t j_tmp = (b * batchSize / 2);
                 j_tmp = (j_tmp & rm) * rb + (j_tmp >> (re - rs));
                 j_tmp = j_tmp % mdiv2;
                 u_int64_t w = root(s + si, j_tmp);
                 u_int64_t wi = root(s + si, 1);

                 u_int64_t mask = mdiv2i - 1;
                 u_int64_t k1 = b * batchSize;
                 u_int64_t k2 = mi / mdiv2i;

                 for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                 {
                     u_int64_t t;
                     u_int64_t u;
                     u_int64_t ki = k1 + k2 * i;
                     u_int64_t ji = i & mask;

                     t = gl_mmul(w, a[ki + ji + mdiv2i]);
                     u = a[ki + ji];
                     a[ki + ji] = gl_add(t, u);
                     a[ki + ji + mdiv2i] = gl_sub(u, t);
                     w = gl_mmul(w, wi);
                 }
                 */
            }
        }

        shuffle(a2, a, n, sInc);
        tmp = a2;
        a2 = a;
        a = tmp;
    }
    if (a != _a)
    {
        // printf("baaaaad!\n");
        shuffle(a, a2, n, 0);
    }
    delete aux_a;
}


void Goldilocks::intt(u_int64_t *a, u_int64_t n)
{
    ntt(a, n);
    u_int64_t domainPow = log2(n);
    u_int64_t nDiv2 = n >> 1;
#pragma omp parallel for
    for (u_int64_t i = 1; i < nDiv2; i++)
    {
        u_int64_t tmp;
        u_int64_t r = n - i;
        tmp = a[i];
        a[i] = gl_mmul2(a[r], powTwoInv[domainPow]);
        a[r] = gl_mmul2(tmp, powTwoInv[domainPow]);
    }
    a[0] = gl_mmul2(a[0], powTwoInv[domainPow]);
    a[n >> 1] = gl_mmul2(a[n >> 1], powTwoInv[domainPow]);
}

void Goldilocks::shuffle(u_int64_t *dst, u_int64_t *src, uint64_t n, uint64_t s)
{
    uint64_t srcRowSize = 1 << s;

    uint64_t srcX = 0;
    uint64_t srcWidth = 1 << s;
    uint64_t srcY = 0;
    uint64_t srcHeight = n / srcRowSize;

    uint64_t dstRowSize = n / srcRowSize;
    uint64_t dstX = 0;
    uint64_t dstY = 0;

#pragma omp parallel
#pragma omp single
    traspose(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight, dstRowSize, dstX, dstY);
#pragma omp taskwait
}

void Goldilocks::traspose(
    u_int64_t *dst,
    u_int64_t *src,
    uint64_t srcRowSize,
    uint64_t srcX,
    uint64_t srcWidth,
    uint64_t srcY,
    uint64_t srcHeight,
    uint64_t dstRowSize,
    uint64_t dstX,
    uint64_t dstY)
{
    if ((srcWidth == 1) || (srcHeight == 1) || (srcWidth * srcHeight < CACHESIZE))
    {
#pragma omp task
        {
            for (uint64_t x = 0; x < srcWidth; x++)
            {
                for (uint64_t y = 0; y < srcHeight; y++)
                {
                    dst[(dstY + +x) * dstRowSize + (dstX + y)] = src[(srcY + +y) * srcRowSize + (srcX + x)];
                }
            }
        }
        return;
    }
    if (srcWidth > srcHeight)
    {
        traspose(dst, src, srcRowSize, srcX, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY);
        traspose(dst, src, srcRowSize, srcX + srcWidth / 2, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY + srcWidth / 2);
    }
    else
    {
        traspose(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight / 2, dstRowSize, dstX, dstY);
        traspose(dst, src, srcRowSize, srcX, srcWidth, srcY + srcHeight / 2, srcHeight / 2, dstRowSize, dstX + srcHeight / 2, dstY);
    }
}

void Goldilocks::reversePermutation(u_int64_t *dst, u_int64_t *a, u_int64_t n)
{
    uint32_t domainSize = log2(n);
#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        r = BR(i, domainSize);
        dst[i] = a[r];
    }
}

u_int32_t Goldilocks::log2(u_int64_t n)
{
    assert(n != 0);
    u_int32_t res = 0;
    while (n != 1)
    {
        n >>= 1;
        res++;
    }
    return res;
}

void Goldilocks::poseidon(uint64_t (&state)[SPONGE_WIDTH])
{
#if ASM == 1
    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = Goldilocks::gl_tom(state[i]);
    }
#endif

    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
#if ASM == 1
        state[i] = Goldilocks::gl_add(state[i], Poseidon_goldilocks_opt_constants::C[i]);
#else
        // state[i] = (state[i] + Poseidon_goldilocks_opt_constants::C[i]) % GOLDILOCKS_PRIME;
        state[i] = add_gl(state[i], Poseidon_goldilocks_opt_constants::C[i]);
#endif
    }

    for (int r = 0; r < HALF_N_FULL_ROUNDS - 1; r++)
    {

        for (int j = 0; j < SPONGE_WIDTH; j++)
        {
            pow7(state[j]);
#if ASM == 1
            state[j] = Goldilocks::gl_add(state[j], Poseidon_goldilocks_opt_constants::C[(r + 1) * SPONGE_WIDTH + j]);
#else
            // state[j] = ((uint128_t)state[j] + (uint128_t)Poseidon_goldilocks_opt_constants::C[(r + 1) * SPONGE_WIDTH + j]) % GOLDILOCKS_PRIME;
            state[j] = add_gl(state[j], Poseidon_goldilocks_opt_constants::C[(r + 1) * SPONGE_WIDTH + j]);
#endif
        }

        uint64_t old_state[SPONGE_WIDTH];
        std::memcpy(&old_state, &state, sizeof(uint64_t) * SPONGE_WIDTH);

        for (int i = 0; i < SPONGE_WIDTH; i++)
        {
            state[i] = 0;
            for (int j = 0; j < SPONGE_WIDTH; j++)
            {
#if ASM == 1
                uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
                mji = Goldilocks::gl_mmul(mji, old_state[j]);
                state[i] = Goldilocks::gl_add(state[i], mji);
#else
                uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
                mji = ((uint128_t)mji * (uint128_t)old_state[j]) % GOLDILOCKS_PRIME;
                state[i] = add_gl(state[i], mji);
                // state[i] = ((uint128_t)state[i] + (uint128_t)mji) % GOLDILOCKS_PRIME;
#endif
            }
        }
    }

    for (int j = 0; j < SPONGE_WIDTH; j++)
    {
        pow7(state[j]);
#if ASM == 1
        state[j] = Goldilocks::gl_add(state[j], Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS * SPONGE_WIDTH)]);
#else
        // state[j] = ((uint128_t)state[j] + (uint128_t)Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS * SPONGE_WIDTH)]) % GOLDILOCKS_PRIME;
        state[j] = add_gl(state[j], Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS * SPONGE_WIDTH)]);

#endif
    }

    uint64_t old_state[SPONGE_WIDTH];
    std::memcpy(&old_state, &state, sizeof(uint64_t) * SPONGE_WIDTH);
    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = 0;
        for (int j = 0; j < SPONGE_WIDTH; j++)
        {
#if ASM == 1
            uint64_t mji = Poseidon_goldilocks_opt_constants::P[j][i];
            mji = Goldilocks::gl_mmul(mji, old_state[j]);
            state[i] = Goldilocks::gl_add(state[i], mji);
#else
            uint64_t mji = Poseidon_goldilocks_opt_constants::P[j][i];
            mji = ((uint128_t)mji * (uint128_t)old_state[j]) % GOLDILOCKS_PRIME;
            state[i] = add_gl(state[i], mji);
#endif
        }
    }

    for (int r = 0; r < N_PARTIAL_ROUNDS; r++)
    {
        pow7(state[0]);
#if ASM == 1
        state[0] = Goldilocks::gl_add(state[0], Poseidon_goldilocks_opt_constants::C[(HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + r]);
#else
        state[0] = add_gl(state[0], Poseidon_goldilocks_opt_constants::C[(HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + r]);
        // state[0] = ((uint128_t)state[0] + (uint128_t)Poseidon_goldilocks_opt_constants::C[(HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + r]) % GOLDILOCKS_PRIME;
#endif
        uint64_t s0 = 0;

        uint64_t accumulator1 = Poseidon_goldilocks_opt_constants::S[(SPONGE_WIDTH * 2 - 1) * r];
        accumulator1 = Goldilocks::gl_mmul(accumulator1, state[0]);
        s0 = Goldilocks::gl_add(s0, accumulator1);

        for (int j = 1; j < SPONGE_WIDTH; j++)
        {
#if ASM == 1
            uint64_t accumulator1 = Poseidon_goldilocks_opt_constants::S[(SPONGE_WIDTH * 2 - 1) * r + j];
            accumulator1 = Goldilocks::gl_mmul(accumulator1, state[j]);
            s0 = Goldilocks::gl_add(s0, accumulator1);

            uint64_t accumulator2 = Poseidon_goldilocks_opt_constants::S[(SPONGE_WIDTH * 2 - 1) * r + SPONGE_WIDTH + j - 1];
            accumulator2 = Goldilocks::gl_mmul(accumulator2, state[0]);
            state[j] = Goldilocks::gl_add(state[j], accumulator2);

#else
            uint64_t accumulator1 = Poseidon_goldilocks_opt_constants::S[(SPONGE_WIDTH * 2 - 1) * r + j];
            accumulator1 = ((uint128_t)accumulator1 * (uint128_t)state[j]) % GOLDILOCKS_PRIME;
            // s0 = ((uint128_t)s0 + (uint128_t)accumulator1) % GOLDILOCKS_PRIME;
            s0 = add_gl(s0, accumulator1);

            if (j > 0)
            {
                uint64_t accumulator2 = Poseidon_goldilocks_opt_constants::S[(SPONGE_WIDTH * 2 - 1) * r + SPONGE_WIDTH + j - 1];
                accumulator2 = ((uint128_t)accumulator2 * (uint128_t)state[0]) % GOLDILOCKS_PRIME;
                // state[j] = ((uint128_t)state[j] + (uint128_t)accumulator2) % GOLDILOCKS_PRIME;
                state[j] = add_gl(state[j], accumulator2);
            }
#endif
        }
        state[0] = s0;
    }
    for (int r = 0; r < HALF_N_FULL_ROUNDS - 1; r++)
    {
        for (int j = 0; j < SPONGE_WIDTH; j++)
        {
            pow7(state[j]);

#if ASM == 1
            state[j] = Goldilocks::gl_add(state[j], Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + N_PARTIAL_ROUNDS + r * SPONGE_WIDTH]);
#else
            // state[j] = ((uint128_t)state[j] + (uint128_t)Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + N_PARTIAL_ROUNDS + r * SPONGE_WIDTH]) % GOLDILOCKS_PRIME;
            state[j] = add_gl(state[j], Poseidon_goldilocks_opt_constants::C[j + (HALF_N_FULL_ROUNDS + 1) * SPONGE_WIDTH + N_PARTIAL_ROUNDS + r * SPONGE_WIDTH]);

#endif
        }

        uint64_t old_state[SPONGE_WIDTH];
        std::memcpy(&old_state, &state, sizeof(uint64_t) * SPONGE_WIDTH);

        for (int i = 0; i < SPONGE_WIDTH; i++)
        {
            state[i] = 0;
            for (int j = 0; j < SPONGE_WIDTH; j++)
            {

#if ASM == 1
                uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
                mji = Goldilocks::gl_mmul(mji, old_state[j]);
                state[i] = Goldilocks::gl_add(state[i], mji);
#else
                uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
                mji = ((uint128_t)mji * (uint128_t)old_state[j]) % GOLDILOCKS_PRIME;
                // state[i] = ((uint128_t)state[i] + (uint128_t)mji) % GOLDILOCKS_PRIME;
                state[i] = add_gl(state[i], mji);
#endif
            }
        }
    }

    for (int j = 0; j < SPONGE_WIDTH; j++)
    {
        pow7(state[j]);
    }

    std::memcpy(&old_state, &state, sizeof(uint64_t) * SPONGE_WIDTH);

    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = 0;
        for (int j = 0; j < SPONGE_WIDTH; j++)
        {
#if ASM == 1
            uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
            mji = Goldilocks::gl_mmul(mji, old_state[j]);
            state[i] = Goldilocks::gl_add(state[i], mji);
#else
            uint64_t mji = Poseidon_goldilocks_opt_constants::M[j][i];
            mji = ((uint128_t)mji * (uint128_t)old_state[j]) % GOLDILOCKS_PRIME;
            // state[i] = ((uint128_t)state[i] + (uint128_t)mji) % GOLDILOCKS_PRIME;
            state[i] = add_gl(state[i], mji);

#endif
        }
    }
#if ASM == 1
    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = Goldilocks::gl_fromm(state[i]);
    }
#endif
}

void Goldilocks::linear_hash(uint64_t *output, uint64_t *input, uint64_t size)
{
    uint64_t remaining = size;
    uint64_t state[SPONGE_WIDTH] = {0};

    while (remaining)
    {
        if (remaining == size)
        {
            memset(state + RATE, 0, CAPACITY * sizeof(uint64_t));
        }
        else
        {
            std::memcpy(state + RATE, state, CAPACITY * sizeof(uint64_t));
        }

        uint64_t n = (remaining < RATE) ? remaining : RATE;

        std::memcpy(state, input + (size - remaining), n * sizeof(uint64_t));

        for (int i = n; i < RATE; i++)
        {
            state[i] = 0;
        }
        std::memcpy(state, input + (size - remaining), n * sizeof(uint64_t));

        poseidon(state);

        remaining -= n;
    }
    std::memcpy(output, state, 4 * sizeof(uint64_t));
}

/// Naive Poseidon Implementation

void Goldilocks::poseidon_naive(uint64_t (&state)[SPONGE_WIDTH])
{
#if ASM == 1
    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = gl_tom(state[i]);
    }
#endif
    uint8_t round_ctr = 0;
    full_rounds_naive(state, round_ctr);
    partial_rounds_naive(state, round_ctr);
    full_rounds_naive(state, round_ctr);
#if ASM == 1
    for (int i = 0; i < SPONGE_WIDTH; i++)
    {
        state[i] = gl_fromm(state[i]);
    }
#endif
}

inline void Goldilocks::full_rounds_naive(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
    for (uint8_t i = 0; i < HALF_N_FULL_ROUNDS; i++)
    {
        constant_layer_naive(state, round_ctr);
        sbox_layer_naive(state);
        mds_layer_naive(state);
        round_ctr += 1;
    }
}

void Goldilocks::constant_layer_naive(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
    for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
    {
#if ASM == 1
        state[i] = gl_add(state[i], Poseidon_goldilocks_constants::ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]);
#else
        state[i] = add_gl(state[i], Poseidon_goldilocks_constants::ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]);
        // state[i] = ((uint128_t)state[i] + (uint128_t)Poseidon_goldilocks_constants::ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]) % GOLDILOCKS_PRIME;
#endif
    }
}

void Goldilocks::sbox_layer_naive(uint64_t (&state)[SPONGE_WIDTH])
{
    for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
    {
        sbox_monomial_naive(state[i]);
    }
}

void Goldilocks::sbox_monomial_naive(uint64_t &x)
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
}

void Goldilocks::mds_layer_naive(uint64_t (&state_)[SPONGE_WIDTH])
{
    uint64_t state[SPONGE_WIDTH] = {0};
    std::memcpy(state, state_, SPONGE_WIDTH * sizeof(uint64_t));

    for (uint8_t r = 0; r < SPONGE_WIDTH; r++)
    {
        state_[r] = mds_row_shf_naive(r, state);
    }
}

uint64_t Goldilocks::mds_row_shf_naive(uint64_t r, uint64_t (&v)[SPONGE_WIDTH])
{

#if ASM == 1
    uint64_t res = 0;
    res = gl_mmul(v[r], Poseidon_goldilocks_constants::MDS_MATRIX_DIAG[r]);

    for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
    {
        res = gl_add(res, gl_mmul(v[(i + r) % SPONGE_WIDTH], Poseidon_goldilocks_constants::MDS_MATRIX_CIRC[i]));
    }
    return res;
#else
    uint128_t res = 0;
    res += (uint128_t)v[r] * (uint128_t)Poseidon_goldilocks_constants::MDS_MATRIX_DIAG[r];

    for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
    {
        res += (((uint128_t)v[(i + r) % SPONGE_WIDTH] * (uint128_t)Poseidon_goldilocks_constants::MDS_MATRIX_CIRC[i]));
    }
    return res % GOLDILOCKS_PRIME;
#endif
}

void Goldilocks::partial_rounds_naive(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
    for (uint8_t i = 0; i < N_PARTIAL_ROUNDS; i++)
    {
        constant_layer_naive(state, round_ctr);
        sbox_monomial_naive(state[0]);
        mds_layer_naive(state);
        round_ctr += 1;
    }
}
