#ifndef POSEIDON_GOLDILOCKS_OPT
#define POSEIDON_GOLDILOCKS_OPT

#include <inttypes.h>
#include <stdint.h> // for uint64_t
#include "poseidon_goldilocks_opt_constants.hpp"
#include "goldilocks/goldilocks.hpp"

#define SPONGE_WIDTH 12
#define MAX_WIDTH 12
#define HALF_N_FULL_ROUNDS 4
#define N_FULL_ROUNDS_TOTAL (2 * HALF_N_FULL_ROUNDS)
#define N_PARTIAL_ROUNDS 22
#define N_ROUNDS (N_FULL_ROUNDS_TOTAL + N_PARTIAL_ROUNDS)
#define RATE 8
#define CAPACITY 4

#define GOLDILOCKS_PRIME 0xFFFFFFFF00000001ULL

#define ASM 1

#define uint128_t __uint128_t

typedef unsigned __int128 uint128_t;

class Poseidon_goldilocks_opt
{
private:
    uint64_t sum_counter = 0;
    uint64_t mul_counter = 0;

    inline void static pow7(uint64_t &x)
    {
#if ASM == 1
        uint64_t x2 = Goldilocks::gl_mmul(x, x);
        uint64_t x3 = Goldilocks::gl_mmul(x, x2);
        uint64_t x4 = Goldilocks::gl_mmul(x2, x2);

        x = Goldilocks::gl_mmul(x3, x4);
#else
        uint128_t x2 = ((uint128_t)x * (uint128_t)x) % GOLDILOCKS_PRIME;
        uint128_t x4 = (x2 * x2) % GOLDILOCKS_PRIME;
        uint128_t x3 = ((uint128_t)x * x2) % GOLDILOCKS_PRIME;
        x = (x3 * x4) % GOLDILOCKS_PRIME;
#endif
    };

public:
    void static hash(uint64_t (&input)[SPONGE_WIDTH]);
    void static linear_hash(uint64_t *output, uint64_t *input, uint64_t size);

    void static print();
};

#endif // POSEIDON_GOLDILOCKS_OPT