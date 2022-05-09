#include "poseidon_goldilocks.hpp"
#include <cstring> //memset
#include <stdio.h>

const uint64_t Poseidon_goldilocks::Q = 0xFFFFFFFF00000001LL;
const uint64_t Poseidon_goldilocks::MM = 0xFFFFFFFeFFFFFFFFLL;
const uint64_t Poseidon_goldilocks::CQ = 0x00000000FFFFFFFFLL;
const uint64_t Poseidon_goldilocks::R2 = 0xFFFFFFFe00000001LL;

void Poseidon_goldilocks::hash(uint64_t (&state)[SPONGE_WIDTH])
{
#if ASM == 1
	for (int i = 0; i < SPONGE_WIDTH; i++)
	{
		state[i] = gl_tom(state[i]);
	}
#endif
	uint8_t round_ctr = 0;
	full_rounds(state, round_ctr);
	partial_rounds_naive(state, round_ctr);
	full_rounds(state, round_ctr);
#if ASM == 1
	for (int i = 0; i < SPONGE_WIDTH; i++)
	{
		state[i] = gl_fromm(state[i]);
	}
#endif
}

inline void Poseidon_goldilocks::full_rounds(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
	for (uint8_t i = 0; i < HALF_N_FULL_ROUNDS; i++)
	{
		constant_layer(state, round_ctr);
		sbox_layer(state);
		mds_layer(state);
		round_ctr += 1;
	}
}

void Poseidon_goldilocks::constant_layer(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
	for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
	{
#if ASM == 1
		state[i] = gl_add(state[i], ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]);
#else
		state[i] = add_gl(state[i], ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]);
		//state[i] = ((uint128_t)state[i] + (uint128_t)ALL_ROUND_CONSTANTS[i + SPONGE_WIDTH * round_ctr]) % GOLDILOCKS_PRIME;
#endif
	}
}

void Poseidon_goldilocks::sbox_layer(uint64_t (&state)[SPONGE_WIDTH])
{
	for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
	{
		sbox_monomial(state[i]);
	}
}

void Poseidon_goldilocks::sbox_layer_new(uint64_t (&state)[SPONGE_WIDTH])
{
	uint64_t s0_2 = gl_mmul(state[0], state[0]);
	uint64_t s1_2 = gl_mmul(state[1], state[1]);
	uint64_t s2_2 = gl_mmul(state[2], state[2]);
	uint64_t s3_2 = gl_mmul(state[3], state[3]);

	uint64_t s0_3 = gl_mmul(state[0], s0_2);
	uint64_t s1_3 = gl_mmul(state[1], s1_2);
	uint64_t s2_3 = gl_mmul(state[2], s2_2);
	uint64_t s3_3 = gl_mmul(state[3], s3_2);

	uint64_t s0_4 = gl_mmul(s0_2, s0_2);
	uint64_t s1_4 = gl_mmul(s1_2, s1_2);
	uint64_t s2_4 = gl_mmul(s2_2, s2_2);
	uint64_t s3_4 = gl_mmul(s3_2, s3_2);

	state[0] = gl_mmul(s0_3, s0_4);
	state[1] = gl_mmul(s1_3, s1_4);
	state[2] = gl_mmul(s2_3, s2_4);
	state[3] = gl_mmul(s3_3, s3_4);

	s0_2 = gl_mmul(state[4], state[4]);
	s1_2 = gl_mmul(state[5], state[5]);
	s2_2 = gl_mmul(state[6], state[6]);
	s3_2 = gl_mmul(state[7], state[7]);

	s0_3 = gl_mmul(state[4], s0_2);
	s1_3 = gl_mmul(state[5], s1_2);
	s2_3 = gl_mmul(state[6], s2_2);
	s3_3 = gl_mmul(state[7], s3_2);

	s0_4 = gl_mmul(s0_2, s0_2);
	s1_4 = gl_mmul(s1_2, s1_2);
	s2_4 = gl_mmul(s2_2, s2_2);
	s3_4 = gl_mmul(s3_2, s3_2);

	state[4] = gl_mmul(s0_3, s0_4);
	state[5] = gl_mmul(s1_3, s1_4);
	state[6] = gl_mmul(s2_3, s2_4);
	state[7] = gl_mmul(s3_3, s3_4);

	s0_2 = gl_mmul(state[8], state[8]);
	s1_2 = gl_mmul(state[9], state[9]);
	s2_2 = gl_mmul(state[10], state[10]);
	s3_2 = gl_mmul(state[11], state[11]);

	s0_3 = gl_mmul(state[8], s0_2);
	s1_3 = gl_mmul(state[9], s1_2);
	s2_3 = gl_mmul(state[10], s2_2);
	s3_3 = gl_mmul(state[11], s3_2);

	s0_4 = gl_mmul(s0_2, s0_2);
	s1_4 = gl_mmul(s1_2, s1_2);
	s2_4 = gl_mmul(s2_2, s2_2);
	s3_4 = gl_mmul(s3_2, s3_2);

	state[0] = gl_mmul(s0_3, s0_4);
	state[1] = gl_mmul(s1_3, s1_4);
	state[2] = gl_mmul(s2_3, s2_4);
	state[3] = gl_mmul(s3_3, s3_4);
}

void Poseidon_goldilocks::sbox_monomial(uint64_t &x)
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

void Poseidon_goldilocks::mds_layer(uint64_t (&state_)[SPONGE_WIDTH])
{
	uint64_t state[SPONGE_WIDTH] = {0};
	std::memcpy(state, state_, SPONGE_WIDTH * sizeof(uint64_t));

	for (uint8_t r = 0; r < SPONGE_WIDTH; r++)
	{
		state_[r] = mds_row_shf(r, state);
	}
}

uint64_t Poseidon_goldilocks::mds_row_shf(uint64_t r, uint64_t (&v)[SPONGE_WIDTH])
{

#if ASM == 1
	uint64_t res = 0;
	res = gl_mmul(v[r], MDS_MATRIX_DIAG[r]);

	for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
	{
		res = gl_add(res, gl_mmul(v[(i + r) % SPONGE_WIDTH], MDS_MATRIX_CIRC[i]));
	}
	return res;
#else
	uint128_t res = 0;
	res += (uint128_t)v[r] * (uint128_t)MDS_MATRIX_DIAG[r];
	
	for (uint8_t i = 0; i < SPONGE_WIDTH; i++)
	{
		res += (((uint128_t)v[(i + r) % SPONGE_WIDTH] * (uint128_t)MDS_MATRIX_CIRC[i]));
	}
	return res % GOLDILOCKS_PRIME;
#endif
}

void Poseidon_goldilocks::partial_rounds_naive(uint64_t (&state)[SPONGE_WIDTH], uint8_t &round_ctr)
{
	for (uint8_t i = 0; i < N_PARTIAL_ROUNDS; i++)
	{
		constant_layer(state, round_ctr);
		sbox_monomial(state[0]);
		mds_layer(state);
		round_ctr += 1;
	}
}
