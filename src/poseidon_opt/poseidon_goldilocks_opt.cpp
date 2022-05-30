#include "poseidon_goldilocks_opt.hpp"
#include <cstring> //memset
#include <stdio.h>

void Poseidon_goldilocks_opt::hash(uint64_t (&state)[SPONGE_WIDTH])
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

void Poseidon_goldilocks_opt::linear_hash(uint64_t *output, uint64_t *input, uint64_t size)
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

		/*
		printf("linear_hash (");
		for (uint64_t j = 0; j < SPONGE_WIDTH; j++)
		{
			printf(" %#lx ", state[j]);
		}*/

		Poseidon_goldilocks_opt::hash(state);
		/*
		printf(") -> ( ");
		for (uint64_t j = 0; j < 4; j++)
		{
			printf(" %#lx ", state[j]);
		}
		printf(")\n");
	*/
		remaining -= n;
	}
	std::memcpy(output, state, 4 * sizeof(uint64_t));
}
