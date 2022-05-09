#include <benchmark/benchmark.h>
#include "goldilocks/ntt_goldilocks.hpp"
#include "poseidon_goldilocks.hpp"
#include "poseidon_goldilocks_opt.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>

#define NUM_HASHES 1

static void FFT_BENCH(benchmark::State &state)
{

    uint64_t length = 1 << state.range(0);
    uint64_t *pol = (uint64_t *)malloc(length * sizeof(uint64_t));
    Goldilocks g(length, 8);

    // Fibonacci
    pol[0] = 0;
    pol[1] = 1;
    for (uint64_t i = 2; i < length; i++)
    {
        pol[i] = g.gl_add(pol[i - 1], pol[i - 2]);
    }

    for (auto _ : state)
    {
#pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < 1; i++)
        {
            g.ntt(pol, length);
            g.intt(pol, length);
        }
    }
}

static void LDE_BENCH(benchmark::State &state)
{

    uint64_t length = 1 << state.range(0);
    uint64_t extensionLength = 1 << (state.range(0) + 1);
    Goldilocks g(length, 8);
    Goldilocks ge(extensionLength, 8);
    uint64_t *pol = (uint64_t *)malloc(length * sizeof(uint64_t));
    uint64_t *pol_ext = (uint64_t *)malloc(extensionLength * sizeof(uint64_t));

    // Fibonacci
    pol[0] = 0;
    pol[1] = 1;
    for (uint64_t i = 2; i < length; i++)
    {
        pol[i] = g.gl_add(pol[i - 1], pol[i - 2]);
    }

    std::memcpy(pol_ext, pol, length * sizeof(uint64_t));

    uint64_t r = 1;
    uint64_t shift = 25;

    for (auto _ : state)
    {
#pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < 100; i++)
        {
            g.intt(pol_ext, length);
            for (uint j = 0; j < length; j++)
            {
                pol_ext[j] = g.gl_mmul2(pol_ext[j], r);
                r = g.gl_mmul2(r, shift);
            }

            ge.ntt(pol_ext, extensionLength);
        }
    }
}

static void main_poseidon(benchmark::State &state)
{
    uint64_t input[SPONGE_WIDTH] = {0};
    for (auto _ : state)
    {
        // memset(input, 0, SPONGE_WIDTH * sizeof(uint64_t));
        Poseidon_goldilocks::hash(input);
    }
    /*
    assert(input[0] == 0x3C18A9786CB0B359);
    assert(input[1] == 0xC4055E3364A246C3);
    assert(input[2] == 0x7953DB0AB48808F4);
    assert(input[3] == 0xC71603F33A1144CA);
    */
}

static void POSEIDON_BENCH(benchmark::State &state)
{
    uint64_t input[NUM_HASHES][SPONGE_WIDTH];

    for (uint64_t i = 0; i < NUM_HASHES; i++)
    {
        for (uint64_t j = 0; j < SPONGE_WIDTH; j++)
        {
            input[i][j] = 0;
        }
    }
    for (auto _ : state)
    {
//#pragma omp parallel for
//        for (uint64_t i = 0; i < NUM_HASHES; i++)
//        {
            memset(input[0], 0, SPONGE_WIDTH * sizeof(uint64_t));
            Poseidon_goldilocks_opt::hash(input[i]);
//        }
    }

    assert(input[0][0] == 0x3C18A9786CB0B359);
    assert(input[0][1] == 0xC4055E3364A246C3);
    assert(input[0][2] == 0x7953DB0AB48808F4);
    assert(input[0][3] == 0xC71603F33A1144CA);
}

// BENCHMARK(FFT_BENCH)->DenseRange(23, 23, 1)->Unit(benchmark::kMillisecond)->Iterations(6);
// BENCHMARK(LDE_BENCH)->DenseRange(23, 23, 1)->Unit(benchmark::kMillisecond)->Iterations(6);
//  BENCHMARK(main_poseidon)->Unit(benchmark::kMicrosecond);
BENCHMARK(POSEIDON_BENCH)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
