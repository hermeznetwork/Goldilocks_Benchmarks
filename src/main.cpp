#include <benchmark/benchmark.h>
#include "ffiasm/fr.hpp"
#include "ffiasm/fft.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <assert.h>

#define CHECK 0
#define NUM_THREADS 128

#define LENGTH (1 << 25)
#define NUM_BLOCKS 1
#define NUM_COLUMNS 64
#define NUM_PHASES 3

#define NUM_COLS_BLOCK NUM_COLUMNS / NUM_BLOCKS

static void FFT_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    uint64_t length = 1 << state.range(0);
    FFT<RawFr> fft(length, NUM_THREADS);
    RawFr::Element *pol, *pol0;

    pol = (RawFr::Element *)malloc(length * sizeof(RawFr::Element));
    pol0 = (RawFr::Element *)malloc(length * sizeof(RawFr::Element));

    // Fibonacci
    field.fromString(pol[0], "0");
    field.fromString(pol[1], "1");
    for (uint64_t i = 2; i < length; i++)
    {
        field.add(pol[i], pol[i - 1], pol[i - 2]);
    }
    std::memcpy(pol0, pol, length * sizeof(RawFr::Element));

    for (auto _ : state)
    {
        fft.fft(pol, length);
        if (CHECK)
        {
            fft.ifft(pol, length);
        }
    }
    if (CHECK)
    {
#pragma omp parallel for
        for (uint64_t i = 0; i < length; i++)
        {
            assert(field.eq(pol[i], pol0[i]));
        }
    }
    free(pol);
}

static void FFT_COLUMNS_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    uint64_t length = 1 << state.range(0);
    FFT<RawFr> fft(length, NUM_THREADS);
    RawFr::Element *pol, *pol0, *J;

    pol = (RawFr::Element *)malloc(length * (uint64_t)NUM_COLUMNS * sizeof(RawFr::Element));
    pol0 = (RawFr::Element *)malloc(length * sizeof(RawFr::Element));
    J = (RawFr::Element *)malloc(NUM_COLUMNS * sizeof(RawFr::Element));

    // Fibonacci
    field.fromString(pol0[0], "0");
    field.fromString(pol0[1], "1");
    for (uint64_t i = 2; i < length; i++)
    {
        field.add(pol0[i], pol0[i - 1], pol0[i - 2]);
    }

    // Column offsets
#pragma omp parallel for
    for (uint64_t j = 0; j < NUM_COLUMNS; ++j)
    {
        field.fromUI(J[j], j);
    }

    // Block initialization fibonacci + column offset
#pragma omp parallel for collapse(2)
    for (uint j = 0; j < NUM_COLUMNS; j++)
    {
        for (uint64_t i = 0; i < length; i++)
        {
            field.add(pol[j * length + i], pol0[i], J[j]);
        }
    }

    // Benchmark
    for (auto _ : state)
    {
        for (int j = 0; j < NUM_COLUMNS; ++j)
        {
            fft.fft(pol + j * NUM_COLUMNS, length);
            if (CHECK)
            {
                fft.ifft(pol + j * NUM_COLUMNS, length);
            }
        }
    }

    // CHECK RESULT
    if (CHECK)
    {
#pragma omp parallel for collapse(2)
        for (uint j = 0; j < NUM_COLUMNS; j++)
        {
            for (uint64_t i = 0; i < length; i++)
            {
                RawFr::Element tmp;
                field.add(tmp, pol0[i], J[j]);
                assert(field.eq(pol[j * length + i], tmp));
            }
        }
    }
    free(pol);
}

static void FFT_BLOCK_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    uint64_t exp = state.range(0);
    uint64_t length = 1 << exp;
    FFT<RawFr> fft(length, NUM_THREADS);
    RawFr::Element *pol, *pol0, *J;

    pol = (RawFr::Element *)malloc(length * (uint64_t)NUM_COLUMNS * sizeof(RawFr::Element));
    pol0 = (RawFr::Element *)malloc(length * sizeof(RawFr::Element));
    J = (RawFr::Element *)malloc(NUM_COLUMNS * sizeof(RawFr::Element));

    // Fibonacci
    field.fromString(pol0[0], "0");
    field.fromString(pol0[1], "1");
    for (uint64_t i = 2; i < length; i++)
    {
        field.add(pol0[i], pol0[i - 1], pol0[i - 2]);
    }
    // Column offsets
#pragma omp parallel for
    for (uint64_t j = 0; j < NUM_COLUMNS; ++j)
    {
        field.fromUI(J[j], j);
    }

    // Block initialization fibonacci + column offset
#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < length; i++)
    {
        for (uint j = 0; j < NUM_COLUMNS; j++)
        {
            field.add(pol[i * NUM_COLUMNS + j], pol0[i], J[j]);
        }
    }
    uint64_t numblocks = NUM_BLOCKS;
    if (exp == 26)
    {
        numblocks = 4;
    }

    // Benchmark
    for (auto _ : state)
    {
        fft.fft_block(pol, pol, length, NUM_COLUMNS, NUM_PHASES, numblocks);
        if (CHECK)
        {
            fft.ifft_block(pol, pol, length, NUM_COLUMNS, NUM_PHASES, numblocks);
        }
    }

    // CHECK RESULT
    if (CHECK)
    {
#pragma omp parallel for collapse(2)
        for (uint64_t i = 0; i < length; i++)
        {
            for (uint j = 0; j < NUM_COLUMNS; j++)
            {
                RawFr::Element tmp;
                field.add(tmp, pol0[i], J[j]);
                assert(field.eq(pol[i * NUM_COLUMNS + j], tmp));
            }
        }
    }
    free(pol);
}

static void FFT_BLOCK_TEMPLATE_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    FFT<RawFr> fft(LENGTH, NUM_THREADS);
    RawFr::Element *pol, *pol0, *J;

    pol = (RawFr::Element *)malloc(LENGTH * (uint64_t)NUM_COLUMNS * sizeof(RawFr::Element));
    pol0 = (RawFr::Element *)malloc(LENGTH * sizeof(RawFr::Element));
    J = (RawFr::Element *)malloc(NUM_COLUMNS * sizeof(RawFr::Element));

    // Fibonacci
    field.fromString(pol0[0], "0");
    field.fromString(pol0[1], "1");
    for (uint64_t i = 2; i < LENGTH; i++)
    {
        field.add(pol0[i], pol0[i - 1], pol0[i - 2]);
    }
    // Column offsets
#pragma omp parallel for
    for (uint64_t j = 0; j < NUM_COLUMNS; ++j)
    {
        field.fromUI(J[j], j);
    }

    // Block initialization fibonacci + column offset
#pragma omp parallel for collapse(2)
    for (uint64_t i = 0; i < LENGTH; i++)
    {
        for (uint j = 0; j < NUM_COLUMNS; j++)
        {
            field.add(pol[i * NUM_COLUMNS + j], pol0[i], J[j]);
        }
    }

    // Benchmark
    for (auto _ : state)
    {
        fft.fft_block_divisible<LENGTH, NUM_COLUMNS, NUM_PHASES, NUM_BLOCKS, NUM_COLS_BLOCK>(pol, pol);
        if (CHECK)
        {
            fft.ifft_block(pol, pol, LENGTH, NUM_COLUMNS, NUM_PHASES, NUM_BLOCKS);
        }
    }

    // CHECK RESULT
    if (CHECK)
    {
#pragma omp parallel for collapse(2)
        for (uint64_t i = 0; i < LENGTH; i++)
        {
            for (uint j = 0; j < NUM_COLUMNS; j++)
            {
                RawFr::Element tmp;
                field.add(tmp, pol0[i], J[j]);
                assert(field.eq(pol[i * NUM_COLUMNS + j], tmp));
            }
        }
    }
    free(pol);
}

BENCHMARK(FFT_BENCHMARK)->DenseRange(23, 27, 1)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(FFT_COLUMNS_BENCHMARK)->DenseRange(23, 26, 1)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(FFT_BLOCK_BENCHMARK)->DenseRange(23, 26, 1)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(FFT_BLOCK_TEMPLATE_BENCHMARK)->DenseRange(25, 25, 1)->Unit(benchmark::kMillisecond)->Iterations(1);

BENCHMARK_MAIN();
