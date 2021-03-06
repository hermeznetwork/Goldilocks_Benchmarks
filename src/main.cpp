#include <benchmark/benchmark.h>
#include "goldilocks/goldilocks.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <math.h> /* ceil */

#define NUM_COLS 100
#define NUM_HASHES 10000
#define FFT_SIZE (1 << 23)
#define NUM_ROWS (1 << 25)
#define BLOWUP_FACTOR 4

#define HASH_SIZE 4
#define RATE 8
#define CAPACITY 4

static void POSEIDON_BENCH(benchmark::State &state)
{
    uint64_t input_size = NUM_HASHES * SPONGE_WIDTH;

    uint64_t fibonacci[input_size];
    uint64_t pol_output[input_size];

    // Test vector: Fibonacci series
    // 0 1 1 2 3 5 8 13 ... NUM_HASHES * SPONGE_WIDTH ...
    fibonacci[0] = 0;
    fibonacci[1] = 1;
    for (uint64_t i = 2; i < input_size; i++)
    {
        fibonacci[i] = Goldilocks::gl_add(fibonacci[i - 1], fibonacci[i - 2]);
    }

    // Benchmark
    for (auto _ : state)
    {
#pragma omp barrier
// Every thread process chunks of SPONGE_WIDTH elements
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_HASHES; i++)
        {
            uint64_t pol_input_t[SPONGE_WIDTH];
            std::memcpy(pol_input_t, &fibonacci[i * SPONGE_WIDTH], SPONGE_WIDTH * sizeof(uint64_t));
            Goldilocks::poseidon(pol_input_t);
            std::memcpy(&pol_output[i * SPONGE_WIDTH], &pol_input_t[0], SPONGE_WIDTH * sizeof(uint64_t));
        }
#pragma omp barrier
    }
    // Check poseidon results poseidon ( 0 1 1 2 3 5 8 13 21 34 55 89 )
    assert(pol_output[0] == 0x3095570037f4605d);
    assert(pol_output[1] == 0x3d561b5ef1bc8b58);
    assert(pol_output[2] == 0x8129db5ec75c3226);
    assert(pol_output[3] == 0x8ec2b67afb6b87ed);

    // Rate = time to process 1 posseidon per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((double)NUM_HASHES / (double)state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void POSEIDON_NAIVE_BENCH(benchmark::State &state)
{
    uint64_t input_size = NUM_HASHES * SPONGE_WIDTH;

    uint64_t fibonacci[input_size];
    uint64_t pol_output[input_size];

    // Test vector: Fibonacci series
    // 0 1 1 2 3 5 8 13 ... NUM_HASHES * SPONGE_WIDTH ...
    fibonacci[0] = 0;
    fibonacci[1] = 1;
    for (uint64_t i = 2; i < input_size; i++)
    {
        fibonacci[i] = Goldilocks::gl_add(fibonacci[i - 1], fibonacci[i - 2]);
    }

    // Benchmark
    for (auto _ : state)
    {
#pragma omp barrier
// Every thread process chunks of SPONGE_WIDTH elements
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_HASHES; i++)
        {
            uint64_t pol_input_t[SPONGE_WIDTH];
            std::memcpy(pol_input_t, &fibonacci[i * SPONGE_WIDTH], SPONGE_WIDTH * sizeof(uint64_t));
            Goldilocks::poseidon_naive(pol_input_t);
            std::memcpy(&pol_output[i * SPONGE_WIDTH], &pol_input_t[0], SPONGE_WIDTH * sizeof(uint64_t));
        }
#pragma omp barrier
    }
    // Check poseidon results poseidon ( 0 1 1 2 3 5 8 13 21 34 55 89 )
    assert(pol_output[0] == 0x3095570037f4605d);
    assert(pol_output[1] == 0x3d561b5ef1bc8b58);
    assert(pol_output[2] == 0x8129db5ec75c3226);
    assert(pol_output[3] == 0x8ec2b67afb6b87ed);

    // Rate = time to process 1 posseidon per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((double)NUM_HASHES / (double)state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_LINEAR_HASH_SINGLE_BENCH(benchmark::State &state)
{
    uint64_t cols[NUM_COLS];

    // Test vector: 1 2 3 4 5 6 ... NUM_COLS
    for (uint64_t i = 0; i < NUM_COLS; i++)
    {
        cols[i] = Goldilocks::gl_add(i, 1);
    }

    uint64_t result[HASH_SIZE];
    for (auto _ : state)
    {
        uint64_t intermediate[NUM_COLS];
        std::memcpy(intermediate, cols, NUM_COLS * sizeof(uint64_t));

        Goldilocks::linear_hash(&result[0], &intermediate[0], NUM_COLS);
    }
    // Check linear_hash results linear_hash ( 1 2 3 4 5 ... 100 )
    if (NUM_COLS == 100)
    {
        assert(result[0] == 0xe82d873de8d5ea68);
        assert(result[1] == 0xe9f839e87c5dc258);
        assert(result[2] == 0x8cfaf75cfe46672d);
        assert(result[3] == 0x8c2147709952c96e);
    }

    // Rate = time to process 1 posseidon per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((double)ceil((float)NUM_COLS / (float)RATE), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(NUM_COLS * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_LINEAR_HASH_BENCH(benchmark::State &state)
{
    uint64_t *cols = (uint64_t *)malloc((uint64_t)NUM_COLS * (uint64_t)NUM_ROWS * sizeof(uint64_t));

    // Test vector: Fibonacci series on the columns and increase the initial values to the right,
    // 1 2 3 4  5  6  ... NUM_COLS
    // 1 2 3 4  5  6  ... NUM_COLS
    // 2 4 6 8  10 12 ... NUM_COLS + NUM_COLS
    // 3 6 9 12 15 18 ... NUM_COLS + NUM_COLS + NUM_COLS
    for (uint64_t i = 0; i < NUM_COLS; i++)
    {
        cols[i] = Goldilocks::gl_add(i, 1);
        cols[i + NUM_COLS] = Goldilocks::gl_add(i, 1);
    }
    for (uint64_t j = 2; j < NUM_ROWS; j++)
    {
        for (uint64_t i = 0; i < NUM_COLS; i++)
        {
            cols[j * NUM_COLS + i] = Goldilocks::gl_add(cols[(j - 2) * NUM_COLS + i], cols[(j - 1) * NUM_COLS + i]);
        }
    }

    uint64_t *result = (uint64_t *)malloc((uint64_t)HASH_SIZE * (uint64_t)NUM_ROWS * sizeof(uint64_t));
    // Benchmark
    for (auto _ : state)
    {
#pragma omp barrier
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_ROWS; i++)
        {
            uint64_t intermediate[NUM_COLS];
            uint64_t temp_result[HASH_SIZE];

            std::memcpy(&intermediate[0], &cols[i * NUM_COLS], NUM_COLS * sizeof(uint64_t));
            Goldilocks::linear_hash(temp_result, intermediate, NUM_COLS);
            std::memcpy(&result[i * HASH_SIZE], &temp_result[0], HASH_SIZE * sizeof(uint64_t));
        }
#pragma omp barrier
    }
    // Check linear_hash results linear_hash ( 1 2 3 4 5 ... 100 )
    if (NUM_COLS == 100)
    {
        assert(result[0] == 0xe82d873de8d5ea68);
        assert(result[1] == 0xe9f839e87c5dc258);
        assert(result[2] == 0x8cfaf75cfe46672d);
        assert(result[3] == 0x8c2147709952c96e);
    }

    free(cols);
    free(result);
    // Rate = time to process 1 linear hash per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((float)NUM_ROWS * (float)ceil((float)NUM_COLS / (float)RATE) / state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter((uint64_t)NUM_ROWS * (uint64_t)NUM_COLS * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_MERKLE_TREE_BENCH(benchmark::State &state)
{
    // Test vector: Fibonacci series on the columns and increase the initial values to the right,
    // 1 2 3 4  5  6  ... NUM_COLS
    // 1 2 3 4  5  6  ... NUM_COLS
    // 2 4 6 8  10 12 ... NUM_COLS + NUM_COLS
    // 3 6 9 12 15 18 ... NUM_COLS + NUM_COLS + NUM_COLS
    uint64_t *cols = (uint64_t *)malloc((uint64_t)NUM_COLS * (uint64_t)NUM_ROWS * sizeof(uint64_t));
    uint64_t *result = (uint64_t *)malloc((uint64_t)HASH_SIZE * (uint64_t)NUM_ROWS * sizeof(uint64_t));

    for (uint64_t i = 0; i < NUM_COLS; i++)
    {
        cols[i] = Goldilocks::gl_add(i, 1);
        cols[i + NUM_COLS] = Goldilocks::gl_add(i, 1);
    }
    for (uint64_t j = 2; j < NUM_ROWS; j++)
    {
        for (uint64_t i = 0; i < NUM_COLS; i++)
        {
            cols[j * NUM_COLS + i] = Goldilocks::gl_add(cols[(j - 2) * NUM_COLS + i], cols[(j - 1) * NUM_COLS + i]);
        }
    }

    // Benchmark
    for (auto _ : state)
    {
        // Linear hash of all rows
#pragma omp barrier
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_ROWS; i++)
        {
            uint64_t intermediate[NUM_COLS];
            uint64_t temp_result[HASH_SIZE];

            std::memcpy(&intermediate[0], &cols[i * NUM_COLS], NUM_COLS * sizeof(uint64_t));
            Goldilocks::linear_hash(temp_result, intermediate, NUM_COLS);
            std::memcpy(&result[i * HASH_SIZE], &temp_result[0], HASH_SIZE * sizeof(uint64_t));
        }
#pragma omp barrier

        // Build the merkle tree
        uint64_t pending = NUM_ROWS;
        while (pending > 1)
        {
#pragma omp barrier
#pragma omp parallel for
            for (uint64_t j = 0; j < NUM_ROWS; j += (2 * NUM_ROWS / pending))
            {
                uint64_t pol_input[SPONGE_WIDTH] = {0};
                std::memcpy(pol_input, &result[j * CAPACITY], CAPACITY * sizeof(uint64_t));
                std::memcpy(&pol_input[CAPACITY], &result[(j + (NUM_ROWS / pending)) * CAPACITY], CAPACITY * sizeof(uint64_t));

                Goldilocks::poseidon(pol_input);
                std::memcpy(&result[j * CAPACITY], pol_input, CAPACITY * sizeof(uint64_t));
            }
            pending = pending / 2;
        }
#pragma omp barrier
    }

    if ((NUM_COLS == 100) && (NUM_ROWS == (1 << 25)))
    {
        assert(result[0] == 0x2c7c38c7291c5340);
        assert(result[1] == 0x98644ad7c8cccab8);
        assert(result[2] == 0xc48264824e14760b);
        assert(result[3] == 0xb049d47bf63da1e4);
    }

    free(cols);
    free(result);
    // Rate = time to process 1 posseidon per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((((float)NUM_ROWS * (float)ceil((float)NUM_COLS / (float)RATE)) + log2(NUM_ROWS)) / state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter((uint64_t)NUM_ROWS * (uint64_t)NUM_COLS * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_iNTT_BENCH(benchmark::State &state)
{
    uint64_t *cols = (uint64_t *)malloc((uint64_t)NUM_COLS * (uint64_t)FFT_SIZE * sizeof(uint64_t));
    Goldilocks g(FFT_SIZE, state.range(0));
    // Test vector: Fibonacci series on the columns and increase the initial values to the right,
    // 1 2 3 4  5  6  ... NUM_COLS
    // 1 2 3 4  5  6  ... NUM_COLS
    // 2 4 6 8  10 12 ... NUM_COLS + NUM_COLS
    // 3 6 9 12 15 18 ... NUM_COLS + NUM_COLS + NUM_COLS

    for (uint64_t j = 0; j < NUM_COLS; j++)
    {
        uint64_t offset = j * FFT_SIZE;
        cols[offset] = j + 1;
        cols[offset + 1] = j + 1;

        for (uint64_t i = 2; i < FFT_SIZE; i++)
        {
            cols[i + offset] = Goldilocks::gl_add(cols[(i - 2) + offset], cols[(i - 1) + offset]);
        }
    }

    // Benchmark
    for (auto _ : state)
    {
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_COLS; i++)
        {
            g.intt(&cols[i * FFT_SIZE], FFT_SIZE);
        }
    }

    free(cols);
    state.counters["Rate"] = benchmark::Counter((float)NUM_COLS / state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter((uint64_t)FFT_SIZE * (uint64_t)NUM_COLS * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_LDE_BENCH(benchmark::State &state)
{
    Goldilocks g(FFT_SIZE, state.range(0));
    Goldilocks ge(FFT_SIZE * BLOWUP_FACTOR, state.range(0));
    uint64_t *pol_ext_intt = (uint64_t *)malloc((uint64_t)FFT_SIZE * (uint64_t)BLOWUP_FACTOR * sizeof(uint64_t));
    uint64_t *pol_ext = (uint64_t *)malloc((uint64_t)FFT_SIZE * (uint64_t)BLOWUP_FACTOR * sizeof(uint64_t));

    // Fibonacci
    pol_ext_intt[0] = 0;
    pol_ext_intt[1] = 1;
    for (uint64_t i = 2; i < FFT_SIZE; i++)
    {
        pol_ext_intt[i] = g.gl_add(pol_ext_intt[i - 1], pol_ext_intt[i - 2]);
    }
    for (uint64_t i = FFT_SIZE; i < FFT_SIZE * BLOWUP_FACTOR; i++)
    {
        pol_ext_intt[i] = 0;
    }

    uint64_t *r = (uint64_t *)malloc((uint64_t)FFT_SIZE * (uint64_t)BLOWUP_FACTOR * sizeof(uint64_t));
    uint64_t shift = 7;
    r[0] = 1;

    for (uint j = 1; j < FFT_SIZE * BLOWUP_FACTOR; j++)
    {
        r[j] = g.gl_mmul(r[j - 1], shift);
    }

    g.intt(pol_ext_intt, FFT_SIZE);

    // Benchmark
    for (auto _ : state)
    {
        std::memcpy(pol_ext, pol_ext_intt, (uint64_t)FFT_SIZE * (uint64_t)BLOWUP_FACTOR * sizeof(uint64_t));
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_COLS; i++)
        {
            for (uint j = 0; j < FFT_SIZE * BLOWUP_FACTOR; j++)
            {
                pol_ext[j] = g.gl_mmul(pol_ext[j], r[j]);
            }

            ge.ntt(pol_ext, FFT_SIZE * BLOWUP_FACTOR);
        }
    }

    free(pol_ext);
    // Rate = time to process 1 posseidon per thread
    // BytesProcessed = total bytes processed per second on every iteration
    state.counters["Rate"] = benchmark::Counter((float)NUM_COLS / state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(((uint64_t)FFT_SIZE + (uint64_t)FFT_SIZE * 2) * (uint64_t)NUM_COLS * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

// DenseRange(1, 1, 1) -> 1 Thread
// RangeMultiplier(2)->Range(2, omp_get_max_threads()) -> From 2 threads to omp_get_max_threads() every two
// DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2) -> +/- 4 around the number of cores every two

BENCHMARK(POSEIDON_BENCH)
    ->Unit(benchmark::kMicrosecond)
    ->DenseRange(1, 1, 1)
    ->RangeMultiplier(2)
    ->Range(2, omp_get_max_threads())
    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->UseRealTime();

BENCHMARK(POSEIDON_NAIVE_BENCH)
    ->Unit(benchmark::kMicrosecond)
    ->DenseRange(1, 1, 1)
    ->RangeMultiplier(2)
    ->Range(2, omp_get_max_threads())
    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->UseRealTime();
BENCHMARK(DISABLED_LINEAR_HASH_SINGLE_BENCH)
    ->Unit(benchmark::kMicrosecond)
    ->UseRealTime();

BENCHMARK(DISABLED_LINEAR_HASH_BENCH)
    ->Unit(benchmark::kSecond)
    //    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->DenseRange(omp_get_max_threads(), omp_get_max_threads(), 1)
    ->UseRealTime();

BENCHMARK(DISABLED_MERKLE_TREE_BENCH)
    ->Unit(benchmark::kSecond)
    //    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->DenseRange(omp_get_max_threads(), omp_get_max_threads(), 1)
    ->UseRealTime();

BENCHMARK(DISABLED_iNTT_BENCH)
    ->Unit(benchmark::kSecond)
    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->DenseRange(omp_get_max_threads(), omp_get_max_threads(), 1)
    ->UseRealTime()
    ->Iterations(1);

BENCHMARK(DISABLED_LDE_BENCH)
    ->Unit(benchmark::kSecond)
    //    ->DenseRange(omp_get_max_threads() / 2 - 4, omp_get_max_threads() / 2 + 4, 2)
    ->DenseRange(omp_get_max_threads(), omp_get_max_threads(), 1)
    ->UseRealTime();

BENCHMARK_MAIN();
