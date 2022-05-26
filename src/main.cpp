#include <benchmark/benchmark.h>
#include "goldilocks/ntt_goldilocks.hpp"
#include "poseidon_goldilocks_opt.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <openssl/md5.h>
#include <sstream>
#include <math.h> /* round, floor, ceil, trunc */

#define NUM_COLS 100
#define RATE 8
#define CAPACITY 4
#define NUM_ROWS (1 << 25)

void linear_hash(uint64_t *output, uint64_t *input, uint64_t size)
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
        }
        printf(") -> ( ");
        */
        Poseidon_goldilocks_opt::hash(state);
        /*
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

void linear_hash_bobbin(uint64_t *output, uint64_t *input, uint64_t size)
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
        }
        printf(") -> ( ");
        */
        Poseidon_goldilocks_opt::hash(state);
        /*
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

static void DISABLED_POSEIDON_SINGLE_BENCH(benchmark::State &state)
{
    uint64_t hash_input_size = NUM_ROWS * (RATE + CAPACITY);
    uint64_t *fibonacci = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));
    uint64_t *pol_input = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));
    uint64_t *pol_output = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));

    Goldilocks g(NUM_ROWS, 8);

    // Fibonacci
    fibonacci[0] = 0;
    fibonacci[1] = 1;
    for (uint64_t i = 2; i < NUM_ROWS * (RATE + CAPACITY); i++)
    {
        fibonacci[i] = g.gl_add(fibonacci[i - 1], fibonacci[i - 2]);
    }

    for (auto _ : state)
    {
#pragma omp parallel for num_threads(state.range(0))
        for (uint64_t i = 0; i < NUM_ROWS; i++)
        {
            uint64_t pol_input_t[12];
            std::memcpy(pol_input_t, &fibonacci[i * (RATE + CAPACITY)], (RATE + CAPACITY) * sizeof(uint64_t));
            Poseidon_goldilocks_opt::hash(pol_input_t);
            std::memcpy(&pol_output[i * (RATE + CAPACITY)], &pol_input_t[0], (RATE + CAPACITY) * sizeof(uint64_t));
        }
    }
    /*
  // FIBONACCI
  assert(pol_input[0] == 0x3095570037f4605d);
  assert(pol_input[1] == 0x3d561b5ef1bc8b58);
  assert(pol_input[2] == 0x8129db5ec75c3226);
  assert(pol_input[3] == 0x8ec2b67afb6b87ed);


  assert(pol_input[0] == 0x3C18A9786CB0B359);
  assert(pol_input[1] == 0xC4055E3364A246C3);
  assert(pol_input[2] == 0x7953DB0AB48808F4);
  assert(pol_input[3] == 0xC71603F33A1144CA);
  */
    unsigned char result[MD5_DIGEST_LENGTH];
    std::string currentHash;

    std::string pol;
    for (uint64_t i = 0; i < NUM_ROWS * (RATE + CAPACITY); i++)
    {

        std::ostringstream os;
        os << pol_output[i];
        pol += os.str();
    }

    MD5((unsigned char *)pol.c_str(), pol.size(), result);

    char tempHash[32];
    for (int i = 0; i < MD5_DIGEST_LENGTH; i++)
    {
        sprintf(tempHash, "%02x", result[i]);
        currentHash.append(tempHash);
    }
    if (NUM_ROWS == 25000)
        assert(currentHash == "894b1d8d817a0bd1ad34b6ac8e6a91e4");
    free(fibonacci);
    free(pol_input);
    state.counters["Rate"] = benchmark::Counter((double)NUM_ROWS / (double)state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(hash_input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_LINEAL_HASH_BENCH(benchmark::State &state)
{
    uint64_t num_batches = ceil((float)NUM_COLS / (float)RATE);
    uint64_t hash_input_size = num_batches * RATE + CAPACITY;

    uint64_t *fibonnaci = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));
    uint64_t *pol_input = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));

    Goldilocks g(NUM_ROWS, 8);

    // Fibonacci
    fibonnaci[0] = 1;
    fibonnaci[1] = 1;
    for (uint64_t i = 2; i < NUM_COLS; i++)
    {
        fibonnaci[i] = g.gl_add(fibonnaci[i - 1], fibonnaci[i - 2]);
    }

    for (uint64_t i = NUM_COLS; i < hash_input_size; i++)
    {
        fibonnaci[i] = 0;
    }

    for (auto _ : state)
    {
        // Reorg
        for (uint64_t i = 0; i < num_batches; i++)
        {
            std::memcpy(&pol_input[hash_input_size - CAPACITY - RATE * (i + 1)], &fibonnaci[i * RATE], std::min((uint64_t)RATE, NUM_COLS - i * RATE) * sizeof(uint64_t));
        }

        for (int i = hash_input_size - SPONGE_WIDTH; i >= 0; i -= RATE)
        {
            Poseidon_goldilocks_opt::hash((uint64_t(&)[SPONGE_WIDTH])pol_input[i]);
        }
    }

    state.counters["Rate"] = benchmark::Counter(num_batches, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(hash_input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);

    free(fibonnaci);
    free(pol_input);
}

static void DISABLED_LINEAL_HASH_BOBIN_BENCH(benchmark::State &state)
{
    uint64_t num_batches = ceil((float)NUM_COLS / (float)RATE);
    uint64_t hash_input_size = num_batches * RATE;

    uint64_t *fibonnaci = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));
    uint64_t *pol_input = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));

    Goldilocks g(NUM_ROWS, 8);

    // Fibonacci
    fibonnaci[0] = 1;
    fibonnaci[1] = 1;
    for (uint64_t i = 2; i < NUM_COLS; i++)
    {
        fibonnaci[i] = g.gl_add(fibonnaci[i - 1], fibonnaci[i - 2]);
    }

    for (uint64_t i = NUM_COLS; i < hash_input_size; i++)
    {
        fibonnaci[i] = 0;
    }

    for (auto _ : state)
    {
        uint64_t pol_input[SPONGE_WIDTH] = {0};
        std::memcpy(&pol_input, fibonnaci, SPONGE_WIDTH * sizeof(uint64_t));

        for (uint64_t i = 0; i < num_batches; i++)
        {
            Poseidon_goldilocks_opt::hash(pol_input);

            for (uint64_t j = 0; j < RATE; j++)
            {
                pol_input[j] = g.gl_add(pol_input[j], fibonnaci[j + i * RATE]);
            }
        }
    }

    state.counters["Rate"] = benchmark::Counter(num_batches, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    state.counters["BytesProcessed"] = benchmark::Counter(hash_input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);

    free(fibonnaci);
    free(pol_input);
}

static void TREE_BENCH(benchmark::State &state)
{
    uint64_t *pol_0 = (uint64_t *)malloc(NUM_ROWS * sizeof(uint64_t));
    uint64_t **cols = (uint64_t **)malloc(NUM_COLS * sizeof(uint64_t **));

    Goldilocks g(NUM_ROWS, 8);

    // Fibonacci
    pol_0[0] = 1;
    pol_0[1] = 1;
    for (uint64_t j = 2; j < NUM_ROWS; j++)
    {
        pol_0[j] = g.gl_add(pol_0[j - 1], pol_0[j - 2]);
    }
    cols[0] = pol_0;

    for (uint64_t i = 1; i < NUM_COLS; i++)
    {
        uint64_t *pol = (uint64_t *)malloc(NUM_ROWS * sizeof(uint64_t));
        pol[0] = g.gl_add(cols[i - 1][0], 1);
        pol[1] = g.gl_add(cols[i - 1][1], 1);
        for (uint64_t j = 2; j < NUM_ROWS; j++)
        {
            pol[j] = g.gl_add(pol[j - 1], pol[j - 2]);
        }
        cols[i] = pol;
    }
    /*
    printf("---------\n");
    for (uint64_t i = 0; i < NUM_ROWS; i++)
    {
        for (uint64_t j = 0; j < NUM_COLS; j++)
        {
            printf("(%lu,%lu): %#lx\t\t", i, j, cols[j][i]);
        }

        printf("\n");
    }
    */
    for (auto _ : state)
    {

        /*
        for (uint64_t i = 0; i < num_batches; i++)
        {
            std::memcpy(&pol_input[hash_input_size - CAPACITY - RATE * (i + 1)], &(rows[j])[i * RATE], std::min((uint64_t)RATE, NUM_COLS - i * RATE) * sizeof(uint64_t));
        }*/
        uint64_t *intermediate = (uint64_t *)malloc((uint64_t)NUM_COLS * (uint64_t)NUM_ROWS * sizeof(uint64_t));

#pragma omp parallel for
        for (uint64_t i = 0; i < NUM_ROWS; i++)
        {
            for (uint64_t j = 0; j < NUM_COLS; j++)
            {
                intermediate[i * NUM_COLS + j] = cols[j][i];
            }
        }
        uint64_t *intermediate_2 = (uint64_t *)malloc(4 * NUM_ROWS * sizeof(uint64_t));

        // printf("---------\n");
#pragma omp parallel for
        for (uint64_t i = 0; i < NUM_ROWS; i++)
        {
            linear_hash(intermediate_2 + i * 4, intermediate + (i * NUM_COLS), NUM_COLS);
            /*
            printf("result: ");
            for (uint64_t j = 0; j < 4; j++)
            {
                printf(" %#lx ", intermediate_2[i * 4 + j]);
            }
            printf("\n---------\n");
            */
        }
        uint64_t pending = NUM_ROWS;
        // printf("Merkle Tree\n");
        while (pending > 1)
        {

            // printf("##########\n");

            for (uint64_t j = 0; j < NUM_ROWS; j += (2 * NUM_ROWS / pending))
            {

                uint64_t pol_input[SPONGE_WIDTH] = {0};
                std::memcpy(pol_input, &intermediate_2[j * CAPACITY], CAPACITY * sizeof(uint64_t));
                std::memcpy(&pol_input[CAPACITY], &intermediate_2[(j + (NUM_ROWS / pending)) * CAPACITY], CAPACITY * sizeof(uint64_t));

                /*
                printf("hash (");
                for (uint64_t k = 0; k < SPONGE_WIDTH; k++)
                {
                    printf(" %#lx ", pol_input[k]);
                }
                printf(") -> \n\t(");
                */

                Poseidon_goldilocks_opt::hash(pol_input);
                std::memcpy(&intermediate_2[j * CAPACITY], pol_input, CAPACITY * sizeof(uint64_t));
                /*
                for (uint64_t j = 0; j < 4; j++)
                {
                    printf(" %#lx ", pol_input[j]);
                }
                printf(")\n");
                */
            }
            pending = pending / 2;
        }
    }
    /*
    printf("________________________________\n");

    for (uint64_t i = 0; i < CAPACITY; i++)
    {
        printf("rows[0][%lu]: %020lu\n", i, (rows[0])[i]);
    }*/

    // state.counters["Rate"] = benchmark::Counter((double)NUM_ROWS * num_batches / (double)state.range(0), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    // state.counters["BytesProcessed"] = benchmark::Counter(NUM_ROWS * hash_input_size * sizeof(uint64_t), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

static void DISABLED_FFT_BENCH(benchmark::State &state)
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

static void DISABLED_LDE_BENCH(benchmark::State &state)
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

BENCHMARK(DISABLED_POSEIDON_SINGLE_BENCH)->Unit(benchmark::kMicrosecond)->DenseRange(1, 1, 1)->RangeMultiplier(2)->Range(2, 2 << 4)->DenseRange(50, 70, 1);
BENCHMARK(DISABLED_LINEAL_HASH_BENCH)->Unit(benchmark::kMicrosecond);
BENCHMARK(DISABLED_LINEAL_HASH_BOBIN_BENCH)->Unit(benchmark::kMicrosecond);

BENCHMARK(TREE_BENCH)->Unit(benchmark::kMicrosecond)->Iterations(1); //->DenseRange(55, 65, 1); //->RangeMultiplier(2)->Range(2, 2 << 4)->DenseRange(50, 70, 1);

BENCHMARK(DISABLED_FFT_BENCH)
    ->DenseRange(23, 23, 1)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(6);
BENCHMARK(DISABLED_LDE_BENCH)->DenseRange(23, 23, 1)->Unit(benchmark::kMillisecond)->Iterations(6);

BENCHMARK_MAIN();
