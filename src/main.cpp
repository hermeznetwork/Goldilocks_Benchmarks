#include <benchmark/benchmark.h>
//#include "ffiasm/fr.hpp"
//#include "ffiasm/fft.hpp"
#include "goldilocks/ntt_goldilocks.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>

#define LENGTH (2 << 27)

extern "C" uint64_t gl_mmul(uint64_t a, uint64_t b);

static void FFT_BENCHMARK(benchmark::State &state)
{

    uint64_t lenght = 1 << state.range(0);

    uint64_t *pol;
    pol = (uint64_t *)malloc(lenght * sizeof(uint64_t));
    Goldilocks g(lenght, 16);

    printf("lenght: %lu\n", lenght);
    // Fibonacci
    pol[0] = 0;
    pol[1] = 1;
    for (uint64_t i = 2; i < lenght; i++)
    {
        pol[i] = g.gl_add(pol[i - 1], pol[i - 2]);
    }

    /*
    uint64_t res = gl_mmul(2930257127136708076LL, 3961981891181468770LL);
    uint64_t res_1 = g.gl_mmul2(2930257127136708076LL, 3961981891181468770LL);
    printf("res: %lu res_1: %lu\n", res, res_1);
    */

    for (auto _ : state)
    {
#pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < 100; i++)
        {
            // g.ntt(pol, lenght);
            g.intt(pol, lenght);
        }
    }
    /*
    for (uint64_t i = 0; i < lenght; i++)
    {
        printf("%lu\n", pol[i]);
    }*/
}

BENCHMARK(FFT_BENCHMARK)->DenseRange(23, 23, 1)->Unit(benchmark::kMillisecond)->Iterations(3);
BENCHMARK_MAIN();
