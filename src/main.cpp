#include <benchmark/benchmark.h>
#include "ffiasm/fr.hpp"
#include "ffiasm/fft.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>

#define LENGTH (2 << 27)

static void FFT_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    uint64_t lenght = 2 << state.range(0);
    FFT<RawFr> fft(lenght);
    RawFr::Element *pol;

    pol = (RawFr::Element *)malloc(lenght * sizeof(RawFr::Element));

    // Fibonacci
    field.fromString(pol[0], "0");
    field.fromString(pol[1], "1");
    for (uint64_t i = 2; i < lenght; i++)
    {
        field.add(pol[i], pol[i - 1], pol[i - 2]);
    }

    for (auto _ : state)
    {
        fft.fft(pol, lenght);
//        fft.ifft(pol, lenght);
    }

/*
    for (uint64_t i = 0; i < 10; i++) {
        printf("%s ", field.toString(pol[i]).c_str());
    }
*/

    free(pol);
}

BENCHMARK(FFT_BENCHMARK)->DenseRange(23, 27, 1)->Unit(benchmark::kMillisecond)->Iterations(5);
BENCHMARK_MAIN();

