#include <benchmark/benchmark.h>
#include "ffiasm/fr.hpp"
#include "ffiasm/fft.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>

#define NUM_ROWS (1 << 24)
#define NUM_COLS 100

static void FFT_BENCHMARK(benchmark::State &state)
{
    RawFr field;
    FFT<RawFr> fft(NUM_ROWS);
    RawFr::Element *pol;
    pol = (RawFr::Element *)malloc((uint64_t)NUM_COLS * (uint64_t)NUM_ROWS * sizeof(RawFr::Element));

#pragma omp parallel for
    for (uint64_t j = 0; j < NUM_COLS; j++)
    {
        uint64_t offset = j * (uint64_t)NUM_ROWS;
	field.fromString(pol[offset], std::to_string((int)(1 + j)));
        field.fromString(pol[offset + 1], std::to_string((int)(1 + j)));
        for (uint64_t i = 2; i < NUM_ROWS; i++)
        {
            field.add(pol[i + offset],pol[(i - 2) + offset], pol[(i - 1) + offset]);
        }
    }
    for (auto _ : state)
    {
#pragma omp parallel for
         for (uint64_t i = 0; i < NUM_COLS; i++)
         {
              fft.fft(&pol[i * NUM_ROWS], NUM_ROWS);
//              fft.ifft(&pol[i * NUM_ROWS], NUM_ROWS);
         }
     }
/*   
 for (uint64_t i = 0; i < 10; i++) {
        printf("%s ", field.toString(pol[i]).c_str());
    }
*/  
  free(pol);
}

BENCHMARK(FFT_BENCHMARK)->DenseRange(11, 1, 1)->Unit(benchmark::kMillisecond)->Iterations(5);
BENCHMARK_MAIN();

