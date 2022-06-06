#include <benchmark/benchmark.h>
#include "goldilocks/goldilocks.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <math.h> /* ceil */

#define NUM_COLS 32
#define FFT_SIZE (1 << 23)

int main()
{
    Goldilocks g(FFT_SIZE, 1);
    uint64_t *pol = (uint64_t *)malloc((uint64_t)NUM_COLS * (uint64_t)FFT_SIZE * sizeof(uint64_t));

    for (uint64_t i = 0; i < NUM_COLS; i++)
    {
        pol[i] = i + 1;
        pol[i + NUM_COLS] = i + 1;
    }

    for (uint64_t i = 2; i < FFT_SIZE; i++)
    {
        uint32_t offset = i * NUM_COLS;
        for (uint64_t k = 0; k < NUM_COLS; k++)
        {
            pol[offset + k] = Goldilocks::gl_add(pol[offset + k - NUM_COLS], pol[offset + k - 2 * NUM_COLS]);
        }
    }

    /*
    for (uint64_t j = 0; j < FFT_SIZE; j++)
    {
        for (uint64_t i = 0; i < NUM_COLS; i++)
        {
            printf("[%lu]: %lu ", j * NUM_COLS + i, pol[j * NUM_COLS + i]);
        }
        printf("\n");
    }*/
    double st = omp_get_wtime();
    g.ntt_block(pol, FFT_SIZE, NUM_COLS);
    double ft = omp_get_wtime();

    printf("runtime: %f\n", ft - st);
}