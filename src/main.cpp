#include <stdio.h>
#include <string.h>
#include <cstring>
#include <openssl/md5.h>
#include <sstream>
#include "mpi.h"
#include "goldilocks/ntt_goldilocks.hpp"
#include "poseidon_goldilocks_opt.hpp"

#define NUM_COLUMNS 1
#define RATE 8
#define CAPACITY 4
#define NUM_ROWS (1 << 25)

int main()
{

    int np, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    Goldilocks g(NUM_COLUMNS, 8);

    uint64_t hash_input_size = (RATE + CAPACITY);
    uint64_t *fibonacci = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));
    uint64_t *pol_output = (uint64_t *)malloc(hash_input_size * sizeof(uint64_t));

    // Fibonacci
    fibonacci[0] = 0;
    fibonacci[1] = 1;
    for (uint64_t i = 2; i < NUM_COLUMNS * (RATE + CAPACITY); i++)
    {
        fibonacci[i] = g.gl_add(fibonacci[i - 1], fibonacci[i - 2]);
    }

    uint64_t pol_input_t[12];
    MPI_Barrier(MPI_COMM_WORLD);
    double st = omp_get_wtime();
    for (uint64_t i = 0; i < 100000; i++)
    {

        // std::memcpy(pol_input_t, fibonacci, (RATE + CAPACITY) * sizeof(uint64_t));
        Poseidon_goldilocks_opt::hash(pol_input_t);
        // std::memcpy(pol_output, &pol_input_t[0], (RATE + CAPACITY) * sizeof(uint64_t));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double runtime = omp_get_wtime() - st;

    printf("rank: %d,np: %d, fibonacci[0]: %lu runtime: %fus\n", rank, np, fibonacci[0], runtime * 1000000 / 100000);
    MPI_Finalize();
    free(fibonacci);
    return 0;
}