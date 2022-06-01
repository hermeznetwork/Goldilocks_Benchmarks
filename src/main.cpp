#include <benchmark/benchmark.h>
#include "goldilocks/goldilocks.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <math.h>
//#include "mpi.h"

#define NUM_HASHES 100
#define SIZE 8388608
#define REPS 1
#define MPI 0


int main(int argc, char **argv ){

    int np=1, rank=0;
    int nThreads;

#if MPI == 1
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
#endif
    #pragma omp parallel
    {
	#pragma omp master
	{
		if(rank==0) printf("Using %d threads\n",omp_get_num_threads());    
		nThreads = omp_get_num_threads();
	}

    }

    uint64_t nhashes = atoi(argv[1]);
    uint64_t length = SIZE;
    uint64_t extensionLength = SIZE*2;
    uint64_t columns=nhashes/np;
    if(rank < (int)(nhashes%np)) columns+=1;
    if(rank==0) printf("Columns per proc=%lu(-1) reps=%d rank=%d\n",columns,REPS,rank);
    nThreads = 1;

    Goldilocks g(length, nThreads);
    Goldilocks ge(extensionLength, nThreads);
    
    uint64_t *pol = (uint64_t *)malloc(length * sizeof(uint64_t));
    uint64_t *pol_ext = (uint64_t *)malloc(extensionLength * sizeof(uint64_t));
    uint64_t *pol_ext_1 = (uint64_t *)malloc(extensionLength * sizeof(uint64_t) *columns);
    uint64_t *pol_ext_2 = (uint64_t *)malloc(extensionLength * sizeof(uint64_t) *columns);

    
    // Fibonacci

    pol[0] = 0;
    pol[1] = 1;
    for (uint64_t i = 2; i < length; i++)
    {
	    pol[i] = g.gl_add(pol[i - 1], pol[i - 2]);
    }
    std::memcpy(pol_ext, pol, length * sizeof(uint64_t));

    for(uint64_t k=0; k< columns; k++){ 
       uint64_t offset = k*extensionLength;;
       std::memcpy(pol_ext_1+offset, pol, length * sizeof(uint64_t));
    }
    for(uint64_t i=0; i<length; ++i){
       uint64_t offset = i*columns;
       for(uint64_t k=0; k< columns; k++){ 
          pol_ext_2[offset+k]=pol[i];
       }
    }

    //uint64_t r = 1;
    //uint64_t shift = 49;

    //
    // Initial casen
    //
    for (int k=0; k<REPS; ++k)
    {
    	g.ntt(pol_ext, length); //ojo que estic fent la ntt no la intt
    }
    /*for (uint j = 0; j < length; j++)
    {
        pol_ext[j] = g.gl_mmul2(pol_ext[j], r);
        r = g.gl_mmul2(r, shift);
    }
    ge.ntt(pol_ext, extensionLength);*/
        

    //
    // NON-vectorized case
    //
    double st = omp_get_wtime();  
    for (int k=0; k<REPS; ++k)
    {
        for (u_int64_t i = 0; i < columns; i++)
        {
            int offset = i*extensionLength;;
            g.ntt(pol_ext_1+offset, length);
            /*for (uint j = 0; j < length; j++)
            {
                pol_ext[j] = g.gl_mmul2(pol_ext[j], r);
                r = g.gl_mmul2(r, shift);
            }
            ge.ntt(pol_ext, extensionLength);*/
        }
    }
    double runtime = omp_get_wtime() - st;
    if(rank==0){
	printf("NON-VECTORIZED// Np: %d, columns: %lu, time:%f\n", np, columns, runtime / REPS);
    }
   

    //
    // VECTORIZED case
    //
    st = omp_get_wtime();  
    for (int k=0; k<REPS; ++k)
    {
        g.ntt_block(pol_ext_2, length, columns);
        /*for (int i = 0; i < columns; i++)
        {
            int offset = i*extensionLength;;
            for (uint j = 0; j < length; j++)
            {
                pol_ext[j] = g.gl_mmul2(pol_ext[j], r);
                r = g.gl_mmul2(r, shift);
            }

            ge.ntt(pol_ext, extensionLength);
        }*/
    }
    runtime = omp_get_wtime() - st;
    if(rank==0){
        printf("VECTORIZED    // Np: %d, columns: %lu, time:%f\n", np, columns, runtime / REPS);
    }
    //
    // CHECK RESULTS
    for (u_int64_t i = 0; i < columns; i++)
    {
        u_int64_t offset = i*extensionLength;;
        for(uint64_t j =0; j<length; ++j){
	   assert(pol_ext_1[offset+j]==pol_ext[j]);
	}
    }

    for(uint64_t j=0; j<length; ++j){
	uint64_t offset = j*columns;
	for (uint64_t i = 0; i < columns; i++)
	{
	    if(pol_ext_2[offset+i] != pol_ext[j]){
	      
	       printf("hola: %lu %lu\n",pol_ext_2[offset+i],pol_ext[j]);
	       assert(pol_ext_2[offset+i]==pol_ext[j]);
	    }
	}
    }
   
    //
    // COUNTING rows carried out
    //
    u_int64_t cont = REPS*columns;
    u_int64_t contg=0;
    contg = cont;
    printf("Np: %d, ffts done: %lu\n", np, contg);

    free(pol);
    free(pol_ext);
    free(pol_ext_1);
    free(pol_ext_2);
}

