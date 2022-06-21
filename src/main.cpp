#include <benchmark/benchmark.h>
#include "goldilocks/goldilocks.hpp"
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <math.h>
#include <random>
#include <iostream>

#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#endif


#define SIZE 1<<4 

int main(int argc, char **argv ){


    uint64_t columns   = atoi(argv[1]);
    uint64_t blocks    = atoi(argv[2]);
    int reps           = atoi(argv[3]);
    uint64_t nThreads  = atoi(argv[4]);
    uint64_t nphase    = atoi(argv[5]);
    if(blocks > columns) blocks=columns;
    if(blocks < 1) blocks=1;
    uint64_t colsBlock = columns/blocks;

    omp_set_num_threads(nThreads);

    printf("Arguments: columns=%lu, blocks=%lu,reps=%d, nThreads=%lu, colsBlock=%lu, nphase=%lu, SIZE=%d\n",columns, blocks, reps, nThreads,colsBlock,nphase,SIZE);

    #pragma omp parallel
    {
	#pragma omp master
	{
           printf("Using %d threads\n",omp_get_num_threads());    
	}
    }

    uint64_t length = SIZE;
    printf("Columns per proc=%lu(-1) reps=%d\n",columns,reps);
    Goldilocks g1(length, nThreads);
    
    uint64_t *pol_ext_1 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);
    uint64_t *pol_ext_2 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);
    uint64_t *pol_ext_3 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);


    // INITIALIZATION NON-VECTORIZED
    #pragma omp parallel for
    for(uint64_t k=0; k< columns; k++){ 
	    uint64_t offset = k*length;;
	    pol_ext_1[offset] = 1+k;
	    pol_ext_1[offset+1] = 1+k;
	    for (uint64_t i = 2; i < length; i++)
	    {
	        pol_ext_1[offset+i] = g1.gl_add(pol_ext_1[offset+i-1], pol_ext_1[offset+i-2]);
	    }
    }
    // INITIALIZATION BLOCK VERSION
    for(uint64_t k=0; k<blocks; k++){
        uint64_t offset0 = colsBlock*length*k;	    
 #pragma omp parallel for
        for(uint64_t j=0; j<length; ++j){
	        uint64_t offset1 = offset0+j*colsBlock;
	        for (uint64_t i = 0; i < colsBlock; i++)
	        {
	            uint64_t offset2 = (i+k*colsBlock)*length;
	            pol_ext_2[offset1+i]=(pol_ext_1[j+offset2]);
	            pol_ext_3[offset1+i]=(pol_ext_1[j+offset2]);
	        }
       }
    }

    //
    // NON-VECTORIZED EXECUTION
    //
    double st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {
        #pragma omp parallel for    
        for (u_int64_t i = 0; i < columns; i++)
        {
            u_int64_t offset = i*length;
            g1.ntt(pol_ext_1+offset, length, nphase);
            g1.intt(pol_ext_1+offset, length,nphase);
        }
        //for (uint64_t i = 0; i < length; i++) pol_ext_1[i] = pol_ext_1[i] % GOLDILOCKS_PRIME;
        for (uint64_t i = 0; i < length; i++){ 
            printf("hola: %d %lu %lu %lu \n", k, pol_ext_1[i], pol_ext_1[i] % GOLDILOCKS_PRIME, pol_ext_1[i]/GOLDILOCKS_PRIME);
        }

    }
    double runtime1 = omp_get_wtime() - st;

    printf("NON-VECTORIZED // Columns: %lu, time:%f\n", columns, runtime1 / reps);
   
    //
    // VECTORIZED EXECUTION
    //
    st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {   
	    for(uint64_t i=0; i<blocks; ++i){
	        u_int64_t offset = i*(columns/blocks)*length;
            g1.ntt_block(&pol_ext_2[offset], length, colsBlock,nphase);
            g1.intt_block(&pol_ext_2[offset], length, colsBlock,nphase);
        }
    }
    double runtime2 = omp_get_wtime() - st;
    printf("VECTORIZED     // Columns: %lu, time:%f\n", columns, runtime2 / reps);
    
    //
    // VECTORIZED EXECUTION 2
    //
    st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {   
	    for(uint64_t i=0; i<blocks; ++i){
	        uint64_t offset = i*(columns/blocks)*length;
            g1.ntt_block_2(&pol_ext_3[offset], length, colsBlock,nphase);
            g1.intt_block_2(&pol_ext_3[offset], length, colsBlock,nphase);
        }
    }
    double runtime3 = omp_get_wtime() - st;
    printf("VECTORIZED 2   // Columns: %lu, time:%f\n", columns, runtime3 / reps);

    //
    // CHECK RESULTS NON-VECT
    //
    #pragma omp parallel for
    for(uint64_t k=0; k< columns; k++){ 
	    uint64_t offset = k*length;
	    assert(pol_ext_1[offset] % GOLDILOCKS_PRIME == 1+k);
        //printf("hola %lu %lu %lu \n", pol_ext_1[offset+1]% GOLDILOCKS_PRIME,1+k, GOLDILOCKS_PRIME);
	    assert(pol_ext_1[offset+1]% GOLDILOCKS_PRIME == 1+k);
	    for (uint64_t i = 2; i < length; i++)
	    {
	        assert(pol_ext_1[offset+i] % GOLDILOCKS_PRIME == g1.gl_add(pol_ext_1[offset+i-1] % GOLDILOCKS_PRIME, pol_ext_1[offset+i-2] % GOLDILOCKS_PRIME));
	    }
    }

    //
    // CHECK RESULTS VECT
    //
    for(uint64_t k=0; k<blocks; k++){
        uint64_t offset0 = colsBlock*length*k;	    
        #pragma omp parallel for
        for(uint64_t j=0; j<length; ++j){
	        uint64_t offset1 = offset0+j*colsBlock;
	        for (uint64_t i = 0; i < colsBlock; i++)
	        {
	            uint64_t offset2 = (i+k*colsBlock)*length;
	            assert(pol_ext_2[offset1+i] % GOLDILOCKS_PRIME  == pol_ext_1[j+offset2] % GOLDILOCKS_PRIME);
                assert(pol_ext_3[offset1+i] % GOLDILOCKS_PRIME  == pol_ext_1[j+offset2]% GOLDILOCKS_PRIME);
	        }
       }
    }
    
    free(pol_ext_1);
    free(pol_ext_2);
    free(pol_ext_3);

    return 0;
}

