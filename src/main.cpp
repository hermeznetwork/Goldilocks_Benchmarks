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


#define SIZE 8388608

int main(int argc, char **argv ){


    uint64_t columns   = atoi(argv[1]);
    uint64_t blocks    = atoi(argv[2]);
    int reps           = atoi(argv[3]);
    uint64_t nThreads  = atoi(argv[4]);
    uint64_t colsBlock = columns/blocks;

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
#endif

    omp_set_num_threads(nThreads);

    printf("Arguments: columns=%lu, blocks=%lu,reps=%d, nThreads=%lu, colsBlock=%lu\n",columns, blocks, reps, nThreads,colsBlock);

    #pragma omp parallel
    {
	#pragma omp master
	{
           printf("Using %d threads\n",omp_get_num_threads());    
	}
#ifdef LIKWID_PERFMON
	LIKWID_MARKER_REGISTER("NONVECT");
	LIKWID_MARKER_REGISTER("VECT");
    LIKWID_MARKER_REGISTER("MEM");
#endif
    }

    uint64_t length = SIZE;
    printf("Columns per proc=%lu(-1) reps=%d\n",columns,reps);
    Goldilocks g1(length, nThreads);
    Goldilocks g2(length, nThreads);
    
    uint64_t *pol_ext_1 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);
    uint64_t *pol_ext_2 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);
    uint64_t *pol_ext_3 = (uint64_t *)malloc(length * sizeof(uint64_t) *columns);

    
    // Fibonacci

#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_START("MEM");
    }
#endif
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
#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_STOP("MEM");
    }
#endif

    // INITIALIZATION VECTORIZED
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
#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_START("NONVECT");
    }
#endif
    double st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {
#pragma omp parallel for    
        for (u_int64_t i = 0; i < columns; i++)
        {
            int offset = i*length;
            g1.ntt(pol_ext_1+offset, length);
        }
    }
    double runtime1 = omp_get_wtime() - st;
#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_STOP("NONVECT");
    }
#endif
    printf("NON-VECTORIZED // Columns: %lu, time:%f\n", columns, runtime1 / reps);
   
    //
    // VECTORIZED EXECUTION
#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_START("VECT");
    }
#endif
    st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {   
	 int i = 0;
	 int offset = i*(columns/blocks)*length;
         g2.ntt_block(&pol_ext_2[offset], length, colsBlock);
    }
    double runtime2 = omp_get_wtime() - st;
#ifdef LIKWID_PERFMON
    #pragma omp parallel
    {
       LIKWID_MARKER_STOP("VECT");
    }
#endif
    printf("VECTORIZED     // Columns: %lu, time:%f\n", columns, runtime2 / reps);
    

    //
    // VECTORIZED EXECUTION 2
    //
    st = omp_get_wtime();  
    for (int k=0; k<reps; ++k)
    {   
	 //#pragma omp task 
	 /*{
	    int i = 0;
	    int offset = i*(columns/blocks)*length;
            g2.ntt_block(&pol_ext_3[offset], length, colsBlock,aux2);
	 }
	 //#pragma omp task 
	 {
	    int i = 1;
	    int offset = i*(columns/blocks)*length;
            g2.ntt_block(&pol_ext_3[offset], length, colsBlock,aux2);
	 }*/
         //#pragma omp taskwait
    }
    double runtime3 = omp_get_wtime() - st;
    printf("VECTORIZED 2   // Columns: %lu, time:%f\n", columns, runtime3 / reps);
    
    //
    // CHECK RESULTS
    //
    for(uint64_t k=0; k<blocks; k++){
       uint64_t offset0 = colsBlock*length*k;	    
 #pragma omp parallel for
       for(uint64_t j=0; j<length; ++j){
	  uint64_t offset1 = offset0+j*colsBlock;
	  for (uint64_t i = 0; i < colsBlock; i++)
	  {
	     uint64_t offset2 = (i+k*colsBlock)*length;
	     assert(pol_ext_2[offset1+i]==(pol_ext_1[j+offset2]));
	     //assert(pol_ext_3[offset1+i]==(pol_ext_1[j+offset2]));
	 }
       }
    }

    u_int64_t cont = reps*columns;
    printf("Ffts done: %lu, length: %lu speedups: %f %f\n", cont,length,runtime1/runtime2,runtime1/runtime3);
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
    
    free(pol_ext_1);
    free(pol_ext_2);

    return 0;
}

