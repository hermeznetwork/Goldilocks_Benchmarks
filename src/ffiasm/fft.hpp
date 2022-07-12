#ifndef FFT_H
#define FFT_H
#define NUM_PHASES_ 2
#define NUM_BLOCKS_ 1

template <typename Field>
class FFT
{
    Field f;
    typedef typename Field::Element Element;

    u_int32_t s;
    Element nqr;
    Element *roots;
    Element *powTwoInv;
    u_int32_t nThreads;

    void reversePermutationInnerLoop(Element *a, u_int64_t from, u_int64_t to, u_int32_t domainPow);
    void reversePermutation(Element *dst, Element *a, u_int64_t n);
    void reversePermutation_inplace(Element *a, u_int64_t n);
    void reversePermutation_block(Element *dst, Element *src, u_int64_t n, u_int64_t offset_cols, u_int64_t ncols, u_int64_t ncols_all);
    void fftInnerLoop(Element *a, u_int64_t from, u_int64_t to, u_int32_t s);
    void finalInverseInner(Element *a, u_int64_t from, u_int64_t to, u_int32_t domainPow);
    void shuffle_old(Element *dst, Element *a, u_int64_t n, u_int64_t s);
    void shuffle(Element *dst, Element *src, u_int64_t n, u_int64_t s);
    void traspose(
        Element *dst,
        Element *src,
        uint64_t srcRowSize,
        uint64_t srcX,
        uint64_t srcWidth,
        uint64_t srcY,
        uint64_t srcHeight,
        uint64_t dstRowSize,
        uint64_t dstX,
        uint64_t dstY);
    void fft_block_iters(Element *dst, Element *src, u_int64_t n, u_int64_t offset_cols, u_int64_t ncols, u_int64_t ncols_all, u_int64_t nphase, Element *aux);

public:
    FFT(u_int64_t maxDomainSize, u_int32_t _nThreads = 0);
    ~FFT();
    void fft(Element *a, u_int64_t n);
    void fft2(Element *a, u_int64_t n);
    void ifft(Element *a, u_int64_t n);
    void fft_block(Element *dst, Element *src, u_int64_t n, u_int64_t ncols, u_int64_t nphase = NUM_PHASES_, u_int64_t nblock = NUM_BLOCKS_);
    void ifft_block(Element *dst, Element *src, u_int64_t n, u_int64_t ncols, u_int64_t nphase = NUM_PHASES_, u_int64_t nblock = NUM_BLOCKS_);
    template <int N_, int NCOLS_, int NBLOCK_, int NCOLS_BLOCK_, int NPHASE_>
    void fft_block_divisible(Element *dst, Element *src);
    template <int N_, int NCOLS_BLOCK_, int NCOLS_, int NPHASE_>
    void fft_block_iters_divisible(Element *dst, Element *src, u_int64_t offset_cols, Element *aux);
    template <int N_, int NCOLS_BLOCK_, int NCOLS_>
    void reversePermutation_block_divisible(Element *dst, Element *src, u_int64_t offset_cols);

    u_int32_t log2(u_int64_t n);
    inline Element &root(u_int32_t domainPow, u_int64_t idx) { return roots[idx << (s - domainPow)]; }

    void printVector(Element *a, u_int64_t n);
};

#include "fft.c.hpp"

#endif // FFT_H