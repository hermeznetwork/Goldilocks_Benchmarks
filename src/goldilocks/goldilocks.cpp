#include "goldilocks.hpp"

const uint64_t Goldilocks::Q = 0xFFFFFFFF00000001LL;
const uint64_t Goldilocks::MM = 0xFFFFFFFeFFFFFFFFLL;
const uint64_t Goldilocks::CQ = 0x00000000FFFFFFFFLL;
const uint64_t Goldilocks::R2 = 0xFFFFFFFe00000001LL;

static inline u_int64_t BR(u_int64_t x, u_int64_t domainPow)
{
    x = (x >> 16) | (x << 16);
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    return (((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1)) >> (32 - domainPow);
}

void Goldilocks::ntt(u_int64_t *_a, u_int64_t n)
{
    u_int64_t *aux_a = new u_int64_t[n];
    u_int64_t *a = _a;
    u_int64_t *a2 = aux_a;
    u_int64_t *tmp;
    reversePermutation(a2, a, n);
    tmp = a2;
    a2 = a;
    a = tmp;

    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    u_int64_t maxBatchPow = s / 4;

    u_int64_t batchSize = 1 << maxBatchPow;
    u_int64_t nBatches = n / batchSize;

    for (u_int64_t s = 1; s <= domainPow; s += maxBatchPow)
    {
        u_int64_t sInc = s + maxBatchPow <= domainPow ? maxBatchPow : domainPow - s + 1;
        omp_set_dynamic(0);
        omp_set_num_threads(nThreads);
#pragma omp parallel for
        for (u_int64_t b = 0; b < nBatches; b++)
        {
            u_int64_t rs = s - 1;
            uint64_t re = domainPow - 1;
            uint64_t rb = 1 << rs;
            uint64_t rm = (1 << (re - rs)) - 1;

            for (u_int64_t si = 0; si < sInc; si++)
            {
                u_int64_t m = 1 << (s + si);
                u_int64_t mdiv2 = m >> 1;
                u_int64_t mdiv2i = 1 << si;
                u_int64_t mi = mdiv2i * 2;
                u_int64_t j_tmp = (b * batchSize / 2);
                j_tmp = (j_tmp & rm) * rb + (j_tmp >> (re - rs));
                j_tmp = j_tmp % mdiv2;
                u_int64_t w = root(s + si, j_tmp);
                u_int64_t wi = root(s + si, 1);
                u_int64_t mask = mdiv2i - 1;
                u_int64_t k1 = b * batchSize;
                u_int64_t k2 = mi / mdiv2i;

                for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                {
                    u_int64_t t;
                    u_int64_t u;
                    u_int64_t ki = k1 + k2 * i;
                    u_int64_t ji = i & mask;

                    t = gl_mmul2(w, a[ki + ji + mdiv2i]);
                    u = a[ki + ji];
                    a[ki + ji] = gl_add(t, u);
                    a[ki + ji + mdiv2i] = gl_sub(u, t);
                    w = gl_mmul2(w, wi);
                }
            }
        }

        shuffle(a2, a, n, sInc);
        tmp = a2;
        a2 = a;
        a = tmp;
    }
    if (a != _a)
    {
        // printf("baaaaad!\n");
        shuffle(a, a2, n, 0);
    }
    delete aux_a;
}

void Goldilocks::ntt_block(u_int64_t *_a, u_int64_t n, u_int64_t ncols)
{
    u_int64_t *aux_a = new u_int64_t[n*ncols];
    u_int64_t *a = _a;
    u_int64_t *a2 = aux_a;
    u_int64_t *tmp;
    reversePermutation_block(a2, a, n,ncols);
    tmp = a2;
    a2 = a;
    a = tmp;

    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    u_int64_t maxBatchPow = s / 4;

    u_int64_t batchSize = 1 << maxBatchPow;
    u_int64_t nBatches = n / batchSize;

    for (u_int64_t s = 1; s <= domainPow; s += maxBatchPow)
    {
        u_int64_t sInc = s + maxBatchPow <= domainPow ? maxBatchPow : domainPow - s + 1;
        omp_set_dynamic(0);
        omp_set_num_threads(nThreads);
#pragma omp parallel for
        for (u_int64_t b = 0; b < nBatches; b++)
        {
            u_int64_t rs = s - 1;
            uint64_t re = domainPow - 1;
            uint64_t rb = 1 << rs;
            uint64_t rm = (1 << (re - rs)) - 1;

            for (u_int64_t si = 0; si < sInc; si++)
            {
                u_int64_t m = 1 << (s + si);
                u_int64_t mdiv2 = m >> 1;
                u_int64_t mdiv2i = 1 << si;
                u_int64_t mi = mdiv2i * 2;
                u_int64_t j_tmp = (b * batchSize / 2);
                j_tmp = (j_tmp & rm) * rb + (j_tmp >> (re - rs));
                j_tmp = j_tmp % mdiv2;
                u_int64_t w = root(s + si, j_tmp);
                u_int64_t wi = root(s + si, 1);
                u_int64_t mask = mdiv2i - 1;
                u_int64_t k1 = b * batchSize;
                u_int64_t k2 = mi / mdiv2i;

                for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                {
                    
                    u_int64_t ki = k1 + k2 * i;
                    u_int64_t ji = i & mask;
                    u_int64_t offset1=(ki + ji + mdiv2i)*ncols;
		            u_int64_t offset2=(ki + ji)*ncols;

                    for(u_int64_t k=0; k<ncols; ++k){ //rick: vectorize
                        u_int64_t t = gl_mmul2(w, a[offset1+k]); 
                        u_int64_t u = a[offset2+k];
                        a[offset2+k] = gl_add(t, u);
                        a[offset1+k] = gl_sub(u, t);
                    }
                    w = gl_mmul2(w, wi);

                }
            }
        }

        shuffle_block(a2, a, n, sInc,ncols);
        tmp = a2;
        a2 = a;
        a = tmp;
    }
    if (a != _a)   //rick: this applyies for al the ntt?
    {
        // printf("baaaaad!\n");
        shuffle_block(a, a2, n, 0,ncols);
    }
    delete aux_a;
}

void Goldilocks::intt(u_int64_t *a, u_int64_t n)
{
    ntt(a, n);
    u_int64_t domainPow = log2(n);
    u_int64_t nDiv2 = n >> 1;
#pragma omp parallel for
    for (u_int64_t i = 1; i < nDiv2; i++)
    {
        u_int64_t tmp;
        u_int64_t r = n - i;
        tmp = a[i];
        a[i] = gl_mmul2(a[r], powTwoInv[domainPow]);
        a[r] = gl_mmul2(tmp, powTwoInv[domainPow]);
    }
    a[0] = gl_mmul2(a[0], powTwoInv[domainPow]);
    a[n >> 1] = gl_mmul2(a[n >> 1], powTwoInv[domainPow]);
}

void Goldilocks::shuffle(u_int64_t *dst, u_int64_t *src, uint64_t n, uint64_t s)
{
    uint64_t srcRowSize = 1 << s;

    uint64_t srcX = 0;
    uint64_t srcWidth = 1 << s;
    uint64_t srcY = 0;
    uint64_t srcHeight = n / srcRowSize;

    uint64_t dstRowSize = n / srcRowSize;
    uint64_t dstX = 0;
    uint64_t dstY = 0;

#pragma omp parallel
#pragma omp single
    traspose(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight, dstRowSize, dstX, dstY);
#pragma omp taskwait
}

void Goldilocks::shuffle_block(u_int64_t *dst, u_int64_t *src, uint64_t n, uint64_t s, uint64_t ncols)
{
    uint64_t srcRowSize = 1 << s;

    uint64_t srcX = 0;
    uint64_t srcWidth = 1 << s;
    uint64_t srcY = 0;
    uint64_t srcHeight = n / srcRowSize;

    uint64_t dstRowSize = n / srcRowSize;
    uint64_t dstX = 0;
    uint64_t dstY = 0;

#pragma omp parallel
#pragma omp single
    traspose_block(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight, dstRowSize, dstX, dstY,ncols);
#pragma omp taskwait
}

void Goldilocks::traspose(
    u_int64_t *dst,
    u_int64_t *src,
    uint64_t srcRowSize,
    uint64_t srcX,
    uint64_t srcWidth,
    uint64_t srcY,
    uint64_t srcHeight,
    uint64_t dstRowSize,
    uint64_t dstX,
    uint64_t dstY)
{
    if ((srcWidth == 1) || (srcHeight == 1) || (srcWidth * srcHeight < CACHESIZE))
    {
#pragma omp task
        {
            for (uint64_t x = 0; x < srcWidth; x++)
            {
                for (uint64_t y = 0; y < srcHeight; y++)
                {
                    dst[(dstY + +x) * dstRowSize + (dstX + y)] = src[(srcY + +y) * srcRowSize + (srcX + x)];
                }
            }
        }
        return;
    }
    if (srcWidth > srcHeight)
    {
        traspose(dst, src, srcRowSize, srcX, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY);
        traspose(dst, src, srcRowSize, srcX + srcWidth / 2, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY + srcWidth / 2);
    }
    else
    {
        traspose(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight / 2, dstRowSize, dstX, dstY);
        traspose(dst, src, srcRowSize, srcX, srcWidth, srcY + srcHeight / 2, srcHeight / 2, dstRowSize, dstX + srcHeight / 2, dstY);
    }
}

void Goldilocks::traspose_block(
    u_int64_t *dst,
    u_int64_t *src,
    uint64_t srcRowSize,
    uint64_t srcX,
    uint64_t srcWidth,
    uint64_t srcY,
    uint64_t srcHeight,
    uint64_t dstRowSize,
    uint64_t dstX,
    uint64_t dstY,
    uint64_t ncols)
{
    if ((srcWidth == 1) || (srcHeight == 1) || (srcWidth * srcHeight < CACHESIZE)) //rick: consider this
    {
#pragma omp task
        {
            for (uint64_t x = 0; x < srcWidth; x++)
            {
                for (uint64_t y = 0; y < srcHeight; y++)
                {
                    uint64_t offset_dstY = ((dstY + +x) * dstRowSize + (dstX + y))*ncols;
		            uint64_t offset_src  = ((srcY + +y) * srcRowSize + (srcX + x))*ncols;
		            for(uint64_t k=0; k<ncols; ++k){
                        	dst[offset_dstY+k] = src[offset_src+k];
		            }
                }
            }
        }
        return;
    }
    if (srcWidth > srcHeight)
    {
        traspose_block(dst, src, srcRowSize, srcX, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY,ncols);
        traspose_block(dst, src, srcRowSize, srcX + srcWidth / 2, srcWidth / 2, srcY, srcHeight, dstRowSize, dstX, dstY + srcWidth / 2,ncols);
    }
    else
    {
        traspose_block(dst, src, srcRowSize, srcX, srcWidth, srcY, srcHeight / 2, dstRowSize, dstX, dstY,ncols);
        traspose_block(dst, src, srcRowSize, srcX, srcWidth, srcY + srcHeight / 2, srcHeight / 2, dstRowSize, dstX + srcHeight / 2, dstY,ncols);
    }
}

void Goldilocks::reversePermutation(u_int64_t *dst, u_int64_t *a, u_int64_t n)
{
    uint32_t domainSize = log2(n);
#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        r = BR(i, domainSize);
        dst[i] = a[r];
    }
}

void Goldilocks::reversePermutation_block(u_int64_t *dst, u_int64_t *a, u_int64_t n, u_int64_t ncols)
{
    uint32_t domainSize = log2(n);
#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        r = BR(i, domainSize);
        u_int64_t offset_i = i*ncols; 
	    u_int64_t offset_r = r*ncols;
	    for(u_int64_t k=0; k< ncols; ++k){
           dst[offset_i+k] = a[offset_r+k];
	    }
    }
}

u_int32_t Goldilocks::log2(u_int64_t n)
{
    assert(n != 0);
    u_int32_t res = 0;
    while (n != 1)
    {
        n >>= 1;
        res++;
    }
    return res;
}
