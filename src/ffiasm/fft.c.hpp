#include <thread>
#include <vector>
#include <omp.h>
#include <cstring>
#include <iostream>
using namespace std;

// The function we want to execute on the new thread.

template <typename Field>
u_int32_t FFT<Field>::log2(u_int64_t n)
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

static inline u_int64_t BR(u_int64_t x, u_int64_t domainPow)
{
    x = (x >> 16) | (x << 16);
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    return (((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1)) >> (32 - domainPow);
}

template <typename Field>
FFT<Field>::FFT(u_int64_t maxDomainSize, uint32_t _nThreads)
{
    nThreads = _nThreads == 0 ? omp_get_max_threads() : _nThreads;
    f = Field::field;
    omp_set_num_threads(nThreads);

    u_int32_t domainPow = log2(maxDomainSize);

    mpz_t m_qm1d2;
    mpz_t m_q;
    mpz_t m_nqr;
    mpz_t m_aux;
    mpz_init(m_qm1d2);
    mpz_init(m_q);
    mpz_init(m_nqr);
    mpz_init(m_aux);

    f.toMpz(m_aux, f.negOne());

    mpz_add_ui(m_q, m_aux, 1);
    mpz_fdiv_q_2exp(m_qm1d2, m_aux, 1);

    mpz_set_ui(m_nqr, 2);
    mpz_powm(m_aux, m_nqr, m_qm1d2, m_q);
    while (mpz_cmp_ui(m_aux, 1) == 0)
    {
        mpz_add_ui(m_nqr, m_nqr, 1);
        mpz_powm(m_aux, m_nqr, m_qm1d2, m_q);
    }

    f.fromMpz(nqr, m_nqr);

    // std::cout << "nqr: " << f.toString(nqr) << std::endl;

    s = 1;
    mpz_set(m_aux, m_qm1d2);
    while ((!mpz_tstbit(m_aux, 0)) && (s < domainPow))
    {
        mpz_fdiv_q_2exp(m_aux, m_aux, 1);
        s++;
    }

    if (s < domainPow)
    {
        throw std::range_error("Domain size too big for the curve");
    }

    uint64_t nRoots = 1LL << s;

    roots = new Element[nRoots];
    powTwoInv = new Element[s + 1];

    f.copy(roots[0], f.one());
    f.copy(powTwoInv[0], f.one());
    if (nRoots > 1)
    {
        mpz_powm(m_aux, m_nqr, m_aux, m_q);
        f.fromMpz(roots[1], m_aux);

        mpz_set_ui(m_aux, 2);
        mpz_invert(m_aux, m_aux, m_q);
        f.fromMpz(powTwoInv[1], m_aux);
    }
#pragma omp parallel
    {
        int idThread = omp_get_thread_num();
        int nThreads = omp_get_num_threads();
        uint64_t increment = nRoots / nThreads;
        uint64_t start = idThread == 0 ? 2 : idThread * increment;
        uint64_t end = idThread == nThreads - 1 ? nRoots : (idThread + 1) * increment;
        if (end > start)
        {
            f.exp(roots[start], roots[1], (uint8_t *)(&start), sizeof(start));
        }
        for (uint64_t i = start + 1; i < end; i++)
        {
            f.mul(roots[i], roots[i - 1], roots[1]);
        }
    }
    Element aux;
    f.mul(aux, roots[nRoots - 1], roots[1]);
    assert(f.eq(aux, f.one()));

    for (uint64_t i = 2; i <= s; i++)
    {
        f.mul(powTwoInv[i], powTwoInv[i - 1], powTwoInv[1]);
    }

    mpz_clear(m_qm1d2);
    mpz_clear(m_q);
    mpz_clear(m_nqr);
    mpz_clear(m_aux);
}

template <typename Field>
FFT<Field>::~FFT()
{
    delete[] roots;
    delete[] powTwoInv;
}

#define PF 5

template <typename Field>
void FFT<Field>::reversePermutation_inplace(Element *a, u_int64_t n)
{
    uint32_t domainSize = log2(n);
#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        Element tmp;
        r = BR(i, domainSize);
        if (i > r)
        {
            f.copy(tmp, a[i]);
            f.copy(a[i], a[r]);
            f.copy(a[r], tmp);
        }
    }
}

template <typename Field>
void FFT<Field>::reversePermutation(Element *dst, Element *a, u_int64_t n)
{
    uint32_t domainSize = log2(n);
#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        r = BR(i, domainSize);
        f.copy(dst[i], a[r]);
    }
}

template <typename Field>
void FFT<Field>::reversePermutation_block(Element *dst, Element *src, u_int64_t n, u_int64_t offset_cols, u_int64_t ncols, u_int64_t ncols_all)
{

    uint32_t domainSize = log2(n);
#pragma omp parallel for schedule(static)
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r = BR(i, domainSize);
        u_int64_t offset_r = r * ncols_all + offset_cols;
        u_int64_t offset_i = i * ncols;
        for (u_int64_t k = 0; k < ncols; ++k)
        {
            f.copy(dst[offset_i + k], src[offset_r + k]);
        }
    }
}

template <typename Field>
void FFT<Field>::fft2(Element *a, u_int64_t n)
{
    omp_set_num_threads(nThreads);
    reversePermutation(a, n);
    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    for (u_int32_t s = 1; s <= domainPow; s++)
    {
        u_int64_t m = 1 << s;
        u_int64_t mdiv2 = m >> 1;
#pragma omp parallel for
        for (u_int64_t i = 0; i < (n >> 1); i++)
        {
            Element t;
            Element u;
            u_int64_t k = (i / mdiv2) * m;
            u_int64_t j = i % mdiv2;

            f.mul(t, root(s, j), a[k + j + mdiv2]);
            f.copy(u, a[k + j]);
            f.add(a[k + j], t, u);
            f.sub(a[k + j + mdiv2], u, t);
        }
    }
}

template <typename Field>
void FFT<Field>::fft(Element *_a, u_int64_t n)
{
    Element *aux_a = new Element[n];
    Element *a = _a;
    Element *a2 = aux_a;
    Element *tmp;

    omp_set_num_threads(nThreads);

    reversePermutation(a2, a, n);
    tmp = a2;
    a2 = a;
    a = tmp;

    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    u_int64_t maxBatchPow = s / 3 + 1;
    // u_int64_t maxBatchPow = s;
    u_int64_t batchSize = 1 << maxBatchPow;
    u_int64_t nBatches = n / batchSize;

    for (u_int64_t s = 1; s <= domainPow; s += maxBatchPow)
    {
        u_int64_t sInc = s + maxBatchPow <= domainPow ? maxBatchPow : domainPow - s + 1;

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
                for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                {
                    Element t;
                    Element u;

                    u_int64_t ki = b * batchSize + (i / mdiv2i) * mi;
                    u_int64_t ji = i % mdiv2i;

                    u_int64_t j = (b * batchSize / 2 + i);
                    j = (j & rm) * rb + (j >> (re - rs));
                    j = j % mdiv2;

                    f.mul(t, root(s + si, j), a[ki + ji + mdiv2i]);
                    f.copy(u, a[ki + ji]);
                    f.add(a[ki + ji], t, u);
                    f.sub(a[ki + ji + mdiv2i], u, t);
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
        shuffle(_a, a, n, 0);
    }
    delete[] aux_a;
}

template <typename Field>
void FFT<Field>::shuffle_old(Element *dst, Element *a, uint64_t n, uint64_t s)
{
    uint64_t e = log2(n);
    uint64_t b = 1 << s;
    uint64_t mask = (1 << (e - s)) - 1;

#pragma omp parallel for
    for (u_int64_t i = 0; i < n; i++)
    {
        u_int64_t r;
        r = (i & mask) * b + (i >> (e - s));
        f.copy(dst[i], a[r]);
    }
}

template <typename Field>
void FFT<Field>::shuffle(Element *dst, Element *src, uint64_t n, uint64_t s)
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

#define CACHESIZE 1 << 18

template <typename Field>
void FFT<Field>::traspose(
    Element *dst,
    Element *src,
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
                    f.copy(dst[(dstY + +x) * dstRowSize + (dstX + y)], src[(srcY + +y) * srcRowSize + (srcX + x)]);
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

template <typename Field>
void FFT<Field>::ifft(Element *a, u_int64_t n)
{
    omp_set_num_threads(nThreads);
    fft(a, n);
    u_int64_t domainPow = log2(n);
    u_int64_t nDiv2 = n >> 1;
#pragma omp parallel for
    for (u_int64_t i = 1; i < nDiv2; i++)
    {
        Element tmp;
        u_int64_t r = n - i;
        f.copy(tmp, a[i]);
        f.mul(a[i], a[r], powTwoInv[domainPow]);
        f.mul(a[r], tmp, powTwoInv[domainPow]);
    }
    f.mul(a[0], a[0], powTwoInv[domainPow]);
    f.mul(a[n >> 1], a[n >> 1], powTwoInv[domainPow]);
}

template <typename Field>
void FFT<Field>::printVector(Element *a, u_int64_t n)
{
    cout << "[" << endl;
    for (u_int64_t i = 0; i < n; i++)
    {
        cout << f.toString(a[i]) << endl;
    }
    cout << "]" << endl;
}

template <typename Field>
void FFT<Field>::fft_block_iters(Element *dst, Element *src, u_int64_t n, u_int64_t offset_cols, u_int64_t ncols, u_int64_t ncols_all, u_int64_t nphase, Element *aux)
{
    Element *dst_;
    if (dst != NULL)
    {
        dst_ = dst;
    }
    else
    {
        dst_ = src;
    }
    Element *a = dst_;
    Element *a2 = aux;
    Element *tmp;

    reversePermutation_block(a2, src, n, offset_cols, ncols, ncols_all);

    tmp = a2;
    a2 = a;
    a = tmp;

    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);
    if (nphase < 1 || domainPow == 0)
    {
        nphase = 1;
    }
    if (nphase > domainPow)
    {
        nphase = domainPow;
    }
    u_int64_t maxBatchPow = s / nphase;
    u_int64_t batchSize = 1 << maxBatchPow;
    u_int64_t nBatches = n / batchSize;

    for (u_int64_t s = 1; s <= domainPow; s += maxBatchPow)
    {
        u_int64_t sInc = s + maxBatchPow <= domainPow ? maxBatchPow : domainPow - s + 1;
#pragma omp parallel for
        for (u_int64_t b = 0; b < nBatches; b++)
        {
            u_int64_t rs = s - 1;
            u_int64_t re = domainPow - 1;
            u_int64_t rb = 1 << rs;
            u_int64_t rm = (1 << (re - rs)) - 1;

            for (u_int64_t si = 0; si < sInc; si++)
            {
                u_int64_t m = 1 << (s + si);
                u_int64_t mdiv2 = m >> 1;
                u_int64_t mdiv2i = 1 << si;
                u_int64_t mi = mdiv2i * 2;
                for (u_int64_t i = 0; i < (batchSize >> 1); i++)
                {
                    u_int64_t ki = b * batchSize + (i / mdiv2i) * mi;
                    u_int64_t ji = i % mdiv2i;

                    u_int64_t offset1 = (ki + ji + mdiv2i) * ncols;
                    u_int64_t offset2 = (ki + ji) * ncols;

                    u_int64_t j = (b * batchSize / 2 + i);
                    j = (j & rm) * rb + (j >> (re - rs));
                    j = j % mdiv2;
                    Element w = root(s + si, j);

                    for (u_int64_t k = 0; k < ncols; ++k)
                    {
                        Element t;
                        Element u;
                        f.mul(t, w, a[offset1 + k]);
                        f.copy(u, a[offset2 + k]);
                        f.add(a[offset2 + k], t, u);
                        f.sub(a[offset1 + k], u, t);
                    }
                }
            }

            u_int64_t srcWidth = 1 << sInc;
            u_int64_t niters = batchSize / srcWidth;
            for (u_int64_t l = 0; l < niters; ++l)
            {
                for (u_int64_t x = 0; x < srcWidth; x++)
                {
                    u_int64_t offset_dstY = (x * (nBatches * niters) + (b * niters + l)) * ncols;
                    u_int64_t offset_src = ((b * niters + l) * srcWidth + x) * ncols;
                    for (u_int64_t k = 0; k < ncols; ++k)
                    {
                        f.copy(a2[offset_dstY + k], a[offset_src + k]);
                    }
                }
            }
        }
        tmp = a2;
        a2 = a;
        a = tmp;
    }
    if (a != dst_)
    {
#pragma omp parallel for schedule(static)
        for (u_int64_t ie = 0; ie < n * ncols; ++ie)
        {
            f.copy(dst_[ie], a[ie]);
        }
    }
}

template <typename Field>
void FFT<Field>::fft_block(Element *dst, Element *src, u_int64_t n, u_int64_t ncols, u_int64_t nphase, u_int64_t nblock)
{
    omp_set_num_threads(nThreads);
    if (nblock < 1)
    {
        nblock = 1;
    }
    if (nblock > ncols)
    {
        nblock = ncols;
    }

    u_int64_t offset_cols = 0;
    u_int64_t ncols_block = ncols / nblock;
    u_int64_t ncols_res = ncols % nblock;
    u_int64_t ncols_alloc = ncols_block;
    if (ncols_res > 0)
    {
        ncols_alloc += 1;
    }
    Element *dst_ = NULL;
    // Element *aux = (Element *)malloc(sizeof(Element) * n * ncols_alloc);
    Element *aux = new Element[n * ncols_alloc];
    if (nblock > 1)
    {
        dst_ = (Element *)malloc(sizeof(Element) * n * ncols_alloc);
    }
    else
    {
        dst_ = dst;
    }
    for (u_int64_t ib = 0; ib < nblock; ++ib)
    {
        uint64_t aux_ncols = ncols_block;
        if (ib < ncols_res)
        {
            aux_ncols += 1;
        }
        fft_block_iters(dst_, src, n, offset_cols, aux_ncols, ncols, nphase, aux);

        if (nblock > 1)
        {
#pragma omp parallel for schedule(static)
            for (u_int64_t ie = 0; ie < n; ++ie)
            {
                u_int64_t offset2 = ie * ncols + offset_cols;
                std::memcpy(&dst[offset2], &dst_[ie * aux_ncols], aux_ncols * sizeof(Element));
            }
        }
        offset_cols += aux_ncols;
    }
    if (nblock > 1)
    {
        free(dst_);
    }
    // free(aux);
    delete[] aux;
}

template <typename Field>
void FFT<Field>::ifft_block(Element *dst, Element *src, u_int64_t n, u_int64_t ncols, u_int64_t nphase, u_int64_t nblock)
{

    Element *dst_;
    if (dst == NULL)
    {
        dst_ = src;
    }
    else
    {
        dst_ = dst;
    }
    omp_set_num_threads(nThreads);
    fft_block(dst_, src, n, ncols, nphase, nblock);

    u_int64_t domainPow = log2(n);
    u_int64_t nDiv2 = n >> 1;
#pragma omp parallel for
    for (u_int64_t i = 1; i < nDiv2; i++)
    {

        u_int64_t r = n - i;
        u_int64_t offset_r = ncols * r;
        u_int64_t offset_i = ncols * i;
        for (uint64_t k = 0; k < ncols; k++)
        {
            Element tmp;
            f.copy(tmp, dst_[offset_i + k]);
            f.mul(dst_[offset_i + k], dst_[offset_r + k], powTwoInv[domainPow]);
            f.mul(dst_[offset_r + k], tmp, powTwoInv[domainPow]);
        }
    }

    u_int64_t offset_n = ncols * (n >> 1);
    for (uint64_t k = 0; k < ncols; k++)
    {
        f.mul(dst_[k], dst_[k], powTwoInv[domainPow]);
        f.mul(dst_[offset_n + k], dst_[offset_n + k], powTwoInv[domainPow]);
    }
}
