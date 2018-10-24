#include <stdio.h>
#include <immintrin.h>

inline void mul_a_(double *Ap,double *Aq,double *G0,double *G1,double *G2,double *G3,int *length)
{
    register int i;
    register __m256d map;
    register __m256d maq;
    register __m256d mg0 = _mm256_set1_pd(*G0);
    register __m256d mg1 = _mm256_set1_pd(*G1);
    register __m256d mg2 = _mm256_set1_pd(*G2);
    register __m256d mg3 = _mm256_set1_pd(*G3);
    register __m256d m1,m2,m3,m4,m5,m6;
    for (i = 0; i < *length; i += 4)
    {
        map = _mm256_load_pd(Ap + i);
        maq = _mm256_load_pd(Aq + i);
        m1 = _mm256_mul_pd(mg0, map);
        m2 = _mm256_mul_pd(mg2, maq);
        m3 = _mm256_mul_pd(mg1, map);
        m4 = _mm256_mul_pd(mg3, maq);
        m5 = _mm256_add_pd(m1, m2);
        m6 = _mm256_add_pd(m3, m4);
        _mm256_store_pd(Ap + i, m5);
        _mm256_store_pd(Aq + i, m6);
    }
}
