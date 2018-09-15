#ifndef SVDCMP_H_INCLUDED
#define SVDCMP_H_INCLUDED

static double PYTHAG(double a, double b);
double *zeros_vector(int n);
int dsvd(float **a, int m, int n, float *w, float **v);

#endif // SVDCMP_H_INCLUDED
