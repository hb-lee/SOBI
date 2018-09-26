#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include "light_matrix.h"
#include "svdcmp.h"

#include "/home/lhb/PAC/software/mkl/include/mkl.h"

#include "schur_lib/rt_nonfinite.h"
#include "schur_lib/schur_decompose.h"
#include "schur_lib/main.h"
#include "schur_lib/schur_decompose_terminate.h"
#include "schur_lib/schur_decompose_initialize.h"

/*
#include "libs/rt_nonfinite.h"
#include "libs/eigvalue.h"
#include "libs/main.h"
#include "libs/eigvalue_terminate.h"
#include "libs/eigvalue_initialize.h"
*/

/*
#include "libs/rt_nonfinite.h"
#include "libs/eigvalue.h"
#include "libs/main.h"
#include "libs/eigvalue_terminate.h"
#include "libs/eigvalue_initialize.h"
*/

static double xarg1,xarg2;
#define MAX(a,b) (xarg1=(a),xarg2=(b),(xarg1) > (xarg2) ?\
(xarg1) : (xarg2))
#define N_128 128
#define N_3 3
#define DEBUG 1


void init_tau(int *tau, int length);
void read_eeg_data(char * file_name, Mat * mat);
void sobi(Mat * eeg_data, int * TAU, int fs, double jthresh, Mat *W);
void store_result(Mat * Res);
void test_func();
int eig_128(Mat * target, Mat * eig_vector, double *eig_value);
int eig_3(Mat * target, Mat * eig_vector, double *eig_value);

int main()
{
    char file_name[20];
    strcpy(file_name, "DATA.bin");
    Mat * eeg_data;
    eeg_data = ConstructMat();
    Mat * W;
    W = ConstructMat();
    int fs = 1000;

    int * TAU;
    TAU = (int *)malloc(42 * sizeof(int));
    double jthresh = 1e-5;
    MatCreate(eeg_data, 128, 81626);
    MatCreate(W, 128, 128);
    printf("begin to read eeg data\n");
    read_eeg_data(file_name, eeg_data);
    printf("begin to init TAU\n");
    init_tau(TAU, 42);

    //test_func();
    printf("begin to execute sobi algorithnm\n");
    sobi(eeg_data, TAU, fs, jthresh, W);
    return 0;
}

void test_func()
{
    int i, j, k;
    Mat * m = NULL;
    m = ConstructMat();
    m = MatCreate(m, 3, 3);
    for (i = 0; i < m->row; i++)
    {
        for (j = 0; j < m->col; j++)
        {
            m->element[i][j] = i + j + 1;
        }
    }
    MatDump(m);
    double * w = zeros_vector(3);
    Mat * v = ConstructMat();
    v = MatCreate(v, 3, 3);
    schur_eig(m, v, w);
    MatDump(v);
    VectorDump(w, 3);
}

void init_tau(int *tau, int length)
{
    int i, j;
    if(tau == NULL){
		printf("tau alloc fail!\n");
		return;
	}
	tau[0] = 0;
    for (i = 1, j = 0; i < length; i++)
    {
        if (j < 10)
        {
            tau[i] = tau[i - 1] + 1;
            j = tau[i];
        } else if (j < 20)
        {
            tau[i] = tau[i - 1] + 2;
            j = tau[i];
        } else if (j < 100)
        {
            tau[i] = tau[i - 1] + 5;
            j = tau[i];
        } else
        {
            tau[i] = tau[i - 1] + 25;
            j = tau[i];
        }
    }
    return;
}

void read_eeg_data(char * file_name, Mat * mat)
{
    FILE *file;
    char * buffer;
    float * data;
    long len;
    size_t result;
    int i,j;
    file = fopen(file_name, "rb");
    if (!file)
    {
        fprintf(stderr, "can't open file %s", file_name);
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    len = ftell(file);
    rewind(file);
    buffer = (char *)malloc(sizeof(char) * len);
    memset(buffer, 0, sizeof(char) * len);
    if (buffer == NULL)
    {
        fprintf(stderr, "Memory error!\n");
        exit(2);
    }
    result = fread(buffer, 1, len, file);
    if (result < len)
    {
        fprintf(stdout, "Not read all data!\n");
        exit(3);
    }
    data = (float *) buffer;
    for (i = 0; i < mat->col; i++)
    {
        for (j = 0; j < mat->row; j++)
        {
            //double tmp = 0.0;
            //tmp = data[i * 128 + j];
            //mat->element[j][i] = floor(tmp * 10000.000f + 0.5) / 10000.000f;
            mat->element[j][i] = data[i * 128 + j];
        }
    }
    fclose(file);
    free(buffer);
    return;
}

void read_eig_vec(Mat *mat)
{
    char * buffer;
    double * data;
    long len;
    int i,j;
    size_t result;
    FILE *file;
    file = fopen("S1.bin", "rb");
    if (!file)
    {
        fprintf(stderr, "can't open file %s", "S1.bin");
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    len = ftell(file);
    rewind(file);
    buffer = (char *)malloc(sizeof(char) * len);
    memset(buffer, 0, sizeof(char) * len);
    if (buffer == NULL)
    {
        fprintf(stderr, "Memory error!\n");
        exit(2);
    }
    result = fread(buffer, 1, len, file);
    if (result < len)
    {
        fprintf(stdout, "Not read all data!\n");
        exit(3);
    }
    data = (double *) buffer;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->col; j++)
        {
            mat->element[i][j] = data[i * 128 + j];
        }
    }
    fclose(file);
    free(buffer);
    return;
}

void read_eig_val(double * val, int n)
{

    char * buffer;
    double * data;
    long len;
    int i;
    size_t result;
    FILE *file;
    file = fopen("lambda.bin", "rb");
    if (!file)
    {
        fprintf(stderr, "can't open file %s", "lambda.bin");
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    len = ftell(file);
    rewind(file);
    buffer = (char *)malloc(sizeof(char) * len);
    memset(buffer, 0, sizeof(char) * len);
    if (buffer == NULL)
    {
        fprintf(stderr, "Memory error!\n");
        exit(2);
    }
    result = fread(buffer, 1, len, file);
    if (result < len)
    {
        fprintf(stdout, "Not read all data!\n");
        exit(3);
    }
    data = (double *) buffer;
    for (i = 0; i < n; i++)
    {
            val[i] = data[i * n + i];
    }
    fclose(file);
    free(buffer);
    return;
}

void read_rtau(Mat *mat)
{
    char * buffer;
    double * data;
    long len;
    int i,j;
    size_t result;
    FILE *file;
    file = fopen("r_all.bin", "rb");
    if (!file)
    {
        fprintf(stderr, "can't open file %s", "S1.bin");
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    len = ftell(file);
    rewind(file);
    buffer = (char *)malloc(sizeof(char) * len);
    memset(buffer, 0, sizeof(char) * len);
    if (buffer == NULL)
    {
        fprintf(stderr, "Memory error!\n");
        exit(2);
    }
    result = fread(buffer, 1, len, file);
    if (result < len)
    {
        fprintf(stdout, "Not read all data!\n");
        exit(3);
    }
    data = (double *) buffer;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->col; j++)
        {
            mat->element[i][j] = data[i * mat->col + j];
        }
    }
    fclose(file);
    free(buffer);
    return;
}

void store_result(Mat * Res)
{
    FILE *fp;
    int i, j;
    if ((fp = fopen("result_W.txt", "w")) == NULL)
    {
        printf("file cann't open...\n");
        exit(1);
    }
    for (i = 0; i < Res->row; i++)
    {
        for (j = 0; j < Res->col; j++)
        {
            fprintf(fp, "%f\t", Res->element[i][j]);
        }
        fprintf(fp, "%c\n", '\n');
    }
    fclose(fp);
}


// real sys matrix
int eig(Mat * target, Mat * eig_vector, double *eig_value)
{
    int i, j;
    int row, col;
    row = target->row;
    col = target->col;
    double *A = NULL;
    A = (double *)malloc(row * row * sizeof(double));
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            A[i * row + j] = target->element[i][j];
        }
    }
    double * wr = NULL;
    //wr = (double *)malloc(col * sizeof(double));
    wr = zeros_vector(col);

    int ret = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', col, A, row, wr);
    for (i = 0; i < col; i++)
    {
        eig_value[i] = wr[i];
    }
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            eig_vector->element[i][j] = A[i * row + j];
        }
    }
    free(A);
    free(wr);
    A = NULL;
    wr = NULL;
    return ret;
}
/*
int eig_128_v2(Mat * target, Mat * eig_vector, double *eig_value)
{
    int i, j;
    double A[N_128 * N_128] = {0};
    creal_T V[N_128 * N_128];
    creal_T D[N_128 * N_128];
    for (i = 0; i < N_128; i++)
    {
        for (j = 0; j < N_128; j++)
        {
            A[i * N_128 + j] = target->element[i][j];
        }
    }

    eigvalue_initialize();
    eigvalue(A, V, D);

    for (i = 0; i < N_128; i++)
    {
        for (j = 0; j < N_128; j++)
        {
            eig_vector->element[i][j] = V[i * N_128 + j].re;
        }
    }
    for (i = 0; i < N_128; i++)
    {
        eig_value[i] = D[i * N_128 + i].re;
    }
    eigvalue_terminate();
    return 0;
}
*/
int eig_128(Mat * target, Mat * eig_vector, double *eig_value)
{
    int i, j;
    double A[N_128 * N_128];
    for (i = 0; i < N_128; i++)
    {
        for (j = 0; j < N_128; j++)
        {
            A[i * N_128 + j] = target->element[i][j];
        }
    }
    double wr[N_128] = {0};
    double wi[N_128] = {0};
    double vl[N_128 * N_128];
    double vr[N_128 * N_128];
    int ret = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', N_128, A, N_128, wr, wi, vl, N_128, vr, N_128);
    for (i = 0; i < N_128; i++)
    {
        eig_value[i] = wr[i];
    }
    for (i = 0; i < N_128; i++)
    {
        for (j = 0; j < N_128; j++)
        {
            eig_vector->element[i][j] = vr[i * N_128 + j];
        }
    }
}

int eig_3(Mat * target, Mat * eig_vector, double *eig_value)
{
    int i, j;
    double A[N_3 * N_3];
    for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            A[i * N_3 + j] = target->element[i][j];
        }
    }
    double wr[N_3] = {0};
    double wi[N_3] = {0};
    double vl[N_3 * N_3];
    double vr[N_3 * N_3];
    int ret = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', N_3, A, N_3, wr, wi, vl, N_3, vr, N_3);
    for (i = 0; i < N_3; i++)
    {
        eig_value[i] = wr[i];
    }
    for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            eig_vector->element[i][j] = vr[i * N_3 + j];
        }
    }
    int max = 0.0;
    int tmp = 0.0;
    max = pow(wr[0],2) + pow(wi[0],2);
    ret = 0;
    tmp = pow(wr[1],2) + pow(wi[1],2);
    if (max < tmp)
    {
        max = tmp;
        ret = 1;
    }
    tmp = pow(wr[2],2) + pow(wi[2],2);
    if (max < tmp)
    {
        max = tmp;
        ret = 2;
    }
    return ret;
}
/*
int eig_3_v2(Mat * target, Mat * eig_vector, double * eig_value)
{
	int i, j;
    double A[N_3 * N_3];
    creal_T V[N_3 * N_3];
    creal_T D[N_3 * N_3];
    for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            A[i * N_3 + j] = target->element[i][j];
        }
    }
    eigvalue_initialize();
    eigvalue(A, V, D);
	for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            eig_vector->element[i][j] = V[i * N_3 + j].re;
        }
    }
    for (i = 0; i < N_3; i++)
    {
        eig_value[i] = D[i * N_3 + i].re;
    }
    eigvalue_terminate();
	return 0;
}
*/
int schur_eig(Mat * target, Mat * eig_vector, double * eig_value)
{
    int i, j;
    int pos;
    double A[N_3 * N_3];
    double U[N_3 * N_3];
    double T[N_3 * N_3];
    for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            A[i * N_3 + j] = target->element[i][j];
        }
    }
    schur_decompose_initialize();
    schur_decompose(A, U, T);
    schur_decompose_terminate();
    for (i = 0; i < N_3; i++)
    {
        for (j = 0; j < N_3; j++)
        {
            eig_vector->element[i][j] = U[i * N_3 + j];
        }
    }
    for (i = 0; i < N_3; i++)
    {
        eig_value[i] = T[i * N_3 + i];
    }
    int max = eig_value[0];
    if (max < eig_value[1])
    {
        max = eig_value[1];
        pos = 1;
    }
    if (max < eig_value[2])
    {
        max = eig_value[2];
        pos = 2;
    }
    return pos;
}

int part_mul(Mat* src1, Mat* src2, Mat* dst, int src1_r_begin, int src2_r_begin, int src1_c_begin, int src2_c_begin, int mid)
{
    // 128 * 1000 * 1000 * 128
    int m = dst->row;
    int n = dst->col;
    int k = mid;
    double * A = (double *)malloc(m * k * sizeof(double));
    double * B = (double *)malloc(n * k * sizeof(double));
    double * C = (double *)malloc(m * n * sizeof(double));
    //double C[N_128][N_128] = {0};
    int i,j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[src1_r_begin + i][src1_c_begin + j];
        }
    }

    // row priority more faster
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            B[i * k + j] = src2->element[src2_r_begin + i][src2_c_begin + j];
        }
    }

    double alpha = 1.0;
    double beta = 0.0;
    int lda = k;
    int ldb = k;
    int ldc = n;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    for (i = 0; i < N_128; i++)
    {
        for (j = 0; j < N_128; j++)
        {
            dst->element[i][j] = C[i * N_128 + j];
        }
    }
    return 0;
}

int matrix_mul_tran(Mat * src1, Mat * src2, Mat * dst)
{
    // 3 * 41 * 41 * 3
    int m = dst->row;
    int n = dst->col;
    int k = src1->col;
    double * A = (double *)malloc(m * k * sizeof(double));
    double * B = (double *)malloc(n * k * sizeof(double));
    double * C = (double *)malloc(m * n * sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    int alpha = 1;
    int beta = 0;
    int lda = k;
    int ldb = k;
    int ldc = m;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            dst->element[i][j] = C[i * n + j];
        }
    }
    return 0;
}

int matrix_mul_notran(Mat * src1, Mat * src2, Mat * dst)
{
    int m = dst->row;
    int n = dst->col;
    int k = src1->col;
    double * A = (double *)malloc(m * k * sizeof(double));
    double * B = (double *)malloc(n * k * sizeof(double));
    double * C = (double *)malloc(m * n * sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    int alpha = 1;
    int beta = 0;
    int lda = k;
    int ldb = k;
    int ldc = m;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            dst->element[i][j] = C[i * n + j];
        }
    }
    return 0;
}

// for complex number
int factor_trans(Mat * src, Mat * dst)
{
    dst->element[0][0] = src->element[0][0];
    dst->element[0][1] = src->element[0][1] + src->element[0][2];
    dst->element[0][2] = 0.0;
    dst->element[1][0] = src->element[1][0] + src->element[2][0];
    dst->element[1][1] = src->element[1][1] + src->element[2][2] + src->element[1][2] + src->element[2][1];
    dst->element[1][2] = 0.0;
    dst->element[2][0] = src->element[2][0];
    dst->element[2][1] = src->element[2][1] + src->element[2][2];
    dst->element[2][2] = src->element[1][2] - src->element[1][1];
}

int set_factor(Mat * mat)
{// for matrix B
    mat->element[0][0] = 1.0;
    mat->element[0][1] = 0.0;
    mat->element[0][2] = 0.0;
    mat->element[1][0] = 0.0;
    mat->element[1][1] = 1.0;
    mat->element[1][2] = 1.0;
    mat->element[2][0] = 0.0;
    mat->element[2][1] = -1.0;
    mat->element[2][2] = 1.0;
}

void sobi(Mat * eeg_data, int * TAU, int fs, double jthresh, Mat *W)
{
    int i, j, k, s, t, taui;
    double temp = 0.0;
    int Ntau = 42;
    int taumax = 350;
    int Nsn = eeg_data->row;
    TriMat * Rtau_all;
    Rtau_all = ConstructTriMat();
    Rtau_all = TriMatCreate(Rtau_all, Nsn, Nsn, Ntau);
    int lblk = 1000;
    int t_end = 1 + eeg_data->col - lblk - taumax;
    int TT = 0;
    Mat * XEEG;
    XEEG = ConstructMat();
    XEEG = MatCreate(XEEG, eeg_data->row, lblk + taumax);
    Mat * X;
    X = ConstructMat();
    X = MatCreate(X, eeg_data->row, lblk + taumax);

    Mat * mult_A;
    mult_A = ConstructMat();
    mult_A = MatCreate(mult_A, eeg_data->row, lblk);
    Mat * mult_B;
    mult_B = ConstructMat();
    mult_B = MatCreate(mult_B, eeg_data->row, lblk);

    Mat * TX;
    TX = ConstructMat();
    TX = MatCreate(TX, lblk + taumax, eeg_data->row);
    Mat * RR;
    RR = ConstructMat();
    RR = MatCreate(RR, eeg_data->row, eeg_data->row);
    Mat * TRR;
    TRR = ConstructMat();
    TRR = MatCreate(TRR, eeg_data->row, eeg_data->row);
    Mat * RR_sub;
    RR_sub = ConstructMat();
    RR_sub = MatCreate(RR_sub, eeg_data->row, eeg_data->row);
    double * avg = (double *)malloc(eeg_data->row * sizeof(double));

    printf("sobi begin to execute step 1\n");
    /*
    // step 1
    for (i = 0; i < t_end; i = i + lblk)
    {
        // Set XEEG
        for (j = 0; j < eeg_data->row; j++)
        {
            for (k = i, s = 0; k < i + lblk + taumax; k++, s++)
            {
                XEEG->element[j][s] = eeg_data->element[j][k];
            }
        }
        printf("calculate mean\n");
        // get mean
        double num = lblk + taumax;
        for (j = 0; j < XEEG->row; j++)
        {
            temp = 0.0;
            for (k = 0; k < XEEG->col; k++)
            {
                temp += XEEG->element[j][k];
            }
            avg[j] = temp / num;
        }
        printf("calculate x - u\n");
        // get x - u
        for (j = 0; j < X->row; j++)
        {
            for (k = 0; k < X->col; k++)
            {
                X->element[j][k] = XEEG->element[j][k] - avg[j];
            }
        }

        // X(:,1:lblk)
        for (j = 0; j < mult_A->row; j++)
        {
            for (k = 0; k < mult_A->col; k++)
            {
                mult_A->element[j][k] = X->element[j][k];
            }
        }

        printf("begin to get Rtau_all\n");
        for (j = 0; j < Ntau; j++)
        {
            taui = TAU[j];

            // X(:,1+taui:lblk+taui)
            for (s =0; s < mult_B->row; s++)
            {
                for (t = 0; t < mult_B->col; t++)
                {
                    mult_B->element[s][t] = X->element[s][taui + t];
                }
            }
            matrix_mul_tran(mult_A, mult_B, RR);
            MatNumMult(RR, 1.0/81000.0);
            //part_mul(mult_A, mult_B, RR, 0, 0, 0, 0, lblk);
            //TX = MatTrans(X, TX);

            TRR = MatTrans(RR, TRR);
            RR_sub = MatAdd(RR, TRR, RR_sub);
            MatNumMult(RR_sub, 0.5);
            for (s = 0; s < Rtau_all->row; s++)
            {
                for (t = 0; t < Rtau_all->col; t++)
                {
                    Rtau_all->element[s][t][j] = Rtau_all->element[s][t][j] + RR_sub->element[s][t];
                }
            }

        }
        TT = TT + lblk;
        printf("index i = %d, TT = %d\n", i, TT);
    }
    */
    //TriMatNumMult(Rtau_all, 1.0 / ((double)TT));

    //debug
    Mat *r_all = ConstructMat();
    r_all = MatCreate(r_all, 128, 5376);
    read_rtau(r_all);

    for (i = 0; i < Rtau_all->h; i++)
    {
        for (j = 0; j < Rtau_all->row; j++)
        {
            for (s = 0; s < Rtau_all->col;s++)
            {
                Rtau_all->element[j][s][i] = r_all->element[j][i * Rtau_all->col + s];
            }
        }
    }

    printf("sobi begin to execute step 2\n");
    // step 2
    Nsn = Rtau_all->row;
    Ntau = Rtau_all->h;
    Ntau = Ntau - 1;
    int N_princComp = Nsn;

    int JOINT_DIAG = 1;
    Mat * R0;
    R0 = ConstructMat();
    R0 = MatCreate(R0, Rtau_all->row, Rtau_all->col);
    Mat * TR0;
    TR0 = ConstructMat();
    TR0 = MatCreate(TR0, Rtau_all->col, Rtau_all->row);
    for (i = 0; i < Rtau_all->row; i++)
    {
        for (j = 0; j < Rtau_all->col; j++)
        {
            R0->element[i][j] = Rtau_all->element[i][j][0];
        }
    }
    TR0 = MatTrans(R0, TR0);
    Nsn = R0->row;

    int ttt = R0->col;
    int N = Nsn;
    Mat * Target;
    Target = ConstructMat();
    Target = MatCreate(Target, Nsn, ttt);
    for (i = 0; i < Nsn; i++)
    {
        for (j = 0; j < ttt; j++)
        {
            Target->element[i][j] = (R0->element[i][j] + TR0->element[i][j]) * 0.5;
        }
    }
    Mat * S;
    S = ConstructMat();
    S = MatCreate(S, Nsn, ttt);

    double * lambda = zeros_vector(Nsn);
    printf("begin to eig 128\n");
    eig(Target, S, lambda);
    //eig_128_v2(Target, S, lambda);

/*
    Nsn = Rtau_all->row;
    int ttt = Rtau_all->row;
    Mat * S = ConstructMat();
    MatCreate(S, Nsn, ttt);
    double *lambda = zeros_vector(Nsn);
    read_eig_vec(S);
    read_eig_val(lambda, Nsn);
    Mat * Target;
    Target = ConstructMat();
    Target = MatCreate(Target, Nsn, ttt);
*/
/*
    Mat *Ts = ConstructMat();
    Ts = MatCreate(Ts, S->col, S->row);
    MatTrans(S, Ts);
    Mat *La = ConstructMat();
    La = MatCreate(La, S->col, S->row);
    for (i = 0; i < S->col; i++)
    {
        La->element[i][i] = lambda[i];
    }
    MatMul(S, La, Target);
    MatMul(Target, Ts, La);
    MatDump(La);
    if(DEBUG)
        return;
    */
    // sort the eigenvalue
    double t0 = 0.0;
    int noswap = 1;
    // sort eigenvalue and eigenvector (big to small sequence >...>)
    for (i = 0; i < Nsn; i++)
    {
        noswap = 1;
        for (j = 0; j < Nsn - 1 - i; j++)
        {
            if (lambda[j] < lambda[j + 1])
            {
                t0 = lambda[j];
                lambda[j] = lambda[j + 1];
                lambda[j + 1] = t0;

                // matrix V in eig decomposition
                for (k = 0; k < S->row; k++)
                {
                    t0 = S->element[k][j];
                    S->element[k][j] = S->element[k][j + 1];
                    S->element[k][j + 1] = t0;
                }
                noswap = 0;
            }
        }
        if (noswap == 1)
                break;
    }


    printf("get sorted eig value done!\n");
    Mat * B;
    B = ConstructMat();
    B = MatCreate(B, Nsn, ttt);
    Mat * TB;
    TB = ConstructMat();
    TB = MatCreate(TB, ttt, Nsn);
    for (i = 0; i < Nsn; i++)
    {
        for (j = 0; j < ttt; j++)
        {
            if (i == j)
                Target->element[i][j] = sqrt(1.0 / (lambda[i] + 1e-20));
            else
                Target->element[i][j] = 0.0;
        }
    }
    TB = MatMul(S, Target, TB);

    B = MatTrans(TB, B);
    TriMat * Rtau_presphered;
    Rtau_presphered = ConstructTriMat();
    Rtau_presphered = TriMatCreate(Rtau_presphered, Rtau_all->row, Rtau_all->col, Ntau);

    //debug
    Mat * input = ConstructMat();
    input = MatCreate(input, 128,128);
    Mat * output = ConstructMat();
    output = MatCreate(output,128,128);
    for (i = 0; i < Ntau; i++)
    {
        for(j=0;j<128;j++)
        {
            for(k=0;k<128;k++)
            {
                input->element[j][k]=Rtau_all->element[j][k][i+1];
            }
        }
        matrix_mul_notran(B,input,output);
        matrix_mul_notran(output,TB,input);
        for(j=0;j<128;j++)
        {
            for(k=0;k<128;k++)
            {
                Rtau_presphered->element[j][k][i] = input->element[j][k];
            }
        }
    }
    /*
    for (i = 0; i < Ntau; i++)
    {
        // B * Rtau_all(:,:,i + 1)
        for (j = 0; j < Target->row; j++)
        {
            for (k = 0; k < Target->col; k++)
            {
                temp = 0.0;
                for (s = 0; s < B->col; s++)
                {
                    temp += B->element[j][s] * Rtau_all->element[s][k][i + 1];
                }
                Target->element[j][k] = temp;
            }
        }
        // Target * B'
        for (j = 0; j < Rtau_presphered->row; j++)
        {
            for (k = 0; k < Rtau_presphered->col; k++)
            {
                temp = 0.0;
                for (s = 0; s < Target->col; s++)
                {
                    temp += Target->element[j][s] * TB->element[s][k];
                }
                Rtau_presphered->element[j][k][i] = temp;
            }
        }
    }
    */
    printf("execute joint diagonalization\n");
    // Joint Diag
    Mat * RRR = NULL;
    RRR = ConstructMat();
    MatZeroConstruct(RRR, N_princComp, N_princComp * Ntau);

    for (i = 0; i < Ntau; i++)
    {
        for (j = 0; j < RRR->row; j++)
        {
            for (k = N_princComp * i, s = 0; s < Rtau_presphered->col; s++, k++)
            {
                RRR->element[j][k] = Rtau_presphered->element[j][s][i]; // 3-D to 2-D
            }
        }
    }
    Mat * B_mult = NULL;
    B_mult = B;
    Mat * A = NULL;
    A = RRR;
    int m = A->row;
    int nm = A->col;
    Mat * B0 = NULL;
    B0 = ConstructMat();
    B0 = MatCreate(B0, 3, 3);
    set_factor(B0);
    Mat * TB0 = NULL;
    TB0 = ConstructMat();
    TB0 = MatCreate(TB0, 3, 3);
    TB0 = MatTrans(B0, TB0);
    Mat * g = NULL;
    g = ConstructMat();
    MatZeroConstruct(g, 3, nm / m); // 41
    Mat * tg = NULL;
    tg = ConstructMat();
    MatZeroConstruct(tg, nm / m, 3);
    Mat * G = NULL;
    G = ConstructMat();
    MatZeroConstruct(G, 2, 2);
    Mat * vcp = NULL;
    vcp = ConstructMat();
    MatZeroConstruct(vcp, 3, 4);
    Mat * D = NULL;
    D = ConstructMat();
    MatZeroConstruct(D, 3, 3);
    double * la = zeros_vector(3); // col vector
    Mat * K = NULL;
    K = ConstructMat();
    MatZeroConstruct(K, 3, 3);
    double * angles = zeros_vector(3); // col vector
    int * pair = NULL;
    pair = (int *) malloc(2 * sizeof(int));
    double c = 0.0;
    double cs = 0.0;
    Mat * V = NULL;
    V = ConstructMat();
    MatEyeConstruct(V, m, m);
    int encore = 1;
    int iter_count = 0;
    double smax = 0.0;
    Mat * W_unscaled = NULL;
    W_unscaled = ConstructMat();
    MatZeroConstruct(W_unscaled, V->row, B_mult->col);
    Mat *W_scaled = NULL;
    W_scaled = ConstructMat();
    MatZeroConstruct(W_scaled, W_unscaled->row, W_unscaled->col);

    int * Ip = NULL;
    int * Iq = NULL;
    Ip = (int *) malloc((nm / m) * sizeof(int));
    Iq = (int *) malloc((nm / m) * sizeof(int));
    Mat * Res_3 = NULL;
    Res_3 = ConstructMat();
    MatZeroConstruct(Res_3, 3, 3);
    Mat * Temp_3 = NULL;
    Temp_3 = ConstructMat();
    MatZeroConstruct(Temp_3, 3, 3);
    Mat * vcp_3 = NULL;
    vcp_3 = ConstructMat();
    MatZeroConstruct(vcp_3, 3, 3);
    double * d_3 = zeros_vector(3);
    double scaling_factor = 0.0;
    double factor_j = 8.0;
    double temp1 = 0.0, temp2 = 0.0;

    Mat * tmp_A = NULL;
    tmp_A = ConstructMat();
    MatCreate(tmp_A, A->row, 2 * (nm / m));
    printf("begin to iteration !\n");
    while (encore == 1)
    {
        encore = 0;
        smax = 0.0;
        for (i = 0; i < m - 1; i++)
        {
            for (k = i, s = 0; k < nm; k = k + m, s++)
            {
                Ip[s] = k;
            }
            for (j = i + 1; j < m; j++)
            {
                for (k = j, s = 0; k < nm; k = k + m, s++)
                {
                    Iq[s] = k;
                }
                for (s = 0; s < (nm / m); s++)
                {
                    g->element[0][s] = A->element[i][Ip[s]] - A->element[j][Iq[s]];
                }
                for (s = 0; s < (nm / m); s++)
                {
                    g->element[1][s] = A->element[i][Iq[s]] + A->element[j][Ip[s]];
                    //g->element[1][s] = A->element[i][Iq[s]];
                }
                for (s = 0; s < (nm / m); s++)
                {
                    g->element[2][s] = A->element[j][Ip[s]] - A->element[i][Iq[s]];
                    //g->element[2][s] = A->element[j][Ip[s]];
                }
                matrix_mul_tran(g, g, Res_3);
                //int tmp_pos = schur_eig(Res_3, vcp_3, d_3);
                //Temp_3 = MatMul(B0, Res_3, Temp_3);
                //Res_3 = MatMul(Temp_3, TB0, Res_3);
                //eig_3(Res_3, vcp_3, d_3);
				//MatDump(Res_3);
				//eig(Res_3, vcp_3, d_3);
				//MatDump(vcp_3);
				//VectorDump(d_3, 3);
               eig(Res_3, vcp_3, d_3);

                //matrix_mul_tran(g, g, Res_3);
                //factor_trans(Res_3, Temp_3);
                //int tmp_pos = eig_3(Temp_3, vcp_3, d_3);

				/*
                double tmp = d_3[0];
                int tmp_pos = 0;
                // sort(diag(D)) K(3)

                if (tmp < d_3[1])
                {
                    tmp_pos = 1;
                    tmp = d_3[1];
                }
                if (tmp < d_3[2])
                {
                    tmp_pos = 2;
                    tmp = d_3[2];
                }
                */
                int tmp_pos = 2;
                for (t = 0; t < vcp_3->row; t++)
                {
                    angles[t] = vcp_3->element[t][tmp_pos];
                }

                if (angles[0] < 0.0)
                {
                    for (t = 0; t < vcp_3->row; t++)
                    {
                        angles[t] = (-1) * angles[t];
                    }
                }
                c = sqrt(0.5 + angles[0] / 2.0);
		//factor_j = 1.0;
                cs = 0.5 * (angles[1] - factor_j * angles[2]) / c;
                //double su = pow(c,2) + pow(cs,2);
                //printf("c = %f, cs = %f, sum = %f\n",c,cs, su);
                if (cs > smax)
                    smax = cs;
                if (fabs(cs) > jthresh)
                {
                    encore = 1;
                    pair[0] = i;
                    pair[1] = j;
                    G->element[0][0] = c;
                    G->element[0][1] = -cs;
                    G->element[1][0] = cs;
                    G->element[1][1] = c;

                    // V(:,pair)       = V(:,pair)*G  in matlab
                    for (t = 0; t < V->row; t++)
                    {
                        temp1 = V->element[t][pair[0]] * G->element[0][0] + V->element[t][pair[1]] * G->element[1][0];
                        temp2 = V->element[t][pair[0]] * G->element[0][1] + V->element[t][pair[1]] * G->element[1][1];
                        V->element[t][pair[0]] = temp1;
                        V->element[t][pair[1]] = temp2;
                    }

					// A(pair,:)       = G' * A(pair,:) ;
                    for (t = 0; t < A->col; t++)
                    {
                        temp1 = G->element[0][0] * A->element[pair[0]][t] + G->element[1][0] * A->element[pair[1]][t];
                        temp2 = G->element[0][1] * A->element[pair[0]][t] + G->element[1][1] * A->element[pair[1]][t];
                        A->element[pair[0]][t] = temp1;
                        A->element[pair[1]][t] = temp2;
                    }

					//  A(:,[Ip Iq])    = [ c*A(:,Ip)+s*A(:,Iq) -conj(s)*A(:,Ip)+c*A(:,Iq) ]
                    for (t = 0; t < A->row; t++)
                    {
                        for (s = 0; s < nm / m; s++)
                        {
                            temp1 = c * A->element[t][Ip[s]] + cs * A->element[t][Iq[s]];
                            tmp_A->element[t][s] = temp1;
                            temp2 = (-1) * cs * A->element[t][Ip[s]] + c * A->element[t][Iq[s]];
                            tmp_A->element[t][s + (nm / m)] = temp2;
                        }
                    }
                    for (t = 0; t < A->row; t++)
                    {
                        for (s = 0; s < nm / m; s++)
                        {
                            A->element[t][Ip[s]] = tmp_A->element[t][s];
                        }
                    }
                    for (t = 0; t < A->row; t++)
                    {
                        for (s = 0; s < nm / m; s++)
                        {
                            A->element[t][Iq[s]] = tmp_A->element[t][s + (nm / m)];
                        }
                    }
                }
            }
        }
        //factor_j = 128.0;
        iter_count = iter_count + 1;
        printf("iter times = %d, accuracy = %f\n", iter_count, smax);
    }
    printf("finish the iteration!\n");
    W_scaled = MatTrans(V, W_scaled);
    W_unscaled = MatMul(W_scaled, B_mult, W_unscaled);
    for (i = 0; i < W_unscaled->row; i++)
        {
            temp = 0.0;
            for (j = 0; j < W_unscaled->col; j++)
            {
                temp1 += pow(W_unscaled->element[i][j], 2);
            }
            scaling_factor = sqrt(temp1);
            for (j = 0; j < W_unscaled->col; j++)
            {
                W_scaled->element[i][j] = W_unscaled->element[i][j] / scaling_factor;
            }
        }
    W = W_scaled;
    store_result(W);
    return;
}


