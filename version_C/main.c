#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "matrix.h"
#include "sse_opt.h"

#include "/home/lhb/PAC/software/mkl/include/mkl.h"
//#include "/home/SystemSoftware/intel2017/mkl/include/mkl.h"

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
    memset(TAU, 0, 42 * sizeof(int));
    double jthresh = 1e-5;
    MatZeroConstruct(eeg_data, 128, 81626);
    MatZeroConstruct(W, 128, 128);
    printf("begin to read eeg data\n");
    read_eeg_data(file_name, eeg_data);
    printf("begin to init TAU\n");
    init_tau(TAU, 42);

    //test_func();
    printf("begin to execute sobi algorithnm\n");
    sobi(eeg_data, TAU, fs, jthresh, W);
    return 0;
}

static inline double myinvsqrt(double x)
{
    double xhalf = 0.5*x;
    float xf = (float)x;
    __asm__
    {
        movss xmm1,xf
        rsqrtss xmm1,xmm1
        movss xf,xmm1
    }
    x = (double)xf;
    x = x*(1.5-xhalf*x*x);
    x = x*(1.5-xhalf*x*x);
    return x;
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
    fclose(file);
    data = (float *) buffer;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->col; j++)
        {
            mat->element[i][j] = data[j * 128 + i];
        }
    }
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
        fprintf(stderr, "can't open file %s", "r_all.bin");
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
    double * result_W = (double*)malloc(128*128*sizeof(double));
    memset(result_W,0,128*128*sizeof(double));
    for (i = 0; i < Res->row; i++)
    {
        for (j = 0; j < Res->col; j++)
        {
            result_W[i * 128 + j] = Res->element[i][j];
        }
    }
    if ((fp = fopen("result_W.bin", "w")) == NULL)
    {
        printf("file cann't open...\n");
        exit(1);
    }
	fwrite(result_W, sizeof(double), 128 * 128, fp);
	fclose(fp);
	return;
}

int eig_128(Mat *target, Mat * eig_vector, double *eig_value)
{
    int i, j;
	MKL_INT N, lda, info;
	MKL_INT isuppz[2 * 128];
    int row, col;
    row = target->row;
    col = target->col;
    N = row;
    lda = N;
    double *A = NULL;
    A = (double *)malloc(row * row * sizeof(double));
    memset(A, 0, row * row * sizeof(double));
    double * lam0 = NULL;
    lam0 = (double *)malloc(col * sizeof(double));
    memset(lam0,0,col*sizeof(double));
    for (i = 0; i < 128; i++)
    {
        for (j = 0; j < 128; j++)
        {
            A[i * col + j] = target->element[i][j];
        }
    }
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', N, A, lda, 1, 1, 1, 1, 2*dlamch("S"), &N, lam0, A, 128, isuppz);
    for (i = 0; i < 128; i++)
    {
        for (j = 0; j < 128; j++)
        {
            eig_vector->element[i][j] = A[i * 128 + 127 - j];
        }
    }
    for (i = 0; i < 128; i++)
    {
        eig_value[i] = lam0[127 - i];
    }
    free(A);
    free(lam0);
    A = NULL;
    lam0 = NULL;
    return info;
}

int eig_3(Mat *target, Mat *eig_vector, double *eig_value)
{
    int i, j;
    MKL_INT N, lda;
    int row, col;
    row = target->row;
    col = target->col;
    N = 3;
    lda = 3;
    double *A = NULL;
    A = (double *)malloc(row * row * sizeof(double));
    memset(A, 0, row * row * sizeof(double));
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            A[i * col + j] = target->element[i][j];
        }
    }
    double * wr = zeros_vector(3);
    int ret = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, A, lda, wr);
    for (i = 0; i < col; i++)
    {
        eig_value[i] = wr[i];
    }
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            eig_vector->element[i][j] = A[i * col + j];
        }
    }
    free(A);
    free(wr);
    A = NULL;
    wr = NULL;
    return ret;
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
    memset(A, 0, row * row * sizeof(double));
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            A[i * row + j] = target->element[i][j];
        }
    }
    double * wr = zeros_vector(col);

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

int matrix_mul_tran(Mat * src1, Mat * src2, Mat * dst)
{
    int m = dst->row;
    int n = dst->col;
    int k = src1->col;
    double * A = (double *)malloc(m * k * sizeof(double));
    double * B = (double *)malloc(n * k * sizeof(double));
    double * C = (double *)malloc(m * n * sizeof(double));
    memset(A, 0, m * k * sizeof(double));
    memset(B, 0, n * k * sizeof(double));
    memset(C, 0, m * n * sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    int alpha = 1.0;
    int beta = 0.0;
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
    free(A);
    free(B);
    free(C);
    A = NULL;
    B = NULL;
    C = NULL;
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
    memset(A, 0, m * k * sizeof(double));
    memset(B, 0, n * k * sizeof(double));
    memset(C, 0, m * n * sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    int alpha = 1.0;
    int beta = 0.0;
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
    free(A);
    free(B);
    free(C);
    A = NULL;
    B = NULL;
    C = NULL;
    return 0;
}

int factor_trans(Mat * src, Mat * dst)
{
    dst->element[0][0] = src->element[0][0];
    dst->element[0][1] = src->element[0][1] + src->element[0][2];
    dst->element[0][2] = 0.0;
    dst->element[1][0] = src->element[1][0] + src->element[2][0];
    dst->element[1][1] = src->element[1][1] + src->element[2][2] + src->element[1][2] + src->element[2][1];
    dst->element[1][2] = 0.0;
    dst->element[2][0] = 0.0;
    dst->element[2][1] = 0.0;
    dst->element[2][2] = src->element[1][2] + src->element[2][1] - src->element[1][1] - src->element[2][2];
    return 0;
}

int set_factor(Mat * mat, bool isTran)
{// for matrix B
    if(isTran)
    {
        mat->element[0][0] = 1.0;
        mat->element[0][1] = 0.0;
        mat->element[0][2] = 0.0;
        mat->element[1][0] = 0.0;
        mat->element[1][1] = 1.0;
        mat->element[1][2] = -41.0;
        mat->element[2][0] = 0.0;
        mat->element[2][1] = 1.0;
        mat->element[2][2] = 41.0;
    } else
    {
        mat->element[0][0] = 1.0;
        mat->element[0][1] = 0.0;
        mat->element[0][2] = 0.0;
        mat->element[1][0] = 0.0;
        mat->element[1][1] = 1.0;
        mat->element[1][2] = 1.0;
        mat->element[2][0] = 0.0;
        mat->element[2][1] = -41.0;
        mat->element[2][2] = 41.0;
    }
    return 0;
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
    //XEEG = MatCreate(XEEG, eeg_data->row, lblk + taumax);
    MatZeroConstruct(XEEG, eeg_data->row, lblk + taumax);
    Mat * X;
    X = ConstructMat();
    MatZeroConstruct(X, eeg_data->row, lblk + taumax);

    Mat * mult_A;
    mult_A = ConstructMat();
    MatZeroConstruct(mult_A, eeg_data->row, lblk);
    Mat * mult_B;
    mult_B = ConstructMat();
    MatZeroConstruct(mult_B, eeg_data->row, lblk);

    Mat * TX;
    TX = ConstructMat();
    MatZeroConstruct(TX, lblk + taumax, eeg_data->row);
    Mat * RR;
    RR = ConstructMat();
    MatZeroConstruct(RR, eeg_data->row, eeg_data->row);
    Mat * TRR;
    TRR = ConstructMat();
    MatZeroConstruct(TRR, eeg_data->row, eeg_data->row);
    Mat * RR_sub;
    RR_sub = ConstructMat();
    MatZeroConstruct(RR_sub, eeg_data->row, eeg_data->row);
    //modify
    float * avg = (float *)malloc(eeg_data->row * sizeof(float));

    printf("sobi begin to execute step 1\n");

    // step 1
    for (i = 0; i < t_end; i = i + lblk)
    {
        // Set XEEG
        for (j = 0; j < XEEG->row; j++)
        {
            for (k = i, s = 0; k < i + lblk + taumax; k++, s++)
            {
                XEEG->element[j][s] = eeg_data->element[j][k];
            }
        }
        printf("calculate mean\n");
        // get mean
        memset(avg,0,eeg_data->row*sizeof(double));
        for (j = 0; j < XEEG->row; j++)
        {
            temp = 0.0;
            for (k = 0; k < XEEG->col; k++)
            {
                temp += XEEG->element[j][k];
            }
            avg[j] = temp / (lblk + taumax);
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

    temp = 1.0 / TT;
    TriMatNumMult(Rtau_all, temp);
/*
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
*/
    printf("sobi begin to execute step 2\n");
    // step 2
    Nsn = Rtau_all->row;
    Ntau = Rtau_all->h;
    Ntau = Ntau - 1;
    int N_princComp = Nsn;

    int JOINT_DIAG = 1;
    Mat * R0;
    R0 = ConstructMat();
    MatZeroConstruct(R0, Rtau_all->row, Rtau_all->col);
    Mat * TR0;
    TR0 = ConstructMat();
    MatZeroConstruct(TR0, Rtau_all->col, Rtau_all->row);
    for (i = 0; i < Rtau_all->row; i++)
    {
        for (j = 0; j < Rtau_all->col; j++)
        {
            R0->element[i][j] = Rtau_all->element[i][j][0];
        }
    }
    TR0 = MatTrans(R0, TR0);
    Nsn = Rtau_all->row;

    int ttt = R0->col;
    int N = Nsn;
    Mat * B;
    B = ConstructMat();
    MatZeroConstruct(B, Nsn, ttt);
    Mat * TB;
    TB = ConstructMat();
    MatZeroConstruct(TB, ttt, Nsn);
    Mat * Target;
    Target = ConstructMat();
    MatZeroConstruct(Target, Nsn, ttt);
    for (i = 0; i < Nsn; i++)
    {
        //for (j = 0; j < ttt; j++)
        for(j = i; j < ttt; j++)
        {
            Target->element[i][j] = (R0->element[i][j] + TR0->element[i][j]) * 0.5;
        }
    }
    Mat * S;
    S = ConstructMat();
    MatZeroConstruct(S, Nsn, ttt);

    double * lambda = zeros_vector(Nsn);
    printf("begin to eig 128\n");
    eig(Target, S, lambda);
    Mat * S1;
    S1 = ConstructMat();
    MatZeroConstruct(S1, ttt, Nsn);
    double * lam1 = zeros_vector(Nsn);
    /*
    for (i = 0; i < ttt; i++)
    {
        for (j = 0; j < Nsn; j++)
        {
            S1->element[i][j] = S->element[i][Nsn - j - 1];
        }
    }*/
    for (i = 0; i < Nsn; i++)
    {
        lam1[i] = sqrt(1.0 / (lambda[Nsn - 1 - i] + 1e-20));
    }

    for (i = 0; i < ttt; i++)
    {
        for (j = 0; j < Nsn; j++)
        {
            TB->element[i][j] = S->element[i][Nsn - 1 -j] * lam1[j];
        }
    }

    B = MatTrans(TB, B);
    TriMat * Rtau_presphered;
    Rtau_presphered = ConstructTriMat();
    Rtau_presphered = TriMatCreate(Rtau_presphered, Rtau_all->row, Rtau_all->col, Ntau);


    Mat * input = ConstructMat();
    MatZeroConstruct(input, B->row,B->row);
    Mat * output = ConstructMat();
    MatZeroConstruct(output,B->row,B->row);
    for (i = 0; i < Ntau; i++)
    {
        for(j=0;j<B->row;j++)
        {
            for(k=0;k<B->row;k++)
            {
                input->element[j][k]=Rtau_all->element[j][k][i+1];
            }
        }
        matrix_mul_notran(B,input,output);
        matrix_mul_notran(output,TB,input);
        for(j=0;j<B->row;j++)
        {
            for(k=0;k<B->row;k++)
            {
                Rtau_presphered->element[j][k][i] = input->element[j][k];
            }
        }
    }
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
    MatZeroConstruct(B0, 3, 3);
    set_factor(B0,false);
    Mat * TB0 = NULL;
    TB0 = ConstructMat();
    MatZeroConstruct(TB0, 3, 3);
    set_factor(TB0,true);

    Mat * g = NULL;
    g = ConstructMat();
    MatZeroConstruct(g, 3, nm / m); // 41
    Mat * G = NULL;
    G = ConstructMat();
    MatZeroConstruct(G, 2, 2);

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
    MatZeroConstruct(W_unscaled, V->row, B->col);
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
    //double scaling_factor = 0.0;
    double factor_j = 8.0;
    double temp1 = 0.0, temp2 = 0.0;

    Mat * tmp_A = NULL;
    tmp_A = ConstructMat();
    MatZeroConstruct(tmp_A, A->row, 2 * (nm / m));
    double givens[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
    printf("begin to iteration !\n");

    //modify
    double * tp_V = (double *)malloc(128*sizeof(double));
    double * tq_V = (double *)malloc(128*sizeof(double));
    double * tp_A = (double *)malloc(5248*sizeof(double));
    double * tq_A = (double *)malloc(5248*sizeof(double));
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
                    //g->element[1][s] = A->element[i][Iq[s]] + A->element[j][Ip[s]];
                    g->element[1][s] = A->element[i][Iq[s]];
                }
                for (s = 0; s < (nm / m); s++)
                {
                    //g->element[2][s] = A->element[j][Ip[s]] - A->element[i][Iq[s]];
                    g->element[2][s] = A->element[j][Ip[s]];
                }

                matrix_mul_tran(g,g,Temp_3);
                matrix_mul_notran(B0, Temp_3, Res_3);
                matrix_mul_notran(Res_3,TB0,Temp_3);
                eig(Temp_3, vcp_3, d_3);

                int tmp_pos = 2;
                for (t = 0; t < vcp_3->row; t++)
                {
                    angles[t] = vcp_3->element[t][tmp_pos];
                }

                if (angles[0] < 0.0)
                {
                    for (t = 0; t < vcp_3->row; t++)
                    {
                        angles[t] = -angles[t];
                    }
                }
                c = sqrt(0.5 + angles[0] / 2.0);
                cs = 0.5 * (angles[1] - 8 * angles[2]) / c;
                //cs = 0.5 * (angles[1]) / c;

                smax = (cs > smax) ? cs : smax;
                if (fabs(cs) > jthresh)
                {
                    encore = 1;
                    pair[0] = i;
                    pair[1] = j;

                    // modify
                    givens[1] = c;
                    givens[2] = -cs;
                    givens[3] = cs;
                    givens[4] = c;
                    memset(tp_V, 0, 128*sizeof(double));
                    memset(tq_V, 0, 128*sizeof(double));
                    memset(tp_A, 0, 5248*sizeof(double));
                    memset(tq_A, 0, 5248*sizeof(double));
                    cblas_drotm(128, &(V->element[i][0]),1,&(V->element[j][0]),1,givens);
                    cblas_drotm(5248,&(A->element[i][0]),1,&(A->element[j][0]),1,givens);
                    for (t = 0; t < A->row; t++)
                    {
                        cblas_drotm(41,&(A->element[t][i]),128,&(A->element[t][j]),128,givens);
                    }
                }
            }
        }
        //factor_j = 128.0;
        iter_count = iter_count + 1;
        printf("iter times = %d, accuracy = %f\n", iter_count, smax);
    }
    printf("finish the iteration!\n");
    //W_scaled = MatTrans(V, W_scaled);
    W_scaled = V;
    matrix_mul_notran(W_scaled,B_mult,W_unscaled);
    double * scaling_factor = (double*)malloc(128*sizeof(double));
    memset(scaling_factor,0,128*sizeof(double));
    for(i = 0; i < W_unscaled->row; i++)
    {
        temp = 0.0;
        for (j = 0; j < W_unscaled->col; j++)
        {
            temp += W_unscaled->element[i][j] * W_unscaled->element[i][j];
        }
        //scaling_factor[i] = sqrt(temp);
        scaling_factor[i] = myinvsqrt(temp);
    }

    for (i = 0; i < W_unscaled->row; i++)
    {
        for (j = 0; j < W_unscaled->col; j++)
        {
            //W_scaled->element[i][j] = W_unscaled->element[i][j] / scaling_factor[i];
            W_scaled->element[i][j] = W_unscaled->element[i][j] * scaling_factor[i];
        }
    }
    //W = W_scaled;
    store_result(W_scaled);
    return;
}


