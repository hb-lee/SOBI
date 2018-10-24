#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "matrix.h"
#include "mpi.h"
#include "omp.h"
#include "/home/pac008/intel/mkl/include/mkl.h"
//#include "/home/SystemSoftware/intel2017/mkl/include/mkl.h"
//#include "mkl.h"

#define DEBUG 1
#define PCOUNT 180


void init_tau(int *tau, int length);
void read_eeg_data(char * file_name, SMat * smat);
void store_result(Mat * Res);
void test_func();

void sobi_prepare(Mat * B, TriMat * Rtau_presphered);
void sobi_iter(Mat *A, double jthresh, Mat * TV, int numofblocks);
void sobi_reduce(Mat * TV, Mat * B_mult);
static int rank,size,stride;

int main(int argc, char * argv[])
{
	MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int i,j,s,k;
    int numofblocks = 0;
    double jthresh = 1e-5;
    double t_start,t_end;
    if(rank == 1){
        numofblocks = 21;
    }
    else{
        numofblocks = 20;
    }
    size = 2;
    stride = 20;    
    t_start = MPI_Wtime();
	    
    Mat * A_block = NULL;
    A_block = ConstructMat();
    MatZeroConstruct(A_block,128,128*numofblocks);
    double * block_tmp = (double*)calloc(128*128*numofblocks,sizeof(double));
    
	TriMat * Rtau_presphered = NULL;
	Mat * B_mult = NULL;
    Mat * TV = NULL;
    Rtau_presphered = ConstructTriMat();
    Rtau_presphered = TriMatOrderCreate(Rtau_presphered,41,128,128);
    B_mult = ConstructMat();
    B_mult = MatCreate(B_mult,128,128);
	TV = ConstructMat();
	MatEyeConstruct(TV,128,128);
        
    sobi_prepare(B_mult,Rtau_presphered);
    
	if (rank == 0)
	{
		for (i = 0; i < numofblocks; i++)
		{
			for (j = 0; j < 128; j++)
			{
				#pragma simd
				for (k = 0; k < 128; k++)
				{
					A_block->element[j][i*128 + k] = Rtau_presphered->element[i][j][k];
				}
			}
		}
	} else
	{
		for (i = 0; i < numofblocks; i++)
		{
			for (j = 0; j < 128; j++)
			{
				#pragma simd
				for (k = 0; k < 128; k++)
				{
					A_block->element[j][i*128 + k] = Rtau_presphered->element[i + 20][j][k];
				}
			}
		}
	}
//    printf("Interation start\n");
	MPI_Barrier(MPI_COMM_WORLD);
    sobi_iter(A_block, jthresh, TV, numofblocks);
    
    if(rank ==0)
    {
        sobi_reduce(TV, B_mult);
    }
    t_end = MPI_Wtime();
    if(rank==0)
    {
        printf("Time consumed by the program is :%f",t_end-t_start);
    }
    MPI_Finalize();
    
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

void read_eeg_data(char * file_name, SMat * smat)
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
    buffer = (char *)calloc(len, sizeof(char));
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
    for (i = 0; i < smat->row; i++)
    {
        for (j = 0; j < smat->col; j++)
        {
            smat->element[i][j] = data[j * 128 + i];
        }
    }
    free(buffer);
    return;
}

void store_result(Mat * Res)
{
    FILE *fp;
    int i, j;
    double * result_W = (double *)calloc(16384,sizeof(double));
    for (i = 0; i < Res->row; i++)
    {
		#pragma simd
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
	fwrite(result_W, sizeof(double), 16384, fp);
	fclose(fp);
	return;
}

// real sys matrix
int eig(Mat * target, Mat * eig_vector, double *eig_value)
{

    int i, j;
    int row, col;
    row = target->row;
    col = target->col;
    double *A = NULL;
    A = (double *)calloc(row * row, sizeof(double));
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
    double * A = (double *)calloc(m * k, sizeof(double));
    double * B = (double *)calloc(n * k, sizeof(double));
    double * C = (double *)calloc(m * n, sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    double alpha = 1.0;
    double beta = 0.0;
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

int matrix_muleig_3(Mat * src1, Mat * src2, Mat * eig_vector, double *eig_value)
{
    // 3 * 41 * 3
    double * A = (double *)calloc(123, sizeof(double));
    double * B = (double *)calloc(123, sizeof(double));
    double * C = (double *)calloc(9, sizeof(double));
    int i, j;
    for (i = 0; i < 3; i++)
    {
		#pragma simd
        for (j = 0; j < 41; j++)
        {
            A[i * 41 + j] = src1->element[i][j];
            B[i * 41 + j] = src2->element[i][j];
        }
    }
    double alpha = 1.0;
    double beta = 0.0;
    int lda = 41;
    int ldb = 41;
    int ldc = 3;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 41, alpha, A, lda, B, ldb, beta, C, ldc);
//    double * wr = zeros_vector(3);
    int ret = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', 3, C, 3, eig_value);
	/*
    for (i = 0; i < 3; i++)
    {
        eig_value[i] = wr[i];
    }
	*/
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            eig_vector->element[i][j] = C[i * 3 + j];
        }
    }
    free(A);
    free(B);
    free(C);
//    free(wr);
    A = NULL;
    B = NULL;
    C = NULL;
 //   wr = NULL;
    return ret;
}

int smatrix_mul_tran(SMat * src1, SMat * src2, SMat * dst)
{
    int m = dst->row;
    int n = dst->col;
    int k = src1->col;
    float * A = (float *)calloc(m * k, sizeof(float));
    float * B = (float *)calloc(n * k, sizeof(float));
    float * C = (float *)calloc(m * n, sizeof(float));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    float alpha = 1.0;
    float beta = 0.0;
    int lda = k;
    int ldb = k;
    int ldc = m;
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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
    double * A = (double *)calloc(m * k, sizeof(double));
    double * B = (double *)calloc(n * k, sizeof(double));
    double * C = (double *)calloc(m * n, sizeof(double));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    double alpha = 1.0;
    double beta = 0.0;
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

int smatrix_mul_notran(SMat * src1, SMat * src2, SMat * dst)
{

    int m = dst->row;
    int n = dst->col;
    int k = src1->col;
    float * A = (float *)calloc(m * k, sizeof(float));
    float * B = (float *)calloc(n * k, sizeof(float));
    float * C = (float *)calloc(m * n, sizeof(float));
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            A[i * k + j] = src1->element[i][j];
            B[i * k + j] = src2->element[i][j];
        }
    }
    float alpha = 1.0;
    float beta = 0.0;
    int lda = k;
    int ldb = k;
    int ldc = m;
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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

void sobi_prepare(Mat * B, TriMat * Rtau_presphered)
{
    char file_name[20];
    strcpy(file_name, "DATA.bin");
    SMat * eeg_data;
    eeg_data = ConstructSMat();
    SMatZeroConstruct(eeg_data, 128, 81626);
    printf("begin to read eeg data\n");
    read_eeg_data(file_name, eeg_data);
    
/*    
    int * TAU;
    TAU = (int *)malloc(42 * sizeof(int));
    printf("begin to init TAU\n");
    init_tau(TAU, 42);
  */
	int TAU[42] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,350};
    int i, j, k, s, t, taui;
    double temp = 0.0;
    int Ntau = 42;
    int taumax = 350;
    int Nsn = eeg_data->row;
    TriMat * Rtau_all = NULL;
    Rtau_all = ConstructTriMat();
    Rtau_all = TriMatOrderCreate(Rtau_all, Ntau, Nsn, Nsn); // 42 * 128 * 128

    int lblk = 1000;
    int t_end = 1 + eeg_data->col - lblk - taumax;
    int TT = 0;
    SMat * XEEG;
    XEEG = ConstructSMat();
    SMatZeroConstruct(XEEG, eeg_data->row, lblk + taumax);
    SMat * X;
    X = ConstructSMat();
    SMatZeroConstruct(X, eeg_data->row, lblk + taumax);

    SMat * mult_A;
    mult_A = ConstructSMat();
    SMatZeroConstruct(mult_A, eeg_data->row, lblk);
    SMat * mult_B;
    mult_B = ConstructSMat();
    SMatZeroConstruct(mult_B, eeg_data->row, lblk);

    SMat * TX;
    TX = ConstructSMat();
    SMatZeroConstruct(TX, lblk + taumax, eeg_data->row);
    SMat * RR;
    RR = ConstructSMat();
    SMatZeroConstruct(RR, eeg_data->row, eeg_data->row);
    SMat * TRR;
    TRR = ConstructSMat();
    SMatZeroConstruct(TRR, eeg_data->row, eeg_data->row);
    SMat * RR_sub;
    RR_sub = ConstructSMat();
    SMatZeroConstruct(RR_sub, eeg_data->row, eeg_data->row);
    float * avg = (float *)malloc(eeg_data->row * sizeof(float));
    printf("sobi begin to execute step 1\n");

    // step 1
    for (i = 0; i < t_end; i = i + lblk)
    {
        // Set XEEG
        for (j = 0; j < XEEG->row; j++)
        {
            //cblas_scopy(XEEG->col,&(eeg_data->element[j][i]),1,&(XEEG->element[j][0]),1);

            for (k = i, s = 0; k < i + lblk + taumax; k++, s++)
            {
                XEEG->element[j][s] = eeg_data->element[j][k];
            }

        }
        //printf("calculate mean\n");
        // get mean
        memset(avg,0,eeg_data->row*sizeof(float));
        for (j = 0; j < XEEG->row; j++)
        {
            //avg[j] = cblas_sasum(1350,&(XEEG->element[j][0]),1);

            for (k = 0; k < XEEG->col; k++)
            {
                avg[j] += XEEG->element[j][k];
            }

            avg[j] /= (lblk + taumax);
        }
        //printf("calculate x - u\n");
        // get x - u
        for (j = 0; j < X->row; j++)
        {
			#pragma simd
            for (k = 0; k < X->col; k++)
            {
                X->element[j][k] = XEEG->element[j][k] - avg[j];
            }
        }

        // X(:,1:lblk)
        for (j = 0; j < mult_A->row; j++)
        {
            //cblas_scopy(mult_A->col,&(X->element[j][0]),1,&(mult_A->element[j][0]),1);
			#pragma simd
            for (k = 0; k < mult_A->col; k++)
            {
                mult_A->element[j][k] = X->element[j][k];
            }

        }

        //printf("begin to get Rtau_all\n");
        for (j = 0; j < Ntau; j++)
        {
            taui = TAU[j];

            // X(:,1+taui:lblk+taui)
            for (s =0; s < mult_B->row; s++)
            {
                //cblas_scopy(mult_B->col,&(X->element[s][taui]),1,&(mult_B->element[s][0]),1);
				#pragma simd
                for (t = 0; t < mult_B->col; t++)
                {
                    mult_B->element[s][t] = X->element[s][taui + t];
                }

            }
            smatrix_mul_tran(mult_A, mult_B, RR);
            TRR = SMatTrans(RR, TRR);
            RR_sub = SMatAdd(RR, TRR, RR_sub);
            SMatNumMult(RR_sub, 0.5);
            for (s = 0; s < Rtau_all->row; s++)
            {
                for (t = 0; t < Rtau_all->col; t++)
                {
                    Rtau_all->element[j][s][t] = Rtau_all->element[j][s][t] + RR_sub->element[s][t];
                    //Rtau_all->element[s][t][j] = Rtau_all->element[s][t][j] + RR_sub->element[s][t];
                }
            }

        }
        TT = TT + lblk;
        //printf("index i = %d, TT = %d\n", i, TT);
    }

    temp = 1.0 / TT;
    TriMatNumMultOrder(Rtau_all, temp);
    //printf("sobi begin to execute step 2\n");
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
		#pragma simd
        for (j = 0; j < Rtau_all->col; j++)
        {
            R0->element[i][j] = Rtau_all->element[0][i][j];
        }
    }
    TR0 = MatTrans(R0, TR0);
    Nsn = Rtau_all->row;

    int ttt = R0->col;
    int N = Nsn;

    Mat * TB;
    TB = ConstructMat();
    MatZeroConstruct(TB, ttt, Nsn);
    Mat * Target;
    Target = ConstructMat();
    MatZeroConstruct(Target, Nsn, ttt);
    for (i = 0; i < Nsn; i++)
    {
        //for (j = 0; j < ttt; j++)
		#pragma simd
        for(j = i; j < ttt; j++)
        {
            Target->element[i][j] = (R0->element[i][j] + TR0->element[i][j]) * 0.5;
        }
    }
    Mat * S;
    S = ConstructMat();
    MatZeroConstruct(S, Nsn, ttt);

    double * lambda = zeros_vector(Nsn);
    //printf("begin to eig 128\n");
    eig(Target, S, lambda);
    double * lam1 = zeros_vector(Nsn);

	#pragma simd
    for (i = 0; i < Nsn; i++)
    {
        //lam1[i] = myinvsqrt(lambda[Nsn -1 - i] + 1e-10);
        //lam1[i] = sqrt(1.0 / (lam1[i] + 1e-20));
        lam1[i] = sqrt(1.0 / (lambda[Nsn - 1 - i] + 1e-20));
    }

    for (i = 0; i < ttt; i++)
    {
		#pragma simd
        for (j = 0; j < Nsn; j++)
        {
            //TB->element[i][j] = S->element[i][j] * lam1[j];
            TB->element[i][j] = S->element[i][Nsn - 1 -j] * lam1[j];
        }
    }

    B = MatTrans(TB, B);

    Mat * input = ConstructMat();
    MatZeroConstruct(input, B->row,B->row);
    Mat * output = ConstructMat();
    MatZeroConstruct(output,B->row,B->row);
    omp_set_num_threads(4);
    for (i = 0; i < Ntau; i++)
    {
        //#pragma omp parallel for private(j, k)
        for(j=0;j<B->row;j++)
        {
			#pragma simd
            for(k=0;k<B->row;k++)
            {
                input->element[j][k]=Rtau_all->element[i+1][j][k];
            }
        }
        matrix_mul_notran(B,input,output);
        matrix_mul_notran(output,TB,input);
        //#pragma omp parallel for private(j, k)
        for(j=0;j<B->row;j++)
        {
			#pragma simd
            for(k=0;k<B->row;k++)
            {
                Rtau_presphered->element[i][j][k] = input->element[j][k];
            }
        }
    }
       
    return;
}

void sobi_iter(Mat * A, double jthresh, Mat * TV, int numofblocks){
    int encore = 1;
    int i,j,k,s,t,ti,tj,tk;
    int iter_count = 0;
    double smax = 0.0;
    double * d_3 = zeros_vector(3);
    double temp1 = 0.0, temp2 = 0.0;
    double factor_j = 8.0;
    int m = A->row; //128
    int nm = A->col; //128*41 now 128
    double * angles = zeros_vector(3); // col vector
    double c = 0.0;
    double cs = 0.0;
    
    int * Ip = NULL;
    int * Iq = NULL;
    Ip = (int *) malloc(numofblocks * sizeof(int));
    Iq = (int *) malloc(numofblocks * sizeof(int));
	int * Ip_sum = NULL;
	int * Iq_sum = NULL;
    Ip_sum = (int *) calloc(41, sizeof(int));
    Iq_sum = (int *) calloc(41, sizeof(int));	
	Mat * A_all = NULL;
	if (rank == 0)
	{
		A_all = ConstructMat();
		MatCreate(A_all, 128, 5248);
	}
    
    double g_sumtmp[123] = {0};
	Mat * g_sum = NULL;
    g_sum = ConstructMat();
    MatZeroConstruct(g_sum, 3, 41);
    
    Mat * g = NULL;
    g = ConstructMat();
    MatCreate(g, 3, numofblocks); // 41 now 1
    
    double * g_block = (double*)calloc(3*numofblocks,sizeof(double));
    
    //用于rank0**********************************
    Mat * Res_3 = NULL;
    Res_3 = ConstructMat();
    MatZeroConstruct(Res_3, 3, 3);
    
    Mat * vcp_3 = NULL;
    vcp_3 = ConstructMat();
    MatZeroConstruct(vcp_3, 3, 3);
    //*******************************************
    
    double ccs[2];
	double givens[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
	int isFirst = 1;
	
    omp_set_num_threads(4);
    while (encore == 1)
    {
        encore = 0;
        smax = 0.0;
		if (iter_count < PCOUNT)
		{
			for (i = 0; i < m - 1; i++)
			{
				#pragma simd
				for (s = 0; s < numofblocks; s++)
				{
					Ip[s] = i + s*128;
				}
				for (j = i + 1; j < m; j++)
				{
					#pragma simd
					for (s = 0; s < numofblocks; s++)
					{
						Iq[s] = j + s*128;
						g_block[3*s] = A->element[i][Ip[s]] - A->element[j][Iq[s]];
						g_block[3*s + 1] = A->element[i][Iq[s]] + A->element[j][Ip[s]];
						g_block[3*s + 2] = A->element[j][Ip[s]] - A->element[i][Iq[s]];
					}
					if(rank == 1)
					{
						MPI_Send(g_block,63,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
					}
						
					if(rank == 0)
					{
						#pragma simd
						for(k = 0;k<60;k++)
						{
							g_sumtmp[k] = g_block[k];
						}
						MPI_Recv(&g_sumtmp[60],63,MPI_DOUBLE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
		   
					if(rank==0)
					{      
						for(k=0;k<3;k++)
							#pragma simd
							for(s=0;s<41;s++)
								g_sum->element[k][s] = g_sumtmp[3*s+k]; 
						
						//matrix_mul_tran(g_sum, g_sum, Res_3);
						//eig(Res_3, vcp_3, d_3);
						matrix_muleig_3(g_sum,g_sum,vcp_3,d_3);
						//int tmp_pos = 2;
						for (t = 0; t < vcp_3->row; t++)
						{
							angles[t] = vcp_3->element[t][2];
						}
						
						if (angles[0] < 0.0)
						{
							for (t = 0; t < vcp_3->row; t++)
							{
								angles[t] = -angles[t];
							}
						}
						c = sqrt(0.5 + angles[0] / 2.0);
						cs = 0.5 * (angles[1] - factor_j * angles[2]) / c;
						ccs[0] = c;
						ccs[1] = cs;
						//
						MPI_Send(ccs,2,MPI_DOUBLE,1,3,MPI_COMM_WORLD);
					}                
					
					if (rank == 1)
						MPI_Recv(ccs,2,MPI_DOUBLE,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					//MPI_Bcast(ccs,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
						
					//-----------------------------------------------------------------------
					//MPI_Barrier(MPI_COMM_WORLD);
					if (ccs[1] > smax)
						smax = ccs[1];
					if (fabs(ccs[1]) > jthresh)
					{
						encore = 1;
						givens[1] = ccs[0];
						givens[2] = -ccs[1];
						givens[3] = ccs[1];
						givens[4] = ccs[0];
						if(rank == 0)
						{
							#pragma omp parallel sections
							{
								#pragma omp section
								{
									cblas_drotm(128, &(TV->element[i][0]),1,&(TV->element[j][0]),1,givens);
								}
								#pragma omp section
								{
									cblas_drotm(1280,&(A->element[i][0]),1,&(A->element[j][0]),1,givens);
								}
								#pragma omp section
								{
									cblas_drotm(1280,&(A->element[i][1280]),1,&(A->element[j][1280]),1,givens);
								}
							}
							
						} else
						{
							#pragma omp parallel sections
							{
								#pragma omp section
								{
									//cblas_drotm(1280,&(A->element[i][0]),1,&(A->element[j][0]),1,givens);
									cblas_drotm(1344,&(A->element[i][0]),1,&(A->element[j][0]),1,givens);
								}
								#pragma omp section
								{
									//cblas_drotm(1408,&(A->element[i][1280]),1,&(A->element[j][1280]),1,givens);
									cblas_drotm(1344,&(A->element[i][1344]),1,&(A->element[j][1344]),1,givens);
								}
							}
						}
						#pragma omp parallel for private(t)
						for (t = 0; t < A->row; t++)
						{
							cblas_drotm(numofblocks,&(A->element[t][i]),128,&(A->element[t][j]),128,givens);
						}
					}
				}
			}
		} else
		{
			if (rank == 1)
			{
				double tmp_2D1[16384*21] = {0};
				#pragma simd
				for (ti = 0; ti < 21; ti++)
					for (tj = 0; tj < 128; tj++)
						for (tk = 0; tk < 128; tk++)
							tmp_2D1[ti*128*128+tj*128+tk] = A->element[tj][ti*128+tk];
				MPI_Send(tmp_2D1,16384*21,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
				return;
			} else 
			{
				if (isFirst == 1)
				{
					#pragma simd
					for (ti = 0; ti < 20; ti++)
						for (tj = 0; tj < 128; tj++)
							for (tk = 0; tk < 128; tk++)
								A_all->element[tj][ti*128+tk] = A->element[tj][ti*128+tk];     
					double tmp_2D0[16384*21] = {0};
					MPI_Recv(tmp_2D0,16384*21,MPI_DOUBLE,1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					#pragma simd
					for (ti = 0; ti < 21; ti++)
						for (tj = 0; tj < 128; tj++)
							for (tk = 0; tk < 128; tk++)
								A_all->element[tj][(20+ti)*128+tk] = tmp_2D0[ti*128*128+tj*128+tk];
					isFirst = 0;
				}
				for(i = 0;i<m-1;i++)
                {
					#pragma simd
					for (s = 0; s < 41; s++)
					{
						Ip_sum[s] = i + s*128;
					}
					for (j = i + 1; j < m; j++)
					{
						#pragma simd
						for (s = 0; s < 41; s++)
						{
							Iq_sum[s] = j + s*128;
							g_sum->element[0][s] = A_all->element[i][Ip_sum[s]] - A_all->element[j][Iq_sum[s]];
							g_sum->element[1][s] = A_all->element[i][Iq_sum[s]] + A_all->element[j][Ip_sum[s]];
							g_sum->element[2][s] = A_all->element[j][Ip_sum[s]] - A_all->element[i][Iq_sum[s]];
						}
						//matrix_mul_tran(g_sum,g_sum,Res_3);
						//eig(Res_3, vcp_3, d_3);
						matrix_muleig_3(g_sum,g_sum,vcp_3,d_3);
						for (t = 0; t < vcp_3->row; t++)
						{
							angles[t] = vcp_3->element[t][2];
						}
						
						if (angles[0] < 0.0)
						{
							for (t = 0; t < vcp_3->row; t++)
							{
								angles[t] = -angles[t];
							}
						}
						c = sqrt(0.5 + angles[0] / 2.0);
						cs = 0.5 * (angles[1] - factor_j * angles[2]) / c;
						if (cs > smax)
							smax = cs;
						if (fabs(cs) > jthresh)
						{
							encore = 1;
							givens[1] = c;
							givens[2] = -cs;
							givens[3] = cs;
							givens[4] = c;
							#pragma omp parallel sections
							{
								#pragma omp section
								{
									cblas_drotm(128, &(TV->element[i][0]),1,&(TV->element[j][0]),1,givens);
								}
								#pragma omp section
								{
									//cblas_drotm(1280,&(A_all->element[i][0]),1,&(A_all->element[j][0]),1,givens);
									cblas_drotm(1312,&(A_all->element[i][0]),1,&(A_all->element[j][0]),1,givens);
								}
								#pragma omp section
								{
									//cblas_drotm(1280,&(A_all->element[i][1280]),1,&(A_all->element[j][1280]),1,givens);
									cblas_drotm(1312,&(A_all->element[i][1312]),1,&(A_all->element[j][1312]),1,givens);
								}
								#pragma omp section
								{
									//cblas_drotm(1280,&(A_all->element[i][2560]),1,&(A_all->element[j][2560]),1,givens);
									cblas_drotm(1312,&(A_all->element[i][2624]),1,&(A_all->element[j][2624]),1,givens);
								}
								#pragma omp section
								{
									//cblas_drotm(1408,&(A_all->element[i][3840]),1,&(A_all->element[j][3840]),1,givens);
									cblas_drotm(1312,&(A_all->element[i][3936]),1,&(A_all->element[j][3936]),1,givens);
								}
							}
							
							//cblas_drotm(5248,&(A_all->element[i][0]),1,&(A_all->element[j][0]),1,givens);
							#pragma omp parallel for private(t)
							for (t = 0; t < A_all->row; t++)
							{
								cblas_drotm(41,&(A_all->element[t][i]),128,&(A_all->element[t][j]),128,givens);
							}
						}
					}
				}
			}
		}

        //factor_j = 128.0;
        iter_count = iter_count + 1;
        if(rank == 0){            
            printf("Rank = %d, iter times = %d, accuracy = %f\n",rank ,iter_count ,smax);      
        }
		//MPI_Barrier(MPI_COMM_WORLD);
    }
    return;
}

void sobi_reduce(Mat * TV, Mat * B_mult)
{
    int i, j;
    double temp = 0.0;
    Mat * W_unscaled = NULL;
    W_unscaled = ConstructMat();
    MatZeroConstruct(W_unscaled, TV->row, 128);
    Mat *W_scaled = NULL;
    W_scaled = ConstructMat();
    MatZeroConstruct(W_scaled, W_unscaled->row, W_unscaled->col);
	
    W_scaled = TV;
    matrix_mul_notran(W_scaled,B_mult,W_unscaled);
    double * scaling_factor = (double*)calloc(128,sizeof(double));
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

    #pragma simd
    for (i = 0; i < W_unscaled->row; i++)
    {
        cblas_dscal(128,scaling_factor[i],&(W_unscaled->element[i][0]),1);
        //for (j = 0; j < W_unscaled->col; j++)
        //{
            //W_scaled->element[i][j] = W_unscaled->element[i][j] / scaling_factor[i];
        //    W_scaled->element[i][j] = W_unscaled->element[i][j] * scaling_factor[i];
            //W_scaled->element[i][j+1] = W_unscaled->element[i][j+1] * scaling_factor[i];
            //W_scaled->element[i][j+2] = W_unscaled->element[i][j+2] * scaling_factor[i];
            //W_scaled->element[i][j+3] = W_unscaled->element[i][j+3] * scaling_factor[i];
        //}
    }
    //W = W_scaled;
    store_result(W_scaled);
    return;
}
