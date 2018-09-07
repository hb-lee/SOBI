#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include "light_matrix.h"
#include "svdcmp.h"

static double xarg1,xarg2;
#define MAX(a,b) (xarg1=(a),xarg2=(b),(xarg1) > (xarg2) ?\
(xarg1) : (xarg2))
#define EOF -1


void init_tau(int *tau, int length);
void read_eeg_data(char * file_name, Mat * mat);
void sobi(Mat * eeg_data, int * TAU, int fs, double jthresh, Mat *W);
void store_result(Mat * Res);

int main()
{
    char file_name[20];
    strcpy(file_name, "DATA.bin");
    Mat eeg_data;
    Mat W;
    int fs = 1000;

    int * TAU;
    double jthresh = 1e-5;
    MatCreate(&eeg_data, 128, 81626);
    MatCreate(&W, 128, 128);
    read_eeg_data(file_name, &eeg_data);
    init_tau(TAU, 42);
    //sobi(&eeg_data, TAU, fs, jthresh, &W);
    return 0;
}

void init_tau(int *tau, int length)
{
    int i, j;
    tau = (int *)malloc(length * sizeof(int));
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
            mat->element[j][i] = data[i * 128 + j];
        }
    }
    /*
    while (fread(buffer, sizeof(float), 128, file) != 0)
    {
        data = (float *) buffer;
        for (i = 0; i < 128; i++)
        {
            mat->element[i][line] = data[i];
            printf("data = %f\n", data[i]);
        }
        break;
        line++;
        fseek(file, line * 128 * sizeof(float), SEEK_SET);
    }
    */
    fclose(file);
    free(buffer);
    /*
    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < 5; j++)
        {
            printf("i = %d, j = %d, data = %lf\n", i, j, mat->element[i][j]);
        }
    }
    */
    return;
}

void store_result(Mat * Res)
{
    FILE *fp;
    int i, j;
    if ((fp = fopen("data.txt", "w")) == NULL)
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

void sobi(Mat * eeg_data, int * TAU, int fs, double jthresh, Mat *W)
{
    int i, j, k, s, t, taui;
    double temp = 0.0f;
    int Ntau = 42;
    int taumax = 350;
    int Nsn = eeg_data->row;
    TriMat * Rtau_all;
    Rtau_all = TriMatCreate(Rtau_all, Nsn, Nsn, Ntau);
    Rtau_all = TriMatZeros(Rtau_all);
    int lblk = 1000;
    int t_end = 1 + eeg_data->col - lblk - taumax;
    int TT = 0;
    Mat * XEEG;
    XEEG = MatCreate(XEEG, eeg_data->row, lblk + taumax);
    Mat * X;
    X = MatCreate(X, eeg_data->row, lblk + taumax);
    Mat * TX;
    TX = MatCreate(TX, lblk + taumax, eeg_data->row);
    Mat * RR;
    RR = MatCreate(RR, eeg_data->row, eeg_data->row);
    Mat * TRR;
    TRR = MatCreate(TRR, eeg_data->row, eeg_data->row);
    Mat * RR_sub;
    RR_sub = MatCreate(RR_sub, eeg_data->row, eeg_data->row);
    double * avg = (double *)malloc(eeg_data->row * sizeof(double));

    // step 1
    for (i = 1; i < t_end; i = i + lblk)
    {
        // 赋值XEEG
        for (j = 0; j < eeg_data->row; j++)
        {
            for (k = i, s = 0; k < i + lblk + taumax - 1; k++, s++)
            {
                XEEG->element[j][s] = eeg_data->element[j][k];
            }
        }
        // 求均值
        for (j = 0; j < XEEG->row; j++)
        {
            temp = 0;
            for (k = 0; k < XEEG->col; k++)
            {
                temp += XEEG->element[j][k];
            }
            avg[j] = temp / XEEG->col;
        }
        // 计算得X矩阵 (中心偏离)
        for (j = 0; j < X->row; j++)
        {
            for (k = 0; k < X->col; k++)
            {
                X->element[j][k] = XEEG->element[j][k] - avg[j];
            }
        }
        for (j = 0; j < Ntau; j++)
        {
            taui = TAU[j];
            TX = MatTrans(X, TX);
            RR = PartMatMul(X, TX, RR, 0, 0, 0, taui, lblk);
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
        TriMatNumMult(Rtau_all, 1 / TT);
    }

    // step 2
    Nsn = Rtau_all->row;
    Ntau = Rtau_all->h;
    Ntau = Ntau - 1;
    int N_princComp = Nsn;
    int JOINT_DIAG = 1;
    Mat * R0;
    R0 = MatCreate(R0, Rtau_all->row, Rtau_all->col);
    for (i = 0; i < Rtau_all->row; i++)
    {
        for (j = 0; j < Rtau_all->col; j++)
        {
            R0->element[i][j] = Rtau_all->element[i][j][0];
        }
    }
    Nsn = R0->row;
    int ttt = R0->col;
    int N = Nsn;
    Mat * Target;
    Target = MatCreate(Target, Nsn, ttt);
    for (i = 0; i < Nsn; i++)
    {
        for (j = 0; j < ttt; j++)
        {
            Target->element[i][j] = R0->element[i][j] + R0->element[j][i];
        }
    }
    MatNumMult(Target, 0.5);
    Mat * S;
    S = MatCreate(S, Nsn, ttt);
    double * lambda = zeros_vector(Nsn);
    svdcmp(Target->element, Nsn, ttt, lambda, S->element);
    // 按特征值排序
    double t0 = 0.0f;
    double * t1 = zeros_vector(Nsn);
    for (i = 0; i < Nsn - 1; i++)
    {
        for (j = i + 1; j < Nsn; j++)
        {
            if (lambda[i] < lambda[j])
            {
                t0 = lambda[i];
                lambda[i] = lambda[j];
                lambda[j] = t0;
            }
            // 矩阵U
            for (k = 0; k < Nsn; k++)
            {
                t1[k] = Target->element[k][i];
            }
            for (k = 0; k < Nsn; k++)
            {
                Target->element[k][i] = Target->element[k][j];
            }
            for (k = 0; k < Nsn; k++)
            {
                Target->element[k][j] = t1[k];
            }
            // 矩阵V
            for (k = 0; k < Nsn; k++)
            {
                t1[k] = S->element[k][i];
            }
            for (k = 0; k < Nsn; k++)
            {
                S->element[k][i] = S->element[k][j];
            }
            for (k = 0; k < Nsn; k++)
            {
                S->element[k][j] = t1[k];
            }
        }
    }

    Mat * B;
    B = MatCreate(B, Nsn, ttt);
    for (i = 0; i < Nsn; i++)
    {
        for (j = 0; j < ttt; j++)
        {
            if (i == j)
                Target->element[i][j] = sqrt(1 / (lambda[i] + 1e-20));
            else
                Target->element[i][j] = 0.0f;
        }
    }
    B = MatMul(S, Target, B);
    TriMat * Rtau_presphered;
    Rtau_presphered = TriMatCreate(Rtau_presphered, Rtau_all->row, Rtau_all->col, Ntau);
    for (i = 0; i < Ntau; i++)
    {
        // B * Rtau_all(:,:,i + 1)
        for (j = 0; j < Target->row; j++)
        {
            for (k = 0; k < Target->col; k++)
            {
                temp = 0.0f;
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
                temp = 0.0f;
                for (s = 0; s < Target->col; s++)
                {
                    temp += Target->element[j][s] * B->element[k][s];
                }
                Rtau_presphered->element[j][k][i] = temp;
            }
        }
    }
    // 联合对角化
    Mat * RRR = NULL;
    MatZeroConstruct(RRR, N_princComp, N_princComp * Ntau);
    for (i = 0; i < Ntau; i++)
    {
        for (j = 0; j < RRR->row; j++)
        {
            for (k = N_princComp * i; k < RRR->col; k++)
            {
                RRR->element[j][k] = Rtau_presphered->element[j][k % N_princComp][i]; // 三维铺成二维
            }
        }
    }
    Mat * B_mult = B;
    Mat * A = RRR;
    int m = A->row;
    int nm = A->col;
    Mat * B0 = NULL;
    B0 = MatCreate(B0, 3, 3);
    Mat * TB0 = NULL;
    TB0 = MatCreate(TB0, 3, 3);
    Mat * g = NULL;
    MatZeroConstruct(g, 3, nm / m); // 41
    Mat * tg = NULL;
    MatZeroConstruct(g, nm / m, 3);
    Mat * G = NULL;
    MatZeroConstruct(G, 2, 2);
    Mat * vcp = NULL;
    MatZeroConstruct(vcp, 3, 4);
    Mat * D = NULL;
    MatZeroConstruct(D, 3, 3);
    double * la = zeros_vector(3); // 列向量
    Mat * K = NULL;
    MatZeroConstruct(K, 3, 3);
    double * angles = zeros_vector(3); // 列向量
    //double * pair = zeros_vector(2);
    int * pair = NULL;
    pair = (int *) malloc(2 * sizeof(int));
    double c = 0.0f;
    double cs = 0.0f;
    Mat * V = NULL;
    MatEyeConstruct(V, m, m);
    int encore = 1;
    int iter_count = 0;
    double smax = 0.0f;
    Mat * W_unscaled = NULL;
    MatZeroConstruct(W_unscaled, V->row, B_mult->col);
    Mat *W_scaled = NULL;
    MatZeroConstruct(W_scaled, W_unscaled->row, W_unscaled->col);

    //double * Ip = zeros_vector(nm / m);
    //double * Iq = zeros_vector(nm / m);
    int * Ip = NULL;
    int * Iq = NULL;
    Ip = (int *) malloc((nm / m) * sizeof(int));
    Iq = (int *) malloc((nm / m) * sizeof(int));
    Mat * Res_3 = NULL;
    MatZeroConstruct(Res_3, 3, 3);
    Mat * Temp_3 = NULL;
    MatZeroConstruct(Temp_3, 3, 3);
    Mat * vcp_3 = NULL;
    MatZeroConstruct(vcp_3, 3, 3);
    double * d_3 = zeros_vector(3);
    double scaling_factor = 0.0f;
    double temp1, temp2;
    while (encore == 1)
    {
        encore = 0;
        smax = 0.0f;
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
                for (s = 0; s < nm / m; s++)
                {
                    g->element[0][s] = A->element[i][Ip[s]] - A->element[j][Iq[s]];
                }
                for (s = 0; s < nm / m; s++)
                {
                    g->element[1][s] = A->element[i][Iq[s]];
                }
                for (s = 0; s < nm / m; s++)
                {
                    g->element[2][s] = A->element[j][Ip[s]];
                }
                tg = MatTrans(g, tg);
                Res_3 = MatMul(g, tg, Res_3);
                Temp_3 = MatMul(B0, Res_3, Temp_3);
                Res_3 = MatMul(Temp_3, TB0, Res_3);
                svdcmp(Res_3->element, Res_3->row, Res_3->col, vcp_3->element, d_3);
                for (t = 0; t < vcp_3->row; t++)
                {
                    angles[t] = vcp->element[t][2];
                }
                if (angles[0] < 0)
                {
                    for (t = 0; t < vcp_3->row; t++)
                    {
                        angles[t] = (-1) * angles[t];
                    }
                }
                c = sqrt(0.5 + angles[0] / 2);
                cs = 0.5 * (angles[1] - 8 * angles[2]) / c;
                smax = MAX(cs, smax);
                if (abs(cs) > jthresh)
                {
                    encore = 1;
                    pair[0] = i;
                    pair[1] = j;
                    G->element[0][0] = c;
                    G->element[0][1] = -cs;
                    G->element[1][0] = s;
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
                            A->element[t][Ip[s]] = temp1;
                            temp2 = (-1) * cs * A->element[t][Ip[s]] + c * A->element[t][Iq[s]];
                            A->element[t][Iq[s]] = temp2;
                        }
                    }
                }
            }
        }
        W_scaled = MatTrans(V, W_scaled);
        W_unscaled = MatMul(W_scaled, B_mult, W_unscaled);
        for (i = 0; i < W_unscaled->row; i++)
        {
            temp = 0.0f;
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
        iter_count = iter_count + 1;
    }
    W = W_scaled;
    return;
}

