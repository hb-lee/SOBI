#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

//#define MAT_LEGAL_CHECKING

#define min(a, b) ((a) > (b) ? (b) : (a))
#define equal(a, b)	((a-b)<1e-7 && (a-b)>-(1e-7))

/************************************************************************/
/*                           Public Function                            */
/************************************************************************/

Mat * ConstructMat()
{
    Mat * mat;
    mat = (Mat *)malloc(sizeof(Mat));
    return mat;
}

SMat * ConstructSMat()
{
    SMat * smat;
    smat = (SMat *)malloc(sizeof(SMat));
    return smat;
}

TriMat * ConstructTriMat()
{
    TriMat * triMat;
    triMat = (TriMat *) malloc(sizeof(TriMat));
    return triMat;
}

TriSMat * ConstructTriSMat()
{
    TriSMat * triSMat;
    triSMat = (TriSMat *) malloc(sizeof(TriSMat));
    return triSMat;
}

TriMat* TriMatCreate(TriMat * triMat, int row, int col, int h)
{
    int i, j;
    triMat->element = (double ***)malloc(row * sizeof(double *));
    for (i = 0; i < row; i++)
    {
        triMat->element[i] = (double **)malloc(col * sizeof(double *));
        for (j = 0; j < col; j++)
        {
			triMat->element[i][j] = (double *) calloc(h, sizeof(double));
            //triMat->element[i][j] = (double *)malloc(h * sizeof(double));
        }
    }
    triMat->row = row;
    triMat->col = col;
    triMat->h = h;
    //TriMatZeros(triMat);
    return triMat;
}

TriSMat* TriSMatCreate(TriSMat * triSMat, int row, int col, int h)
{
    int i, j;
    triSMat->element = (float ***)malloc(row * sizeof(float *));
    for (i = 0; i < row; i++)
    {
        triSMat->element[i] = (float **)malloc(col * sizeof(float *));
        for (j = 0; j < col; j++)
        {
            triSMat->element[i][j] = (float *)malloc(h * sizeof(float));
        }
    }
    triSMat->row = row;
    triSMat->col = col;
    triSMat->h = h;
    TriSMatZeros(triSMat);
    return triSMat;
}

TriMat* TriMatOrderCreate(TriMat * triMat, int h, int row, int col)
{
    int i, j;
    triMat->element = (double ***)malloc(h * sizeof(double *));
    for (i = 0; i < h; i++)
    {
        triMat->element[i] = (double **)malloc(row * sizeof(double *));
        for (j = 0; j < row; j++)
        {
            //triMat->element[i][j] = (double *)malloc(col * sizeof(double));
			triMat->element[i][j] = (double *) calloc(col, sizeof(double));
        }
    }
    triMat->row = row;
    triMat->col = col;
    triMat->h = h;
    //TriMatZerosOrder(triMat);
    return triMat;
}

TriSMat* TriSMatOrderCreate(TriSMat * triSMat, int h, int row, int col)
{
    int i, j;
    triSMat->element = (float ***)malloc(h * sizeof(float *));
    for (i = 0; i < h; i++)
    {
        triSMat->element[i] = (float **)malloc(row * sizeof(float *));
        for (j = 0; j < row; j++)
        {
			triSMat->element[i][j] = (float *) calloc(col, sizeof(float));
            //triSMat->element[i][j] = (float *)malloc(col * sizeof(float));
        }
    }
    triSMat->row = row;
    triSMat->col = col;
    triSMat->h = h;
    //TriSMatZerosOrder(triSMat);
    return triSMat;
}

Mat* MatCreate(Mat* mat, int row, int col)
{
    Mat * p;
	int i;

	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return NULL;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double*) calloc(col, sizeof(double));
		//mat->element[i] = (double*)malloc(col * sizeof(double));
		if(mat->element[i] == NULL){
			int j;
			printf("mat create fail!\n");
			for(j = 0 ; j < i ; j++)
				free(mat->element[j]);
			free(mat->element);
			return NULL;
		}
	}

	mat->row = row;
	mat->col = col;
    //mat = MatZeros(mat);
	return mat;
}

SMat* SMatCreate(SMat* smat, int row, int col)
{
    SMat * p;
	int i;

	smat->element = (float**)malloc(row * sizeof(float*));
	if(smat->element == NULL){
		printf("smat create fail!\n");
		return NULL;
	}
	for(i = 0 ; i < row ; i++){
		smat->element[i] = (float*) calloc(col, sizeof(float));
		//smat->element[i] = (float*)malloc(col * sizeof(float));
		if(smat->element[i] == NULL){
			int j;
			printf("smat create fail!\n");
			for(j = 0 ; j < i ; j++)
				free(smat->element[j]);
			free(smat->element);
			return NULL;
		}
	}

	smat->row = row;
	smat->col = col;
    //smat = SMatZeros(smat);
	return smat;
}

void MatZeroConstruct(Mat * mat, int row, int col)
{
	int i, j;

	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double *) calloc(col, sizeof(double));
		//mat->element[i] = (double*)malloc(col * sizeof(double));
		//for (j = 0; j < col; j++)
        //    mat->element[i][j] = 0.0;
	}

	mat->row = row;
	mat->col = col;
}

void SMatZeroConstruct(SMat * smat, int row, int col)
{
	int i, j;

	smat->element = (float**)malloc(row * sizeof(float*));
	if(smat->element == NULL){
		printf("smat create fail!\n");
		return;
	}
	for(i = 0 ; i < row ; i++){
		smat->element[i] = (float *)calloc(col, sizeof(float));
		//smat->element[i] = (float*)malloc(col * sizeof(float));
		//for (j = 0; j < col; j++)
        //    smat->element[i][j] = 0.0;
	}

	smat->row = row;
	smat->col = col;
}

void MatEyeConstruct(Mat * mat, int row, int col)
{
    int i, j;

	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double*)malloc(col * sizeof(double));
		for (j = 0; j < col; j++)
            if (i == j)
                mat->element[i][j] = 1.0;
            else
                mat->element[i][j] = 0.0;
	}

	mat->row = row;
	mat->col = col;
}

void SMatEyeConstruct(SMat * smat, int row, int col)
{
    int i, j;

	smat->element = (float**)malloc(row * sizeof(float*));
	if(smat->element == NULL){
		printf("smat create fail!\n");
		return;
	}
	for(i = 0 ; i < row ; i++){
		smat->element[i] = (float*)malloc(col * sizeof(float));
		for (j = 0; j < col; j++)
            if (i == j)
                smat->element[i][j] = 1.0;
            else
                smat->element[i][j] = 0.0;
	}

	smat->row = row;
	smat->col = col;
}

void MatDelete(Mat* mat)
{
	int i;

	for(i = 0 ; i<mat->row ; i++)
		free(mat->element[i]);
	free(mat->element);
}

void SMatDelete(SMat* smat)
{
	int i;

	for(i = 0 ; i<smat->row ; i++)
		free(smat->element[i]);
	free(smat->element);
}

void MatDump(const Mat* mat)
{
    FILE * fp = fopen("matrix_debug.txt", "a+");
	int row,col;

#ifdef MAT_LEGAL_CHECKING
	if(mat == NULL){
		return ;
	}
#endif

	fprintf(fp, "Mat %dx%d:\n", mat->row, mat->col);
	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
            fprintf(fp, "%f\t", mat->element[row][col]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void SMatDump(const SMat* smat)
{
    FILE * fp = fopen("matrix_debug.txt", "a+");
	int row,col;

#ifdef MAT_LEGAL_CHECKING
	if(smat == NULL){
		return ;
	}
#endif

	fprintf(fp, "SMat %dx%d:\n", smat->row, smat->col);
	for(row = 0 ; row < smat->row ; row++){
		for(col = 0 ; col < smat->col ; col++){
            fprintf(fp, "%f\t", smat->element[row][col]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void TriMatDump(const TriMat * triMat, int h)
{
    FILE * fp = fopen("matrix_debug.txt", "a+");
    int row, col;
    #ifdef MAT_LEGAL_CHECKING
	if(triMat == NULL){
		return ;
	}
    #endif
    fprintf(fp, "TriMat %dx%d:\n", triMat->row, triMat->col);
    for (row = 0; row < triMat->row; row++)
    {
        for (col = 0; col < triMat->col; col++)
        {
            fprintf(fp, "%f\t", triMat->element[row][col][h]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void TriSMatDump(const TriSMat * triSMat, int h)
{
    FILE * fp = fopen("matrix_debug.txt", "a+");
    int row, col;
    #ifdef MAT_LEGAL_CHECKING
	if(triSMat == NULL){
		return ;
	}
    #endif
    fprintf(fp, "TriSMat %dx%d:\n", triSMat->row, triSMat->col);
    for (row = 0; row < triSMat->row; row++)
    {
        for (col = 0; col < triSMat->col; col++)
        {
            fprintf(fp, "%f\t", triSMat->element[row][col][h]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void VectorDump (const double * v, int n)
{
    int i = 0;
    FILE * fp = fopen("matrix_debug.txt", "a+");
    fprintf(fp, "Vector %d:\n", n);
    for (i = 0; i < n; i++)
    {
        fprintf(fp, "%f\t", v[i]);
    }
    fclose(fp);
}

Mat* MatZeros(Mat* mat)
{
	int row,col;

	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			mat->element[row][col] = 0.0;
		}
	}

	return mat;
}

SMat* SMatZeros(SMat* smat)
{
	int row,col;

	for(row = 0 ; row < smat->row ; row++){
		for(col = 0 ; col < smat->col ; col++){
			smat->element[row][col] = 0.0;
		}
	}

	return smat;
}

TriMat* TriMatZeros(TriMat *triMat)
{
    int row, col, h;
    for (row = 0; row < triMat->row; row++)
    {
        for (col = 0; col < triMat->col; col++)
        {
            for (h = 0; h < triMat->h; h++)
            {
                triMat->element[row][col][h] = 0.0;
            }
        }
    }
    return triMat;
}

TriSMat* TriSMatZeros(TriSMat *triSMat)
{
    int row, col, h;
    for (row = 0; row < triSMat->row; row++)
    {
        for (col = 0; col < triSMat->col; col++)
        {
            for (h = 0; h < triSMat->h; h++)
            {
                triSMat->element[row][col][h] = 0.0;
            }
        }
    }
    return triSMat;
}

TriMat* TriMatZerosOrder(TriMat *triMat)
{
    int row, col, h;
    for (row = 0; h < triMat->h; h++)
    {
        for (col = 0; row < triMat->row; row++)
        {
            for (h = 0; col < triMat->col; col++)
            {
                triMat->element[h][row][col] = 0.0;
            }
        }
    }
    return triMat;
}

TriSMat* TriSMatZerosOrder(TriSMat *triSMat)
{
    int row, col, h;
    for (row = 0; h < triSMat->h; h++)
    {
        for (col = 0; row < triSMat->row; row++)
        {
            for (h = 0; col < triSMat->col; col++)
            {
                triSMat->element[h][row][col] = 0.0;
            }
        }
    }
    return triSMat;
}

double *zeros_vector(int n)
{
    double *v;
    int i = 0;
    v = (double *)malloc(n * sizeof(double));
    memset(v, 0, n*sizeof(double));
    return v;
}

Mat* MatEye(Mat* mat)
{
	int i;

	MatZeros(mat);
	for(i = 0 ; i < min(mat->row, mat->col) ; i++){
		mat->element[i][i] = 1.0;
	}

	return mat;
}

SMat* SMatEye(SMat* smat)
{
	int i;

	SMatZeros(smat);
	for(i = 0 ; i < min(smat->row, smat->col) ; i++){
		smat->element[i][i] = 1.0;
	}

	return smat;
}


Mat* MatAdd(Mat* src1, Mat* src2, Mat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst->row && src1->col == src2->col && src2->col == dst->col) ){
		printf("err check, unmatch matrix for MatAdd\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		#pragma simd
		for(col = 0 ; col < src1->col ; col++){
			dst->element[row][col] = src1->element[row][col] + src2->element[row][col];
		}
	}

	return dst;
}

SMat* SMatAdd(SMat* src1, SMat* src2, SMat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst->row && src1->col == src2->col && src2->col == dst->col) ){
		printf("err check, unmatch matrix for SMatAdd\n");
		SMatDump(src1);
		SMatDump(src2);
		SMatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		#pragma simd
		for(col = 0 ; col < src1->col ; col++){
			dst->element[row][col] = src1->element[row][col] + src2->element[row][col];
		}
	}

	return dst;
}

// dst = src1 - src2
Mat* MatSub(Mat* src1, Mat* src2, Mat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst->row && src1->col == src2->col && src2->col == dst->col) ){
		printf("err check, unmatch matrix for MatSub\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		for(col = 0 ; col < src1->col ; col++){
			dst->element[row][col] = src1->element[row][col] - src2->element[row][col];
		}
	}

	return dst;
}

SMat* SMatSub(SMat* src1, SMat* src2, SMat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst->row && src1->col == src2->col && src2->col == dst->col) ){
		printf("err check, unmatch matrix for SMatSub\n");
		SMatDump(src1);
		SMatDump(src2);
		SMatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		for(col = 0 ; col < src1->col ; col++){
			dst->element[row][col] = src1->element[row][col] - src2->element[row][col];
		}
	}

	return dst;
}

/* dst = src1 * src2 */
Mat* MatMul(Mat* src1, Mat* src2, Mat* dst)
{
	int row, col;
	int i;
	double temp;

#ifdef MAT_LEGAL_CHECKING
	if( src1->col != src2->row || src1->row != dst->row || src2->col != dst->col ){
		printf("err check, unmatch matrix for MatMul\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < dst->row ; row++){
		for(col = 0 ; col < dst->col ; col++){
			temp = 0.0;
			for(i = 0 ; i < src1->col ; i++){
				temp += src1->element[row][i] * src2->element[i][col];
			}
			dst->element[row][col] = temp;
		}
	}

	return dst;
}

SMat* SMatMul(SMat* src1, SMat* src2, SMat* dst)
{
	int row, col;
	int i;
	float temp;

#ifdef MAT_LEGAL_CHECKING
	if( src1->col != src2->row || src1->row != dst->row || src2->col != dst->col ){
		printf("err check, unmatch matrix for SMatMul\n");
		SMatDump(src1);
		SMatDump(src2);
		SMatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < dst->row ; row++){
		for(col = 0 ; col < dst->col ; col++){
			temp = 0.0;
			for(i = 0 ; i < src1->col ; i++){
				temp += src1->element[row][i] * src2->element[i][col];
			}
			dst->element[row][col] = temp;
		}
	}

	return dst;
}

/* dst = src' */
Mat* MatTrans(Mat* src, Mat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != dst->col || src->col != dst->row ){
		printf("err check, unmatch matrix for MatTranspose\n");
		MatDump(src);
		MatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < dst->row ; row++){
		for(col = 0 ; col < dst->col ; col++){
			dst->element[row][col] = src->element[col][row];
		}
	}

	return dst;
}

SMat* SMatTrans(SMat* src, SMat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != dst->col || src->col != dst->row ){
		printf("err check, unmatch matrix for SMatTranspose\n");
		SMatDump(src);
		SMatDump(dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < dst->row ; row++){
		for(col = 0 ; col < dst->col ; col++){
			dst->element[row][col] = src->element[col][row];
		}
	}

	return dst;
}

void MatCopy(Mat* src, Mat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != dst->row || src->col != dst->col){
		printf("err check, unmathed matrix for MatCopy\n");
		MatDump(src);
		MatDump(dst);
		return ;
	}
#endif

	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++)
			dst->element[row][col] = src->element[row][col];
	}
}

void SMatCopy(SMat* src, SMat* dst)
{
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != dst->row || src->col != dst->col){
		printf("err check, unmathed matrix for SMatCopy\n");
		SMatDump(src);
		SMatDump(dst);
		return ;
	}
#endif

	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++)
			dst->element[row][col] = src->element[row][col];
	}
}

void MatNumMult(Mat * src, double k)
{
    int row, col;
    for (row = 0;  row < src->row; row++)
    {
		cblas_dscal(src->col,k,&(src->element[row][0]),1);
        //for (col = 0; col < src->col; col++)
            //src->element[row][col] = k * src->element[row][col];
    }
}

void SMatNumMult(SMat * src, float k)
{
    int row, col;
    for (row = 0;  row < src->row; row++)
    {
		cblas_sscal(src->col,k,&(src->element[row][0]),1);
        //for (col = 0; col < src->col; col++)
        //    src->element[row][col] = k * src->element[row][col];
    }
}

void TriMatNumMult(TriMat * src, double k)
{
	int row, col, h;
    for (row = 0; row < src->row; row++)
    {
        for (col = 0; col < src->col; col++)
        {
            for (h = 0; h < src->h; h++)
                src->element[row][col][h] = k * src->element[row][col][h];
        }
    }
}

void TriSMatNumMult(TriSMat * src, float k)
{
    int row, col, h;
    for (row = 0; row < src->row; row++)
    {
        for (col = 0; col < src->col; col++)
        {
            for (h = 0; h < src->h; h++)
                src->element[row][col][h] = k * src->element[row][col][h];
        }
    }
}

void TriMatNumMultOrder(TriMat * src, double k)
{
	int row, col, h;
	for (h = 0; h < src->h; h++)
	{
		for (row = 0; row < src->row; row++)
		{
			cblas_dscal(src->col,k,&(src->element[h][row][0]),1);
			//for (col = 0; col < src->col; col++)
			//{
			//		src->element[h][row][col] = k * src->element[h][row][col];
			//}
		}
	}
}

void TriSMatNumMultOrder(TriSMat * src, float k)
{
	int row, col, h;
	for (h = 0; h < src->h; h++)
	{
		for (row = 0; row < src->row; row++)
		{
			cblas_sscal(src->col,k,&(src->element[h][row][0]),1);
			//for (col = 0; col < src->col; col++)
			//{
			//		src->element[h][row][col] = k * src->element[h][row][col];
			//}
		}
	}
}