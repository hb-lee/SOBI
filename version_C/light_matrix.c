#include "light_matrix.h"
#include <stdio.h>
#include <stdlib.h>

#define MAT_LEGAL_CHECKING

#define min(a, b) ((a) > (b) ? (b) : (a))
#define equal(a, b)	((a-b)<1e-7 && (a-b)>-(1e-7))

/************************************************************************/
/*                          Private Function                            */
/************************************************************************/

void swap(int *a, int *b)
{
	int m;
	m = *a;
	*a = *b;
	*b = m;
}

void perm(int list[], int k, int m, int* p, Mat* mat, double* det)
{
	int i;

	if(k > m){
		double res = mat->element[0][list[0]];

		for(i = 1; i < mat->row ; i++){
			res *= mat->element[i][list[i]];
		}

		if(*p%2){
			//odd is negative
			*det -= res;
		}else{
			//even is positive
			*det += res;
		}
	}
	else{
		// if the element is 0, we don't need to calculate the value for this permutation
		if(!equal(mat->element[k][list[k]], (double)0.0))
			perm(list, k + 1, m, p, mat, det);
		for(i = k+1; i <= m; i++)
		{
			if(equal(mat->element[k][list[i]], (double)0.0))
				continue;
			swap(&list[k], &list[i]);
			*p += 1;
			perm(list, k + 1, m, p, mat, det);
			swap(&list[k], &list[i]);
			*p -= 1;
		}
	}
}

/************************************************************************/
/*                           Public Function                            */
/************************************************************************/

Mat * ConstructMat()
{
    Mat * mat;
    mat = (Mat *)malloc(sizeof(Mat));
    return mat;
}

TriMat * ConstructTriMat()
{
    TriMat * triMat;
    triMat = (TriMat *) malloc(sizeof(TriMat));
    return triMat;
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
            triMat->element[i][j] = (double *)malloc(h * sizeof(double));
        }
    }
    triMat->row = row;
    triMat->col = col;
    triMat->h = h;
    TriMatZeros(triMat);
    return triMat;
}

// ´´½¨¾ØÕó
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
		mat->element[i] = (double*)malloc(col * sizeof(double));
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
    mat = MatZeros(mat);
	return mat;
}

// ¾ØÕóÏú»Ù
void MatDelete(Mat* mat)
{
	int i;

	for(i = 0 ; i<mat->row ; i++)
		free(mat->element[i]);
	free(mat->element);
}
/*
// ¸³Öµ¾ØÕóµÄÖµ
Mat* MatSetVal(Mat* mat, double* val)
{
	int row,col;

	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			mat->element[row][col] = val[col + row * mat->col];
		}
	}

	return mat;
}
*/
// Êä³ö¾ØÕóµÄÖµ
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
            fprintf(fp, "%0.10lf\t", mat->element[row][col]);
			//printf("%0.10lf\t", mat->element[row][col]);
		}
		fprintf(fp, "\n");
		//printf("\n");
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
            fprintf(fp, "%0.10lf\t", triMat->element[row][col][h]);
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
        fprintf(fp, "%0.10lf\t", v[i]);
    }
    fclose(fp);
}

// ½«¾ØÕóÖÃÁã
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

// ½«¾ØÕó¶Ô½ÇÔªËØÖÃ1
Mat* MatEye(Mat* mat)
{
	int i;

	MatZeros(mat);
	for(i = 0 ; i < min(mat->row, mat->col) ; i++){
		mat->element[i][i] = 1.0;
	}

	return mat;
}

// ¾ØÕóÏà¼Ó
/* dst = src1 + src2 */
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
		for(col = 0 ; col < src1->col ; col++){
			dst->element[row][col] = src1->element[row][col] + src2->element[row][col];
		}
	}

	return dst;
}
/*
// ¾ØÕó¼õ·¨
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
*/
// ¾ØÕó³Ë·¨ M * N
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
/*
Mat* PartMatMul(Mat* src1, Mat* src2, Mat* dst, int src1_r_begin, int src2_r_begin, int src1_c_begin, int src2_c_begin, int mid)
{
    int row, col;
	int i;
	double temp;

	for(row = 0 ; row < dst->row ; row++){
		for(col = 0 ; col < dst->col ; col++){
			temp = 0.0;
			for(i = 0 ; i < mid ; i++){
				temp = temp + src1->element[src1_r_begin + row][src1_c_begin + i] * src2->element[src2_r_begin + i][col + src2_c_begin];
			}
			dst->element[row][col] = temp;
		}
	}

	return dst;
}
*/
// ¾ØÕóµÄ×ªÖÃ
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
/*
// ¾ØÕóÐÐÁÐÊ½
// return det(mat)
double MatDet(Mat* mat)
{
	double det = 0.0;
	int plarity = 0;
	int *list;
	int i;

#ifdef MAT_LEGAL_CHECKING
	if( mat->row != mat->col){
		printf("err check, not a square matrix for MatDetermine\n");
		MatDump(mat);
		return 0.0;
	}
#endif

	list = (int*)malloc(sizeof(int)*mat->col);
	if(list == NULL){
		printf("malloc list fail\n");
		return 0.0;
	}
	for(i = 0 ; i < mat->col ; i++)
		list[i] = i;

	perm(list, 0, mat->row-1, &plarity, mat, &det);
	free(list);

	return det;
}
*/
/*
// °éËæ¾ØÕó
// dst = adj(src)
Mat* MatAdj(Mat* src, Mat* dst)
{
	Mat smat;
	int row, col;
	int i,j,r,c;
	double det;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != src->col || src->row != dst->row || src->col != dst->col){
		printf("err check, not a square matrix for MatAdj\n");
		MatDump(src);
		MatDump(dst);
		return NULL;
	}
#endif

	MatCreate(&smat, src->row-1, src->col-1);

	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++){
			r = 0;
			for(i = 0 ; i < src->row ; i++){
				if(i == row)
					continue;
				c = 0;
				for(j = 0; j < src->col ; j++){
					if(j == col)
						continue;
					smat.element[r][c] = src->element[i][j];
					c++;
				}
				r++;
			}
			det = MatDet(&smat);
			if((row+col)%2)
				det = -det;
			dst->element[col][row] = det;
		}
	}

	MatDelete(&smat);

	return dst;
}
*/
/*
// ¾ØÕóµÄÄæ
// dst = src^(-1)
Mat* MatInv(Mat* src, Mat* dst)
{
	Mat adj_mat;
	double det;
	int row, col;

#ifdef MAT_LEGAL_CHECKING
	if( src->row != src->col || src->row != dst->row || src->col != dst->col){
		printf("err check, not a square matrix for MatInv\n");
		MatDump(src);
		MatDump(dst);
		return NULL;
	}
#endif
	MatCreate(&adj_mat, src->row, src->col);
	MatAdj(src, &adj_mat);
	det = MatDet(src);

	if(equal(det, (double)0.0)){
		printf("err, determinate is 0 for MatInv\n");
		return NULL;
	}

	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++)
			dst->element[row][col] = adj_mat.element[row][col]/det;
	}

	MatDelete(&adj_mat);

	return dst;
}
*/
// ¾ØÕó¿½±´
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

void MatNumMult(Mat * src, double k)
{
    int row, col;
    for (row = 0;  row < src->row; row++)
    {
        for (col = 0; col < src->col; col++)
            src->element[row][col] = k * src->element[row][col];
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

void MatZeroConstruct(Mat * mat, int row, int col)
{
	int i, j;

	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return NULL;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double*)malloc(col * sizeof(double));
		for (j = 0; j < col; j++)
            mat->element[i][j] = 0.0;
	}

	mat->row = row;
	mat->col = col;
}

void MatEyeConstruct(Mat * mat, int row, int col)
{
    int i, j;

	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return NULL;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double*)malloc(col * sizeof(double));
		for (j = 0; j < col; j++)
            if (i == j)
                mat->element[i][j] = 1;
            else
                mat->element[i][j] = 0.0;
	}

	mat->row = row;
	mat->col = col;
}
