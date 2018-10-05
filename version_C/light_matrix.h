#ifndef LIGHT_MATRIX_H_INCLUDED
#define LIGHT_MATRIX_H_INCLUDED
typedef struct  {
	int row, col;
	double **element;
}Mat;

typedef struct {
    int row, col, h;
    double ***element;
}TriMat;

Mat * ConstructMat();
TriMat * ConstructTriMat();

TriMat* TriMatCreate(TriMat * triMat, int row, int col, int h);
Mat* MatCreate(Mat* mat, int row, int col);
void MatDelete(Mat* mat);
Mat* MatSetVal(Mat* mat, double* val);
void MatDump(const Mat* mat);

Mat* MatZeros(Mat* mat);
TriMat* TriMatZeros(TriMat *triMat);
double *zeros_vector(int n);
Mat* MatEye(Mat* mat);

Mat* MatAdd(Mat* src1, Mat* src2, Mat* dst);
Mat* MatSub(Mat* src1, Mat* src2, Mat* dst);
Mat* MatMul(Mat* src1, Mat* src2, Mat* dst);
Mat* PartMatMul(Mat* src1, Mat* src2, Mat* dst, int src1_r_begin, int src2_r_begin, int src1_c_begin, int src2_c_begin, int mid);
Mat* MatTrans(Mat* src, Mat* dst);
double MatDet(Mat* mat);
Mat* MatAdj(Mat* src, Mat* dst);
Mat* MatInv(Mat* src, Mat* dst);

void MatCopy(Mat* src, Mat* dst);
void MatNumMult(Mat * src, double k);
void TriMatNumMult(TriMat * src, double k);
void MatZeroConstruct(Mat * mat, int row, int col);
void MatEyeConstruct(Mat * mat, int row, int col);


#endif // LIGHT_MATRIX_H_INCLUDED
