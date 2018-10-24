#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
typedef struct  {
	int row, col;
	float **element;
}SMat;

typedef struct  {
	int row, col;
	double **element;
}Mat;

typedef struct {
    int row, col, h;
    float ***element;
}TriSMat;

typedef struct {
    int row, col, h;
    double ***element;
}TriMat;

Mat * ConstructMat();
SMat * ConstructSMat();
TriMat * ConstructTriMat();
TriSMat * ConstructTriSMat();

TriMat* TriMatCreate(TriMat * triMat, int row, int col, int h);
TriSMat* TriSMatCreate(TriSMat * triSMat, int row, int col, int h);
TriMat* TriMatOrderCreate(TriMat * triMat, int h, int row, int col);
TriSMat* TriSMatOrderCreate(TriSMat * triSMat, int h, int row, int col);
Mat* MatCreate(Mat* mat, int row, int col);
SMat* SMatCreate(SMat* smat, int row, int col);
void MatZeroConstruct(Mat * mat, int row, int col);
void SMatZeroConstruct(SMat * smat, int row, int col);
void MatEyeConstruct(Mat * mat, int row, int col);
void SMatEyeConstruct(SMat * smat, int row, int col);
void MatDelete(Mat* mat);
void SMatDelete(SMat* smat);
void MatDump(const Mat* mat);
void SMatDump(const SMat* smat);
void TriMatDump(const TriMat * triMat, int h);
void TriSMatDump(const TriSMat * triSMat, int h);
void VectorDump (const double * v, int n);
Mat* MatZeros(Mat* mat);
SMat* SMatZeros(SMat* smat);
TriMat* TriMatZeros(TriMat *triMat);
TriMat* TriMatZerosOrder(TriMat *triMat);
TriSMat* TriSMatZerosOrder(TriSMat *triSMat);
TriSMat* TriSMatZeros(TriSMat *triSMat);
double *zeros_vector(int n);
Mat* MatEye(Mat* mat);
SMat* SMatEye(SMat* smat);
Mat* MatAdd(Mat* src1, Mat* src2, Mat* dst);
SMat* SMatAdd(SMat* src1, SMat* src2, SMat* dst);
Mat* MatMul(Mat* src1, Mat* src2, Mat* dst);
SMat* SMatMul(SMat* src1, SMat* src2, SMat* dst);
Mat* MatTrans(Mat* src, Mat* dst);
SMat* SMatTrans(SMat* src, SMat* dst);
void MatCopy(Mat* src, Mat* dst);
void SMatCopy(SMat* src, SMat* dst);
void MatNumMult(Mat * src, double k);
void SMatNumMult(SMat * src, float k);
void TriMatNumMult(TriMat * src, double k);
void TriSMatNumMult(TriSMat * src, float k);
void TriMatNumMultOrder(TriMat * src, double k);
void TriSMatNumMultOrder(TriSMat * src, float k);

#endif // MATRIX_H_INCLUDED
