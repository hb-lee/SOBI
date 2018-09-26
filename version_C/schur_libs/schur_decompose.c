/*
 * File: schur_decompose.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 26-Sep-2018 10:17:51
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "schur_decompose.h"
#include "xhseqr.h"
#include "xungorghr.h"
#include "xgehrd.h"

/* Function Definitions */

/*
 * Arguments    : const double W[9]
 *                double U[9]
 *                double T[9]
 * Return Type  : void
 */
void schur_decompose(const double W[9], double U[9], double T[9])
{
  double tau[2];
  memcpy(&T[0], &W[0], 9U * sizeof(double));
  xgehrd(T, tau);
  memcpy(&U[0], &T[0], 9U * sizeof(double));
  xungorghr(3, 1, 3, U, 1, tau, 1);
  xhseqr(T, U);
}

/*
 * File trailer for schur_decompose.c
 *
 * [EOF]
 */
