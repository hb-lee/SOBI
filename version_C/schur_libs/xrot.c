/*
 * File: xrot.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 26-Sep-2018 10:17:51
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "schur_decompose.h"
#include "xrot.h"

/* Function Definitions */

/*
 * Arguments    : int n
 *                double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
void b_xrot(int n, double x[9], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  if (n < 1) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 1; k <= n; k++) {
      temp = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = temp;
      iy++;
      ix++;
    }
  }
}

/*
 * Arguments    : int n
 *                double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
void xrot(int n, double x[9], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  if (n < 1) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 1; k <= n; k++) {
      temp = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = temp;
      iy += 3;
      ix += 3;
    }
  }
}

/*
 * File trailer for xrot.c
 *
 * [EOF]
 */
