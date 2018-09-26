/*
 * File: xrot.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 26-Sep-2018 10:17:51
 */

#ifndef XROT_H
#define XROT_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "schur_decompose_types.h"

/* Function Declarations */
extern void b_xrot(int n, double x[9], int ix0, int iy0, double c, double s);
extern void xrot(int n, double x[9], int ix0, int iy0, double c, double s);

#endif

/*
 * File trailer for xrot.h
 *
 * [EOF]
 */
