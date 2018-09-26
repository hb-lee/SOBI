/*
 * File: xungorghr.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 26-Sep-2018 10:17:51
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "schur_decompose.h"
#include "xungorghr.h"
#include "xgerc.h"

/* Function Definitions */

/*
 * Arguments    : int n
 *                int ilo
 *                int ihi
 *                double A[9]
 *                int ia0
 *                const double tau[2]
 *                int itau0
 * Return Type  : void
 */
void xungorghr(int n, int ilo, int ihi, double A[9], int ia0, const double tau[2],
               int itau0)
{
  int nh;
  int i;
  int ia;
  int b_i;
  int itau;
  int b_ia;
  double work[3];
  int iaii;
  int i3;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int iac;
  int exitg1;
  int ix;
  double c;
  int i4;
  if (n == 0) {
  } else {
    nh = ihi - ilo;
    for (i = ihi; i >= ilo + 1; i--) {
      ia = (ia0 + (i - 1) * 3) - 2;
      for (b_i = 1; b_i < i; b_i++) {
        A[ia + b_i] = 0.0;
      }

      for (b_i = i + 1; b_i <= ihi; b_i++) {
        A[ia + b_i] = A[(ia + b_i) - 3];
      }

      for (b_i = ihi + 1; b_i <= n; b_i++) {
        A[ia + b_i] = 0.0;
      }
    }

    for (i = 0; i + 1 <= ilo; i++) {
      ia = (ia0 + i * 3) - 1;
      for (b_i = 1; b_i <= n; b_i++) {
        A[(ia + b_i) - 1] = 0.0;
      }

      A[ia + i] = 1.0;
    }

    for (i = ihi; i + 1 <= n; i++) {
      ia = (ia0 + i * 3) - 1;
      for (b_i = 1; b_i <= n; b_i++) {
        A[(ia + b_i) - 1] = 0.0;
      }

      A[ia + i] = 1.0;
    }

    ia = (ia0 + ilo) + ilo * 3;
    if (nh < 1) {
    } else {
      for (i = nh; i < nh; i++) {
        b_ia = ia + i * 3;
        for (b_i = 0; b_i < nh; b_i++) {
          A[(b_ia + b_i) - 1] = 0.0;
        }

        A[(b_ia + i) - 1] = 1.0;
      }

      itau = ((itau0 + ilo) + nh) - 3;
      for (b_i = 0; b_i < 3; b_i++) {
        work[b_i] = 0.0;
      }

      for (b_i = nh; b_i >= 1; b_i--) {
        iaii = ((ia + b_i) + (b_i - 1) * 3) - 1;
        if (b_i < nh) {
          A[iaii - 1] = 1.0;
          i = nh - b_i;
          if (tau[itau] != 0.0) {
            lastv = i + 1;
            i += iaii;
            while ((lastv > 0) && (A[i - 1] == 0.0)) {
              lastv--;
              i--;
            }

            lastc = nh - b_i;
            exitg2 = false;
            while ((!exitg2) && (lastc > 0)) {
              i = (iaii + (lastc - 1) * 3) + 3;
              b_ia = i;
              do {
                exitg1 = 0;
                if (b_ia <= (i + lastv) - 1) {
                  if (A[b_ia - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    b_ia++;
                  }
                } else {
                  lastc--;
                  exitg1 = 2;
                }
              } while (exitg1 == 0);

              if (exitg1 == 1) {
                exitg2 = true;
              }
            }
          } else {
            lastv = 0;
            lastc = 0;
          }

          if (lastv > 0) {
            if (lastc == 0) {
            } else {
              for (i = 1; i <= lastc; i++) {
                work[i - 1] = 0.0;
              }

              i = 0;
              i3 = (iaii + 3 * (lastc - 1)) + 3;
              for (iac = iaii + 3; iac <= i3; iac += 3) {
                ix = iaii;
                c = 0.0;
                i4 = (iac + lastv) - 1;
                for (b_ia = iac; b_ia <= i4; b_ia++) {
                  c += A[b_ia - 1] * A[ix - 1];
                  ix++;
                }

                work[i] += c;
                i++;
              }
            }

            xgerc(lastv, lastc, -tau[itau], iaii, work, A, iaii + 3);
          }
        }

        if (b_i < nh) {
          i3 = (iaii + nh) - b_i;
          for (i = iaii; i + 1 <= i3; i++) {
            A[i] *= -tau[itau];
          }
        }

        A[iaii - 1] = 1.0 - tau[itau];
        for (i = 1; i < b_i; i++) {
          A[(iaii - i) - 1] = 0.0;
        }

        itau--;
      }
    }
  }
}

/*
 * File trailer for xungorghr.c
 *
 * [EOF]
 */
