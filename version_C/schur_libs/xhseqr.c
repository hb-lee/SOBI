/*
 * File: xhseqr.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 26-Sep-2018 10:17:51
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "schur_decompose.h"
#include "xhseqr.h"
#include "xrot.h"
#include "xdlanv2.h"
#include "xzlarfg.h"

/* Function Definitions */

/*
 * Arguments    : double h[9]
 *                double z[9]
 * Return Type  : int
 */
int xhseqr(double h[9], double z[9])
{
  int info;
  double v[3];
  int i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  double tst;
  double htmp1;
  double aa;
  double ab;
  double ba;
  double h22;
  double s;
  double cs;
  double sn;
  boolean_T guard1 = false;
  int hoffset;
  int nr;
  double d1;
  int j;
  info = -1;
  h[2] = 0.0;
  i = 2;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 31)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && ((k + 1 > L) && (!(fabs(h[k + 3 * (k - 1)]) <=
                3.0062525400134592E-292)))) {
        tst = fabs(h[(k + 3 * (k - 1)) - 1]) + fabs(h[k + 3 * k]);
        if (tst == 0.0) {
          if (k - 1 >= 1) {
            tst = fabs(h[(k + 3 * (k - 2)) - 1]);
          }

          if (k + 2 <= 3) {
            tst += fabs(h[(k + 3 * k) + 1]);
          }
        }

        guard1 = false;
        if (fabs(h[k + 3 * (k - 1)]) <= 2.2204460492503131E-16 * tst) {
          htmp1 = fabs(h[k + 3 * (k - 1)]);
          tst = fabs(h[(k + 3 * k) - 1]);
          if (htmp1 > tst) {
            ab = htmp1;
            ba = tst;
          } else {
            ab = tst;
            ba = htmp1;
          }

          htmp1 = fabs(h[k + 3 * k]);
          tst = fabs(h[(k + 3 * (k - 1)) - 1] - h[k + 3 * k]);
          if (htmp1 > tst) {
            aa = htmp1;
            htmp1 = tst;
          } else {
            aa = tst;
          }

          s = aa + ab;
          tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
          if ((3.0062525400134592E-292 >= tst) || rtIsNaN(tst)) {
            d1 = 3.0062525400134592E-292;
          } else {
            d1 = tst;
          }

          if (ba * (ab / s) <= d1) {
            exitg3 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          k--;
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + 3 * (k - 1)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          s = fabs(h[1]) + fabs(h[5]);
          tst = 0.75 * s + h[0];
          aa = -0.4375 * s;
          htmp1 = s;
          h22 = tst;
        } else if (its == 20) {
          s = fabs(h[i + 3 * (i - 1)]) + fabs(h[(i + 3 * (i - 2)) - 1]);
          tst = 0.75 * s + h[i + 3 * i];
          aa = -0.4375 * s;
          htmp1 = s;
          h22 = tst;
        } else {
          tst = h[(i + 3 * (i - 1)) - 1];
          htmp1 = h[i + 3 * (i - 1)];
          aa = h[(i + 3 * i) - 1];
          h22 = h[i + 3 * i];
        }

        s = ((fabs(tst) + fabs(aa)) + fabs(htmp1)) + fabs(h22);
        if (s == 0.0) {
          ba = 0.0;
          tst = 0.0;
          ab = 0.0;
          htmp1 = 0.0;
        } else {
          tst /= s;
          htmp1 /= s;
          aa /= s;
          h22 /= s;
          ab = (tst + h22) / 2.0;
          tst = (tst - ab) * (h22 - ab) - aa * htmp1;
          htmp1 = sqrt(fabs(tst));
          if (tst >= 0.0) {
            ba = ab * s;
            ab = ba;
            tst = htmp1 * s;
            htmp1 = -tst;
          } else {
            ba = ab + htmp1;
            ab -= htmp1;
            if (fabs(ba - h22) <= fabs(ab - h22)) {
              ba *= s;
              ab = ba;
            } else {
              ab *= s;
              ba = ab;
            }

            tst = 0.0;
            htmp1 = 0.0;
          }
        }

        if (i - 1 >= 1) {
          s = (fabs(h[0] - ab) + fabs(htmp1)) + fabs(h[1]);
          aa = h[1] / s;
          v[0] = (aa * h[3] + (h[0] - ba) * ((h[0] - ab) / s)) - tst * (htmp1 /
            s);
          v[1] = aa * (((h[0] + h[4]) - ba) - ab);
          v[2] = aa * h[5];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
        }

        for (k = i - 1; k <= i; k++) {
          nr = (i - k) + 2;
          if (3 <= nr) {
            nr = 3;
          }

          if (k > i - 1) {
            hoffset = k + 3 * (k - 2);
            for (j = 1; j <= nr; j++) {
              v[j - 1] = h[(j + hoffset) - 2];
            }
          }

          tst = v[0];
          h22 = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (k > i - 1) {
            h[k - 1] = tst;
            h[k] = 0.0;
            if (k < i) {
              h[2] = 0.0;
            }
          }

          tst = v[1];
          htmp1 = h22 * v[1];
          if (nr == 3) {
            ab = v[2];
            ba = h22 * v[2];
            for (j = k - 1; j + 1 < 4; j++) {
              aa = (h[(k + 3 * j) - 1] + tst * h[k + 3 * j]) + ab * h[(k + 3 * j)
                + 1];
              h[(k + 3 * j) - 1] -= aa * h22;
              h[k + 3 * j] -= aa * htmp1;
              h[2 + 3 * j] -= aa * ba;
            }

            if (k + 3 <= i + 1) {
              nr = k;
            } else {
              nr = i - 2;
            }

            for (j = 0; j + 1 <= nr + 3; j++) {
              aa = (h[j + 3 * (k - 1)] + tst * h[j + 3 * k]) + ab * h[j + 3 * (k
                + 1)];
              h[j + 3 * (k - 1)] -= aa * h22;
              h[j + 3 * k] -= aa * htmp1;
              h[6 + j] -= aa * ba;
            }

            for (j = 0; j < 3; j++) {
              aa = (z[j + 3 * (k - 1)] + tst * z[j + 3 * k]) + ab * z[j + 3 * (k
                + 1)];
              z[j + 3 * (k - 1)] -= aa * h22;
              z[j + 3 * k] -= aa * htmp1;
              z[6 + j] -= aa * ba;
            }
          } else {
            if (nr == 2) {
              for (j = k - 1; j + 1 < 4; j++) {
                aa = h[(k + 3 * j) - 1] + tst * h[k + 3 * j];
                h[(k + 3 * j) - 1] -= aa * h22;
                h[k + 3 * j] -= aa * htmp1;
              }

              for (j = 0; j + 1 <= i + 1; j++) {
                aa = h[j + 3 * (k - 1)] + tst * h[j + 3 * k];
                h[j + 3 * (k - 1)] -= aa * h22;
                h[j + 3 * k] -= aa * htmp1;
              }

              for (j = 0; j < 3; j++) {
                aa = z[j + 3 * (k - 1)] + tst * z[j + 3 * k];
                z[j + 3 * (k - 1)] -= aa * h22;
                z[j + 3 * k] -= aa * htmp1;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = i;
      exitg1 = true;
    } else {
      if ((L == i + 1) || (!(L == i))) {
      } else {
        tst = h[(i + 3 * i) - 1];
        htmp1 = h[i + 3 * (i - 1)];
        aa = h[i + 3 * i];
        xdlanv2(&h[(i + 3 * (i - 1)) - 1], &tst, &htmp1, &aa, &ab, &ba, &h22, &s,
                &cs, &sn);
        h[(i + 3 * i) - 1] = tst;
        h[i + 3 * (i - 1)] = htmp1;
        h[i + 3 * i] = aa;
        if (3 > i + 1) {
          xrot(2 - i, h, i + (i + 1) * 3, (i + (i + 1) * 3) + 1, cs, sn);
        }

        b_xrot(i - 1, h, 1 + (i - 1) * 3, 1 + i * 3, cs, sn);
        hoffset = (i - 1) * 3;
        nr = i * 3;
        for (k = 0; k < 3; k++) {
          tst = cs * z[hoffset] + sn * z[nr];
          z[nr] = cs * z[nr] - sn * z[hoffset];
          z[hoffset] = tst;
          nr++;
          hoffset++;
        }
      }

      i = L - 2;
    }
  }

  info++;
  return info;
}

/*
 * File trailer for xhseqr.c
 *
 * [EOF]
 */
