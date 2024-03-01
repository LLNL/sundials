#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "test_utilities.hpp"

#define ZERO SUN_RCONST(0.0) /* real 0.0     */
#define ONE  SUN_RCONST(1.0) /* real 1.0     */
#define PI   3.14159265358979323846

/* -----------------------------------------------------------------------------
 * Build Adams Nordsieck array from f(t,y) history
 * ---------------------------------------------------------------------------*/

int BuildNordsieckArrayAdams(sunrealtype* t, N_Vector y, N_Vector* f,
                             N_Vector* wrk, int order, sunrealtype hscale,
                             N_Vector* zn)
{
  /* Check for valid inputs */
  if (!t || !y || !f || !wrk || order < 1 || !zn) { return CV_ILL_INPUT; }

  for (int i = 0; i < order; i++)
  {
    if (!f[i]) { return CV_ILL_INPUT; }
    if (!wrk[i]) { return CV_ILL_INPUT; }
  }

  /* Compute Nordsieck array */
  if (order > 1)
  {
    /* Compute Newton polynomial coefficients interpolating f history */
    for (int i = 0; i < order; i++) { N_VScale(ONE, f[i], wrk[i]); }

    for (int i = 1; i < order; i++)
    {
      for (int j = order - 1; j >= i; j--)
      {
        /* Divided difference */
        sunrealtype delta_t = ONE / (t[j - i] - t[j]);
        N_VLinearSum(delta_t, wrk[j - 1], -delta_t, wrk[j], wrk[j]);
      }
    }

    for (int i = 0; i < order; i++)
    {
      fprintf(stdout, "wrk[%d]\n", i);
      N_VPrintFile(wrk[i], stdout);
    }

    /* Compute derivatives of Newton polynomial of f history */
    N_VScale(ONE, wrk[order - 1], zn[1]);
    for (int i = 2; i <= order; i++) { N_VConst(ZERO, zn[i]); }

    for (int i = order - 2; i >= 0; i--)
    {
      for (int j = order - 1; j > 0; j--)
      {
        N_VLinearSum(t[0] - t[i], zn[j + 1], j, zn[j], zn[j + 1]);
      }
      N_VLinearSum(t[0] - t[i], zn[1], ONE, wrk[i], zn[1]);
    }
  }

  /* Overwrite first two columns with input values */
  N_VScale(ONE, y, zn[0]);
  N_VScale(ONE, f[0], zn[1]);

  /* Scale entries */
  sunrealtype scale = ONE;
  for (int i = 1; i <= order; i++)
  {
    scale *= hscale / ((sunrealtype)i);
    N_VScale(scale, zn[i], zn[i]);
  }

  fprintf(stdout, "\n");
  for (int i = 0; i <= order; i++)
  {
    fprintf(stdout, "zn[%d]\n", i);
    N_VPrintFile(zn[i], stdout);
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Function to build BDF Nordsieck array
 * ---------------------------------------------------------------------------*/

int BuildNordsieckArrayBDF(sunrealtype* t, N_Vector* y, N_Vector f, N_Vector* wrk,
                           int order, sunrealtype hscale, N_Vector* zn)
{
  /* Check for valid inputs */
  if (!t || !y || !f || !wrk || order < 1 || !zn) { return CV_ILL_INPUT; }

  for (int i = 0; i < order; i++)
  {
    if (!y[i]) { return CV_ILL_INPUT; }
  }

  for (int i = 0; i < order + 1; i++)
  {
    if (!wrk[i]) { return CV_ILL_INPUT; }
  }

  if (order > 1)
  {
    /* Setup extended array of times to incorporate derivative value */
    sunrealtype t_ext[6];

    t_ext[0] = t[0];
    for (int i = 1; i <= order; i++) { t_ext[i] = t[i - 1]; }

    /* Compute Hermite polynomial coefficients interpolating y history and f */
    N_VScale(ONE, y[0], wrk[0]);
    for (int i = 1; i <= order; i++) { N_VScale(ONE, y[i - 1], wrk[i]); }

    for (int i = 1; i <= order; i++)
    {
      for (int j = order; j > i - 1; j--)
      {
        if (i == 1 && j == 1)
        {
          /* Replace with actual derivative value */
          N_VScale(ONE, f, wrk[j]);
        }
        else
        {
          /* Divided difference */
          sunrealtype delta_t = ONE / (t_ext[j - i] - t_ext[j]);
          N_VLinearSum(delta_t, wrk[j - 1], -delta_t, wrk[j], wrk[j]);
        }
        fprintf(stdout, "i = %d, j = %d, wrk[%d]\n", i, j, j);
        N_VPrintFile(wrk[j], stdout);
      }
    }

    for (int i = 0; i <= order; i++)
    {
      fprintf(stdout, "wrk[%d]\n", i);
      N_VPrintFile(wrk[i], stdout);
    }

    /* Compute derivatives of Hermite polynomial */
    N_VScale(ONE, wrk[order], zn[0]);
    for (int i = 1; i <= order; i++) { N_VConst(ZERO, zn[i]); }

    for (int i = order - 1; i >= 0; i--)
    {
      for (int j = order; j > 0; j--)
      {
        N_VLinearSum(t_ext[0] - t_ext[i], zn[j], j, zn[j - 1], zn[j]);
      }
      N_VLinearSum(t_ext[0] - t_ext[i], zn[0], ONE, wrk[i], zn[0]);
    }
  }

  /* Overwrite first two columns with input values */
  N_VScale(ONE, y[0], zn[0]);
  N_VScale(ONE, f, zn[1]);

  /* Scale entries */
  sunrealtype scale = ONE;
  for (int i = 1; i <= order; i++)
  {
    scale *= hscale / ((sunrealtype)i);
    N_VScale(scale, zn[i], zn[i]);
  }

  fprintf(stdout, "\n");
  for (int i = 0; i <= order; i++)
  {
    fprintf(stdout, "zn[%d]\n", i);
    N_VPrintFile(zn[i], stdout);
  }

  return CV_SUCCESS;
}

// Test main
int main(int argc, char* argv[])
{
  // Set output formatting
  std::cout << std::scientific;
  // std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::setprecision(16);
  std::cout << std::endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create template vector
  N_Vector tmp = N_VNew_Serial(1, sunctx);
  if (check_ptr(tmp, "N_VNew_Serial")) return 1;

  // -------------------
  // Create history data
  // -------------------

  const int n_hist = 12;
  sunrealtype t_hist[n_hist];
  sunrealtype dt = 1;

  for (int i = 0; i < n_hist; i++) { t_hist[i] = 1 - i * dt; }

  N_Vector* y = N_VCloneVectorArray(n_hist, tmp);
  if (check_ptr(y, "N_VCloneVectorArray")) return 1;

  N_Vector* yp = N_VCloneVectorArray(n_hist, tmp);
  if (check_ptr(yp, "N_VClone")) return 1;

  std::cout << "y_hist:" << std::endl;
  for (int i = 0; i < n_hist; i++)
  {
    sunrealtype* data = N_VGetArrayPointer(y[i]);
    data[0]           = std::sin(t_hist[i]);
    std::cout << i << ": t = " << t_hist[i] << " y = " << std::endl;
    N_VPrint(y[i]);
  }
  std::cout << std::endl;
  for (int i = 0; i < n_hist; i++)
  {
    sunrealtype* ypdata = N_VGetArrayPointer(yp[i]);
    ypdata[0]           = std::cos(t_hist[i]);
    std::cout << i << ": t = " << t_hist[i] << " y' = " << std::endl;
    N_VPrint(yp[i]);
  }
  std::cout << std::endl;

  // -------------------------------
  // Create interpolating polynomial
  // -------------------------------

  N_Vector* wrk = N_VCloneVectorArray(13, tmp);
  if (check_ptr(wrk, "N_VCloneVectorArray")) return 1;

  N_Vector* zn = N_VCloneVectorArray(13, tmp);
  if (check_ptr(zn, "N_VCloneVectorArray")) return 1;

  // Generate coefficients
  std::cout << std::endl << "BDF" << std::endl;
  int order = 5;
  int flag  = BuildNordsieckArrayBDF(t_hist, y, yp[0], wrk, order, dt, zn);
  if (flag) return 1;

  std::cout << std::endl << "ADAMS" << std::endl;
  order = 12;
  flag  = BuildNordsieckArrayAdams(t_hist, y[0], yp, wrk, order, dt, zn);
  if (flag) return 1;

  N_VDestroy(tmp);
  N_VDestroyVectorArray(y, n_hist);
  N_VDestroyVectorArray(yp, n_hist);
  N_VDestroyVectorArray(wrk, 13);
  N_VDestroyVectorArray(zn, 13);

  return 0;
}
