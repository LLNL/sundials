/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Newton Polynomial Test
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include <cvode/cvode_impl.h>

#include "sundials/sundials_nvector.h"
#include "test_pr.hpp"
#include "test_utilities.hpp"

#define PI 3.14159265358979323846

/*
 * Compute Hermite interpolating polynomial coefficients
 *
 * Inputs:
 *   t  -- array (length M) of time values t_i
 *   y  -- array (length M) of state vectors y_i
 *   yp -- derivative of y_0 = f(t_0, y_0)
 *   M  -- number of interpolation times
 *
 * Outputs:
 *   c -- array (length M + 1) of coefficient vectors c_i
 */
int HermitePolyCoef(sunrealtype* t, N_Vector* y, N_Vector f, int n_times,
                    N_Vector* coeff)
{
  int i, j;
  int n_coeff = n_times + 1; // n_times state values + 1 derivative
  sunrealtype* t_ext = NULL;

  /* Check for valid inputs */
  if (!t || !y || !f || n_times < 1 || !coeff) { return CV_ILL_INPUT; }

  /* Check for valid solution vectors */
  for (i = 0; i < n_times; i++)
  {
    if (!y[i]) { return CV_ILL_INPUT; }
  }

  /* Check for valid coefficient vector */
  for (i = 0; i < n_coeff; i++)
  {
    if (!coeff[i]) { return CV_ILL_INPUT; }
  }

  /* Setup extended array of times to incorporate derivative value */
  t_ext = (sunrealtype*) malloc(sizeof(sunrealtype) * (n_coeff));

  t_ext[0] = t[0];
  t_ext[1] = t[0];

  for (i = 1; i < n_times; i++)
  {
    t_ext[i + 1] = t[i];
  }

  std::cout << "t_ext" << std::endl;
  for (i = 0; i < n_coeff; i++)
  {
    std::cout << t_ext[i] << std::endl;
  }

  /* Initialize coefficient arrays with values to interpolate */
  N_VScale(ONE, y[0], coeff[0]);
  N_VScale(ONE, y[0], coeff[1]);

  for (i = 1; i < n_times; i++)
  {
    N_VScale(ONE, y[i], coeff[i + 1]);
  }

  std::cout << "coeff" << std::endl;
  for (i = 0; i < n_coeff; i++)
  {
    N_VPrint(coeff[i]);
  }
  std::cout << std::endl;

  /* Compute coefficients from bottom up to write in place */
  for (i = 1; i < n_coeff; i++)
  {
    for (j = n_coeff - 1; j > i - 1; j--)
    {
      /* Replace with actual derivative value */
      if (i == 1 && j == 1)
      {
        N_VScale(ONE, f, coeff[j]);
      }
      else
      {
        sunrealtype denom = ONE / (t_ext[j - 1] - t_ext[j]);
        N_VLinearSum(denom, coeff[j - 1], -denom, coeff[j], coeff[j]);
      }
      std::cout << "i, j " << i << "," << j << std::endl;
      N_VPrint(coeff[j]);
    }
  }

  free(t_ext);
  t_ext = NULL;

  return CV_SUCCESS;
}

/*
 * Evaluate the interpolating polynomial and its derivatives up to order d at a
 * given time
 *
 * Inputs:
 *
 * t -- array (length M) of interpolated time values, t_i
 * c -- array (length M + 1) of polynomial coefficient vectors (length N), c_i
 * s -- time at which to evaluate the polynomial, s
 *
 * Output:
 *
 * p -- array (length d + 1) of values p[0] = p(s), p[1] = p'[s], etc.
 *
 * Example:
 *
 * Consider the M = 3 case, we have the interpolation times t0, t1, t2 and
 * the coefficients c0, c1, c2, c3. The Hermite interpolating polynomial is
 *
 * p(s) = c0 + c1 (s - t0) + c2 (s - t0)(s - t0) + c3 (s - t0)(s - t0)(s - t1)
 *
 * This can be rewritten as
 *
 * p(s) = a0 + (s - t0) [ a1 + (s - t1) [ a2 + (s - t2) a3 ] ]
 *
 * We can then write this recursively as P{k} = P{k + 1} (s - t{k}) + a{k} for
 * k = M-1,...,0 where P{M - 1} = c{M - 1} and P0 = p
 *
 * Using this recursive definition we can also compute derivatives of the
 * polynomial as P^(d){k} = P^(d){k + 1} (t - t{k}) + d P^(d - 1){k + 1} where
 * P^(d){M - 1} = 0 for d > 0 and p^(d) = P^(d)0. For example, the first three
 * derivatives are given by
 *
 * P'{k}   = P'{k + 1}   (t - t{k}) +   P{k + 1}
 * P''{k}  = P''{k + 1}  (t - t{k}) + 2 P'{k + 1}
 * P'''{k} = P'''{k + 1} (t - t{k}) + 3 P''{k + 1}
 */

int HermitePolyMultiDerEval(sunrealtype* t, N_Vector* coeff, int n_times,
                            sunrealtype s, int derv, N_Vector* p)
{
  int i, j;
  int n_coeff = n_times + 1; // n_times state values + 1 derivative
  sunrealtype* t_ext = NULL;

  /* Check for valid inputs */
  if (!t || !coeff || n_times < 1 || derv < 0 || !p) { return CV_ILL_INPUT; }

  for (i = 0; i < n_coeff; i++)
  {
    if (!coeff[i]) return CV_ILL_INPUT;
  }

  for (i = 0; i < derv + 1; i++)
  {
    if (!p[i]) { return CV_ILL_INPUT; }
  }

  t_ext = (sunrealtype*) malloc(sizeof(sunrealtype) * (n_coeff));

  t_ext[0] = t[0];
  t_ext[1] = t[0];

  for (i = 1; i < n_times; i++)
  {
    t_ext[i + 1] = t[i];
  }

  std::cout << "t_ext" << std::endl;
  for (i = 0; i < n_coeff; i++)
  {
    std::cout << t_ext[i] << std::endl;
  }

  /* Initialize interpolation output to P_{n_times-1} = c_{n_times-1} and derivative output
     to zero */
  N_VScale(ONE, coeff[n_coeff - 1], p[0]);
  for (i = 1; i < derv + 1; i++)
  {
    N_VConst(ZERO, p[i]);
  }

  /* Accumulate polynomial terms in-place i.e., P_{i+1} is stored in p[0] and
     overwritten each iteration by P_{i} and similarly for P', P'', etc. */
  for (i = n_coeff - 2; i >= 0; i--)
  {
    for (j = derv; j > 0; j--)
    {
      // P^(j)_i = P_{j + 1} * (s - t_i) + j * P^(j - 1)_i
      N_VLinearSum(s - t_ext[i], p[j], j, p[j-1], p[j]);
    }
    // P_i = P_{i + 1} * (s - t_i) + c_i
    N_VLinearSum(s - t_ext[i], p[0], ONE, coeff[i], p[0]);
  }

  free(t_ext);
  t_ext = NULL;

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

  const int n_hist = 5;
  sunrealtype t_hist[n_hist];

  for (int i = 0; i < n_hist; i++)
  {
    t_hist[i] = (5.0 * PI / 6.0) - i * PI / 6.0;
  }

  N_Vector* y = N_VCloneVectorArray(n_hist, tmp);
  if (check_ptr(y, "N_VCloneVectorArray")) return 1;

  N_Vector yp = N_VClone(tmp);
  if (check_ptr(yp, "N_VClone")) return 1;


  std::cout << "y_hist:" << std::endl;
  for (int i = 0; i < n_hist; i++)
  {
    sunrealtype* data = N_VGetArrayPointer(y[i]);
    data[0] = std::sin(t_hist[i]);
    std::cout << i << ": t = " << t_hist[i] << std::endl;
    N_VPrint(y[i]);
  }
  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  ypdata[0] = std::cos(t_hist[0]);
  std::cout << "yp = " << std::endl;
  N_VPrint(yp);
  std::cout << std::endl;

  // -------------------------------
  // Create interpolating polynomial
  // -------------------------------

  N_Vector* coeff = N_VCloneVectorArray(n_hist + 1, tmp);
  if (check_ptr(coeff, "N_VCloneVectorArray")) return 1;

  // Generate coefficients
  int flag = HermitePolyCoef(t_hist, y, yp, n_hist, coeff);
  if (flag) return 1;

  for (int i = 0; i < n_hist + 1; i++)
  {
    std::cout << "coeff[" << i << "] = " << std::endl;
    N_VPrint(coeff[i]);
  }
  std::cout << std::endl;

  // ---------------------------------
  // Evaluate interpolating polynomial
  // ---------------------------------

  N_Vector* p = N_VCloneVectorArray(n_hist + 1, tmp);
  if (check_ptr(p, "N_VCloneVectorArray")) return 1;

  flag = HermitePolyMultiDerEval(t_hist, coeff, n_hist, t_hist[0], n_hist, p);
  if (flag) return 1;

  for (int i = 0; i < n_hist + 1; i++)
  {
    sunrealtype p_true;
    if (i % 2 == 0)
    {
      p_true = std::sin(t_hist[0]);
    }
    else
    {
      p_true = std::cos(t_hist[0]);
    }
    p_true *= std::pow(-1, (i * (i - 1)) / 2);

    sunrealtype* p_data = N_VGetArrayPointer(p[i]);

    std::cout << i << std::endl
              << " p_data: " << std::setw(22) << p_data[0] << std::endl
              << " p_true: " << std::setw(22) << p_true  << std::endl
              << " error:  " << std::setw(22) << p_data[0] - p_true << std::endl
              << std::endl;
  }
  std::cout << std::endl;

  N_VDestroy(tmp);
  N_VDestroyVectorArray(y, n_hist);
  N_VDestroy(yp);
  N_VDestroyVectorArray(coeff, n_hist + 1);
  N_VDestroyVectorArray(p, n_hist + 1);

  return 0;
}
