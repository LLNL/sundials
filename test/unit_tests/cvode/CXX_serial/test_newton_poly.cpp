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

/*
 * Compute Newton interpolating polynomial coefficients
 *
 * Inputs:
 *   t -- array (length M) of time values, t_i
 *   y -- array (length M) of state vectors (length N), y_i
 *   N -- number of interpolation times
 *
 * Outputs:
 *   c -- array (length M) of coefficient vectors (length N)
 */
int NewtonPolyCoef(sunrealtype* t, N_Vector* y, int M, N_Vector* c)
{
  int i, j;

  /* Check for valid inputs */
  if (!t || !y || M < 1 || !c) return CV_ILL_INPUT;

  for (i = 0; i < M; i++)
  {
    if (!y[i]) return CV_ILL_INPUT;
    if (!c[i]) return CV_ILL_INPUT;
  }

  /* Initialize coefficient arrays with values to interpolate */
  for (i = 0; i < M; i++)
  {
    N_VScale(ONE, y[i], c[i]);
  }

  /* printf("y hist\n"); */
  /* for (int k = 0; k < M; k++) */
  /* { */
  /*   printf("y[%d]\n", k); */
  /*   N_VPrint(y[k]); */
  /* } */

  /* printf("Initial values\n"); */
  /* for (int k = 0; k < M; k++) */
  /* { */
  /*   printf("c[%d]\n", k); */
  /*   N_VPrint(c[k]); */
  /* } */

  if (M == 1) return CV_SUCCESS;

  /* Compute coefficients from bottom up to write in place */
  /* printf("Iteration Newton Coef\n"); */
  for (i = 1; i < M; i++)
  {
    for (j = M - 1; j > i - 1; j--)
    {
      printf("%d %d c[%d]\n", i, j, j);
      printf("t[%d] = %.16e\n", j, t[j]);
      printf("t[%d] = %.16e\n", j-1, t[j-1]);
      printf("c[%d]\n", j);
      N_VPrint(c[j]);
      printf("c[%d]\n", j - 1);
      N_VPrint(c[j - 1]);

      /* c_j = (c_j - c_{j - 1}) / (t_j - t_{j - i}) */
      N_VLinearSum(ONE / (t[j] - t[j - i]), c[j],
                   -ONE / (t[j] - t[j - i]), c[j - 1], c[j]);
      printf("result c[%d]\n", j);
      N_VPrint(c[j]);
    }
  }

  /* printf("Newton Coef\n"); */
  /* for (int k = 0; k < M; k++) */
  /* { */
  /*   printf("c[%d]\n", k); */
  /*   N_VPrint(c[k]); */
  /* } */

  return CV_SUCCESS;
}

/*
 * Evaluate the Newton interpolating polynomial and its derivatives up to
 * order d at a given time
 *
 * Inputs:
 *
 * t -- array (length M) of interpolated time values, t_i
 * c -- array (length M) of polynomial coefficient vectors (length N), c_i
 * s -- time at which to evaluate the polynomial, s
 *
 * Output:
 *
 * p -- array (length d + 1) of values p[0] = p(s), p[1] = p'[s], etc.
 *
 * Example:
 *
 * Consider the M = 4 case, we have the interpolation times t0, t1, t2, t3 and
 * the coefficients c0, c1, c2, c3. The Newton interpolating polynomial is
 *
 * p(s) = c0 + c1 (s - t0) + c2 (s - t0)(s - t1) + c3 (s - t0)(s - t1)(s - t2)
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

int NewtonPolyMultiDerEval(sunrealtype* t, N_Vector* c, int M, sunrealtype s,
                           int d, N_Vector* p)
{
  int i, j;

  /* Check for valid inputs */
  if (!t || !c || M < 1 || d < 0 || !p) return CV_ILL_INPUT;

  for (i = 0; i < M; i++)
  {
    if (!c[i]) return CV_ILL_INPUT;
  }

  for (i = 0; i < d + 1; i++)
  {
    if (!p[i]) return CV_ILL_INPUT;
  }

  /* Initialize interpolation output to P_{M-1} = c_{M-1} and derivative output
     to zero */
  N_VScale(ONE, c[M - 1], p[0]);
  for (i = 1; i < d + 1; i++)
  {
    N_VConst(ZERO, p[i]);
  }

  /* Accumulate polynomial terms in-place i.e., P_{i+1} is stored in p[0] and
     overwritten each iteration by P_{i} and similarly for P', P'', etc. */
  for (i = M - 2; i >= 0; i--)
  {
    for (j = d; j > 0; j--)
    {
      // P^(j)_i = P_{j + 1} * (s - t_i) + j * P^(j - 1)_i
      N_VLinearSum(s - t[i], p[j], j, p[j-1], p[j]);
    }
    // P_i = P_{i + 1} * (s - t_i) + c_i
    N_VLinearSum(s - t[i], p[0], ONE, c[i], p[0]);
  }

  /* printf("Newton Eval\n"); */
  /* for (int k = 0; k < d + 1; k++) */
  /* { */
  /*   printf("p[%d]\n", k); */
  /*   N_VPrint(p[k]); */
  /* } */

  return CV_SUCCESS;
}


// Test main
int main(int argc, char* argv[])
{
  // Set output formatting
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create data
  N_Vector tmp = N_VNew_Serial(2, sunctx);
  if (check_ptr(tmp, "N_VNew_Serial")) return 1;

  N_Vector* y = N_VCloneVectorArray(5, tmp);
  if (check_ptr(y, "N_VCloneVectorArray")) return 1;

  N_Vector* c = N_VCloneVectorArray(5, tmp);
  if (check_ptr(c, "N_VCloneVectorArray")) return 1;

  N_Vector* p = N_VCloneVectorArray(5, tmp);
  if (check_ptr(p, "N_VCloneVectorArray")) return 1;

  // (0, 9)
  // (3, 24)
  // (4, 15)
  // (8, -27)
  // (8, 66)
  sunrealtype sample_times[5] = {0.0, 3.0, 4.0, 8.0, 10.0};

  for (int i = 0; i < 5; i++)
  {
    sunrealtype t = sample_times[i];
    sunrealtype* ydata = N_VGetArrayPointer(y[i]);
    // (t + 2)(t - 5)(t - 9)(t + 1) / 10
    //  = 1/10 (t^4 - 11 t^3 + 5 t^2 + 107 t + 90)
    ydata[0] = (std::pow(t, 4) - 11 * std::pow(t, 3) + 5 * std::pow(t, 2) + 107 * t + 90) / 10.0;
    ydata[1] = std::sin(t);
    N_VPrint(y[i]);
  }
  std::cout << std::endl;

  // Generate coefficients
  int flag = NewtonPolyCoef(sample_times, y, 5, c);
  if (flag) return 1;

  for (int i = 0; i < 5; i++)
  {
    N_VPrint(c[i]);
  }
  std::cout << std::endl;

  // Evaluate polynomial
  sunrealtype eval_times[5] = {0.0, 3.0, 4.0, 8.0, 10.0};
  //sunrealtype eval_times[4] = {1.0, 3.0, 5.0, 7.0};

  for (int i = 0; i < 4; i++)
  {
    sunrealtype t = eval_times[i];

    flag = NewtonPolyMultiDerEval(sample_times, c, 5, t, 0, p);
    if (flag) return 1;

    sunrealtype p0 = (std::pow(t, 4)
                      - 11 * std::pow(t, 3)
                      + 5 * std::pow(t, 2)
                      + 107 * t
                      + 90) / 10.0;

    // sunrealtype p0p1 = (4.0 * std::pow(t, 3)
    //                     - 33.0 * std::pow(t, 2)
    //                     + 10.0 * t
    //                     + 107.0) / 10.0;

    // sunrealtype p0p2 = 3.0 * (4.0 * t - 11.0) / 5.0;

    // sunrealtype p0p3 = (3.0 * (5.0 * t * t - 4.0 * t - 47.0)) / 2500.0;

    // sunrealtype p0p4 = 12.0 / 5.0;

    // sunrealtype p1   =  std::sin(t);
    // sunrealtype p1p1 =  std::cos(t);
    // sunrealtype p1p2 = -std::sin(t);
    // sunrealtype p1p3 = -std::cos(t);
    // sunrealtype p1p4 =  std::sin(t);

    sunrealtype* pdata = N_VGetArrayPointer(p[0]);

    std::cout << "True" << std::setw(25) << p0 << std::endl;
    std::cout << "Poly" << std::setw(25) << pdata[0] << std::endl;

    // std::cout << "True" << std::setw(25) << p1 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[1] << std::endl;

    // pdata = N_VGetArrayPointer(p[1]);

    // std::cout << "True" << std::setw(25) << p0p1 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[0] << std::endl;

    // std::cout << "True"<< std::setw(25) << p1p1 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[1] << std::endl;

    // pdata = N_VGetArrayPointer(p[2]);

    // std::cout << "True" << std::setw(25) << p0p2 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[0] << std::endl;

    // std::cout << "True" << std::setw(25) << p1p2 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[1] << std::endl;

    // pdata = N_VGetArrayPointer(p[3]);

    // std::cout << "True" << std::setw(25) << p0p3 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[0] << std::endl;

    // std::cout << "True" << std::setw(25) << p1p3 << std::endl;
    // std::cout << "Poly" << std::setw(25) << pdata[1] << std::endl;
  }

  return 0;
}
