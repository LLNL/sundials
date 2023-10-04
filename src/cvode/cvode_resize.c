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
 * Build Nordsieck array from solution history
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "cvode_impl.h"
#include <sundials/sundials_types.h>

#define ZERO    RCONST(0.0)     /* real 0.0     */
#define ONE     RCONST(1.0)     /* real 1.0     */

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
  int n_coeff = n_times + 1; /* n_times state values + 1 derivative */
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

  /* Initialize coefficient arrays with values to interpolate */
  N_VScale(ONE, y[0], coeff[0]);
  N_VScale(ONE, y[0], coeff[1]);

  for (i = 1; i < n_times; i++)
  {
    N_VScale(ONE, y[i], coeff[i + 1]);
  }

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
        sunrealtype denom = ONE / (t_ext[j - i] - t_ext[j]);
        N_VLinearSum(denom, coeff[j - 1], -denom, coeff[j], coeff[j]);
      }
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
  int n_coeff = n_times + 1; /* n_times state values + 1 derivative */
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
      /* P^(j)_i = P_{j + 1} * (s - t_i) + j * P^(j - 1)_i */
      N_VLinearSum(s - t_ext[i], p[j], j, p[j-1], p[j]);
    }
    /* P_i = P_{i + 1} * (s - t_i) + c_i */
    N_VLinearSum(s - t_ext[i], p[0], ONE, coeff[i], p[0]);
  }

  free(t_ext);
  t_ext = NULL;

  return CV_SUCCESS;
}


/*
 * t_hist = [ t_{n}, t_{n - 1}, ..., t_{n - n_hist - 1} ]
 * y_hist = [ y_{n}, y_{n - 1}, ..., y_{n - n_hist - 1} ]
 * f_hist = [ f_{n}, f_{n - 1}, ..., f_{n - n_hist - 1} ]
 */

int CVodeResizeHistory(void *cvode_mem, sunrealtype* t_hist, N_Vector* y_hist,
                       N_Vector* f_hist, int n_hist, CVResizeVecFn resize_fn,
                       FILE* debug_file)
{
  CVodeMem cv_mem = NULL;
  int retval = 0;
  int i, j, ithist, ord;
  N_Vector tmpl;
  int maxord;
  sunrealtype scale;

  if (!cvode_mem)
  {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeResizeHistory",
                   MSGCV_NO_MEM);
    return CV_MEM_NULL;
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (!t_hist)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "Time history array is NULL");
    return CV_ILL_INPUT;
  }

  if (!y_hist)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "State history array is NULL");
    return CV_ILL_INPUT;
  }

  if (n_hist < 1)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "Invalid history size value");
    return CV_ILL_INPUT;
  }

  if (!resize_fn)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "Resize function is NULL");
    return CV_ILL_INPUT;
  }

  if (debug_file)
  {
    fprintf(debug_file, "---------------\n");
    fprintf(debug_file, "Start Resize\n");
    fprintf(debug_file, "Input values\n");
    fprintf(debug_file, "n_hist         = %d\n", n_hist);
    for (ithist = 0; ithist < n_hist; ithist++)
    {
      fprintf(debug_file, "t_hist[%d]      = %g\n", ithist, t_hist[ithist]);
    }
    for (ithist = 0; ithist < n_hist; ithist++)
    {
      fprintf(debug_file, "y_hist[%d]\n", ithist);
      N_VPrintFile(y_hist[ithist], debug_file);
    }
    fprintf(debug_file, "CVODE values\n");
    fprintf(debug_file, "tn             = %g\n", cv_mem->cv_tn);
    fprintf(debug_file, "current h      = %g\n", cv_mem->cv_h);
    fprintf(debug_file, "next h         = %g\n", cv_mem->cv_hprime);
    fprintf(debug_file, "next h (?)     = %g\n", cv_mem->cv_next_h);
    fprintf(debug_file, "h scale        = %g\n", cv_mem->cv_hscale);
    fprintf(debug_file, "current order  = %d\n", cv_mem->cv_q);
    fprintf(debug_file, "next order     = %d\n", cv_mem->cv_qprime);
    fprintf(debug_file, "next order (?) = %d\n", cv_mem->cv_next_q);
    for (ord = 0; ord <= cv_mem->cv_qmax_alloc; ord++)
    {
      fprintf(debug_file, "zn[%d]\n", ord);
      N_VPrintFile(cv_mem->cv_zn[ord], debug_file);
    }
  }

  /* Make sure number of inputs is sufficient for the current (next) order */
  /* Make sure times[0] == tn */

  tmpl = y_hist[0];
  maxord = cv_mem->cv_qmax_alloc;

  /* -------------- *
   * Resize vectors *
   * -------------- */

  N_VDestroy(cv_mem->cv_ewt);
  cv_mem->cv_ewt = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_acor);
  cv_mem->cv_acor = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_tempv);
  cv_mem->cv_tempv = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_ftemp);
  cv_mem->cv_ftemp = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_vtemp1);
  cv_mem->cv_vtemp1 = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_vtemp2);
  cv_mem->cv_vtemp2 = N_VClone(tmpl);

  N_VDestroy(cv_mem->cv_vtemp3);
  cv_mem->cv_vtemp3 = N_VClone(tmpl);

  for (j = 0; j <= maxord; j++)
  {
    N_VDestroy(cv_mem->cv_zn[j]);
    cv_mem->cv_zn[j] = N_VClone(tmpl);
    N_VConst(NAN, cv_mem->cv_zn[j]);
  }

  for (j = 0; j <= maxord; j++)
  {
    N_VDestroy(cv_mem->resize_wrk[j]);
    cv_mem->resize_wrk[j] = N_VClone(tmpl);
    N_VConst(NAN, cv_mem->resize_wrk[j]);
  }

  /* ------------------------------------------------------------------------ *
   * Construct Nordsieck array at the old time but with the new size to
   * compute correction vector at the new state size.
   * ------------------------------------------------------------------------ */

  if (cv_mem->cv_q < cv_mem->cv_qmax)
  {
    if (!f_hist)
    {
      retval = cv_mem->cv_f(t_hist[1], y_hist[1],
                            cv_mem->cv_vtemp2, cv_mem->cv_user_data);
      cv_mem->cv_nfe++;
      if (retval)
      {
        cvProcessError(cv_mem, CV_RHSFUNC_FAIL, "CVODE", "CVode",
                       MSGCV_RHSFUNC_FAILED, cv_mem->cv_tn);
        return CV_RHSFUNC_FAIL;
      }
    }
    else
    {
      N_VScale(ONE, f_hist[1], cv_mem->cv_vtemp2);
    }

    /* Compute interpolation coefficients */
    retval = HermitePolyCoef(t_hist + 1, y_hist + 1, cv_mem->cv_vtemp2, cv_mem->cv_q,
                             cv_mem->resize_wrk);
    if (retval)
    {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                     "HermitePolyCoef failed");
      return CV_ILL_INPUT;
    }

    /* Get predicted value */
    retval = HermitePolyMultiDerEval(t_hist + 1, cv_mem->resize_wrk,
                                     cv_mem->cv_q,
                                     cv_mem->cv_tn, 0,
                                     &(cv_mem->cv_vtemp2));
    if (retval)
    {
      cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                     "NewtonPolyMultiDerEval failed");
      return CV_ILL_INPUT;
    }

    N_VLinearSum(ONE, y_hist[0], -ONE, cv_mem->cv_vtemp2, cv_mem->cv_vtemp2);
    N_VScale(ONE, cv_mem->cv_vtemp2, cv_mem->cv_zn[cv_mem->cv_qmax]);
  }

  /* ----------------------------- *
   * Construct new Nordsieck Array *
   * ----------------------------- */

  if (!f_hist)
  {
    retval = cv_mem->cv_f(cv_mem->cv_tn, y_hist[0],
                          cv_mem->cv_vtemp2, cv_mem->cv_user_data);
    cv_mem->cv_nfe++;
    if (retval)
    {
      cvProcessError(cv_mem, CV_RHSFUNC_FAIL, "CVODE", "CVode",
                     MSGCV_RHSFUNC_FAILED, cv_mem->cv_tn);
      return CV_RHSFUNC_FAIL;
    }
  }
  else
  {
    N_VScale(ONE, f_hist[0], cv_mem->cv_vtemp2);
  }

  /* Compute interpolation coefficients */
  /* >>> TODO(DJG): use q' for BDF and q for ADAMS <<< */
  retval = HermitePolyCoef(t_hist, y_hist, cv_mem->cv_vtemp2, cv_mem->cv_qprime,
                           cv_mem->resize_wrk);
  if (retval)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "HermitePolyCoef failed");
    return CV_ILL_INPUT;
  }

  /*  >>> TODO(DJG): use q' for BDF and q for ADAMS <<< */
  retval = HermitePolyMultiDerEval(t_hist, cv_mem->resize_wrk,
                                   cv_mem->cv_qprime,
                                   cv_mem->cv_tn, cv_mem->cv_qprime,
                                   cv_mem->cv_zn);
  if (retval)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeResizeHistory",
                   "NewtonPolyMultiDerEval failed");
    return CV_ILL_INPUT;
  }

  N_VScale(ONE, y_hist[0], cv_mem->cv_zn[0]);
  N_VScale(ONE, cv_mem->cv_vtemp2, cv_mem->cv_zn[1]);

  /* >>> TODO(DJG): use q' for BDF and q for ADAMS <<< */
  scale = ONE;
  for (i = 1; i < cv_mem->cv_qprime + 1; i++)
  {
    scale *= cv_mem->cv_hscale / ((sunrealtype) i);
    N_VScale(scale, cv_mem->cv_zn[i], cv_mem->cv_zn[i]);
  }

  if (debug_file)
  {
    fprintf(debug_file, "Finish Resize\n");
    fprintf(debug_file, "tn = %g\n", cv_mem->cv_tn);
    for (ord = 0; ord <= cv_mem->cv_qmax_alloc; ord++)
    {
      fprintf(debug_file, "zn[%d]\n", ord);
      N_VPrintFile(cv_mem->cv_zn[ord], debug_file);
    }
    fprintf(debug_file, "---------------\n");
  }

  if (cv_mem->cv_VabstolMallocDone)
  {
    N_VDestroy(cv_mem->cv_Vabstol);
    cv_mem->cv_Vabstol = N_VClone(tmpl);
  }

  /* User will need to set a new constraints vector */
  if (cv_mem->cv_constraintsMallocDone)
  {
    N_VDestroy(cv_mem->cv_constraints);
    cv_mem->cv_constraintsMallocDone = SUNFALSE;
    cv_mem->cv_constraintsSet = SUNFALSE;
  }

  cv_mem->first_step_after_resize = SUNTRUE;

  return CV_SUCCESS;
}
