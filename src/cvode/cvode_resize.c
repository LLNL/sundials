/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>

#include "cvode/cvode.h"
#include "cvode_impl.h"

#define ZERO SUN_RCONST(0.0) /* real 0.0     */
#define ONE  SUN_RCONST(1.0) /* real 1.0     */

/* -----------------------------------------------------------------------------
 * Build Adams Nordsieck array from f(t,y) history and y(t) value
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

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Build BDF Nordsieck array from y(t) history and f(t,y) value
 * ---------------------------------------------------------------------------*/

int BuildNordsieckArrayBDF(sunrealtype* t, N_Vector* y, N_Vector f,
                           N_Vector* wrk, int order, sunrealtype hscale,
                           N_Vector* zn)
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

  /* Compute Nordsieck array */
  if (order > 1)
  {
    /* Setup extended array of times to incorporate derivative value */
    sunrealtype t_ext[BDF_Q_MAX + 1];

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
      }
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

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Function to compute predicted new state (simplified cvPredict)
 * ---------------------------------------------------------------------------*/

int PredictY(int order, N_Vector* zn, N_Vector ypred)
{
  N_VScale(ONE, zn[0], ypred);
  for (int j = 1; j <= order; j++)
  {
    N_VLinearSum(ONE, zn[j], ONE, ypred, ypred);
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Function to resize CVODE and initialize with new history
 *
 * t_hist = [ t_{n}, t_{n - 1}, ..., t_{n - (n_hist - 1)} ]
 * y_hist = [ y_{n}, y_{n - 1}, ..., y_{n - (n_hist - 1)} ]
 * f_hist = [ f_{n}, f_{n - 1}, ..., f_{n - (n_hist - 1)} ]
 * ---------------------------------------------------------------------------*/

int CVodeResizeHistory(void* cvode_mem, sunrealtype* t_hist, N_Vector* y_hist,
                       N_Vector* f_hist, int n_hist)
{
  CVodeMem cv_mem = NULL;
  int retval      = 0;

  if (!cvode_mem)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__,
                   MSGCV_NO_MEM);
    return CV_MEM_NULL;
  }
  cv_mem = (CVodeMem)cvode_mem;

  if (!t_hist)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "Time history array is NULL");
    return CV_ILL_INPUT;
  }

  if (!y_hist)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "State history array is NULL");
    return CV_ILL_INPUT;
  }

  if (!f_hist)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "RHS history array is NULL");
    return CV_ILL_INPUT;
  }

  if (n_hist < 1)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "Invalid history size value");
    return CV_ILL_INPUT;
  }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "---------------\n");
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "Start Resize\n");
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "Input values\n");
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "n_hist         = %d\n", n_hist);
  for (int ithist = 0; ithist < n_hist; ithist++)
  {
    fprintf(cv_mem->cv_sunctx->logger->debug_fp, "t_hist[%d]      = %g\n",
            ithist, t_hist[ithist]);
  }
  for (int ithist = 0; ithist < n_hist; ithist++) /* n_hist is for BDF */
  {
    fprintf(cv_mem->cv_sunctx->logger->debug_fp, "y_hist[%d]\n", ithist);
    N_VPrintFile(y_hist[ithist], cv_mem->cv_sunctx->logger->debug_fp);
  }
  for (int ithist = 0; ithist < 2; ithist++) /* 2 is for BDF */
  {
    fprintf(cv_mem->cv_sunctx->logger->debug_fp, "f_hist[%d]\n", ithist);
    N_VPrintFile(f_hist[ithist], cv_mem->cv_sunctx->logger->debug_fp);
  }
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "CVODE values\n");
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "tn             = %g\n",
          cv_mem->cv_tn);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "current h      = %g\n",
          cv_mem->cv_h);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "next h         = %g\n",
          cv_mem->cv_hprime);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "next h (?)     = %g\n",
          cv_mem->cv_next_h);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "h scale        = %g\n",
          cv_mem->cv_hscale);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "current order  = %d\n",
          cv_mem->cv_q);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "next order     = %d\n",
          cv_mem->cv_qprime);
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "next order (?) = %d\n",
          cv_mem->cv_next_q);
  for (int ord = 0; ord <= cv_mem->cv_qmax_alloc; ord++)
  {
    fprintf(cv_mem->cv_sunctx->logger->debug_fp, "zn[%d]\n", ord);
    N_VPrintFile(cv_mem->cv_zn[ord], cv_mem->cv_sunctx->logger->debug_fp);
  }
#endif

  /* Make sure number of inputs is sufficient for the current (next) order */
  /* Make sure times[0] == tn */

  N_Vector tmpl = y_hist[0];
  int maxord    = cv_mem->cv_qmax_alloc;

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

  for (int j = 0; j <= maxord; j++)
  {
    N_VDestroy(cv_mem->cv_zn[j]);
    cv_mem->cv_zn[j] = N_VClone(tmpl);
    N_VConst(NAN, cv_mem->cv_zn[j]);
  }

  for (int j = 0; j <= maxord; j++)
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
    /* Compute z_{n-1} with new history size */
    if (cv_mem->cv_lmm == CV_ADAMS)
    {
      retval = BuildNordsieckArrayAdams(t_hist + 1, y_hist[1], f_hist + 1,
                                        cv_mem->resize_wrk, cv_mem->cv_q,
                                        cv_mem->cv_hscale, cv_mem->cv_zn);
    }
    else
    {
      retval = BuildNordsieckArrayBDF(t_hist + 1, y_hist + 1, f_hist[1],
                                      cv_mem->resize_wrk, cv_mem->cv_q,
                                      cv_mem->cv_hscale, cv_mem->cv_zn);
    }

    if (retval)
    {
      cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "BuildNordsieckArray failed");
      return CV_ILL_INPUT;
    }

    /* Get predicted value */
    retval = PredictY(cv_mem->cv_q, cv_mem->cv_zn, cv_mem->cv_vtemp1);

    if (retval)
    {
      cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "NewtonPolyMultiDerEval failed");
      return CV_ILL_INPUT;
    }

    N_VLinearSum(ONE, y_hist[0], -ONE, cv_mem->cv_vtemp1,
                 cv_mem->cv_zn[cv_mem->cv_qmax]);
  }

  /* ----------------------------- *
   * Construct new Nordsieck Array *
   * ----------------------------- */

  if (cv_mem->cv_lmm == CV_ADAMS)
  {
    retval = BuildNordsieckArrayAdams(t_hist, y_hist[0], f_hist,
                                      cv_mem->resize_wrk, cv_mem->cv_qprime,
                                      cv_mem->cv_hscale, cv_mem->cv_zn);
  }
  else
  {
    retval = BuildNordsieckArrayBDF(t_hist, y_hist, f_hist[0],
                                    cv_mem->resize_wrk, cv_mem->cv_qprime,
                                    cv_mem->cv_hscale, cv_mem->cv_zn);
  }

  if (retval)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "BuildNordsieckArray failed");
    return CV_ILL_INPUT;
  }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "Finish Resize\n");
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "tn = %g\n", cv_mem->cv_tn);
  for (ord = 0; ord <= cv_mem->cv_qmax_alloc; ord++)
  {
    fprintf(cv_mem->cv_sunctx->logger->debug_fp, "zn[%d]\n", ord);
    N_VPrintFile(cv_mem->cv_zn[ord], cv_mem->cv_sunctx->logger->debug_fp);
  }
  fprintf(cv_mem->cv_sunctx->logger->debug_fp, "---------------\n");
#endif

  /*
   * Why are these at the end?
   */

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
    cv_mem->cv_constraintsSet        = SUNFALSE;
  }

  cv_mem->first_step_after_resize = SUNTRUE;

  return CV_SUCCESS;
}
