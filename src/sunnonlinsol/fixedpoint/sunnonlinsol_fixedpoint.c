/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the SUNNonlinearSolver module
 * implementation of the Anderson-accelerated Fixed-Point method.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_nvector_senswrapper.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Internal utility routines */
static SUNErrCode AndersonAccelerate(SUNNonlinearSolver NLS, N_Vector gval,
                                     N_Vector x, N_Vector xold, int iter);

static SUNErrCode AllocateContent(SUNNonlinearSolver NLS, N_Vector tmpl);
static void FreeContent(SUNNonlinearSolver NLS);

/* Content structure accessibility macros */
#define FP_CONTENT(S) ((SUNNonlinearSolverContent_FixedPoint)(S->content))

/* Constant macros */
#define ONE  SUN_RCONST(1.0)
#define ZERO SUN_RCONST(0.0)

/*==============================================================================
  Constructor to create a new fixed point solver
  ============================================================================*/

SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS                       = NULL;
  SUNNonlinearSolverContent_FixedPoint content = NULL;

  /* Check that the supplied N_Vector supports all required operations */
  SUNAssertNull(y->ops->nvclone && y->ops->nvdestroy && y->ops->nvscale &&
                  y->ops->nvlinearsum && y->ops->nvdotprod,
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create nonlinear linear solver */
  NLS = SUNNonlinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  NLS->ops->gettype         = SUNNonlinSolGetType_FixedPoint;
  NLS->ops->initialize      = SUNNonlinSolInitialize_FixedPoint;
  NLS->ops->solve           = SUNNonlinSolSolve_FixedPoint;
  NLS->ops->free            = SUNNonlinSolFree_FixedPoint;
  NLS->ops->setsysfn        = SUNNonlinSolSetSysFn_FixedPoint;
  NLS->ops->setctestfn      = SUNNonlinSolSetConvTestFn_FixedPoint;
  NLS->ops->setmaxiters     = SUNNonlinSolSetMaxIters_FixedPoint;
  NLS->ops->getnumiters     = SUNNonlinSolGetNumIters_FixedPoint;
  NLS->ops->getcuriter      = SUNNonlinSolGetCurIter_FixedPoint;
  NLS->ops->getnumconvfails = SUNNonlinSolGetNumConvFails_FixedPoint;

  /* Create nonlinear solver content structure */
  content = NULL;
  content = (SUNNonlinearSolverContent_FixedPoint)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Initialize all components of content to 0/NULL */
  memset(content, 0, sizeof(struct _SUNNonlinearSolverContent_FixedPoint));

  /* Attach content */
  NLS->content = content;

  /* Fill general content */
  content->Sys        = NULL;
  content->CTest      = NULL;
  content->m          = m;
  content->damping    = SUNFALSE;
  content->beta       = ONE;
  content->curiter    = 0;
  content->maxiters   = 3;
  content->niters     = 0;
  content->nconvfails = 0;
  content->ctest_data = NULL;

  /* Fill allocatable content */
  SUNCheckCallNull(AllocateContent(NLS, y));

  return (NLS);
}

/*==============================================================================
  Constructor wrapper to create a new fixed point solver for sensitivity solvers
  ============================================================================*/

SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m,
                                               SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS = NULL;
  N_Vector w             = NULL;

  /* create sensitivity vector wrapper */
  w = N_VNew_SensWrapper(count, y);
  SUNCheckLastErrNull();

  /* create nonlinear solver using sensitivity vector wrapper */
  NLS = SUNNonlinSol_FixedPoint(w, m, sunctx);
  SUNCheckLastErrNull();

  /* free sensitivity vector wrapper */
  N_VDestroy(w);
  SUNCheckLastErrNull();

  /* return NLS object */
  return (NLS);
}

/*==============================================================================
  GetType, Initialize, Setup, Solve, and Free operations
  ============================================================================*/

SUNNonlinearSolver_Type SUNNonlinSolGetType_FixedPoint(
  SUNDIALS_MAYBE_UNUSED SUNNonlinearSolver NLS)
{
  return (SUNNONLINEARSOLVER_FIXEDPOINT);
}

SUNErrCode SUNNonlinSolInitialize_FixedPoint(SUNNonlinearSolver NLS)
{
  SUNFunctionBegin(NLS->sunctx);
  /* check that all required function pointers have been set */
  SUNAssert(FP_CONTENT(NLS)->Sys && FP_CONTENT(NLS)->CTest, SUN_ERR_ARG_CORRUPT);

  /* reset the total number of iterations and convergence failures */
  FP_CONTENT(NLS)->niters     = 0;
  FP_CONTENT(NLS)->nconvfails = 0;

  return SUN_SUCCESS;
}

/*-----------------------------------------------------------------------------
  SUNNonlinSolSolve_FixedPoint: Performs the fixed-point solve g(y) = y

  Successful solve return code:
   SUN_SUCCESS = 0

  Recoverable failure return codes (positive):
    SUN_NLS_CONV_RECVR
    *_RHSFUNC_RECVR (ODEs) or *_RES_RECVR (DAEs)

  Unrecoverable failure return codes (negative):
    *_MEM_NULL
    *_RHSFUNC_FAIL (ODEs) or *_RES_FAIL (DAEs)

  Note that return values beginning with * are package specific values returned
  by the Sys function provided to the nonlinear solver.
  ---------------------------------------------------------------------------*/

int SUNNonlinSolSolve_FixedPoint(SUNNonlinearSolver NLS,
                                 SUNDIALS_MAYBE_UNUSED N_Vector y0,
                                 N_Vector ycor, N_Vector w, sunrealtype tol,
                                 SUNDIALS_MAYBE_UNUSED sunbooleantype callSetup,
                                 void* mem)
{
  SUNFunctionBegin(NLS->sunctx);
  /* local variables */
  int retval;
  N_Vector yprev, gy, delta;

  /* check that all required function pointers have been set */
  SUNAssert(FP_CONTENT(NLS)->Sys && FP_CONTENT(NLS)->CTest, SUN_ERR_ARG_CORRUPT);

  SUNLogInfo(NLS->sunctx->logger, "nonlinear-solver", "solver = Fixed-Point");

  /* set local shortcut variables */
  yprev = FP_CONTENT(NLS)->yprev;
  gy    = FP_CONTENT(NLS)->gy;
  delta = FP_CONTENT(NLS)->delta;

  /* initialize iteration and convergence fail counters for this solve */
  FP_CONTENT(NLS)->niters     = 0;
  FP_CONTENT(NLS)->nconvfails = 0;

  /* Looping point for attempts at solution of the nonlinear system:
       Evaluate fixed-point function (store in gy).
       Performs the accelerated fixed-point iteration.
       Performs stopping tests. */
  for (FP_CONTENT(NLS)->curiter = 0;
       FP_CONTENT(NLS)->curiter < FP_CONTENT(NLS)->maxiters;
       FP_CONTENT(NLS)->curiter++)
  {
    SUNLogInfo(NLS->sunctx->logger, "begin-nonlinear-iterate", "");

    /* update previous solution guess */
    N_VScale(ONE, ycor, yprev);
    SUNCheckLastErr();

    /* Compute fixed-point iteration function, store in gy.
       We do not use SUNCheck macros here because Sys is an integrator-provided
       callback and returns integrator specific error values  where 0 == success,
       < 0 is a failure, > 0 is recoverable error. */
    retval = FP_CONTENT(NLS)->Sys(ycor, gy, mem);
    if (retval != 0)
    {
      SUNLogInfo(NLS->sunctx->logger, "end-nonlinear-iterate",
                 "status = failed nonlinear system evaluation, retval = %d",
                 retval);
      return retval;
    }

    /* perform fixed point update, based on choice of acceleration or not */
    if (FP_CONTENT(NLS)->m == 0)
    { /* basic fixed-point solver */
      N_VScale(ONE, gy, ycor);
      SUNCheckLastErr();
    }
    else
    { /* Anderson-accelerated solver */
      SUNCheckCall(
        AndersonAccelerate(NLS, gy, ycor, yprev, FP_CONTENT(NLS)->curiter));
    }

    /* increment nonlinear solver iteration counter */
    FP_CONTENT(NLS)->niters++;

    /* compute change in solution, and call the convergence test function */
    N_VLinearSum(ONE, ycor, -ONE, yprev, delta);
    SUNCheckLastErr();

    /* test for convergence */
    retval = FP_CONTENT(NLS)->CTest(NLS, ycor, delta, tol, w,
                                    FP_CONTENT(NLS)->ctest_data);

    SUNLogInfo(NLS->sunctx->logger, "nonlinear-iterate",
               "cur-iter = %d, update-norm = %.16g", FP_CONTENT(NLS)->niters,
               N_VWrmsNorm(delta, w));

    /* return if successful */
    if (retval == 0)
    {
      SUNLogInfo(NLS->sunctx->logger, "end-nonlinear-iterate",
                 "status = success");
      return SUN_SUCCESS;
    }

    /* check if the iterations should continue; otherwise increment the
       convergence failure count and return error flag */
    if (retval != SUN_NLS_CONTINUE)
    {
      SUNLogInfo(NLS->sunctx->logger, "end-nonlinear-iterate",
                 "status = failed, retval = %i", retval);
      FP_CONTENT(NLS)->nconvfails++;
      return (retval);
    }

    SUNLogInfoIf(FP_CONTENT(NLS)->curiter < FP_CONTENT(NLS)->maxiters - 1,
                 NLS->sunctx->logger, "end-nonlinear-iterate",
                 "status = continue");
  }

  SUNLogInfo(NLS->sunctx->logger, "end-nonlinear-iterate",
             "status = failed max iterations");

  /* if we've reached this point, then we exhausted the iteration limit;
     increment the convergence failure count and return */
  FP_CONTENT(NLS)->nconvfails++;
  return SUN_NLS_CONV_RECVR;
}

SUNErrCode SUNNonlinSolFree_FixedPoint(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL) { return SUN_SUCCESS; }

  /* free items from content structure, then the structure itself */
  if (NLS->content)
  {
    FreeContent(NLS);
    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops)
  {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  /* free the overall NLS structure */
  free(NLS);

  return SUN_SUCCESS;
}

/*==============================================================================
  Set functions
  ============================================================================*/

SUNErrCode SUNNonlinSolSetSysFn_FixedPoint(SUNNonlinearSolver NLS,
                                           SUNNonlinSolSysFn SysFn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(SysFn, SUN_ERR_ARG_CORRUPT);
  FP_CONTENT(NLS)->Sys = SysFn;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetConvTestFn_FixedPoint(SUNNonlinearSolver NLS,
                                                SUNNonlinSolConvTestFn CTestFn,
                                                void* ctest_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(CTestFn, SUN_ERR_ARG_CORRUPT);

  FP_CONTENT(NLS)->CTest = CTestFn;

  /* attach convergence test data */
  FP_CONTENT(NLS)->ctest_data = ctest_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetMaxIters_FixedPoint(SUNNonlinearSolver NLS, int maxiters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(maxiters >= 1, SUN_ERR_ARG_OUTOFRANGE);
  FP_CONTENT(NLS)->maxiters = maxiters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetDamping_FixedPoint(SUNNonlinearSolver NLS,
                                             sunrealtype beta)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(beta > 0, SUN_ERR_ARG_OUTOFRANGE);

  if (beta < ONE)
  {
    /* enable damping */
    FP_CONTENT(NLS)->beta    = beta;
    FP_CONTENT(NLS)->damping = SUNTRUE;
  }
  else
  {
    /* disable damping */
    FP_CONTENT(NLS)->beta    = ONE;
    FP_CONTENT(NLS)->damping = SUNFALSE;
  }

  return SUN_SUCCESS;
}

/*==============================================================================
  Get functions
  ============================================================================*/

SUNErrCode SUNNonlinSolGetNumIters_FixedPoint(SUNNonlinearSolver NLS,
                                              long int* niters)
{
  /* return number of nonlinear iterations in the last solve */
  *niters = FP_CONTENT(NLS)->niters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetCurIter_FixedPoint(SUNNonlinearSolver NLS, int* iter)
{
  /* return the current nonlinear solver iteration count */
  *iter = FP_CONTENT(NLS)->curiter;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNumConvFails_FixedPoint(SUNNonlinearSolver NLS,
                                                  long int* nconvfails)
{
  /* return the total number of nonlinear convergence failures */
  *nconvfails = FP_CONTENT(NLS)->nconvfails;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetSysFn_FixedPoint(SUNNonlinearSolver NLS,
                                           SUNNonlinSolSysFn* SysFn)
{
  /* return the nonlinear system defining function */
  *SysFn = FP_CONTENT(NLS)->Sys;
  return SUN_SUCCESS;
}

/*=============================================================================
  Utility routines
  ===========================================================================*/

/*---------------------------------------------------------------
  AndersonAccelerate

  This routine computes the Anderson-accelerated fixed point
  iterate.  Upon entry, the predicted solution is held in xold;
  this array is never changed throughout this routine.

  The result of the routine is held in x.
  -------------------------------------------------------------*/
static SUNErrCode AndersonAccelerate(SUNNonlinearSolver NLS, N_Vector gval,
                                     N_Vector x, N_Vector xold, int iter)
{
  SUNFunctionBegin(NLS->sunctx);
  /* local variables */
  int nvec, i_pt, i, j, lAA, maa, *ipt_map;
  sunrealtype a, b, rtemp, c, s, beta, onembeta, *cvals, *R, *gamma;
  N_Vector fv, vtemp, gold, fold, *df, *dg, *Q, *Xvecs;
  sunbooleantype damping;

  /* local shortcut variables */
  vtemp   = x; /* use result as temporary vector */
  ipt_map = FP_CONTENT(NLS)->imap;
  maa     = FP_CONTENT(NLS)->m;
  gold    = FP_CONTENT(NLS)->gold;
  fold    = FP_CONTENT(NLS)->fold;
  df      = FP_CONTENT(NLS)->df;
  dg      = FP_CONTENT(NLS)->dg;
  Q       = FP_CONTENT(NLS)->q;
  cvals   = FP_CONTENT(NLS)->cvals;
  Xvecs   = FP_CONTENT(NLS)->Xvecs;
  R       = FP_CONTENT(NLS)->R;
  gamma   = FP_CONTENT(NLS)->gamma;
  fv      = FP_CONTENT(NLS)->delta;
  damping = FP_CONTENT(NLS)->damping;
  beta    = FP_CONTENT(NLS)->beta;

  /* reset ipt_map, i_pt */
  for (i = 0; i < maa; i++) { ipt_map[i] = 0; }
  i_pt = iter - 1 - ((iter - 1) / maa) * maa;

  /* update dg[i_pt], df[i_pt], fv, gold and fold*/
  N_VLinearSum(ONE, gval, -ONE, xold, fv);
  SUNCheckLastErr();
  if (iter > 0)
  {
    N_VLinearSum(ONE, gval, -ONE, gold, dg[i_pt]);
    SUNCheckLastErr(); /* dg_new = gval - gold */
    N_VLinearSum(ONE, fv, -ONE, fold, df[i_pt]);
    SUNCheckLastErr(); /* df_new = fv - fold */
  }
  N_VScale(ONE, gval, gold);
  SUNCheckLastErr();
  N_VScale(ONE, fv, fold);
  SUNCheckLastErr();

  /* on first iteration, just do basic fixed-point update */
  if (iter == 0)
  {
    N_VScale(ONE, gval, x);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* update data structures based on current iteration index */

  if (iter == 1)
  { /* second iteration */

    R[0] = N_VDotProd(df[i_pt], df[i_pt]);
    SUNCheckLastErr();
    R[0] = SUNRsqrt(R[0]);
    N_VScale(ONE / R[0], df[i_pt], Q[i_pt]);
    SUNCheckLastErr();
    ipt_map[0] = 0;
  }
  else if (iter <= maa)
  { /* another iteration before we've reached maa */

    N_VScale(ONE, df[i_pt], vtemp);
    SUNCheckLastErr();
    for (j = 0; j < iter - 1; j++)
    {
      ipt_map[j]              = j;
      R[(iter - 1) * maa + j] = N_VDotProd(Q[j], vtemp);
      SUNCheckLastErr();
      N_VLinearSum(ONE, vtemp, -R[(iter - 1) * maa + j], Q[j], vtemp);
      SUNCheckLastErr();
    }
    R[(iter - 1) * maa + iter - 1] = N_VDotProd(vtemp, vtemp);
    SUNCheckLastErr();
    R[(iter - 1) * maa + iter - 1] = SUNRsqrt(R[(iter - 1) * maa + iter - 1]);
    if (R[(iter - 1) * maa + iter - 1] == ZERO)
    {
      N_VScale(ZERO, vtemp, Q[i_pt]);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale((ONE / R[(iter - 1) * maa + iter - 1]), vtemp, Q[i_pt]);
      SUNCheckLastErr();
    }
    ipt_map[iter - 1] = iter - 1;
  }
  else
  { /* we've filled the acceleration subspace, so start recycling */

    /* delete left-most column vector from QR factorization */
    for (i = 0; i < maa - 1; i++)
    {
      a                        = R[(i + 1) * maa + i];
      b                        = R[(i + 1) * maa + i + 1];
      rtemp                    = SUNRsqrt(a * a + b * b);
      c                        = a / rtemp;
      s                        = b / rtemp;
      R[(i + 1) * maa + i]     = rtemp;
      R[(i + 1) * maa + i + 1] = ZERO;
      if (i < maa - 1)
      {
        for (j = i + 2; j < maa; j++)
        {
          a                  = R[j * maa + i];
          b                  = R[j * maa + i + 1];
          rtemp              = c * a + s * b;
          R[j * maa + i + 1] = -s * a + c * b;
          R[j * maa + i]     = rtemp;
        }
      }
      N_VLinearSum(c, Q[i], s, Q[i + 1], vtemp);
      SUNCheckLastErr();
      N_VLinearSum(-s, Q[i], c, Q[i + 1], Q[i + 1]);
      SUNCheckLastErr();
      N_VScale(ONE, vtemp, Q[i]);
      SUNCheckLastErr();
    }

    /* ahift R to the left by one */
    for (i = 1; i < maa; i++)
    {
      for (j = 0; j < maa - 1; j++) { R[(i - 1) * maa + j] = R[i * maa + j]; }
    }

    /* add the new df vector */
    N_VScale(ONE, df[i_pt], vtemp);
    SUNCheckLastErr();
    for (j = 0; j < maa - 1; j++)
    {
      R[(maa - 1) * maa + j] = N_VDotProd(Q[j], vtemp);
      SUNCheckLastErr();
      N_VLinearSum(ONE, vtemp, -R[(maa - 1) * maa + j], Q[j], vtemp);
      SUNCheckLastErr();
    }
    R[(maa - 1) * maa + maa - 1] = N_VDotProd(vtemp, vtemp);
    SUNCheckLastErr();
    R[(maa - 1) * maa + maa - 1] = SUNRsqrt(R[(maa - 1) * maa + maa - 1]);
    N_VScale((ONE / R[(maa - 1) * maa + maa - 1]), vtemp, Q[maa - 1]);
    SUNCheckLastErr();

    /* update the iteration map */
    j = 0;
    for (i = i_pt + 1; i < maa; i++) { ipt_map[j++] = i; }
    for (i = 0; i < i_pt + 1; i++) { ipt_map[j++] = i; }
  }

  /* solve least squares problem and update solution */
  lAA = iter;
  if (maa < iter) { lAA = maa; }
  SUNCheckCall(N_VDotProdMulti(lAA, fv, Q, gamma));

  /* set arrays for fused vector operation */
  cvals[0] = ONE;
  Xvecs[0] = gval;
  nvec     = 1;
  for (i = lAA - 1; i > -1; i--)
  {
    for (j = i + 1; j < lAA; j++) { gamma[i] -= R[j * maa + i] * gamma[j]; }
    if (gamma[i] == ZERO) { gamma[i] = ZERO; }
    else { gamma[i] /= R[i * maa + i]; }
    cvals[nvec] = -gamma[i];
    Xvecs[nvec] = dg[ipt_map[i]];
    nvec += 1;
  }

  /* if enabled, apply damping */
  if (damping)
  {
    onembeta    = (ONE - beta);
    cvals[nvec] = -onembeta;
    Xvecs[nvec] = fv;
    nvec += 1;
    for (i = lAA - 1; i > -1; i--)
    {
      cvals[nvec] = onembeta * gamma[i];
      Xvecs[nvec] = df[ipt_map[i]];
      nvec += 1;
    }
  }

  /* update solution */
  SUNCheckCall(N_VLinearCombination(nvec, cvals, Xvecs, x));

  return SUN_SUCCESS;
}

static SUNErrCode AllocateContent(SUNNonlinearSolver NLS, N_Vector y)
{
  SUNFunctionBegin(NLS->sunctx);
  int m = FP_CONTENT(NLS)->m;

  FP_CONTENT(NLS)->yprev = N_VClone(y);
  SUNCheckLastErr();

  FP_CONTENT(NLS)->gy = N_VClone(y);
  SUNCheckLastErr();

  FP_CONTENT(NLS)->delta = N_VClone(y);
  SUNCheckLastErr();

  /* Allocate all m-dependent content */
  if (m > 0)
  {
    FP_CONTENT(NLS)->fold = N_VClone(y);
    SUNCheckLastErr();

    FP_CONTENT(NLS)->gold = N_VClone(y);
    SUNCheckLastErr();

    FP_CONTENT(NLS)->imap = (int*)malloc(m * sizeof(int));
    SUNAssert(FP_CONTENT(NLS)->imap, SUN_ERR_MALLOC_FAIL);

    FP_CONTENT(NLS)->R = (sunrealtype*)malloc((m * m) * sizeof(sunrealtype));
    SUNAssert(FP_CONTENT(NLS)->R, SUN_ERR_MALLOC_FAIL);

    FP_CONTENT(NLS)->gamma = (sunrealtype*)malloc(m * sizeof(sunrealtype));
    SUNAssert(FP_CONTENT(NLS)->gamma, SUN_ERR_MALLOC_FAIL);

    FP_CONTENT(NLS)->cvals =
      (sunrealtype*)malloc(2 * (m + 1) * sizeof(sunrealtype));
    SUNAssert(FP_CONTENT(NLS)->cvals, SUN_ERR_MALLOC_FAIL);

    FP_CONTENT(NLS)->df = N_VCloneVectorArray(m, y);
    SUNCheckLastErr();

    FP_CONTENT(NLS)->dg = N_VCloneVectorArray(m, y);
    SUNCheckLastErr();

    FP_CONTENT(NLS)->q = N_VCloneVectorArray(m, y);
    SUNCheckLastErr();

    FP_CONTENT(NLS)->Xvecs = (N_Vector*)malloc(2 * (m + 1) * sizeof(N_Vector));
    SUNAssert(FP_CONTENT(NLS)->Xvecs, SUN_ERR_MALLOC_FAIL);
  }

  return SUN_SUCCESS;
}

static void FreeContent(SUNNonlinearSolver NLS)
{
  if (FP_CONTENT(NLS)->yprev)
  {
    N_VDestroy(FP_CONTENT(NLS)->yprev);
    FP_CONTENT(NLS)->yprev = NULL;
  }

  if (FP_CONTENT(NLS)->gy)
  {
    N_VDestroy(FP_CONTENT(NLS)->gy);
    FP_CONTENT(NLS)->gy = NULL;
  }

  if (FP_CONTENT(NLS)->fold)
  {
    N_VDestroy(FP_CONTENT(NLS)->fold);
    FP_CONTENT(NLS)->fold = NULL;
  }

  if (FP_CONTENT(NLS)->gold)
  {
    N_VDestroy(FP_CONTENT(NLS)->gold);
    FP_CONTENT(NLS)->gold = NULL;
  }

  if (FP_CONTENT(NLS)->delta)
  {
    N_VDestroy(FP_CONTENT(NLS)->delta);
    FP_CONTENT(NLS)->delta = NULL;
  }

  if (FP_CONTENT(NLS)->imap)
  {
    free(FP_CONTENT(NLS)->imap);
    FP_CONTENT(NLS)->imap = NULL;
  }

  if (FP_CONTENT(NLS)->R)
  {
    free(FP_CONTENT(NLS)->R);
    FP_CONTENT(NLS)->R = NULL;
  }

  if (FP_CONTENT(NLS)->gamma)
  {
    free(FP_CONTENT(NLS)->gamma);
    FP_CONTENT(NLS)->gamma = NULL;
  }

  if (FP_CONTENT(NLS)->cvals)
  {
    free(FP_CONTENT(NLS)->cvals);
    FP_CONTENT(NLS)->cvals = NULL;
  }

  if (FP_CONTENT(NLS)->df)
  {
    N_VDestroyVectorArray(FP_CONTENT(NLS)->df, FP_CONTENT(NLS)->m);
    FP_CONTENT(NLS)->df = NULL;
  }

  if (FP_CONTENT(NLS)->dg)
  {
    N_VDestroyVectorArray(FP_CONTENT(NLS)->dg, FP_CONTENT(NLS)->m);
    FP_CONTENT(NLS)->dg = NULL;
  }

  if (FP_CONTENT(NLS)->q)
  {
    N_VDestroyVectorArray(FP_CONTENT(NLS)->q, FP_CONTENT(NLS)->m);
    FP_CONTENT(NLS)->q = NULL;
  }

  if (FP_CONTENT(NLS)->Xvecs)
  {
    free(FP_CONTENT(NLS)->Xvecs);
    FP_CONTENT(NLS)->Xvecs = NULL;
  }

  return;
}
