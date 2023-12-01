/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNLINEARSOLVER
 * package.  It contains the implementation of the SUNLinearSolver
 * operations listed in sundials_linearsolver.h
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_iterative_impl.h"
#include "sundials_logger_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNLinearSolver S)
{
  return (S->sunctx->profiler);
}
#endif

/* -----------------------------------------------------------------
 * Create a new empty SUNLinearSolver object
 * ----------------------------------------------------------------- */

SUNLinearSolver SUNLinSolNewEmpty(SUNContext sunctx)
{
  SUNLinearSolver LS;
  SUNLinearSolver_Ops ops;

  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);

  /* create linear solver object */
  LS = NULL;
  LS = (SUNLinearSolver)malloc(sizeof *LS);
  SUNAssert(LS, SUN_ERR_MALLOC_FAIL);

  /* create linear solver ops structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops)malloc(sizeof *ops);
  SUNAssert(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype           = NULL;
  ops->getid             = NULL;
  ops->setatimes         = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->setzeroguess      = NULL;
  ops->initialize        = NULL;
  ops->setup             = NULL;
  ops->solve             = NULL;
  ops->numiters          = NULL;
  ops->resnorm           = NULL;
  ops->resid             = NULL;
  ops->lastflag          = NULL;
  ops->space             = NULL;
  ops->free              = NULL;

  /* attach ops and initialize content and context to NULL */
  LS->ops     = ops;
  LS->content = NULL;
  LS->sunctx  = sunctx;

  return (LS);
}

/* -----------------------------------------------------------------
 * Free a generic SUNLinearSolver (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNLinSolFreeEmpty(SUNLinearSolver S)
{
  if (S == NULL) { return; }

  /* free non-NULL ops structure */
  if (S->ops) { free(S->ops); }
  S->ops = NULL;

  /* free overall N_Vector object and return */
  free(S);
  return;
}

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver S)
{
  return (S->ops->gettype(S));
}

SUNLinearSolver_ID SUNLinSolGetID(SUNLinearSolver S)
{
  if (S->ops->getid) { return (S->ops->getid(S)); }
  else { return (SUNLINEARSOLVER_CUSTOM); }
}

SUNErrCode SUNLinSolSetATimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setatimes) { ier = S->ops->setatimes(S, A_data, ATimes); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                      SUNPSetupFn Pset, SUNPSolveFn Psol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setpreconditioner)
  {
    ier = S->ops->setpreconditioner(S, P_data, Pset, Psol);
  }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setscalingvectors) { ier = S->ops->setscalingvectors(S, s1, s2); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetZeroGuess(SUNLinearSolver S, sunbooleantype onoff)
{
  if (S->ops->setzeroguess) { return (S->ops->setzeroguess(S, onoff)); }
  else { return SUN_SUCCESS; }
}

SUNErrCode SUNLinSolInitialize(SUNLinearSolver S)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->initialize) { ier = S->ops->initialize(S); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setup) { ier = S->ops->setup(S, A); }
  else { ier = SUNLS_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                   sunrealtype tol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  ier = S->ops->solve(S, A, x, b, tol);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolNumIters(SUNLinearSolver S)
{
  int result;
  if (S->ops->numiters) { result = S->ops->numiters(S); }
  else { result = 0; }
  return (result);
}

sunrealtype SUNLinSolResNorm(SUNLinearSolver S)
{
  sunrealtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->resnorm) { result = S->ops->resnorm(S); }
  else { result = SUN_RCONST(0.0); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (result);
}

N_Vector SUNLinSolResid(SUNLinearSolver S)
{
  N_Vector resid;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->resid) { resid = S->ops->resid(S); }
  else { resid = NULL; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (resid);
}

sunindextype SUNLinSolLastFlag(SUNLinearSolver S)
{
  if (S->ops->lastflag) { return ((sunindextype)S->ops->lastflag(S)); }
  else { return 0; }
}

SUNErrCode SUNLinSolSpace(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS)
{
  if (S->ops->space) { return (S->ops->space(S, lenrwLS, leniwLS)); }
  else
  {
    *lenrwLS = 0;
    *leniwLS = 0;
    return SUNLS_SUCCESS;
  }
}

SUNErrCode SUNLinSolFree(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  /* if the free operation exists use it */
  if (S->ops)
  {
    if (S->ops->free) { return (S->ops->free(S)); }
  }

  /* if we reach this point, either ops == NULL or free == NULL,
     try to cleanup by freeing the content, ops, and solver */
  if (S->content)
  {
    free(S->content);
    S->content = NULL;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;

  return SUN_SUCCESS;
}
