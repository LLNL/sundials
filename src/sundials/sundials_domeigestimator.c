/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNDomEigEst
 * package. It contains the implementation of the SUNDomEigEstimator
 * operations listed in sundials_domeigestimator.h
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include <sundials/priv/sundials_domeigestimator_impl.h>
#include <sundials/sundials_domeigestimator.h>

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNDomEigEstimator DEE)
{
  return (DEE->sunctx->profiler);
}
#endif

/* -----------------------------------------------------------------
 * Create a new empty SUNDomEigEstimator object
 * ----------------------------------------------------------------- */

SUNDomEigEstimator SUNDomEigEst_NewEmpty(SUNContext sunctx)
{
  SUNDomEigEstimator DEE;
  SUNDomEigEstimator_Ops ops;

  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);

  /* create dominant eigenvalue estimator object */
  DEE = NULL;
  DEE = (SUNDomEigEstimator)malloc(sizeof *DEE);
  SUNAssertNull(DEE, SUN_ERR_MALLOC_FAIL);

  /* create dominant eigenvalue estimator ops structure */
  ops = NULL;
  ops = (SUNDomEigEstimator_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->setatimes        = NULL;
  ops->setmaxiters      = NULL;
  ops->settol           = NULL;
  ops->setnumpreprocess = NULL;
  ops->initialize       = NULL;
  ops->preprocess       = NULL;
  ops->computehess      = NULL;
  ops->estimate         = NULL;
  ops->getcurniters     = NULL;
  ops->getmaxniters     = NULL;
  ops->getminniters     = NULL;
  ops->getcurres        = NULL;
  ops->printstats       = NULL;
  ops->free             = NULL;

  /* attach ops and initialize content and context to NULL */
  DEE->ops     = ops;
  DEE->content = NULL;
  DEE->sunctx  = sunctx;

  return (DEE);
}

/* -----------------------------------------------------------------
 * Free a generic SUNDomEigEstimator (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNDomEigEst_FreeEmpty(SUNDomEigEstimator DEE)
{
  if (DEE == NULL) { return; }

  /* free non-NULL ops structure */
  if (DEE->ops) { free(DEE->ops); }
  DEE->ops = NULL;

  /* free overall N_Vector object and return */
  free(DEE);
  return;
}

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------*/

SUNErrCode SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* A_data,
                                  SUNATimesFn ATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setatimes) { ier = DEE->ops->setatimes(DEE, A_data, ATimes); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, int max_iters)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setmaxiters) { ier = DEE->ops->setmaxiters(DEE, max_iters); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE, int numpreprocess)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setnumpreprocess)
  {
    ier = DEE->ops->setnumpreprocess(DEE, numpreprocess);
  }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_SetTol(SUNDomEigEstimator DEE, sunrealtype tol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->settol) { ier = DEE->ops->settol(DEE, tol); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_Initialize(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->initialize) { ier = DEE->ops->initialize(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_PreProcess(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->preprocess) { ier = DEE->ops->preprocess(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_ComputeHess(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->computehess) { ier = DEE->ops->computehess(DEE); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                              sunrealtype* lambdaI)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->estimate) { ier = DEE->ops->estimate(DEE, lambdaR, lambdaI); }
  else { ier = SUN_ERR_DEE_NULL_ESTIMATE; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetCurRes(SUNDomEigEstimator DEE, sunrealtype* curres)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getcurres) { ier = DEE->ops->getcurres(DEE, curres); }
  else
  {
    *curres = SUN_RCONST(0.0);
    ier     = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetCurNumIters(SUNDomEigEstimator DEE, int* curniter)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getcurniters) { ier = DEE->ops->getcurniters(DEE, curniter); }
  else
  {
    *curniter = 0;
    ier       = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetMaxNumIters(SUNDomEigEstimator DEE, int* maxniter)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getmaxniters) { ier = DEE->ops->getmaxniters(DEE, maxniter); }
  else
  {
    *maxniter = 0;
    ier       = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetMinNumIters(SUNDomEigEstimator DEE, int* minniter)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getminniters) { ier = DEE->ops->getminniters(DEE, minniter); }
  else
  {
    *minniter = 0;
    ier       = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                          long int* nATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getnumatimescalls)
  {
    ier = DEE->ops->getnumatimescalls(DEE, nATimes);
  }
  else
  {
    *nATimes = 0;
    ier      = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_PrintStats(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->printstats) { ier = DEE->ops->printstats(DEE, outfile); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEstFree(SUNDomEigEstimator DEE)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->free) { ier = DEE->ops->free(DEE); }
  else { ier = SUN_ERR_DEE_NULL_FREE; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}
