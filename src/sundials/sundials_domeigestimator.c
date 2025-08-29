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
  if (sunctx == NULL) { return NULL; }
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimator_Ops ops;

  /* create dominant eigenvalue estimator object */
  DEE = NULL;
  DEE = (SUNDomEigEstimator)malloc(sizeof *DEE);
  SUNAssertNull(DEE, SUN_ERR_MALLOC_FAIL);

  /* create dominant eigenvalue estimator ops structure */
  ops = NULL;
  ops = (SUNDomEigEstimator_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->setatimes             = NULL;
  ops->setmaxiters           = NULL;
  ops->setreltol             = NULL;
  ops->setnumpreprocessiters = NULL;
  ops->initialize            = NULL;
  ops->estimate              = NULL;
  ops->getnumiters           = NULL;
  ops->getres                = NULL;
  ops->write                 = NULL;
  ops->destroy               = NULL;

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

  /* free overall SUNDomEigEstimator object and return */
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

SUNErrCode SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, long int max_iters)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setmaxiters) { ier = DEE->ops->setmaxiters(DEE, max_iters); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_SetNumPreprocessIters(SUNDomEigEstimator DEE,
                                              int num_iters)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setnumpreprocessiters)
  {
    ier = DEE->ops->setnumpreprocessiters(DEE, num_iters);
  }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_SetRelTol(SUNDomEigEstimator DEE, sunrealtype tol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->setreltol) { ier = DEE->ops->setreltol(DEE, tol); }
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

SUNErrCode SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                              sunrealtype* lambdaI)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->estimate) { ier = DEE->ops->estimate(DEE, lambdaR, lambdaI); }
  else { ier = SUN_ERR_NOT_IMPLEMENTED; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetRes(SUNDomEigEstimator DEE, sunrealtype* cur_res)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getres) { ier = DEE->ops->getres(DEE, cur_res); }
  else
  {
    *cur_res = SUN_RCONST(0.0);
    ier      = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetNumIters(SUNDomEigEstimator DEE, long int* curniter)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getnumiters) { ier = DEE->ops->getnumiters(DEE, curniter); }
  else
  {
    *curniter = 0;
    ier       = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                          long int* num_ATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->getnumatimescalls)
  {
    ier = DEE->ops->getnumatimescalls(DEE, num_ATimes);
  }
  else
  {
    *num_ATimes = 0;
    ier         = SUN_SUCCESS;
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_Write(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(DEE));
  if (DEE->ops->write) { ier = DEE->ops->write(DEE, outfile); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(DEE));
  return (ier);
}

SUNErrCode SUNDomEigEst_Destroy(SUNDomEigEstimator* DEEptr)
{
  SUNDomEigEstimator DEE = *DEEptr;
  SUNErrCode ier;
  if (DEE->ops->destroy) { ier = DEE->ops->destroy(&DEE); }
  else
  {
    SUNDomEigEst_FreeEmpty(DEE);
    *DEEptr = NULL;
    ier     = SUN_SUCCESS;
  }
  return (ier);
}
