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
 * This is the implementation file for the Power Iteration (PI)
 * implementation of the SUNDomEigEst package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_domeigestimator_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sundomeigest/sundomeigest_pi.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Power itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define PI_CONTENT(DEE) ((SUNDomEigEstimatorContent_PI)(DEE->content))

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new PI estimator
 */

SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, int max_powiter, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_PI content;

  /* check for legal q; if illegal return NULL */
  // TO DO: check required vector operations
  SUNAssertNull(!((q->ops->nvclone == NULL) || (q->ops->nvdestroy == NULL) ||
                  (q->ops->nvdotprod == NULL) || (q->ops->nvscale == NULL) ||
                  (q->ops->nvgetlength == NULL) || (q->ops->nvspace == NULL)),
                SUN_ERR_DEE_BAD_NVECTOR);

  /* check for max_powiter values; if illegal use defaults */
  if (max_powiter <= 0) { max_powiter = SUNDOMEIGEST_MAX_PI_DEFAULT; }

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->getid              = SUNDomEigEst_PIGetID;
  DEE->ops->setatimes          = SUNDomEigEstSetATimes_PI;
  DEE->ops->setmaxpoweriter    = SUNDomEigEst_PISetMaxPowerIter;
  DEE->ops->settol             = SUNDomEigEstSetTol_PI;
  DEE->ops->setnumofperprocess = SUNDomEigEstSetNumPreProcess_PI;
  DEE->ops->initialize         = SUNDomEigEstInitialize_PI;
  DEE->ops->preprocess         = SUNDomEigEstPreProcess_PI;
  DEE->ops->computehess        = NULL;
  DEE->ops->estimate           = SUNDomEigEstimate_PI;
  DEE->ops->getnumofiters      = SUNDomEigEstNumIters_PI;
  DEE->ops->getres             = SUNDomEigEstRes_PI;
  DEE->ops->free               = SUNDomEigEstFree_PI;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_PI)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes      = NULL;
  content->ATdata      = NULL;
  content->V           = NULL;
  content->q           = NULL;
  content->numwarmups  = 0;
  content->powiter_tol = ZERO;
  content->res         = ZERO;
  content->max_powiter = max_powiter;
  content->numiters    = 0;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  content->V = N_VClone(q);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNDomEigEstimator_ID SUNDomEigEst_PIGetID(
  SUNDIALS_MAYBE_UNUSED SUNDomEigEstimator DEE)
{
  return (SUNDSOMEIGESTIMATOR_POWER);
}

SUNErrCode SUNDomEigEstInitialize_PI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (PI_CONTENT(DEE)->powiter_tol <= 0)
  {
    PI_CONTENT(DEE)->powiter_tol = SUNDOMEIGEST_PI_TOL_DEFAULT;
  }
  if (PI_CONTENT(DEE)->numwarmups <= 0)
  {
    PI_CONTENT(DEE)->numwarmups = DEE_NUM_OF_WARMUPS_DEFAULT;
  }
  if (PI_CONTENT(DEE)->max_powiter <= 0)
  {
    PI_CONTENT(DEE)->max_powiter = SUNDOMEIGEST_MAX_PI_DEFAULT;
  }

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  N_VRandom(PI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  /* Initialize the vector V */
  sunrealtype normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstSetNumPreProcess_PI(SUNDomEigEstimator DEE,
                                           int numofperprocess)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the number of warmups */
  PI_CONTENT(DEE)->numwarmups = numofperprocess;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstSetTol_PI(SUNDomEigEstimator DEE, sunrealtype tol)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the tolerance */
  PI_CONTENT(DEE)->powiter_tol = tol;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstSetATimes_PI(SUNDomEigEstimator DEE, void* A_data,
                                    SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  PI_CONTENT(DEE)->ATimes = ATimes;
  PI_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_PISetMaxPowerIter(SUNDomEigEstimator DEE, int max_powiter)
{
  SUNFunctionBegin(DEE->sunctx);

  /* Check for legal number of iters */
  if (max_powiter <= 0) { max_powiter = SUNDOMEIGEST_MAX_PI_DEFAULT; }

  /* Set max iters */
  PI_CONTENT(DEE)->max_powiter = max_powiter;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstPreProcess_PI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  sunrealtype normq;
  sunindextype i, retval;

  /* Set the initial q = A^{numwarmups}q/||A^{numwarmups}q|| */
  for (i = 0; i < PI_CONTENT(DEE)->numwarmups; i++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DEE_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DEE_ATIMES_FAIL_REC; }
    }
    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimate_PI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT); // Maybe NULL, should we allow?
  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert((PI_CONTENT(DEE)->max_powiter >= 0), SUN_ERR_ARG_CORRUPT);

  sunrealtype newlambdaR = ZERO;
  sunrealtype oldlambdaR = ZERO;

  sunindextype retval, k;
  sunrealtype normq;

  for (k = 0; k < PI_CONTENT(DEE)->max_powiter; k++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DEE_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DEE_ATIMES_FAIL_REC; }
    }

    newlambdaR = N_VDotProd(PI_CONTENT(DEE)->V,
                            PI_CONTENT(DEE)->q); //Rayleigh quotient
    SUNCheckLastErr();

    PI_CONTENT(DEE)->res = fabs(newlambdaR - oldlambdaR);

    if (PI_CONTENT(DEE)->res < PI_CONTENT(DEE)->powiter_tol) { break; }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();

    oldlambdaR = newlambdaR;
  }

  *lambdaI                  = ZERO;
  *lambdaR                  = newlambdaR;
  PI_CONTENT(DEE)->numiters = k;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstNumIters_PI(SUNDomEigEstimator DEE, int* niter)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *niter = PI_CONTENT(DEE)->numiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstRes_PI(SUNDomEigEstimator DEE, sunrealtype* res)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *res = PI_CONTENT(DEE)->res;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstFree_PI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (DEE == NULL) { return SUN_SUCCESS; }

  if (DEE->content)
  {
    /* delete items from within the content structure */
    if (PI_CONTENT(DEE)->q)
    {
      N_VDestroy(PI_CONTENT(DEE)->q);
      PI_CONTENT(DEE)->q = NULL;
    }
    if (PI_CONTENT(DEE)->V)
    {
      N_VDestroy(PI_CONTENT(DEE)->V);
      PI_CONTENT(DEE)->V = NULL;
    }

    free(DEE->content);
    DEE->content = NULL;
  }
  if (DEE->ops)
  {
    free(DEE->ops);
    DEE->ops = NULL;
  }
  free(DEE);
  DEE = NULL;
  return SUN_SUCCESS;
}
