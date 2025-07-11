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

SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, int max_iters, SUNContext sunctx)
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

  /* check for max_iters values; if illegal use defaults */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEst_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes          = SUNDomEigEst_SetATimes_PI;
  DEE->ops->setmaxiters        = SUNDomEigEst_SetMaxIters_PI;
  DEE->ops->settol             = SUNDomEigEst_SetTol_PI;
  DEE->ops->setnumpreprocess   = SUNDomEigEst_SetNumPreProcess_PI;
  DEE->ops->initialize         = SUNDomEigEst_Initialize_PI;
  DEE->ops->preprocess         = SUNDomEigEst_PreProcess_PI;
  DEE->ops->computehess        = NULL;
  DEE->ops->estimate           = SUNDomEig_Estimate_PI;
  DEE->ops->getnumofiters      = SUNDomEigEst_GetNumIters_PI;
  DEE->ops->getmaxnumofiters   = SUNDomEigEstGetMaxNumIters_PI;
  DEE->ops->getminnumofiters   = SUNDomEigEstGetMinNumIters_PI;
  DEE->ops->getres             = SUNDomEigEst_GetRes_PI;
  DEE->ops->printstats         = SUNDomEigEst_PrintStats_PI;
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
  content->max_iters = max_iters;
  content->curnumiters    = 0;
  content->maxnumiters = 0;
  content->minnumiters = 0;

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

SUNErrCode SUNDomEigEst_SetATimes_PI(SUNDomEigEstimator DEE, void* A_data,
                                    SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  PI_CONTENT(DEE)->ATimes = ATimes;
  PI_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Initialize_PI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (PI_CONTENT(DEE)->powiter_tol <= 0)
  {
    PI_CONTENT(DEE)->powiter_tol = DEE_TOL_DEFAULT;
  }
  if (PI_CONTENT(DEE)->numwarmups <= 0)
  {
    PI_CONTENT(DEE)->numwarmups = DEE_NUM_OF_WARMUPS_DEFAULT;
  }
  if (PI_CONTENT(DEE)->max_iters <= 0)
  {
    PI_CONTENT(DEE)->max_iters = DEE_MAX_ITER_DEFAULT;
  }

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  /* Note: this previously called N_VRandom (now removed) to initialize q, so PI_CONTENT(DEE)->q needs to be initialized in some other manner. */


  /* Initialize the vector V */
  sunrealtype normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetNumPreProcess_PI(SUNDomEigEstimator DEE,
                                            int numpreprocess)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the number of warmups */
  PI_CONTENT(DEE)->numwarmups = numpreprocess;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetTol_PI(SUNDomEigEstimator DEE, sunrealtype tol)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the tolerance */
  PI_CONTENT(DEE)->powiter_tol = tol;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetMaxIters_PI(SUNDomEigEstimator DEE, int max_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  /* Check for legal number of iters */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Set max iters */
  PI_CONTENT(DEE)->max_iters = max_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_PreProcess_PI(SUNDomEigEstimator DEE)
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

SUNErrCode SUNDomEig_Estimate_PI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT); // Maybe NULL, should we allow?
  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert((PI_CONTENT(DEE)->max_iters >= 0), SUN_ERR_ARG_CORRUPT);

  sunrealtype newlambdaR = ZERO;
  sunrealtype oldlambdaR = ZERO;

  sunindextype retval, k;
  sunrealtype normq;

  for (k = 0; k < PI_CONTENT(DEE)->max_iters; k++)
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

    PI_CONTENT(DEE)->res = SUNRabs(newlambdaR - oldlambdaR);

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
  PI_CONTENT(DEE)->curnumiters = k;
  PI_CONTENT(DEE)->maxnumiters =
    (k > PI_CONTENT(DEE)->maxnumiters) ? k : PI_CONTENT(DEE)->maxnumiters;
  PI_CONTENT(DEE)->minnumiters =
    (k < PI_CONTENT(DEE)->minnumiters) ? k : PI_CONTENT(DEE)->minnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetNumIters_PI(SUNDomEigEstimator DEE, int* niter)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *niter = PI_CONTENT(DEE)->curnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstGetMaxNumIters_PI(SUNDomEigEstimator DEE, int* maxniter)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *maxniter = PI_CONTENT(DEE)->maxnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstGetMinNumIters_PI(SUNDomEigEstimator DEE, int* minniter)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *minniter = PI_CONTENT(DEE)->minnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetRes_PI(SUNDomEigEstimator DEE, sunrealtype* res)
{
  SUNFunctionBegin(DEE->sunctx);
  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);

  *res = PI_CONTENT(DEE)->res;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_PrintStats_PI(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  if (DEE == NULL || outfile == NULL) { return SUN_ERR_ARG_CORRUPT; }

  // TODO: update each stat in the corresponding function
  fprintf(outfile, "Power Iteration DEE Statistics:\n");
  fprintf(outfile, "Max. num. of allowed iterations: %d\n", PI_CONTENT(DEE)->max_iters);
  fprintf(outfile, "Number of warmups: %d\n", PI_CONTENT(DEE)->numwarmups);
  fprintf(outfile, "Power iteration tolerance: %g\n", PI_CONTENT(DEE)->powiter_tol);
  fprintf(outfile, "Current num. of power iterations: %d\n", PI_CONTENT(DEE)->curnumiters);
  fprintf(outfile, "Current residual: %g\n", PI_CONTENT(DEE)->res);
  fprintf(outfile, "Max. Num. of power iterations: %d\n", PI_CONTENT(DEE)->maxnumiters);
  fprintf(outfile, "Min. Num. of power iterations: %d\n", PI_CONTENT(DEE)->minnumiters);

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
