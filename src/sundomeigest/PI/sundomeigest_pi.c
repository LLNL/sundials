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
  DEE->ops->getcurres          = SUNDomEigEst_GetCurRes_PI;
  DEE->ops->getcurniters       = SUNDomEigEst_GetCurNumIters_PI;
  DEE->ops->getmaxniters       = SUNDomEigEst_GetMaxNumIters_PI;
  DEE->ops->getminniters       = SUNDomEigEst_GetMinNumIters_PI;
  DEE->ops->getnumatimescalls  = SUNDomEigEst_GetNumATimesCalls_PI;
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
  content->numwarmups  = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  content->max_iters   = max_iters;
  content->powiter_tol = ZERO;
  content->curres      = ZERO;
  content->curnumiters = 0;
  content->maxnumiters = 0;
  content->minnumiters = max_iters;
  content->nATimes     = 0;
  
  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  N_VScale(ONE, q, content->q);
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

  if (PI_CONTENT(DEE)->powiter_tol < SUN_SMALL_REAL)
  {
    PI_CONTENT(DEE)->powiter_tol = DEE_TOL_DEFAULT;
  }
  if (PI_CONTENT(DEE)->numwarmups < 0)
  {
    PI_CONTENT(DEE)->numwarmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  }
  if (PI_CONTENT(DEE)->max_iters <= 0)
  {
    PI_CONTENT(DEE)->max_iters = DEE_MAX_ITER_DEFAULT;
  }

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

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

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  sunrealtype normq;
  int i, retval;

  /* Set the initial q = A^{numwarmups}q/||A^{numwarmups}q|| */
  for (i = 0; i < PI_CONTENT(DEE)->numwarmups; i++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->nATimes++;
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

  int retval, k;
  sunrealtype normq;

  for (k = 0; k < PI_CONTENT(DEE)->max_iters; k++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->nATimes++;
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DEE_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DEE_ATIMES_FAIL_REC; }
    }

    newlambdaR = N_VDotProd(PI_CONTENT(DEE)->V,
                            PI_CONTENT(DEE)->q); //Rayleigh quotient
    SUNCheckLastErr();

    PI_CONTENT(DEE)->curres = SUNRabs(newlambdaR - oldlambdaR)/SUNRabs(newlambdaR);

    if (PI_CONTENT(DEE)->curres < PI_CONTENT(DEE)->powiter_tol) { break; }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();

    oldlambdaR = newlambdaR;
  }

  k++;
  *lambdaI                  = ZERO;
  *lambdaR                  = newlambdaR;
  PI_CONTENT(DEE)->curnumiters = k;
  PI_CONTENT(DEE)->maxnumiters =
    (k > PI_CONTENT(DEE)->maxnumiters) ? k : PI_CONTENT(DEE)->maxnumiters;
  PI_CONTENT(DEE)->minnumiters =
    (k < PI_CONTENT(DEE)->minnumiters) ? k : PI_CONTENT(DEE)->minnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetCurRes_PI(SUNDomEigEstimator DEE, sunrealtype* curres)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);
  SUNAssert(curres, SUN_ERR_ARG_CORRUPT);

  *curres = PI_CONTENT(DEE)->curres;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetCurNumIters_PI(SUNDomEigEstimator DEE, int* curniter)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);
  SUNAssert(curniter, SUN_ERR_ARG_CORRUPT);

  *curniter = PI_CONTENT(DEE)->curnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetMaxNumIters_PI(SUNDomEigEstimator DEE, int* maxniter)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);
  SUNAssert(maxniter, SUN_ERR_ARG_CORRUPT);

  *maxniter = PI_CONTENT(DEE)->maxnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetMinNumIters_PI(SUNDomEigEstimator DEE, int* minniter)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);
  SUNAssert(minniter, SUN_ERR_ARG_CORRUPT);

  *minniter = PI_CONTENT(DEE)->minnumiters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetNumATimesCalls_PI(SUNDomEigEstimator DEE,
                                            long int* nATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_DEE_NULL_MEM);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_DEE_NULL_CONTENT);
  SUNAssert(nATimes, SUN_ERR_ARG_CORRUPT);

  *nATimes = PI_CONTENT(DEE)->nATimes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_PrintStats_PI(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  if (DEE == NULL || outfile == NULL) { return SUN_ERR_ARG_CORRUPT; }

  fprintf(outfile, "\nPower Iteration DEE Statistics:");
  fprintf(outfile, "\n------------------------------------------------\n");
  fprintf(outfile, "Max. num. of iters allowed    = %d\n", PI_CONTENT(DEE)->max_iters);
  fprintf(outfile, "Num. of warmups               = %d\n", PI_CONTENT(DEE)->numwarmups);
  fprintf(outfile, "Power iteration tolerance     = " SUN_FORMAT_G "\n", PI_CONTENT(DEE)->powiter_tol);
  fprintf(outfile, "Cur. residual                 = " SUN_FORMAT_G "\n", PI_CONTENT(DEE)->curres);
  fprintf(outfile, "Cur. num. of power iters      = %d\n", PI_CONTENT(DEE)->curnumiters);
  fprintf(outfile, "Max. num. of power iters      = %d\n", PI_CONTENT(DEE)->maxnumiters);
  fprintf(outfile, "Min. num. of power iters      = %d\n", PI_CONTENT(DEE)->minnumiters);
  fprintf(outfile, "Num. of ATimes calls          = %ld\n", PI_CONTENT(DEE)->nATimes);

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
