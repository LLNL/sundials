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

#include <sundomeigest/sundomeigest_power.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Default estimator parameters */
#define DEE_NUM_OF_WARMUPS_PI_DEFAULT 0

/* Default Power Iteration parameters */
#define DEE_TOL_DEFAULT      SUN_RCONST(0.01)
#define DEE_MAX_ITER_DEFAULT 100

/*
 * -----------------------------------------------------------------
 * Power itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define PI_CONTENT(DEE) ((SUNDomEigEstimatorContent_Power)(DEE->content))

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new PI estimator
 */

SUNDomEigEstimator SUNDomEigEst_Power(N_Vector q, long int max_iters,
                                      sunrealtype rel_tol, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_Power content;

  /* check for legal q; if illegal return NULL */
  SUNAssertNull(!((q->ops->nvclone == NULL) || (q->ops->nvdestroy == NULL) ||
                  (q->ops->nvdotprod == NULL) || (q->ops->nvscale == NULL)),
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Check if q != 0 vector */
  SUNAssertNull(N_VDotProd(q, q) > SUN_SMALL_REAL, SUN_ERR_ARG_INCOMPATIBLE);

  /* check for max_iters values; if illegal use defaults */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Check if rel_tol > 0 */
  if (rel_tol < SUN_SMALL_REAL) { rel_tol = DEE_TOL_DEFAULT; }

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEst_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes         = SUNDomEigEst_SetATimes_Power;
  DEE->ops->setmaxiters       = SUNDomEigEst_SetMaxIters_Power;
  DEE->ops->settol            = SUNDomEigEst_SetRelTol_Power;
  DEE->ops->setnumpreprocess  = SUNDomEigEst_SetNumPreProcess_Power;
  DEE->ops->initialize        = SUNDomEigEst_Initialize_Power;
  DEE->ops->estimate          = SUNDomEig_Estimate_Power;
  DEE->ops->getcurres         = SUNDomEigEst_GetRes_Power;
  DEE->ops->getcurniters      = SUNDomEigEst_GetNumIters_Power;
  DEE->ops->getnumatimescalls = SUNDomEigEst_GetNumATimesCalls_Power;
  DEE->ops->write             = SUNDomEigEst_Write_Power;
  DEE->ops->destroy           = SUNDomEigEst_Destroy_Power;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_Power)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes        = NULL;
  content->ATdata        = NULL;
  content->V             = NULL;
  content->q             = NULL;
  content->max_iters     = max_iters;
  content->num_warmups   = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  content->powiter_tol   = rel_tol;
  content->cur_res       = ZERO;
  content->cur_num_iters = 0;
  content->num_ATimes    = 0;

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

SUNErrCode SUNDomEigEst_SetATimes_Power(SUNDomEigEstimator DEE, void* A_data,
                                        SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);  
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT); 
  SUNAssert(ATimes, SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  PI_CONTENT(DEE)->ATimes = ATimes;
  PI_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Initialize_Power(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->powiter_tol < SUN_SMALL_REAL)
  {
    PI_CONTENT(DEE)->powiter_tol = DEE_TOL_DEFAULT;
  }
  if (PI_CONTENT(DEE)->num_warmups < 0)
  {
    PI_CONTENT(DEE)->num_warmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  }
  if (PI_CONTENT(DEE)->max_iters <= 0)
  {
    PI_CONTENT(DEE)->max_iters = DEE_MAX_ITER_DEFAULT;
  }

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
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

SUNErrCode SUNDomEigEst_SetNumPreProcess_Power(SUNDomEigEstimator DEE,
                                               int num_warmups)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if num_warmups >= 0 */
  if (num_warmups < 0) { num_warmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;}

  /* set the number of warmups */
  PI_CONTENT(DEE)->num_warmups = num_warmups;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetRelTol_Power(SUNDomEigEstimator DEE, sunrealtype tol)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if rel_tol > 0 */
  if (tol < SUN_SMALL_REAL) { tol = DEE_TOL_DEFAULT; }

  /* set the tolerance */
  PI_CONTENT(DEE)->powiter_tol = tol;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetMaxIters_Power(SUNDomEigEstimator DEE,
                                          long int max_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check for legal number of iters */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Set max iters */
  PI_CONTENT(DEE)->max_iters = max_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEig_Estimate_Power(SUNDomEigEstimator DEE,
                                    sunrealtype* lambdaR, sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert((PI_CONTENT(DEE)->max_iters >= 0), SUN_ERR_ARG_CORRUPT);

  sunrealtype newlambdaR = ZERO;
  sunrealtype oldlambdaR = ZERO;

  int retval;
  sunrealtype normq;
  PI_CONTENT(DEE)->num_ATimes = 0;

  for (int i = 0; i < PI_CONTENT(DEE)->num_warmups; i++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    if (retval != 0)
    {
      return SUN_ERR_USER_FCN_FAIL;
    }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();
  }

  int k; // k will be used out of the loop as a counter
  for (k = 0; k < PI_CONTENT(DEE)->max_iters; k++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    if (retval != 0)
    {
      return SUN_ERR_USER_FCN_FAIL;
    }

    newlambdaR = N_VDotProd(PI_CONTENT(DEE)->V,
                            PI_CONTENT(DEE)->q); //Rayleigh quotient
    SUNCheckLastErr();

    PI_CONTENT(DEE)->cur_res = SUNRabs(newlambdaR - oldlambdaR) /
                               SUNRabs(newlambdaR);

    if (PI_CONTENT(DEE)->cur_res < PI_CONTENT(DEE)->powiter_tol) { break; }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();

    oldlambdaR = newlambdaR;
  }

  *lambdaI                       = ZERO;
  *lambdaR                       = newlambdaR;

  /* Set the current number of iterations */
  PI_CONTENT(DEE)->cur_num_iters = k + PI_CONTENT(DEE)->num_warmups + 1;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetRes_Power(SUNDomEigEstimator DEE,
                                        sunrealtype* cur_res)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(cur_res, SUN_ERR_ARG_CORRUPT);

  *cur_res = PI_CONTENT(DEE)->cur_res;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetNumIters_Power(SUNDomEigEstimator DEE,
                                             long int* curniter)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(curniter, SUN_ERR_ARG_CORRUPT);

  *curniter = PI_CONTENT(DEE)->cur_num_iters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetNumATimesCalls_Power(SUNDomEigEstimator DEE,
                                                long int* num_ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_ATimes, SUN_ERR_ARG_CORRUPT);

  *num_ATimes = PI_CONTENT(DEE)->num_ATimes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Write_Power(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(outfile, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (DEE == NULL || outfile == NULL) { return SUN_ERR_ARG_CORRUPT; }

  fprintf(outfile, "\nPower Iteration DEE Statistics:");
  fprintf(outfile, "\n------------------------------------------------\n");
  fprintf(outfile, "Max. num. of iters allowed    = %ld\n",
          PI_CONTENT(DEE)->max_iters);
  fprintf(outfile, "Num. of warmups               = %d\n",
          PI_CONTENT(DEE)->num_warmups);
  fprintf(outfile, "Power iteration tolerance     = " SUN_FORMAT_G "\n",
          PI_CONTENT(DEE)->powiter_tol);
  fprintf(outfile, "Cur. residual                 = " SUN_FORMAT_G "\n",
          PI_CONTENT(DEE)->cur_res);
  fprintf(outfile, "Cur. num. of power iters      = %ld\n",
          PI_CONTENT(DEE)->cur_num_iters);
  fprintf(outfile, "Num. of ATimes calls          = %ld\n\n",
          PI_CONTENT(DEE)->num_ATimes);

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Destroy_Power(SUNDomEigEstimator* DEEptr)
{
  SUNFunctionBegin((*DEEptr)->sunctx);

  if ((*DEEptr) == NULL) { return SUN_SUCCESS; }

  SUNDomEigEstimator DEE = *DEEptr;

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
  *DEEptr = NULL;
  return SUN_SUCCESS;
}
