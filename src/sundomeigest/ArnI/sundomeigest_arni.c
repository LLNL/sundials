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
 * This is the implementation file for the Arnoldi Iteration (ArnI)
 * implementation of the SUNDomEigEst package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_domeigestimator_impl.h>
#include <sundomeigest/sundomeigest_arni.h>

#include "sundials_lapack_defs.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgeev_f77 dgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgeev_f77 sgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Arnoldi itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define ArnI_CONTENT(DEE) ((SUNDomEigEstimatorContent_ArnI)(DEE->content))

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new ArnI estimator
 */

SUNDomEigEstimator SUNDomEigEst_ArnI(N_Vector q, int krydim, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_ArnI content;

  /* Check if krydim >= 2 */
  if (krydim < 3) { krydim = DEE_KRYLOV_DIM_DEFAULT; }

  /* check for legal q; if illegal return NULL */
  // TO DO: check required vector operations
  SUNAssertNull(!((q->ops->nvclone == NULL) || (q->ops->nvdestroy == NULL) ||
                  (q->ops->nvdotprod == NULL) || (q->ops->nvscale == NULL) ||
                  (q->ops->nvgetlength == NULL) || (q->ops->nvspace == NULL)),
                SUN_ERR_DEE_BAD_NVECTOR);

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->getid              = SUNDomEigEst_ArnIGetID;
  DEE->ops->setatimes          = SUNDomEigEstSetATimes_ArnI;
  DEE->ops->setmaxpoweriter    = NULL;
  DEE->ops->settol             = NULL;
  DEE->ops->setnumofperprocess = SUNDomEigEstSetNumPreProcess_ArnI;
  DEE->ops->initialize         = SUNDomEigEstInitialize_ArnI;
  DEE->ops->preprocess         = SUNDomEigEstPreProcess_ArnI;
  DEE->ops->computehess        = SUNDomEigEstComputeHess_ArnI;
  DEE->ops->estimate           = SUNDomEigEstimate_ArnI;
  DEE->ops->getnumofiters      = NULL;
  DEE->ops->getres             = NULL;
  DEE->ops->free               = SUNDomEigEstFree_ArnI;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_ArnI)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes      = NULL;
  content->ATdata      = NULL;
  content->V           = NULL;
  content->q           = NULL;
  content->krydim      = krydim;
  content->numwarmups  = 0;
  content->LAPACK_A    = NULL;
  content->LAPACK_wr   = NULL;
  content->LAPACK_wi   = NULL;
  content->LAPACK_work = NULL;
  content->LAPACK_arr  = NULL;
  content->Hes         = NULL;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  content->V = N_VCloneVectorArray(krydim + 1, q);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNDomEigEstimator_ID SUNDomEigEst_ArnIGetID(
  SUNDIALS_MAYBE_UNUSED SUNDomEigEstimator DEE)
{
  return (SUNDSOMEIGESTIMATOR_ARNOLDI);
}

SUNErrCode SUNDomEigEstSetATimes_ArnI(SUNDomEigEstimator DEE, void* A_data,
                                      SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  ArnI_CONTENT(DEE)->ATimes = ATimes;
  ArnI_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstInitialize_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (ArnI_CONTENT(DEE)->krydim < 2)
  {
    ArnI_CONTENT(DEE)->krydim = DEE_KRYLOV_DIM_DEFAULT;
  }
  if (ArnI_CONTENT(DEE)->numwarmups < 0)
  {
    ArnI_CONTENT(DEE)->numwarmups = DEE_NUM_OF_WARMUPS_DEFAULT;
  }

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  ArnI_CONTENT(DEE)->LAPACK_A = (sunrealtype*)malloc(
    (ArnI_CONTENT(DEE)->krydim * ArnI_CONTENT(DEE)->krydim) * sizeof(sunrealtype));
  if (ArnI_CONTENT(DEE)->LAPACK_A == NULL) { return SUN_ERR_MALLOC_FAIL; }
  ArnI_CONTENT(DEE)->LAPACK_wr =
    malloc(ArnI_CONTENT(DEE)->krydim * sizeof(sunrealtype));
  if (ArnI_CONTENT(DEE)->LAPACK_wr == NULL) { return SUN_ERR_MALLOC_FAIL; }
  ArnI_CONTENT(DEE)->LAPACK_wi =
    malloc(ArnI_CONTENT(DEE)->krydim * sizeof(sunrealtype));
  if (ArnI_CONTENT(DEE)->LAPACK_wi == NULL) { return SUN_ERR_MALLOC_FAIL; }

  /* query the workspace size (call with lwork = -1) */
  char jobvl = 'N';
  char jobvr = 'N';
  sunindextype N    = ArnI_CONTENT(DEE)->krydim;
  sunindextype lda  = ArnI_CONTENT(DEE)->krydim;
  sunindextype ldvl = 1;
  sunindextype ldvr = 1;
  sunindextype info = 0;
  sunindextype lwork = -1;
  sunrealtype work = SUN_RCONST(0.0);

  xgeev_f77(&jobvl, &jobvr, &N, ArnI_CONTENT(DEE)->LAPACK_A, &lda,
            ArnI_CONTENT(DEE)->LAPACK_wr, ArnI_CONTENT(DEE)->LAPACK_wi, NULL,
            &ldvl, NULL, &ldvr, &work, &lwork, &info);

  ArnI_CONTENT(DEE)->LAPACK_work =
    (sunrealtype*)malloc(((sunindextype)work) * sizeof(sunrealtype));
  if (ArnI_CONTENT(DEE)->LAPACK_work == NULL) { return SUN_ERR_MALLOC_FAIL; }

  int k;

  /* LAPACK array */
  if (ArnI_CONTENT(DEE)->LAPACK_arr == NULL)
  {
    ArnI_CONTENT(DEE)->LAPACK_arr =
      (sunrealtype**)malloc(ArnI_CONTENT(DEE)->krydim * sizeof(sunrealtype*));

    for (k = 0; k < ArnI_CONTENT(DEE)->krydim; k++)
    {
      ArnI_CONTENT(DEE)->LAPACK_arr[k] =
        (sunrealtype*)malloc(2 * sizeof(sunrealtype));
    }
  }

  /* Hessenberg matrix Hes */
  if (ArnI_CONTENT(DEE)->Hes == NULL)
  {
    ArnI_CONTENT(DEE)->Hes = (sunrealtype**)malloc(
      (ArnI_CONTENT(DEE)->krydim + 1) * sizeof(sunrealtype*));

    for (k = 0; k <= ArnI_CONTENT(DEE)->krydim; k++)
    {
      ArnI_CONTENT(DEE)->Hes[k] =
        (sunrealtype*)malloc(ArnI_CONTENT(DEE)->krydim * sizeof(sunrealtype));
    }
  }

  /* Initialize the vector V[0] */
  N_VRandom(ArnI_CONTENT(DEE)->q);
  SUNCheckLastErr();
  sunrealtype normq = N_VDotProd(ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->V[0]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstSetNumPreProcess_ArnI(SUNDomEigEstimator DEE,
                                             int numofperprocess)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the number of warmups */
  ArnI_CONTENT(DEE)->numwarmups = numofperprocess;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstPreProcess_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  sunrealtype normq;
  sunindextype i, retval;

  /* Set the initial q = A^{numwarmups}q/||A^{numwarmups}q|| */
  for (i = 0; i < ArnI_CONTENT(DEE)->numwarmups; i++)
  {
    retval = ArnI_CONTENT(DEE)->ATimes(ArnI_CONTENT(DEE)->ATdata,
                                       ArnI_CONTENT(DEE)->V[0],
                                       ArnI_CONTENT(DEE)->q);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DEE_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DEE_ATIMES_FAIL_REC; }
    }
    normq = N_VDotProd(ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->V[0]);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstComputeHess_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->Hes, SUN_ERR_DEE_NULL_HES);

  sunindextype retval, i, j;
  /* Initialize the Hessenberg matrix Hes with zeros */
  for (i = 0; i < ArnI_CONTENT(DEE)->krydim; i++)
  {
    for (j = 0; j < ArnI_CONTENT(DEE)->krydim; j++)
    {
      ArnI_CONTENT(DEE)->Hes[i][j] = ZERO;
    }
  }

  for (i = 0; i < ArnI_CONTENT(DEE)->krydim; i++)
  {
    /* Compute the next Krylov vector */
    retval = ArnI_CONTENT(DEE)->ATimes(ArnI_CONTENT(DEE)->ATdata,
                                       ArnI_CONTENT(DEE)->V[i],
                                       ArnI_CONTENT(DEE)->V[i + 1]);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DEE_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DEE_ATIMES_FAIL_REC; }
    }

    SUNCheckCall(SUNModifiedGS(ArnI_CONTENT(DEE)->V, ArnI_CONTENT(DEE)->Hes,
                               i + 1, ArnI_CONTENT(DEE)->krydim,
                               &(ArnI_CONTENT(DEE)->Hes[i + 1][i])));

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0) / ArnI_CONTENT(DEE)->Hes[i + 1][i],
             ArnI_CONTENT(DEE)->V[i + 1], ArnI_CONTENT(DEE)->V[i + 1]);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimate_ArnI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                  sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_DEE_NULL_ATIMES);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->Hes, SUN_ERR_DEE_NULL_HES);

  sunindextype n = ArnI_CONTENT(DEE)->krydim;

  /* Reshape the Hessenberg matrix as an input vector for the LAPACK dgeev_ function */
  sunindextype i, j, k = 0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      ArnI_CONTENT(DEE)->LAPACK_A[k] = ArnI_CONTENT(DEE)->Hes[i][j];
      k++;
    }
  }

  char jobvl = 'N'; // Do not compute left eigenvectors
  char jobvr = 'N'; // Do not compute right eigenvectors

  /* Call LAPACK's dgeev function
      return info values refer to
    = 0:  successful exit
    < 0:  if info = -i, the i-th argument had an illegal value.
    > 0:  if info = i, the QR algorithm failed to compute all the
          eigenvalues, and no eigenvectors have been computed;
          elements i+1:N of LAPACK_wr and LAPACK_wi contain
          eigenvalues which have converged.
  */
  sunindextype lda  = n;
  sunindextype ldvl = n;
  sunindextype ldvr = n;
  sunindextype info;
  sunindextype lwork = 4 * n;
  xgeev_f77(&jobvl, &jobvr, &n, ArnI_CONTENT(DEE)->LAPACK_A, &lda,
            ArnI_CONTENT(DEE)->LAPACK_wr, ArnI_CONTENT(DEE)->LAPACK_wi, NULL,
            &ldvl, NULL, &ldvr, ArnI_CONTENT(DEE)->LAPACK_work, &lwork, &info);

  if (info != 0)
  {
    // printf(DEE_LAPACK_FAIL, info); // Need to use the error handler

    return SUN_ERR_DEE_LAPACK_FAIL;
  }

  /* order the eigenvalues by their magnitude */
  for (i = 0; i < n; i++)
  {
    ArnI_CONTENT(DEE)->LAPACK_arr[i][0] = ArnI_CONTENT(DEE)->LAPACK_wr[i];
    ArnI_CONTENT(DEE)->LAPACK_arr[i][1] = ArnI_CONTENT(DEE)->LAPACK_wi[i];
  }

  /* Sort the array using qsort */
  qsort(ArnI_CONTENT(DEE)->LAPACK_arr, n,
        sizeof(ArnI_CONTENT(DEE)->LAPACK_arr[0]), sundomeigest_Compare);

  // alternatively we can return a vector of all computed dom_eigs (up to krydim)
  // TODO: Get opinions

  /* Copy the dominant eigenvalue */
  *lambdaR = ArnI_CONTENT(DEE)->LAPACK_arr[0][0];
  *lambdaI = ArnI_CONTENT(DEE)->LAPACK_arr[0][1];

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstFree_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (DEE == NULL) { return SUN_SUCCESS; }

  if (DEE->content)
  {
    /* delete items from within the content structure */
    if (ArnI_CONTENT(DEE)->q)
    {
      N_VDestroy(ArnI_CONTENT(DEE)->q);
      ArnI_CONTENT(DEE)->q = NULL;
    }
    if (ArnI_CONTENT(DEE)->V)
    {
      N_VDestroyVectorArray(ArnI_CONTENT(DEE)->V, ArnI_CONTENT(DEE)->krydim + 1);
      ArnI_CONTENT(DEE)->V = NULL;
    }
    if (ArnI_CONTENT(DEE)->LAPACK_A != NULL)
    {
      free(ArnI_CONTENT(DEE)->LAPACK_A);
      ArnI_CONTENT(DEE)->LAPACK_A = NULL;
    }
    if (ArnI_CONTENT(DEE)->LAPACK_wr != NULL)
    {
      free(ArnI_CONTENT(DEE)->LAPACK_wr);
      ArnI_CONTENT(DEE)->LAPACK_wr = NULL;
    }
    if (ArnI_CONTENT(DEE)->LAPACK_wi != NULL)
    {
      free(ArnI_CONTENT(DEE)->LAPACK_wi);
      ArnI_CONTENT(DEE)->LAPACK_wi = NULL;
    }
    if (ArnI_CONTENT(DEE)->LAPACK_work != NULL)
    {
      free(ArnI_CONTENT(DEE)->LAPACK_work);
      ArnI_CONTENT(DEE)->LAPACK_work = NULL;
    }

    int k;
    /* free LAPACK_arr */
    if (ArnI_CONTENT(DEE)->LAPACK_arr != NULL)
    {
      for (k = 0; k < ArnI_CONTENT(DEE)->krydim; k++)
      {
        free(ArnI_CONTENT(DEE)->LAPACK_arr[k]);
      }
      free(ArnI_CONTENT(DEE)->LAPACK_arr);
      ArnI_CONTENT(DEE)->LAPACK_arr = NULL;
    }
    /* free Hes */
    if (ArnI_CONTENT(DEE)->Hes != NULL)
    {
      for (k = 0; k <= ArnI_CONTENT(DEE)->krydim; k++)
      {
        free(ArnI_CONTENT(DEE)->Hes[k]);
      }
      free(ArnI_CONTENT(DEE)->Hes);
      ArnI_CONTENT(DEE)->Hes = NULL;
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

// Comparison function for qsort
int sundomeigest_Compare(const void* a, const void* b)
{
  const sunrealtype* c1 = *(const sunrealtype* const*)a;
  const sunrealtype* c2 = *(const sunrealtype* const*)b;

  sunrealtype mag1 = SUNRsqrt(c1[0] * c1[0] + c1[1] * c1[1]);
  sunrealtype mag2 = SUNRsqrt(c2[0] * c2[0] + c2[1] * c2[1]);
  return (mag2 > mag1) - (mag2 < mag1); // Descending order
}
