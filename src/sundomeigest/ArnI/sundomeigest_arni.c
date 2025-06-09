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
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sundomeigest/sundomeigest_arni.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

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

SUNDomEigEstimator SUNDomEigEst_ArnI(N_Vector q, sunindextype maxl,
                                     SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_ArnI content;

  /* Check if maxl >= 2 */
  SUNAssert(maxl >= 2, SUN_ERR_DOMEIG_NOT_ENOUGH_ITER);

  /* check for legal q; if illegal return NULL */
  // TO DO: check required vector operations
  SUNAssert(!((q->ops->nvclone == NULL) || (q->ops->nvdestroy == NULL) ||
              (q->ops->nvdotprod == NULL) || (q->ops->nvscale == NULL) ||
              (q->ops->nvgetlength == NULL) || (q->ops->nvspace == NULL)),
            SUN_ERR_DOMEIG_BAD_NVECTOR);

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->gettype            = SUNDomEigEst_ArnIGetType;
  DEE->ops->setatimes          = SUNDomEigEstSetATimes_ArnI;
  DEE->ops->setmaxpoweriter    = NULL;
  DEE->ops->setnumofperprocess = SUNDomEigEstSetNumofPreProcess_ArnI;
  DEE->ops->initialize         = SUNDomEigEstInitialize_ArnI;
  DEE->ops->preprocess         = SUNDomEigEstPreProcess_ArnI;
  DEE->ops->computehess        = SUNDomEigEstComputeHess_ArnI;
  DEE->ops->estimate           = SUNDomEigEstimate_ArnI;
  DEE->ops->getnumofiters      = NULL;
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
  content->maxl        = maxl;
  content->power_of_A  = 0;
  content->LAPACK_A    = NULL;
  content->LAPACK_wr   = NULL;
  content->LAPACK_wi   = NULL;
  content->LAPACK_work = NULL;
  content->LAPACK_arr  = NULL;
  content->Hes         = NULL;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  content->V = N_VCloneVectorArray(maxl + 1, q);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNDomEigEstimator_Type SUNDomEigEst_ArnIGetType(
  SUNDIALS_MAYBE_UNUSED SUNDomEigEstimator DEE)
{
  return (SUNDOMEIG_ARNOLDI);
}

SUNErrCode SUNDomEigEstInitialize_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  if (ArnI_CONTENT(DEE)->maxl < 2)
  {
    ArnI_CONTENT(DEE)->maxl = SUNDOMEIGEST_ARN_MAXL_DEFAULT;
  }
  if (ArnI_CONTENT(DEE)->power_of_A < 0)
  {
    ArnI_CONTENT(DEE)->power_of_A = SUNDOMEIGEST_PI_POWER_OF_A_DEFAULT;
  }

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  ArnI_CONTENT(DEE)->LAPACK_A = (sunrealtype*)malloc(
    (ArnI_CONTENT(DEE)->maxl * ArnI_CONTENT(DEE)->maxl) * sizeof(sunrealtype));
  ArnI_CONTENT(DEE)->LAPACK_wr =
    malloc(ArnI_CONTENT(DEE)->maxl * sizeof(sunrealtype));
  ArnI_CONTENT(DEE)->LAPACK_wi =
    malloc(ArnI_CONTENT(DEE)->maxl * sizeof(sunrealtype));
  ArnI_CONTENT(DEE)->LAPACK_work =
    malloc((4 * ArnI_CONTENT(DEE)->maxl) * sizeof(sunrealtype));
  ArnI_CONTENT(DEE)->LAPACK_arr =
    (suncomplextype*)malloc(ArnI_CONTENT(DEE)->maxl * sizeof(suncomplextype));

  N_VRandom(ArnI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  /* Hessenberg matrix Hes */
  if (ArnI_CONTENT(DEE)->Hes == NULL)
  {
    sunindextype k;
    ArnI_CONTENT(DEE)->Hes = (sunrealtype**)malloc(
      (ArnI_CONTENT(DEE)->maxl + 1) * sizeof(sunrealtype*));

    for (k = 0; k <= ArnI_CONTENT(DEE)->maxl; k++)
    {
      ArnI_CONTENT(DEE)->Hes[k] = NULL;
      ArnI_CONTENT(DEE)->Hes[k] =
        (sunrealtype*)malloc(ArnI_CONTENT(DEE)->maxl * sizeof(sunrealtype));
    }
  }

  /* Initialize the vector V[0] */
  sunrealtype normq = N_VDotProd(ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, ArnI_CONTENT(DEE)->q, ArnI_CONTENT(DEE)->V[0]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
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

SUNErrCode SUNDomEigEstSetNumofPreProcess_ArnI(SUNDomEigEstimator DEE,
                                               sunindextype numofperprocess)
{
  SUNFunctionBegin(DEE->sunctx);

  /* set the number of warmups */
  ArnI_CONTENT(DEE)->power_of_A = numofperprocess;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstPreProcess_ArnI(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  sunrealtype normq;
  sunindextype i, retval;

  /* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
  for (i = 0; i < ArnI_CONTENT(DEE)->power_of_A; i++)
  {
    retval = ArnI_CONTENT(DEE)->ATimes(ArnI_CONTENT(DEE)->ATdata,
                                       ArnI_CONTENT(DEE)->V[0],
                                       ArnI_CONTENT(DEE)->q);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DOMEIG_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DOMEIG_ATIMES_FAIL_REC; }
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

  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->Hes, SUN_ERR_ARG_CORRUPT);

  sunindextype retval, i, j;
  /* Initialize the Hessenberg matrix Hes with zeros */
  for (i = 0; i < ArnI_CONTENT(DEE)->maxl; i++)
  {
    for (j = 0; j < ArnI_CONTENT(DEE)->maxl; j++)
    {
      ArnI_CONTENT(DEE)->Hes[i][j] = ZERO;
    }
  }

  for (i = 0; i < ArnI_CONTENT(DEE)->maxl; i++)
  {
    /* Compute the next Krylov vector */
    retval = ArnI_CONTENT(DEE)->ATimes(ArnI_CONTENT(DEE)->ATdata,
                                       ArnI_CONTENT(DEE)->V[i],
                                       ArnI_CONTENT(DEE)->V[i + 1]);
    if (retval != 0)
    {
      if (retval < 0) { return SUN_ERR_DOMEIG_ATIMES_FAIL_UNREC; }
      else { return SUN_ERR_DOMEIG_ATIMES_FAIL_REC; }
    }

    SUNCheckCall(SUNModifiedGS(ArnI_CONTENT(DEE)->V, ArnI_CONTENT(DEE)->Hes,
                               i + 1, ArnI_CONTENT(DEE)->maxl,
                               &(ArnI_CONTENT(DEE)->Hes[i + 1][i])));

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0) / ArnI_CONTENT(DEE)->Hes[i + 1][i],
             ArnI_CONTENT(DEE)->V[i + 1], ArnI_CONTENT(DEE)->V[i + 1]);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimate_ArnI(SUNDomEigEstimator DEE, suncomplextype* dom_eig)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(dom_eig, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(ArnI_CONTENT(DEE)->Hes, SUN_ERR_ARG_CORRUPT);

  suncomplextype dom_eig_new;
  dom_eig_new.real = ZERO;
  dom_eig_new.imag = ZERO;

  int n = ArnI_CONTENT(DEE)->maxl;

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

  int lda = n, ldvl = n, ldvr = n;
  int info, lwork = 4 * n;

  char jobvl = 'N'; // Do not compute left eigenvectors
  char jobvr = 'N'; // Do not compute right eigenvectors

  /* Call LAPACK's dgeev function */
  dgeev_(&jobvl, &jobvr, &n, ArnI_CONTENT(DEE)->LAPACK_A, &lda,
         ArnI_CONTENT(DEE)->LAPACK_wr, ArnI_CONTENT(DEE)->LAPACK_wi, NULL,
         &ldvl, NULL, &ldvr, ArnI_CONTENT(DEE)->LAPACK_work, &lwork, &info);

  if (info != 0)
  {
    printf(SUNDOMEIGEST_LAPACK_FAIL, info);

    return SUN_ERR_DOMEIG_LAPACK_FAIL;
  }

  /* order the eigenvalues by their magnitude */
  for (i = 0; i < n; i++)
  {
    ArnI_CONTENT(DEE)->LAPACK_arr[i].real = ArnI_CONTENT(DEE)->LAPACK_wr[i];
    ArnI_CONTENT(DEE)->LAPACK_arr[i].imag = ArnI_CONTENT(DEE)->LAPACK_wi[i];
  }

  /* Sort the array using qsort */
  qsort(ArnI_CONTENT(DEE)->LAPACK_arr, n, sizeof(suncomplextype), domeig_Compare);

  /* Update the original arrays */
  for (i = 0; i < n; i++)
  {
    ArnI_CONTENT(DEE)->LAPACK_wr[i] = ArnI_CONTENT(DEE)->LAPACK_arr[i].real;
    ArnI_CONTENT(DEE)->LAPACK_wi[i] = ArnI_CONTENT(DEE)->LAPACK_arr[i].imag;
  }

  // alternatively we can return a vector of all computed dom_eigs (up to maxl)
  // TODO: Get opinions

  /* Copy the dominant eigenvalue */
  dom_eig_new.real = ArnI_CONTENT(DEE)->LAPACK_wr[0];
  dom_eig_new.imag = ArnI_CONTENT(DEE)->LAPACK_wi[0];

  *dom_eig = dom_eig_new;

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
      N_VDestroyVectorArray(ArnI_CONTENT(DEE)->V, 1);
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
    if (ArnI_CONTENT(DEE)->LAPACK_arr != NULL)
    {
      free(ArnI_CONTENT(DEE)->LAPACK_arr);
      ArnI_CONTENT(DEE)->LAPACK_arr = NULL;
    }
    if (ArnI_CONTENT(DEE)->Hes != NULL)
    {
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

// Function to calculate the magnitude of a suncomplextype number
sunrealtype domeig_Magnitude(const suncomplextype* c)
{
  return sqrt(c->real * c->real + c->imag * c->imag);
}

// Comparison function for qsort
sunindextype domeig_Compare(const void* a, const void* b)
{
  const suncomplextype* c1 = (const suncomplextype*)a;
  const suncomplextype* c2 = (const suncomplextype*)b;
  sunrealtype mag1         = domeig_Magnitude(c1);
  sunrealtype mag2         = domeig_Magnitude(c2);
  return (mag2 > mag1) - (mag2 < mag1); // Descending order
}