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
 * This is the implementation file for the Arnoldi Iteration
 * implementation of the SUNDomEigEst package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundomeigest/sundomeigest_arnoldi.h>

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

#define ONE SUN_RCONST(1.0)

/* Default estimator parameters */
#define DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT 100

/* Default Arnoldi Iteration parameters */
#define DEE_KRYLOV_DIM_DEFAULT 3

/*
 * -----------------------------------------------------------------
 * Arnoldi itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define Arnoldi_CONTENT(DEE) ((SUNDomEigEstimatorContent_Arnoldi)(DEE->content))

/*
 * -----------------------------------------------------------------
 * Arnoldi module private function prototypes
 * -----------------------------------------------------------------
 */

int sundomeigest_Compare(const void* a, const void* b);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new Arnoldi estimator
 */

SUNDomEigEstimator SUNDomEigEst_Arnoldi(N_Vector q, int kry_dim,
                                        int num_warmups, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_Arnoldi content;

  /* Check if kry_dim >= 2 */
  if (kry_dim < 3) { kry_dim = DEE_KRYLOV_DIM_DEFAULT; }

  /* Check if num_warmups >= 0 */
  if (num_warmups < 0) { kry_dim = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT; }

  /* check for legal q; if illegal return NULL */
  SUNAssertNull(!((q->ops->nvclone == NULL) || (q->ops->nvdestroy == NULL) ||
                  (q->ops->nvdotprod == NULL) || (q->ops->nvscale == NULL)),
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Check if q != 0 vector */
  SUNAssertNull(N_VDotProd(q, q) > SUN_SMALL_REAL, SUN_ERR_ARG_INCOMPATIBLE);

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEst_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes         = SUNDomEigEst_SetATimes_Arnoldi;
  DEE->ops->setnumpreprocess  = SUNDomEigEst_SetNumPreProcess_Arnoldi;
  DEE->ops->initialize        = SUNDomEigEst_Initialize_Arnoldi;
  DEE->ops->estimate          = SUNDomEig_Estimate_Arnoldi;
  DEE->ops->getnumatimescalls = SUNDomEigEst_GetNumATimesCalls_Arnoldi;
  DEE->ops->write             = SUNDomEigEst_Write_Arnoldi;
  DEE->ops->destroy           = SUNDomEigEst_Destroy_Arnoldi;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_Arnoldi)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes      = NULL;
  content->ATdata      = NULL;
  content->V           = NULL;
  content->q           = NULL;
  content->kry_dim     = kry_dim;
  content->num_warmups = num_warmups;
  content->LAPACK_A    = NULL;
  content->LAPACK_wr   = NULL;
  content->LAPACK_wi   = NULL;
  content->LAPACK_work = NULL;
  content->LAPACK_arr  = NULL;
  content->Hes         = NULL;
  content->num_ATimes  = 0;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  N_VScale(ONE, q, content->q);
  SUNCheckLastErrNull();

  content->V = N_VCloneVectorArray(kry_dim + 1, q);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNErrCode SUNDomEigEst_SetATimes_Arnoldi(SUNDomEigEstimator DEE, void* A_data,
                                          SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT); 

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  Arnoldi_CONTENT(DEE)->ATimes = ATimes;
  Arnoldi_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Initialize_Arnoldi(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT); 

  if (Arnoldi_CONTENT(DEE)->kry_dim < 2)
  {
    Arnoldi_CONTENT(DEE)->kry_dim = DEE_KRYLOV_DIM_DEFAULT;
  }
  if (Arnoldi_CONTENT(DEE)->num_warmups < 0)
  {
    Arnoldi_CONTENT(DEE)->num_warmups = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;
  }

  SUNAssert(Arnoldi_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  if (Arnoldi_CONTENT(DEE)->LAPACK_A == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_A = (sunrealtype*)malloc(
      (Arnoldi_CONTENT(DEE)->kry_dim * Arnoldi_CONTENT(DEE)->kry_dim) *
      sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_A, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wr == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wr =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wr, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wi == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wi =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wi, SUN_ERR_MALLOC_FAIL);
  }

  /* query the workspace size (call with lwork = -1) */
  char jobvl         = 'N';
  char jobvr         = 'N';
  sunindextype N     = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype lda   = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype ldvl  = 1;
  sunindextype ldvr  = 1;
  sunindextype info  = 0;
  sunindextype lwork = -1;
  sunrealtype work   = SUN_RCONST(0.0);

  xgeev_f77(&jobvl, &jobvr, &N, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr, Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr, &work, &lwork, &info);

  Arnoldi_CONTENT(DEE)->LAPACK_work =
    (sunrealtype*)malloc(((sunindextype)work) * sizeof(sunrealtype));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_work, SUN_ERR_MALLOC_FAIL);

  /* LAPACK array */
  Arnoldi_CONTENT(DEE)->LAPACK_arr =
    (sunrealtype**)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype*));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr, SUN_ERR_MALLOC_FAIL); 

  for (int k = 0; k < Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_arr[k] =
      (sunrealtype*)malloc(2 * sizeof(sunrealtype));
      SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr[k], SUN_ERR_MALLOC_FAIL);
  }

  /* Hessenberg matrix Hes */
  Arnoldi_CONTENT(DEE)->Hes = (sunrealtype**)malloc(
    (Arnoldi_CONTENT(DEE)->kry_dim + 1) * sizeof(sunrealtype*));
    SUNAssert(Arnoldi_CONTENT(DEE)->Hes, SUN_ERR_MALLOC_FAIL); 

  for (int k = 0; k <= Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->Hes[k] =
      (sunrealtype*)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
      SUNAssert(Arnoldi_CONTENT(DEE)->Hes[k], SUN_ERR_MALLOC_FAIL);
  }

  sunrealtype normq = N_VDotProd(Arnoldi_CONTENT(DEE)->q,
                                 Arnoldi_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->V[0]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_SetNumPreProcess_Arnoldi(SUNDomEigEstimator DEE,
                                                 int numpreprocess)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if numpreprocess >= 0 */
  if (numpreprocess < 0) { numpreprocess = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;}

  /* set the number of warmups */
  Arnoldi_CONTENT(DEE)->num_warmups = numpreprocess;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEig_Estimate_Arnoldi(SUNDomEigEstimator DEE,
                                      sunrealtype* lambdaR, sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->Hes, SUN_ERR_ARG_CORRUPT);  

  int retval;
  sunindextype n = Arnoldi_CONTENT(DEE)->kry_dim;
  sunrealtype normq;

  /* Set the initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
  for (int i = 0; i < Arnoldi_CONTENT(DEE)->num_warmups; i++)
  {
    retval = Arnoldi_CONTENT(DEE)->ATimes(Arnoldi_CONTENT(DEE)->ATdata,
                                          Arnoldi_CONTENT(DEE)->V[0],
                                          Arnoldi_CONTENT(DEE)->q);
    Arnoldi_CONTENT(DEE)->num_ATimes++;
    if (retval != 0)
    {
      return SUN_ERR_USER_FCN_FAIL; 
    }

    normq = N_VDotProd(Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->V[0]);
    SUNCheckLastErr();
  }

  for (int i = 0; i < n; i++)
  {
    /* Compute the next Krylov vector */
    retval = Arnoldi_CONTENT(DEE)->ATimes(Arnoldi_CONTENT(DEE)->ATdata,
                                          Arnoldi_CONTENT(DEE)->V[i],
                                          Arnoldi_CONTENT(DEE)->V[i + 1]);
    Arnoldi_CONTENT(DEE)->num_ATimes++;
    if (retval != 0)
    {
      return SUN_ERR_USER_FCN_FAIL;
    }

    SUNCheckCall(SUNModifiedGS(Arnoldi_CONTENT(DEE)->V,
                               Arnoldi_CONTENT(DEE)->Hes, i + 1, (int)n,
                               &(Arnoldi_CONTENT(DEE)->Hes[i + 1][i])));

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0) / Arnoldi_CONTENT(DEE)->Hes[i + 1][i],
             Arnoldi_CONTENT(DEE)->V[i + 1], Arnoldi_CONTENT(DEE)->V[i + 1]);
    SUNCheckLastErr();
  }

  /* Reshape the Hessenberg matrix as an input vector for the LAPACK dgeev_ function */
  int k = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      Arnoldi_CONTENT(DEE)->LAPACK_A[k] = Arnoldi_CONTENT(DEE)->Hes[i][j];
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
  xgeev_f77(&jobvl, &jobvr, &n, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr, Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr, Arnoldi_CONTENT(DEE)->LAPACK_work, &lwork,
            &info);

  if (info != 0)
  {
    return SUN_ERR_EXT_FAIL;
  }

  /* order the eigenvalues by their magnitude */
  for (int i = 0; i < n; i++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_arr[i][0] = Arnoldi_CONTENT(DEE)->LAPACK_wr[i];
    Arnoldi_CONTENT(DEE)->LAPACK_arr[i][1] = Arnoldi_CONTENT(DEE)->LAPACK_wi[i];
  }

  /* Sort the array using qsort */
  qsort(Arnoldi_CONTENT(DEE)->LAPACK_arr, n,
        sizeof(Arnoldi_CONTENT(DEE)->LAPACK_arr[0]), sundomeigest_Compare);

  /* Substitute the ordered eigenvalues back in LAPACK_w* */
  for (int i = 0; i < n; i++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wr[i] = Arnoldi_CONTENT(DEE)->LAPACK_arr[i][0];
    Arnoldi_CONTENT(DEE)->LAPACK_wi[i] = Arnoldi_CONTENT(DEE)->LAPACK_arr[i][1];
  }

  /* Copy the dominant eigenvalue */
  *lambdaR = Arnoldi_CONTENT(DEE)->LAPACK_wr[0];
  *lambdaI = Arnoldi_CONTENT(DEE)->LAPACK_wi[0];

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_GetNumATimesCalls_Arnoldi(SUNDomEigEstimator DEE,
                                                  long int* num_ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);  
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT); 
  SUNAssert(num_ATimes, SUN_ERR_ARG_CORRUPT);

  *num_ATimes = Arnoldi_CONTENT(DEE)->num_ATimes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Write_Arnoldi(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(outfile, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (DEE == NULL || outfile == NULL) { return SUN_ERR_ARG_CORRUPT; }

  fprintf(outfile, "\nArnoldi Iteration DEE Statistics:");
  fprintf(outfile, "\n------------------------------------------------\n");
  fprintf(outfile, "Krylov dimensions             = %d\n",
          Arnoldi_CONTENT(DEE)->kry_dim);
  fprintf(outfile, "Num. of warmups               = %d\n",
          Arnoldi_CONTENT(DEE)->num_warmups);
  fprintf(outfile, "Num. of ATimes calls          = %ld\n\n",
          Arnoldi_CONTENT(DEE)->num_ATimes);

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEst_Destroy_Arnoldi(SUNDomEigEstimator* DEEptr)
{
  SUNFunctionBegin((*DEEptr)->sunctx);

  if ((*DEEptr) == NULL) { return SUN_SUCCESS; }

  SUNDomEigEstimator DEE = *DEEptr;

  if (DEE->content)
  {
    /* delete items from within the content structure */
    if (Arnoldi_CONTENT(DEE)->q)
    {
      N_VDestroy(Arnoldi_CONTENT(DEE)->q);
      Arnoldi_CONTENT(DEE)->q = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->V)
    {
      N_VDestroyVectorArray(Arnoldi_CONTENT(DEE)->V,
                            Arnoldi_CONTENT(DEE)->kry_dim + 1);
      Arnoldi_CONTENT(DEE)->V = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_A != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_A);
      Arnoldi_CONTENT(DEE)->LAPACK_A = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_wr != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_wr);
      Arnoldi_CONTENT(DEE)->LAPACK_wr = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_wi != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_wi);
      Arnoldi_CONTENT(DEE)->LAPACK_wi = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_work != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_work);
      Arnoldi_CONTENT(DEE)->LAPACK_work = NULL;
    }

    /* free LAPACK_arr */
    if (Arnoldi_CONTENT(DEE)->LAPACK_arr != NULL)
    {
      for (int k = 0; k < Arnoldi_CONTENT(DEE)->kry_dim; k++)
      {
        free(Arnoldi_CONTENT(DEE)->LAPACK_arr[k]);
      }
      free(Arnoldi_CONTENT(DEE)->LAPACK_arr);
      Arnoldi_CONTENT(DEE)->LAPACK_arr = NULL;
    }
    /* free Hes */
    if (Arnoldi_CONTENT(DEE)->Hes != NULL)
    {
      for (int k = 0; k <= Arnoldi_CONTENT(DEE)->kry_dim; k++)
      {
        free(Arnoldi_CONTENT(DEE)->Hes[k]);
      }
      free(Arnoldi_CONTENT(DEE)->Hes);
      Arnoldi_CONTENT(DEE)->Hes = NULL;
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

// Comparison function for qsort
int sundomeigest_Compare(const void* a, const void* b)
{
  const sunrealtype* c1 = *(const sunrealtype* const*)a;
  const sunrealtype* c2 = *(const sunrealtype* const*)b;

  sunrealtype mag1 = SUNRsqrt(c1[0] * c1[0] + c1[1] * c1[1]);
  sunrealtype mag2 = SUNRsqrt(c2[0] * c2[0] + c2[1] * c2[1]);
  return (mag2 > mag1) - (mag2 < mag1); // Descending order
}
