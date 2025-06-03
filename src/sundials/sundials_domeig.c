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
 * This is the implementation file for the eigenvalue estimation of
 * the SUNDOMEIG package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_domeig_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>

#define ZERO SUN_RCONST(0.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

SUNErrCode DomEigCreate(SUNATimesFn ATimes, void* Adata,
                N_Vector q, int maxl, SUNContext sunctx, void** domeig_mem_out)
{
  SUNFunctionBegin(sunctx); // is this correct?
  DOMEIGMem domeig_mem;

  /* Test if Atimes and q are provided */
  if (ATimes == NULL || q == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_ATIMES;
  }

  /* Check if maxl >= 2 */
  if (maxl < 2)
  {
    return SUN_ERR_DOMEIG_NOT_ENOUGH_ITER;
  }

  /* Test if sunctx is provided */
  if (sunctx == NULL)
  {
    return SUN_ERR_SUNCTX_CORRUPT;
  }

  /* Test if all required vector operations are implemented */
  SUNCheckCall(domeig_CheckNVector(q));

  /* Allocate DOMEIGMem structure, and initialize to zero */
  domeig_mem = (DOMEIGMem)calloc(1, sizeof(*domeig_mem));
  SUNAssert(domeig_mem, SUN_ERR_MALLOC_FAIL);

  /* Copy the inputs into DOMEIG memory */
  domeig_mem->sunctx      = q->sunctx;
  domeig_mem->ATimes      = ATimes;
  domeig_mem->Adata       = Adata;
  domeig_mem->q           = q;
  domeig_mem->maxl        = maxl;
  domeig_mem->length      = q->ops->nvgetlength(q);

  /* Set the default power of A to start with (# of warm-ups) */
  domeig_mem->power_of_A  = DEFAULT_POWER_OF_A;

  /* Set the default tolerance of the power iteration */
  domeig_mem->powiter_tol = DEFAULT_POWER_ITER_TOL;

  /* Set the default max number of the power iteration */
  domeig_mem->max_powiter = DEFAULT_MAX_POWER_ITER;

  if (domeig_mem->length > 2)
  {
    domeig_mem->LAPACK_A = (sunrealtype*)malloc((maxl*maxl) * sizeof(sunrealtype));
    domeig_mem->LAPACK_wr = malloc(maxl * sizeof(sunrealtype));
    domeig_mem->LAPACK_wi = malloc(maxl * sizeof(sunrealtype));
    domeig_mem->LAPACK_work = malloc((4 * maxl) * sizeof(sunrealtype));
    domeig_mem->LAPACK_arr = (suncomplextype *)malloc(maxl * sizeof(suncomplextype));
  }

  /* Hessenberg matrix Hes */
  if (domeig_mem->Hes == NULL && domeig_mem->length > 2)
  {
    int k;
    domeig_mem->Hes =
      (sunrealtype**)malloc((maxl + 1) * sizeof(sunrealtype*));

    for (k = 0; k <= maxl; k++)
    {
      domeig_mem->Hes[k] = NULL;
      domeig_mem->Hes[k] = (sunrealtype*)malloc(maxl * sizeof(sunrealtype));
    }
  }

  /* Krylov subspace vectors */
  if (domeig_mem->V == NULL)
  {
    if(domeig_mem->length > 2)
    {
      domeig_mem->V = N_VCloneVectorArray(maxl + 1, q);
    }
    else
    {
      domeig_mem->V = N_VCloneVectorArray(1, q);
    }
  }

  /* Unitize the initial vector V[0] */
  sunrealtype normq = N_VDotProd(q, q);
  normq = SUNRsqrt(normq);
  N_VScale(SUN_RCONST(1.0)/normq, domeig_mem->q, domeig_mem->V[0]);

  *domeig_mem_out = (void*)domeig_mem;

  return SUN_SUCCESS;
}

/* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
SUNErrCode DomEigPreProcess(DOMEIGMem domeig_mem, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx); // is this correct?

  int retval;

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    return SUN_ERR_DOMEIG_MEM_FAIL;
  }

  /* Check if ATimes is provided */
  if(domeig_mem->ATimes == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_ATIMES;
  }

  sunrealtype normq;
  int i;

  /* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
  for(i = 0; i < domeig_mem->power_of_A; i++) {
    retval = domeig_mem->ATimes(domeig_mem->Adata, domeig_mem->V[0], domeig_mem->q);
    if (retval != 0)
    {
      DomEigDestroy(&domeig_mem);

      if(retval < 0)
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_UNREC;
      }
      else
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_REC;
      }
    }
    normq = N_VDotProd(domeig_mem->q, domeig_mem->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(SUN_RCONST(1.0)/normq, domeig_mem->q, domeig_mem->V[0]);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

/* Compute the Hessenberg matrix DomEig_mem->Hes*/
SUNErrCode DomEigComputeHess(DOMEIGMem domeig_mem, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx); // is this correct?

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    return SUN_ERR_DOMEIG_MEM_FAIL;
  }

  /* Check if the dim of the matrix is less than 3.
     Return immediately if dim <= 2;
     no need to compute Hessenberg matrix
     since the default is power iteration for dim <= 2 */
  if(domeig_mem->length < 3)
  {
    return SUN_SUCCESS;
  }

  /* Check if ATimes is provided */
  if(domeig_mem->ATimes == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_ATIMES;
  }

  /* Check if Hes is allocated */
  if(domeig_mem->Hes == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_HES;
  }

  int retval, i, j;
  /* Initialize the Hessenberg matrix Hes with zeros */
  for (i = 0; i < domeig_mem->maxl; i++)
  {
    for (j = 0; j < domeig_mem->maxl; j++) { domeig_mem->Hes[i][j] = ZERO; }
  }

  for (i = 0; i < domeig_mem->maxl; i++)
  {
    /* Compute the next Krylov vector */
    retval = domeig_mem->ATimes(domeig_mem->Adata, domeig_mem->V[i], domeig_mem->V[i+1]);
    if (retval != 0)
    {
      DomEigDestroy(&domeig_mem);

      if(retval < 0)
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_UNREC;
      }
      else
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_REC;
      }
    }

    SUNCheckCall(SUNModifiedGS(domeig_mem->V, domeig_mem->Hes, i + 1, domeig_mem->maxl, &(domeig_mem->Hes[i + 1][i])));

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0)/domeig_mem->Hes[i + 1][i], domeig_mem->V[i+1], domeig_mem->V[i+1]);
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

SUNErrCode DomEigPowerIteration(DOMEIGMem domeig_mem, suncomplextype* dom_eig, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx); // is this correct?

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    return SUN_ERR_DOMEIG_MEM_FAIL;
  }

  /* Check if ATimes is provided */
  if(domeig_mem->ATimes == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_ATIMES;
  }

  suncomplextype dom_eig_new = *dom_eig;
  suncomplextype dom_eig_old;

  dom_eig_new.real = ZERO;
  dom_eig_new.imag = ZERO;
  dom_eig_old.real = ZERO;
  dom_eig_old.imag = ZERO;

  int retval, k;
  sunrealtype normq;

  for (k = 0; k < domeig_mem->max_powiter; k++)
  {
    retval = domeig_mem->ATimes(domeig_mem->Adata, domeig_mem->V[0], domeig_mem->q);
    if (retval != 0)
    {
      DomEigDestroy(&domeig_mem);

      if(retval < 0)
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_UNREC;
      }
      else
      {
        return SUN_ERR_DOMEIG_ATIMES_FAIL_REC;
      }
    }

    dom_eig_new.real = N_VDotProd(domeig_mem->V[0], domeig_mem->q); //Rayleigh quotient
    SUNCheckLastErr();

    if(fabs(dom_eig_new.real - dom_eig_old.real) < domeig_mem->powiter_tol)
    {
      break;
    }

    normq = N_VDotProd(domeig_mem->q, domeig_mem->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(SUN_RCONST(1.0)/normq, domeig_mem->q, domeig_mem->V[0]);
    SUNCheckLastErr();

    dom_eig_old.real = dom_eig_new.real;
  }

  *dom_eig = dom_eig_new;

  return SUN_SUCCESS;
}

/* Estimate the dominant eigvalues of the Hessenberg matrix or
   run power iterations for problem size less than 3 */
SUNErrCode DomEigEstimate(DOMEIGMem domeig_mem, suncomplextype* dom_eig, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx); // is this correct?

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    return SUN_ERR_DOMEIG_MEM_FAIL;
  }

  suncomplextype dom_eig_new;
  dom_eig_new.real = ZERO;
  dom_eig_new.imag = ZERO;

  /* run power iterations for problem size less than 3 */
  if(domeig_mem->length < 3)
  {
    SUNCheckCall(DomEigPowerIteration(domeig_mem, &dom_eig_new, domeig_mem->sunctx));

    *dom_eig = dom_eig_new;

    return SUN_SUCCESS;
  }

  /* Check if Hes is allocated */
  if(domeig_mem->Hes == NULL)
  {
    return SUN_ERR_DOMEIG_NULL_HES;
  }

  int n = domeig_mem->maxl;

  /* Reshape the Hessenberg matrix as an input vector for the LAPACK dgeev_ function */
  int i, j, k = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
        domeig_mem->LAPACK_A[k] = domeig_mem->Hes[i][j];
        k++;
    }
  }

  int lda = n, ldvl = n, ldvr = n;
  int info, lwork = 4 * n;

  char jobvl = 'N'; // Do not compute left eigenvectors
  char jobvr = 'N'; // Do not compute right eigenvectors

  /* Call LAPACK's dgeev function */
  dgeev_(&jobvl, &jobvr, &n, domeig_mem->LAPACK_A, &lda, domeig_mem->LAPACK_wr, domeig_mem->LAPACK_wi, NULL, &ldvl, NULL, &ldvr, domeig_mem->LAPACK_work, &lwork, &info);

  if (info != 0) {
      printf(MSG_DOMEIG_LAPACK_FAIL, info);
      DomEigDestroy(&domeig_mem);

      return SUN_ERR_DOMEIG_LAPACK_FAIL;
  }

  /* order the eigenvalues by their magnitude */
  for (i = 0; i < n; i++) {
      domeig_mem->LAPACK_arr[i].real = domeig_mem->LAPACK_wr[i];
      domeig_mem->LAPACK_arr[i].imag = domeig_mem->LAPACK_wi[i];
  }

  /* Sort the array using qsort */
  qsort(domeig_mem->LAPACK_arr, n, sizeof(suncomplextype), domeig_Compare);

  /* Update the original arrays */
  for (i = 0; i < n; i++) {
      domeig_mem->LAPACK_wr[i] = domeig_mem->LAPACK_arr[i].real;
      domeig_mem->LAPACK_wi[i] = domeig_mem->LAPACK_arr[i].imag;
  }

  // alternatively we can return a vector of all computed dom_eigs (up to maxl)
  // TODO: Get opinions

  /* Copy the dominant eigenvalue */
  dom_eig_new.real = domeig_mem->LAPACK_wr[0];
  dom_eig_new.imag = domeig_mem->LAPACK_wi[0];

  *dom_eig = dom_eig_new;

  return SUN_SUCCESS;
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  domeig_CheckNVector:

  This routine checks if all required vector operations are
  present. If any of them is missing it returns the corresponding
  SUNErrCode.
  ---------------------------------------------------------------*/
SUNErrCode domeig_CheckNVector(N_Vector tmpl)
{ // TO DO: check required vector operations
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvdotprod == NULL) || (tmpl->ops->nvscale == NULL) ||
      (tmpl->ops->nvgetlength == NULL) || (tmpl->ops->nvspace == NULL))
  {
    return SUN_ERR_DOMEIG_BAD_NVECTOR;
  }
  return SUN_SUCCESS;
}

// Function to calculate the magnitude of a suncomplextype number
sunrealtype domeig_Magnitude(const suncomplextype *c) {
    return sqrt(c->real * c->real + c->imag * c->imag);
}

// Comparison function for qsort
int domeig_Compare(const void *a, const void *b) {
    const suncomplextype *c1 = (const suncomplextype *)a;
    const suncomplextype *c2 = (const suncomplextype *)b;
    sunrealtype mag1 = domeig_Magnitude(c1);
    sunrealtype mag2 = domeig_Magnitude(c2);
    return (mag2 > mag1) - (mag2 < mag1); // Descending order
}

/*---------------------------------------------------------------
  DomEigDestroy frees all DomEig memory.
  ---------------------------------------------------------------*/
void DomEigDestroy(void** domeig_mem)
{
  DOMEIGMem dom_eig_mem;

  /* nothing to do if domeig_mem is already NULL */
  if (*domeig_mem == NULL) { return; }

  dom_eig_mem = (DOMEIGMem)(*domeig_mem);

  if (dom_eig_mem->q != NULL)
  {
    N_VDestroy(dom_eig_mem->q);
    dom_eig_mem->q = NULL;
  }
  if (dom_eig_mem->V != NULL)
  {
    N_VDestroyVectorArray(dom_eig_mem->V, dom_eig_mem->maxl + 1);
    dom_eig_mem->V = NULL;
  }

  if (dom_eig_mem->LAPACK_A != NULL)
  {
    free(dom_eig_mem->LAPACK_A);
    dom_eig_mem->LAPACK_A = NULL;
  }

  if (dom_eig_mem->LAPACK_wr != NULL)
  {
    free(dom_eig_mem->LAPACK_wr);
    dom_eig_mem->LAPACK_wr = NULL;
  }

  if (dom_eig_mem->LAPACK_wi != NULL)
  {
    free(dom_eig_mem->LAPACK_wi);
    dom_eig_mem->LAPACK_wi = NULL;
  }

  if (dom_eig_mem->LAPACK_arr != NULL)
  {
    free(dom_eig_mem->LAPACK_arr);
    dom_eig_mem->LAPACK_arr = NULL;
  }

  if (dom_eig_mem->Hes != NULL)
  {
    free(dom_eig_mem->Hes);
    dom_eig_mem->Hes = NULL;
  }

  free(*domeig_mem);
  *domeig_mem = NULL;
}