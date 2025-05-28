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

#define ZERO SUN_RCONST(0.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

void* DomEigCreate(SUNATimesFn ATimes, void* Adata,
                N_Vector q, int maxl, SUNContext sunctx)
{
  DOMEIGMem domeig_mem;

  /* Test if Atimes is provided */
  if (ATimes == NULL)
  {
    printf(MSG_DOMEIG_NULL_ATIMES);

    return NULL;
  }

  /* Test if q is provided */
  if (q == NULL)
  {
    printf(MSG_DOMEIG_NULL_q);

    return NULL;
  }

  /* Check if maxl >= 2 */
  if (maxl < 2)
  {
    printf(MSG_DOMEIG_NOT_ENOUGH_ITER);

    return NULL;
  }

  /* Test if sunctx is provided */
  if (sunctx == NULL)
  {
    printf(MSG_DOMEIG_NULL_SUNCTX);

    return NULL;
  }

  /* Test if all required vector operations are implemented */
  if (!domeig_CheckNVector(q))
  {
    printf(MSG_DOMEIG_BAD_NVECTOR);

    return NULL;
  }

  /* Allocate DOMEIGMem structure, and initialize to zero */
  domeig_mem = (DOMEIGMem)calloc(1, sizeof(*domeig_mem));

  /* Copy the inputs into DOMEIG memory */
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

  return (void*)domeig_mem;
}

/* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
int DomEigPreProcess(DOMEIGMem domeig_mem)
{
  int retval;

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    printf(MSG_DOMEIG_MEM_FAIL);

    return -1;
  }

  /* Check if ATimes is provided */
  if(domeig_mem->ATimes == NULL)
  {
    printf(MSG_DOMEIG_NULL_ATIMES);
    DomEigFree(&domeig_mem);

    return -1;
  }

  sunrealtype normq;
  int i;

  /* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
  for(i = 0; i < domeig_mem->power_of_A; i++) {
    retval = domeig_mem->ATimes(domeig_mem->Adata, domeig_mem->V[0], domeig_mem->q);
    if (retval != 0)
    {
      (retval < 0) ?
      printf(MSG_DOMEIG_ATIMES_FAIL_UNREC) :
      printf(MSG_DOMEIG_ATIMES_FAIL_REC);
      DomEigFree(&domeig_mem);

      return retval;
    }
    normq = N_VDotProd(domeig_mem->q, domeig_mem->q);
    normq = SUNRsqrt(normq);
    N_VScale(SUN_RCONST(1.0)/normq, domeig_mem->q, domeig_mem->V[0]);
  }

  return 0;
}

/* Compute the Hessenberg matrix DomEig_mem->Hes*/
int DomEigComputeHess(DOMEIGMem domeig_mem)
{
  int retval;

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    printf(MSG_DOMEIG_MEM_FAIL);

    return -1;
  }

  /* Check if the dim of the matrix is less than 3
     Return immediately if dim <= 2 */
  if(domeig_mem->length < 3)
  {
    return 0;
  }

  int i, j;
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
      (retval < 0) ?
      printf(MSG_DOMEIG_ATIMES_FAIL_UNREC) :
      printf(MSG_DOMEIG_ATIMES_FAIL_REC);
      DomEigFree(&domeig_mem);

      return retval;
    }

    if(SUNModifiedGS(domeig_mem->V, domeig_mem->Hes, i + 1, domeig_mem->maxl, &(domeig_mem->Hes[i + 1][i])) != SUN_SUCCESS)
    {
      printf(MSG_DOMEIG_GS_FAIL);
      DomEigFree(&domeig_mem);

      return -1;
    }

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0)/domeig_mem->Hes[i + 1][i], domeig_mem->V[i+1], domeig_mem->V[i+1]);
  }

  return 0;
}

suncomplextype DomEigPowerIteration(DOMEIGMem domeig_mem)
{
  int retval;

  suncomplextype dom_eig, dom_eig_old;
  dom_eig.real = ZERO;
  dom_eig.imag = ZERO;
  dom_eig_old.real = ZERO;
  dom_eig_old.imag = ZERO;

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    printf(MSG_DOMEIG_MEM_FAIL);

    return dom_eig;
  }

  /* Check if ATimes is provided */
  if(domeig_mem->ATimes == NULL)
  {
    printf(MSG_DOMEIG_NULL_ATIMES);
    DomEigFree(&domeig_mem);

    return dom_eig;
  }

  int k;
  sunrealtype normq;

  for (k = 0; k < domeig_mem->max_powiter; k++)
  {
    retval = domeig_mem->ATimes(domeig_mem->Adata, domeig_mem->V[0], domeig_mem->q);
    if (retval != 0)
    {
      (retval < 0) ?
      printf(MSG_DOMEIG_ATIMES_FAIL_UNREC) :
      printf(MSG_DOMEIG_ATIMES_FAIL_REC);
      DomEigFree(&domeig_mem);

      dom_eig.real = ZERO;
      dom_eig.imag = ZERO;

      return dom_eig;
    }

    dom_eig.real = N_VDotProd(domeig_mem->V[0], domeig_mem->q); //Rayleigh quotient

    if(fabs(dom_eig.real - dom_eig_old.real) < domeig_mem->powiter_tol)
    {
      break;
    }

    normq = N_VDotProd(domeig_mem->q, domeig_mem->q);
    normq = SUNRsqrt(normq);
    N_VScale(SUN_RCONST(1.0)/normq, domeig_mem->q, domeig_mem->V[0]);

    dom_eig_old.real = dom_eig.real;
  }

  return dom_eig;
}

/* Estimate the dominant eigvalues of the Hessenberg matrix */
suncomplextype DomEigEstimate(DOMEIGMem domeig_mem) {
  suncomplextype dom_eig;
  dom_eig.real = ZERO;
  dom_eig.imag = ZERO;

  /* Check if DomEig memory is allocated */
  if(domeig_mem == NULL)
  {
    printf(MSG_DOMEIG_MEM_FAIL);

    return dom_eig;
  }

  if(domeig_mem->length < 3)
  {
    dom_eig = DomEigPowerIteration(domeig_mem);

    return dom_eig;
  }

  int n = domeig_mem->maxl;

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

  // Call LAPACK's dgeev function
  dgeev_(&jobvl, &jobvr, &n, domeig_mem->LAPACK_A, &lda, domeig_mem->LAPACK_wr, domeig_mem->LAPACK_wi, NULL, &ldvl, NULL, &ldvr, domeig_mem->LAPACK_work, &lwork, &info);

  if (info != 0) {
      printf(MSG_DOMEIG_LAPACK_FAIL, info);
      DomEigFree(&domeig_mem);

      return dom_eig;
  }

  //Following part will order the eigenvalues by their magnitude
  for (i = 0; i < n; i++) {
      domeig_mem->LAPACK_arr[i].real = domeig_mem->LAPACK_wr[i];
      domeig_mem->LAPACK_arr[i].imag = domeig_mem->LAPACK_wi[i];
  }

  // Sort the array using qsort
  qsort(domeig_mem->LAPACK_arr, n, sizeof(suncomplextype), domeig_Compare);

  // Update the original arrays
  for (i = 0; i < n; i++) {
      domeig_mem->LAPACK_wr[i] = domeig_mem->LAPACK_arr[i].real;
      domeig_mem->LAPACK_wi[i] = domeig_mem->LAPACK_arr[i].imag;
  }

  dom_eig.real = domeig_mem->LAPACK_wr[0];
  dom_eig.imag = domeig_mem->LAPACK_wi[0];

  return dom_eig;
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  domeig_CheckNVector:

  This routine checks if all required vector operations are
  present. If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype domeig_CheckNVector(N_Vector tmpl)
{ // TO DO: check required vector operations
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvdotprod == NULL) || (tmpl->ops->nvscale == NULL) ||
      (tmpl->ops->nvgetlength == NULL) || (tmpl->ops->nvspace == NULL))
  {
    return SUNFALSE;
  }
  return SUNTRUE;
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
  DomEigFree frees all DomEig memory.
  ---------------------------------------------------------------*/
void DomEigFree(DOMEIGMem* domeig_mem)
{
  DOMEIGMem arn_mem;

  /* nothing to do if domeig_mem is already NULL */
  if (domeig_mem == NULL) { return; }

  arn_mem = (DOMEIGMem)(*domeig_mem);

  if (arn_mem->q != NULL)
  {
    N_VDestroy(arn_mem->q);
    arn_mem->q = NULL;
  }
  if (arn_mem->V != NULL)
  {
    N_VDestroyVectorArray(arn_mem->V, arn_mem->maxl + 1);
    arn_mem->V = NULL;
  }

  if (arn_mem->LAPACK_A != NULL)
  {
    free(arn_mem->LAPACK_A);
    arn_mem->LAPACK_A = NULL;
  }

  if (arn_mem->LAPACK_wr != NULL)
  {
    free(arn_mem->LAPACK_wr);
    arn_mem->LAPACK_wr = NULL;
  }

  if (arn_mem->LAPACK_wi != NULL)
  {
    free(arn_mem->LAPACK_wi);
    arn_mem->LAPACK_wi = NULL;
  }

  if (arn_mem->LAPACK_arr != NULL)
  {
    free(arn_mem->LAPACK_arr);
    arn_mem->LAPACK_arr = NULL;
  }

  if (arn_mem->Hes != NULL)
  {
    free(arn_mem->Hes);
    arn_mem->Hes = NULL;
  }

  free(*domeig_mem);
  *domeig_mem = NULL;
}