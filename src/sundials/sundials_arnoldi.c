/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the eigenvalue estimation of
 * the SUNARNOLDI package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_arnoldi_impl.h>

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

void* ArnoldiCreate(SUNATimesFn ATimes, void* Adata,
                N_Vector q, int maxl, SUNContext sunctx)
{
  ARNOLDIMem arnoldi_mem;

  /* Test if Atimes is provided */
  if (ATimes == NULL)
  {
    printf(MSG_ARNOLDI_NULL_ATIMES);

    return NULL;
  }

  /* Test if q is provided */
  if (q == NULL)
  {
    printf(MSG_ARNOLDI_NULL_q);

    return NULL;
  }

  /* Check if maxl > 2 */
  if (maxl < 3)
  {
    printf(MSG_ARNOLDI_NOT_ENOUGH_ITER);

    return NULL;
  }

  if (sunctx == NULL)
  {
    printf(MSG_ARNOLDI_NULL_SUNCTX);

    return NULL;
  }

  /* Test if all required vector operations are implemented */
  if (!arnoldi_CheckNVector(q))
  {
    printf(MSG_ARNOLDI_BAD_NVECTOR);

    return NULL;
  }

  /* Allocate ARNOLDIMem structure, and initialize to zero */
  arnoldi_mem = (ARNOLDIMem)calloc(1, sizeof(*arnoldi_mem));

  /* Copy the inputs into ARNOLDI memory */
  arnoldi_mem->ATimes      = ATimes;
  arnoldi_mem->Adata       = Adata;
  arnoldi_mem->q           = q;
  arnoldi_mem->maxl        = maxl;

  /* Set the default power of A to start with (# of warm-ups) */
  arnoldi_mem->power_of_A  = DEFAULT_POWER_OF_A;

  /* Hessenberg matrix Hes */
  if (arnoldi_mem->Hes == NULL)
  {
    arnoldi_mem->Hes =
      (sunrealtype**)malloc((maxl + 1) * sizeof(sunrealtype*));

    for (int k = 0; k <= maxl; k++)
    {
      arnoldi_mem->Hes[k] = NULL;
      arnoldi_mem->Hes[k] = (sunrealtype*)malloc(maxl * sizeof(sunrealtype));
    }
  }

  /* Krylov subspace vectors */
  if (arnoldi_mem->V == NULL)
  {
    arnoldi_mem->V = N_VCloneVectorArray(maxl + 1, q);
  }

  /* Unitize the initial vector V[0] */
  sunrealtype normq = N_VDotProd(q, q);
  normq = SUNRsqrt(normq);
  N_VScale(SUN_RCONST(1.0)/normq, arnoldi_mem->q, arnoldi_mem->V[0]);

  return (void*)arnoldi_mem;
}

/* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
int ArnoldiPreProcess(ARNOLDIMem arnoldi_mem)
{
  int retval;

  /* Check if Arnoldi memory is allocated */
  if(arnoldi_mem == NULL)
  {
    printf(MSG_ARNOLDI_MEM_FAIL);

    return -1;
  }

  /* Check if ATimes is provided */
  if(arnoldi_mem->ATimes == NULL)
  {
    printf(MSG_ARNOLDI_NULL_ATIMES);
    ArnoldiFree(&arnoldi_mem);

    return -1;
  }

  sunrealtype normq;

  /* Set the initial q = A^{power_of_A}q/||A^{power_of_A}q|| */
  for(int i = 0; i < arnoldi_mem->power_of_A; i++) {
    retval = arnoldi_mem->ATimes(arnoldi_mem->Adata, arnoldi_mem->V[0], arnoldi_mem->q);
    if (retval != 0)
    {
      (retval < 0) ?
      printf(MSG_ARNOLDI_ATIMES_FAIL_UNREC) :
      printf(MSG_ARNOLDI_ATIMES_FAIL_REC);
      ArnoldiFree(&arnoldi_mem);

      return retval;
    }
    normq = N_VDotProd(arnoldi_mem->q, arnoldi_mem->q);
    normq = SUNRsqrt(normq);
    N_VScale(SUN_RCONST(1.0)/normq, arnoldi_mem->q, arnoldi_mem->V[0]);
  }

  return 0;
}

/* Compute the Hessenberg matrix Arnoldi_mem->Hes*/
int ArnoldiComputeHess(ARNOLDIMem arnoldi_mem)
{
  int retval;

  /* Check if Arnoldi memory is allocated */
  if(arnoldi_mem == NULL)
  {
    printf(MSG_ARNOLDI_MEM_FAIL);

    return -1;
  }

  /* Initialize the Hessenberg matrix Hes with zeros */
  for (int i = 0; i < arnoldi_mem->maxl; i++)
  {
    for (int j = 0; j < arnoldi_mem->maxl; j++) { arnoldi_mem->Hes[i][j] = ZERO; }
  }

  for (int i = 0; i < arnoldi_mem->maxl; i++)
  {
    /* Compute the next Krylov vector */
    retval = arnoldi_mem->ATimes(arnoldi_mem->Adata, arnoldi_mem->V[i], arnoldi_mem->V[i+1]);
    if (retval != 0)
    {
      (retval < 0) ?
      printf(MSG_ARNOLDI_ATIMES_FAIL_UNREC) :
      printf(MSG_ARNOLDI_ATIMES_FAIL_REC);
      ArnoldiFree(&arnoldi_mem);

      return retval;
    }

    if(SUNModifiedGS(arnoldi_mem->V, arnoldi_mem->Hes, i + 1, arnoldi_mem->maxl, &(arnoldi_mem->Hes[i + 1][i])) != SUN_SUCCESS)
    {
      printf(MSG_ARNOLDI_GS_FAIL);
      ArnoldiFree(&arnoldi_mem);

      return -1;
    }

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0)/arnoldi_mem->Hes[i + 1][i], arnoldi_mem->V[i+1], arnoldi_mem->V[i+1]);
  }

  return 0;
}

/* Estimate the dominant eigvalues of the Hessenberg matrix */
suncomplextype ArnoldiEstimate(ARNOLDIMem arnoldi_mem) {
  suncomplextype dom_eig;
  dom_eig.real = ZERO;
  dom_eig.imag = ZERO;

  /* Check if Arnoldi memory is allocated */
  if(arnoldi_mem == NULL)
  {
    printf(MSG_ARNOLDI_MEM_FAIL);

    return dom_eig;
  }

  int n = arnoldi_mem->maxl;

  /* Create the vector A which holds rows of the Hessenberg matrix in the given order */
  sunrealtype* A;
  A = (sunrealtype*)malloc((n*n) * sizeof(sunrealtype));
  int k = 0;
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          A[k] = arnoldi_mem->Hes[i][j];
          k++;
      }
  }

  sunrealtype wr[n], wi[n];     // Real and imaginary parts of eigenvalues
  sunrealtype vl[n*n], vr[n*n]; // Left and right eigenvectors (optional, can be NULL)
  int lda = n, ldvl = n, ldvr = n;
  int info, lwork = 4 * n;
  sunrealtype work[4 * n]; // Workspace array

  char jobvl = 'N'; // Do not compute left eigenvectors
  char jobvr = 'N'; // Do not compute right eigenvectors

  // Call LAPACK's dgeev function
  dgeev_(&jobvl, &jobvr, &n, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  if (info != 0) {
      printf(MSG_ARNOLDI_LAPACK_FAIL, info);
      ArnoldiFree(&arnoldi_mem);

      return dom_eig;
  }

  //Following part will order the eigenvalues by their magnitude
  {
    // Create an array of suncomplextype structs
    suncomplextype *arr = (suncomplextype *)malloc(n * sizeof(suncomplextype));
    for (int i = 0; i < n; i++) {
        arr[i].real = wr[i];
        arr[i].imag = wi[i];
    }

    // Sort the array using qsort
    qsort(arr, n, sizeof(suncomplextype), arnoldi_Compare);

    // Update the original arrays
    for (int i = 0; i < n; i++) {
        wr[i] = arr[i].real;
        wi[i] = arr[i].imag;
    }

    // Cleanup
    free(arr);
  }

  // // Print eigenvalues
  // printf("\nEigenvalues:\n");
  // for (int i = 0; i < n; i++) {
  //     if (wi[i] == 0.0) {
  //         printf("%f\n", wr[i]); // sunrealtype eigenvalue
  //     } else {
  //         printf("%f + %fi\n", wr[i], wi[i]); // suncomplextype eigenvalue
  //     }
  // }

  dom_eig.real = wr[0];
  dom_eig.imag = wi[0];

  return dom_eig;
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  arnoldi_CheckNVector:

  This routine checks if all required vector operations are
  present. If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arnoldi_CheckNVector(N_Vector tmpl)
{ // TO DO: check required vector operations
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvscale == NULL) || (tmpl->ops->nvwrmsnorm == NULL) ||
      (tmpl->ops->nvspace == NULL))
  {
    return SUNFALSE;
  }
  return SUNTRUE;
}

// Function to calculate the magnitude of a suncomplextype number
sunrealtype arnoldi_Magnitude(const suncomplextype *c) {
    return sqrt(c->real * c->real + c->imag * c->imag);
}

// Comparison function for qsort
int arnoldi_Compare(const void *a, const void *b) {
    const suncomplextype *c1 = (const suncomplextype *)a;
    const suncomplextype *c2 = (const suncomplextype *)b;
    sunrealtype mag1 = arnoldi_Magnitude(c1);
    sunrealtype mag2 = arnoldi_Magnitude(c2);
    return (mag2 > mag1) - (mag2 < mag1); // Descending order
}

/*---------------------------------------------------------------
  ArnoldiFree frees all Arnoldi memory.
  ---------------------------------------------------------------*/
void ArnoldiFree(ARNOLDIMem* arnoldi_mem)
{
  ARNOLDIMem arn_mem;

  /* nothing to do if arnoldi_mem is already NULL */
  if (arnoldi_mem == NULL) { return; }

  arn_mem = (ARNOLDIMem)(*arnoldi_mem);

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

  if (arn_mem->Hes != NULL)
  {
    free(arn_mem->Hes);
    arn_mem->Hes = NULL;
  }

  free(*arnoldi_mem);
  *arnoldi_mem = NULL;
}