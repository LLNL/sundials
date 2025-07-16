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
 * This is the header file for the Arnoldi Iteration (ArnI)
 * implementation of the SUNDomEigEst package.
 *
 * Note:
 *   - The definition of the generic SUNDomEigEstimator structure can
 *     be found in the header file sundials_domeigestimator.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDOMEIGEST_ARNOLDI_H
#define _SUNDOMEIGEST_ARNOLDI_H

#include <sundials/sundials_domeigestimator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------
 * Arnoldi Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct _SUNDomEigEstimatorContent_ArnI
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector *V, q; /* Krylov subspace vectors */

  int krydim;     /* Krylov subspace dimension */
  int numwarmups; /* Power of A in the preprocessing; initial q = A^{numwarmups}q/||A^{numwarmups}q|| */

  long int nATimes; /* Number of ATimes calls */

  sunrealtype* LAPACK_A; /* The vector which holds rows of the Hessenberg matrix in the given order */
  sunrealtype* LAPACK_wr;   /* Real parts of eigenvalues */
  sunrealtype* LAPACK_wi;   /* Imaginary parts of eigenvalues */
  sunrealtype* LAPACK_work; /* Workspace array */
  sunrealtype** LAPACK_arr; /* an array to sort eigenvalues*/

  sunrealtype** Hes; /* Hessenberg matrix Hes */
};

typedef struct _SUNDomEigEstimatorContent_ArnI* SUNDomEigEstimatorContent_ArnI;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_ArnI
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_ArnI(N_Vector q, int krydim, SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetATimes_ArnI(SUNDomEigEstimator DEE, void* A_data,
                                       SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetNumPreProcess_ArnI(SUNDomEigEstimator DEE,
                                              int numpreprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Initialize_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PreProcess_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_ComputeHess_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEig_Estimate_ArnI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                   sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumATimesCalls_ArnI(SUNDomEigEstimator DEE,
                                               long int* nATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PrintStats_ArnI(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Destroy_ArnI(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
