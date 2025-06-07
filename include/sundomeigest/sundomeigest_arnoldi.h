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

#ifndef _DOMEIGEST_ARNI_H
#define _DOMEIGEST_ARNI_H

#include <sundials/sundials_domeigestimator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default Arnoldi Iteration parameters */
#define SUNDOMEIGEST_ARN_MAXL_DEFAULT        3
#define SUNDOMEIGEST_PI_POWER_OF_A_DEFAULT  10
#define SUNDOMEIGEST_LAPACK_FAIL "Error: LAPACK dgeev failed with info = %d\n"

/* -----------------------------------------------------
 * Arnoldi Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct _SUNDomEigEstimatorContent_ArnI
{
  SUNATimesFn ATimes;   /* User provided ATimes function */
  void* ATdata;         /* ATimes function data*/

  N_Vector *V, q;       /* Krylov subspace vectors */

  int maxl;             /* Krylov subspace dimension */
  int power_of_A;       /* Power of A in the preprocessing; initial q = A^{power_of_A}q/||A^{power_of_A}q|| */

  sunrealtype* LAPACK_A;      /* The vector which holds rows of the Hessenberg matrix in the given order */
  sunrealtype* LAPACK_wr;     /* Real parts of eigenvalues */
  sunrealtype* LAPACK_wi;     /* Imaginary parts of eigenvalues */
  sunrealtype* LAPACK_work;   /* Workspace array */
  suncomplextype* LAPACK_arr; /* an array to sort eigenvalues*/

  sunrealtype **Hes;    /* Hessenberg matrix Hes */
};

typedef struct _SUNDomEigEstimatorContent_ArnI* SUNDomEigEstimatorContent_ArnI;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_ArnI
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_ArnI(N_Vector q, int maxl, SUNContext sunctx);

SUNDIALS_EXPORT
SUNDomEigEstimator_Type SUNDomEigEst_ArnIGetType(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetATimes_ArnI(SUNDomEigEstimator DEE, void* A_data,
                                    SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstInitialize_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstPreProcess_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstComputeHess_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimate_ArnI(SUNDomEigEstimator DEE, suncomplextype* dom_eig);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstFree_ArnI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
int SUNDomEigEstNumIters_ArnI(SUNDomEigEstimator DEE);

#ifdef __cplusplus
}
#endif

#endif
