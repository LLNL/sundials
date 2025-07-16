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
 * This is the header file for the Power Iteration (PI) implementation
 * of the SUNDomEigEst package.
 *
 * Note:
 *   - The definition of the generic SUNDomEigEstimator structure can
 *     be found in the header file sundials_domeigestimator.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDOMEIGEST_POWER_H
#define _SUNDOMEIGEST_POWER_H

#include <sundials/sundials_domeigestimator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------
 * Power Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct _SUNDomEigEstimatorContent_PI
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector V, q; /* workspace vectors */

  int numwarmups; /* Power of A in the preprocessing; initial q = A^{numwarmups}q/||A^{numwarmups}q|| */
  int max_iters;  /* Maximum number of power iterations */
  int curnumiters; /* Current number of power iterations */

  int maxnumiters; /* Maximum number of power iterations so far */
  int minnumiters; /* Minimum number of power iterations so far */

  long int nATimes; /* Number of ATimes calls */

  sunrealtype powiter_tol; /* Convergence criteria for the power iteration */
  sunrealtype curres;      /* Current residual of power iterations */
};

typedef struct _SUNDomEigEstimatorContent_PI* SUNDomEigEstimatorContent_PI;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_PI
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, int max_iters, SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetATimes_PI(SUNDomEigEstimator DEE, void* A_data,
                                     SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetMaxIters_PI(SUNDomEigEstimator DEE, int max_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetNumPreProcess_PI(SUNDomEigEstimator DEE,
                                            int numpreprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetTol_PI(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Initialize_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PreProcess_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEig_Estimate_PI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                 sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurRes_PI(SUNDomEigEstimator DEE, sunrealtype* curres);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurNumIters_PI(SUNDomEigEstimator DEE, int* curniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMaxNumIters_PI(SUNDomEigEstimator DEE, int* maxniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMinNumIters_PI(SUNDomEigEstimator DEE, int* minniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumATimesCalls_PI(SUNDomEigEstimator DEE,
                                             long int* nATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PrintStats_PI(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Destroy_PI(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
