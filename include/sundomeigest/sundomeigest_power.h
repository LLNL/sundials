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

struct _SUNDomEigEstimatorContent_Power
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector V, q; /* workspace vectors */

  int num_warmups; /* Power of A in the preprocessing; initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
  int max_iters;   /* Maximum number of power iterations */
  int cur_num_iters; /* Current number of power iterations */

  int max_num_iters; /* Maximum number of power iterations so far */
  int min_num_iters; /* Minimum number of power iterations so far */

  long int num_ATimes; /* Number of ATimes calls */

  sunrealtype powiter_tol; /* Convergence criteria for the power iteration */
  sunrealtype cur_res;     /* Current residual of power iterations */
};

typedef struct _SUNDomEigEstimatorContent_Power* SUNDomEigEstimatorContent_Power;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_Power
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_Power(N_Vector q, int max_iters,
                                      SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetATimes_Power(SUNDomEigEstimator DEE, void* A_data,
                                        SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetMaxIters_Power(SUNDomEigEstimator DEE, int max_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetNumPreProcess_Power(SUNDomEigEstimator DEE,
                                               int numpreprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetTol_Power(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Initialize_Power(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEig_Estimate_Power(SUNDomEigEstimator DEE,
                                    sunrealtype* lambdaR, sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurRes_Power(SUNDomEigEstimator DEE,
                                        sunrealtype* cur_res);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurNumIters_Power(SUNDomEigEstimator DEE,
                                             int* curniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMaxNumIters_Power(SUNDomEigEstimator DEE,
                                             int* max_niter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMinNumIters_Power(SUNDomEigEstimator DEE,
                                             int* min_niter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumATimesCalls_Power(SUNDomEigEstimator DEE,
                                                long int* num_ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PrintStats_Power(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Destroy_Power(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
