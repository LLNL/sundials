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

struct SUNDomEigEstimatorContent_Power_
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector V, q; /* workspace vectors */

  int num_warmups;     /* Number of preprocessing iterations */
  long int max_iters;  /* Maximum number of power iterations */
  long int num_iters;  /* Number of iterations in last Estimate call */

  long int num_ATimes; /* Number of ATimes calls */

  sunrealtype rel_tol; /* Convergence criteria for the power iteration */
  sunrealtype res;     /* Residual from the last Estimate call */
};

typedef struct SUNDomEigEstimatorContent_Power_* SUNDomEigEstimatorContent_Power;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_Power
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_Power(N_Vector q, long int max_iters,
                                      sunrealtype rel_tol, SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetATimes_Power(SUNDomEigEstimator DEE, void* A_data,
                                        SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetMaxIters_Power(SUNDomEigEstimator DEE,
                                          long int max_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetNumPreprocessIters_Power(SUNDomEigEstimator DEE,
                                                    int num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetRelTol_Power(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Initialize_Power(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEig_Estimate_Power(SUNDomEigEstimator DEE,
                                    sunrealtype* lambdaR, sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetRes_Power(SUNDomEigEstimator DEE,
                                     sunrealtype* res);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumIters_Power(SUNDomEigEstimator DEE,
                                          long int* num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumATimesCalls_Power(SUNDomEigEstimator DEE,
                                                long int* num_ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Write_Power(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Destroy_Power(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
