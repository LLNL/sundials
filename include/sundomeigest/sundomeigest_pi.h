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

#ifndef _DOMEIGEST_PI_H
#define _DOMEIGEST_PI_H

#include <sundials/sundials_domeigestimator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default Power Iteration parameters */
#define SUNDOMEIGEST_PI_TOL_DEFAULT SUN_RCONST(0.01)
#define SUNDOMEIGEST_MAX_PI_DEFAULT 100

/* -----------------------------------------------------
 * Power Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct _SUNDomEigEstimatorContent_PI
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector V, q; /* workspace vectors */

  int numwarmups; /* Power of A in the preprocessing; initial q = A^{numwarmups}q/||A^{numwarmups}q|| */

  sunrealtype powiter_tol; /* Convergence criteria for the power iteration */
  sunrealtype res;         /* Current residual of power iterations */
  int max_powiter;         /* Maximum number of power iterations */
  int numiters;            /* Current number of power iterations */
};

typedef struct _SUNDomEigEstimatorContent_PI* SUNDomEigEstimatorContent_PI;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_PI
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, int max_powiter,
                                   SUNContext sunctx);

SUNDIALS_EXPORT
SUNDomEigEstimator_ID SUNDomEigEst_PIGetID(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetATimes_PI(SUNDomEigEstimator DEE, void* A_data,
                                    SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstInitialize_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetNumPreProcess_PI(SUNDomEigEstimator DEE,
                                           int numofperprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetTol_PI(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetMaxPowerIter_PI(SUNDomEigEstimator DEE,
                                          int max_powiter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstPreProcess_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimate_PI(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstNumIters_PI(SUNDomEigEstimator DEE, int* niter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstRes_PI(SUNDomEigEstimator DEE, sunrealtype* res);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstFree_PI(SUNDomEigEstimator DEE);

#ifdef __cplusplus
}
#endif

#endif
