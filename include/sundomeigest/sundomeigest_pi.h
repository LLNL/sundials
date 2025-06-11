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
#define SUNDOMEIGEST_PI_TOL_DEFAULT        SUN_RCONST(0.01)
#define SUNDOMEIGEST_MAX_PI_DEFAULT        100
#define SUNDOMEIGEST_PI_POWER_OF_A_DEFAULT 10

/* -----------------------------------------------------
 * Power Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct _SUNDomEigEstimatorContent_PI
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  N_Vector V, q; /* workspace vectors */

  sunindextype power_of_A; /* Power of A in the preprocessing; initial q = A^{power_of_A}q/||A^{power_of_A}q|| */

  sunrealtype powiter_tol;  /* Convergence criteria for the power iteration */
  sunrealtype resnorm;      /* Current residual of power iterations */
  sunindextype max_powiter; /* Maximum number of power iterations */
  sunindextype numiters;    /* Number of power iterations */
};

typedef struct _SUNDomEigEstimatorContent_PI* SUNDomEigEstimatorContent_PI;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_PI
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEst_PI(N_Vector q, sunindextype max_powiter,
                                   SUNContext sunctx);

SUNDIALS_EXPORT
SUNDomEigEstimator_Type SUNDomEigEst_PIGetType(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstInitialize_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetNumPreProcess_PI(SUNDomEigEstimator DEE,
                                             sunindextype numofperprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetATimes_PI(SUNDomEigEstimator DEE, void* A_data,
                                    SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PISetMaxPowerIter(SUNDomEigEstimator DEE,
                                          sunindextype max_powiter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstPreProcess_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimate_PI(SUNDomEigEstimator DEE, suncomplextype* dom_eig);

SUNDIALS_EXPORT
sunindextype SUNDomEigEstNumIters_PI(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstFree_PI(SUNDomEigEstimator DEE);

#ifdef __cplusplus
}
#endif

#endif
