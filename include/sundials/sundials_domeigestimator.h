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
 * This is the header file for a generic SUNDomEigEst package.
 * -----------------------------------------------------------------*/

#ifndef _DOMEIGEST_H
#define _DOMEIGEST_H

/* TODO: Check to see if they are all required */
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default estimator parameters */
#define SUNDOMEIGEST_NUM_OF_WARMUPS_DEFAULT 0

/* -----------------------------------------------------------------
 * Implemented SUNDomEigEstimator types
 * ----------------------------------------------------------------- */

typedef enum
{
  SUNDSOMEIGESTIMATOR_POWER,
  SUNDSOMEIGESTIMATOR_ARNOLDI
} SUNDomEigEstimator_ID;

/* -----------------------------------------------------------------
 * Generic definition of SUNDomEigEstimator
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNDomEigEstimator_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNDomEigEstimator_Ops* SUNDomEigEstimator_Ops;

/* Forward reference for pointer to SUNDomEigEstimator object */
typedef _SUNDIALS_STRUCT_ _generic_SUNDomEigEstimator* SUNDomEigEstimator;

/* Structure containing function pointers to estimator operations */
struct _generic_SUNDomEigEstimator_Ops
{
  SUNDomEigEstimator_ID (*getid)(SUNDomEigEstimator);
  SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn);
  SUNErrCode (*setmaxpoweriter)(SUNDomEigEstimator, int);
  SUNErrCode (*setnumofperprocess)(SUNDomEigEstimator, int);
  SUNErrCode (*settol)(SUNDomEigEstimator, sunrealtype);
  SUNErrCode (*initialize)(SUNDomEigEstimator);
  SUNErrCode (*preprocess)(SUNDomEigEstimator);
  SUNErrCode (*computehess)(SUNDomEigEstimator);
  SUNErrCode (*estimate)(SUNDomEigEstimator, sunrealtype*, sunrealtype*);
  SUNErrCode (*getnumofiters)(SUNDomEigEstimator, int*);
  SUNErrCode (*getres)(SUNDomEigEstimator, sunrealtype*);
  SUNErrCode (*free)(SUNDomEigEstimator);
};

/* A estimator is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of estimator
   operations corresponding to that implementation. */
struct _generic_SUNDomEigEstimator
{
  void* content;
  SUNDomEigEstimator_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNDomEigEstimator module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEstNewEmpty(SUNContext sunctx);

SUNDIALS_EXPORT
void SUNDomEigEstFreeEmpty(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNDomEigEstimator_ID SUNDomEigEstGetID(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetATimes(SUNDomEigEstimator DEE, void* A_data,
                                 SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetMaxPowerIter(SUNDomEigEstimator DEE, int max_powiter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetNumPreProcess(SUNDomEigEstimator DEE,
                                        int numofperprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetTol(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstInitialize(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstPreProcess(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstComputeHess(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                             sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstNumIters(SUNDomEigEstimator DEE, int* niter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstRes(SUNDomEigEstimator DEE, sunrealtype* res);

#ifdef __cplusplus
}
#endif

#endif
