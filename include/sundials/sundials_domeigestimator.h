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

#ifndef _SUNDOMEIGEST_H
#define _SUNDOMEIGEST_H

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
#define DEE_NUM_OF_WARMUPS_PI_DEFAULT   0
#define DEE_NUM_OF_WARMUPS_ARNI_DEFAULT 100

/* Default Power Iteration parameters */
#define DEE_TOL_DEFAULT      SUN_RCONST(0.01)
#define DEE_MAX_ITER_DEFAULT 100

/* Default Arnoldi Iteration parameters */
#define DEE_KRYLOV_DIM_DEFAULT 3

//#define DEE_LAPACK_FAIL        "Error: LAPACK dgeev failed with info = %d\n"

/* -----------------------------------------------------------------
 * Generic definition of SUNDomEigEstimator (DEE)
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNDomEigEstimator_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNDomEigEstimator_Ops* SUNDomEigEstimator_Ops;

/* Forward reference for pointer to SUNDomEigEstimator object */
typedef _SUNDIALS_STRUCT_ _generic_SUNDomEigEstimator* SUNDomEigEstimator;

/* Structure containing function pointers to estimator operations */
struct _generic_SUNDomEigEstimator_Ops
{
  SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn);
  SUNErrCode (*setmaxiters)(SUNDomEigEstimator, int);
  SUNErrCode (*setnumpreprocess)(SUNDomEigEstimator, int);
  SUNErrCode (*settol)(SUNDomEigEstimator, sunrealtype);
  SUNErrCode (*initialize)(SUNDomEigEstimator);
  SUNErrCode (*preprocess)(SUNDomEigEstimator);
  SUNErrCode (*computehess)(SUNDomEigEstimator);
  SUNErrCode (*estimate)(SUNDomEigEstimator, sunrealtype*, sunrealtype*);
  SUNErrCode (*getcurres)(SUNDomEigEstimator, sunrealtype*);
  SUNErrCode (*getcurniters)(SUNDomEigEstimator, int*);
  SUNErrCode (*getmaxniters)(SUNDomEigEstimator, int*);
  SUNErrCode (*getminniters)(SUNDomEigEstimator, int*);
  SUNErrCode (*getnumatimescalls)(SUNDomEigEstimator, long int*);
  SUNErrCode (*printstats)(SUNDomEigEstimator, FILE*);
  SUNErrCode (*free)(SUNDomEigEstimator);
};

/* An estimator is a structure with an implementation-dependent
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
SUNDomEigEstimator SUNDomEigEst_NewEmpty(SUNContext sunctx);

SUNDIALS_EXPORT
void SUNDomEigEst_FreeEmpty(SUNDomEigEstimator DEE); 

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* A_data,
                                  SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, int max_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE,
                                         int numpreprocess);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_SetTol(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_Initialize(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PreProcess(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_ComputeHess(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                              sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurRes(SUNDomEigEstimator DEE, sunrealtype* curres);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetCurNumIters(SUNDomEigEstimator DEE, int* curniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMaxNumIters(SUNDomEigEstimator DEE, int* maxniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetMinNumIters(SUNDomEigEstimator DEE, int* minniter);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE, long int* nATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEst_PrintStats(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstFree(SUNDomEigEstimator DEE);

#ifdef __cplusplus
}
#endif

#endif
