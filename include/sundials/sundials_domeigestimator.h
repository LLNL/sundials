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

#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include "sundials/sundials_errors.h"
#include <sundials/sundials_iterative.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

// Struct to hold the real and imaginary parts
typedef struct {
    sunrealtype real;
    sunrealtype imag;
} suncomplextype;

/* -----------------------------------------------------------------
 * Implemented SUNDomEigEstimator types
 * ----------------------------------------------------------------- */

typedef enum
{
  SUNDOMEIG_POWER,
  SUNDOMEIG_ARNOLDI
} SUNDomEigEstimator_Type;

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
  SUNDomEigEstimator_Type (*gettype)(SUNDomEigEstimator);
  SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn);
  SUNErrCode (*setmaxpoweriter)(SUNDomEigEstimator, int);
  SUNErrCode (*initialize)(SUNDomEigEstimator);
  SUNErrCode (*preprocess)(SUNDomEigEstimator);
  SUNErrCode (*computehess)(SUNDomEigEstimator);
  SUNErrCode (*estimate)(SUNDomEigEstimator, suncomplextype*);
  int (*getnumofiters)(SUNDomEigEstimator);
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
void SUNDomEigEstFreeEmpty(SUNDomEigEstimator D);

SUNDIALS_EXPORT
SUNDomEigEstimator_Type SUNDomEigEstGetType(SUNDomEigEstimator D);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstSetATimes(SUNDomEigEstimator D, void* A_data,
                              SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstInitialize(SUNDomEigEstimator D, void* A_data,
                              SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstPreProcess(SUNDomEigEstimator D);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstComputeHess(SUNDomEigEstimator D);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimate(SUNDomEigEstimator D, suncomplextype* dom_eig);

/* -----------------------------------------------------------------
 * SUNDomEigEstimator return values
 * ----------------------------------------------------------------- */

#define SUNDEE_LAPACK_FAIL        "Error: LAPACK dgeev failed with info = %d\n"

#ifdef __cplusplus
}
#endif

#endif