/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS accuracy-based adaptivity controller class. These
 * objects estimate step sizes for time integration methods such
 * that the next step solution should satisfy a desired temporal
 * accuracy.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_ADAPTCONTROLLER_H
#define _SUNDIALS_ADAPTCONTROLLER_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * SUNAdaptController types (currently, only "H" is implemented;
 * others are planned):
 *    NONE - empty controller (does nothing)
 *    H    - controls a single-rate step size
 * ----------------------------------------------------------------- */

typedef enum
{
  SUN_ADAPTCONTROLLER_NONE,
  SUN_ADAPTCONTROLLER_H
} SUNAdaptController_Type;

/* -----------------------------------------------------------------
 * Generic definition of SUNAdaptController
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNAdaptController_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNAdaptController_Ops* SUNAdaptController_Ops;

/* Forward reference for pointer to SUNAdaptController object */
typedef _SUNDIALS_STRUCT_ _generic_SUNAdaptController* SUNAdaptController;

/* Structure containing function pointers to controller operations  */
struct _generic_SUNAdaptController_Ops
{
  /* REQUIRED of all controller implementations. */
  SUNAdaptController_Type (*gettype)(SUNAdaptController C);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_H type. */
  int (*estimatestep)(SUNAdaptController C, sunrealtype h, int p,
                      sunrealtype dsm, sunrealtype* hnew);

  /* OPTIONAL for all SUNAdaptController implementations. */
  int (*destroy)(SUNAdaptController C);
  int (*reset)(SUNAdaptController C);
  int (*setdefaults)(SUNAdaptController C);
  int (*write)(SUNAdaptController C, FILE* fptr);
  int (*seterrorbias)(SUNAdaptController C, sunrealtype bias);
  int (*updateh)(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
  int (*space)(SUNAdaptController C, long int *lenrw, long int *leniw);
};

/* A SUNAdaptController is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct _generic_SUNAdaptController
{
  void* content;
  SUNAdaptController_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNAdaptController module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNAdaptController data structure. */
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx);

/* Function to report the type of a SUNAdaptController object. */
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C);

/* Function to deallocate a SUNAdaptController object.

   Any return value other than SUNADAPTCONTROLLER_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNAdaptController_Destroy(SUNAdaptController C);

/* Main step size controller function.  This is called following
   a time step with size 'h' and local error factor 'dsm', and the
   controller should estimate 'hnew' so that the ensuing step
   will have 'dsm' value JUST BELOW 1.

   Any return value other than SUNADAPTCONTROLLER_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h,
                                    int p, sunrealtype dsm, sunrealtype* hnew);

/* Function to reset the controller to its initial state, e.g., if
   it stores a small number of previous dsm or step size values. */
SUNDIALS_EXPORT
int SUNAdaptController_Reset(SUNAdaptController C);

/* Function to set the controller parameters to their default values. */
SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults(SUNAdaptController C);

/* Function to write all controller parameters to the indicated file
   pointer. */
SUNDIALS_EXPORT
int SUNAdaptController_Write(SUNAdaptController C, FILE* fptr);

/* Function to set an error bias factor to use for scaling the local error
   'dsm' factors above. */
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias);

/* Function to notify a controller of type SUN_ADAPTCONTROLLER_H that
   a successful time step was taken with stepsize h and local error factor
   dsm, indicating that these can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNAdaptController_UpdateH(SUNAdaptController C, sunrealtype h, sunrealtype dsm);

/* Function to return the memory requirements of the controller object. */
SUNDIALS_EXPORT
int SUNAdaptController_Space(SUNAdaptController C, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNAdaptController error codes
 * ----------------------------------------------------------------- */

#define SUNADAPTCONTROLLER_SUCCESS           0     /* function successfull        */
#define SUNADAPTCONTROLLER_ILL_INPUT         -1001 /* illegal function input      */
#define SUNADAPTCONTROLLER_MEM_FAIL          -1002 /* failed memory access/alloc  */
#define SUNADAPTCONTROLLER_USER_FCN_FAIL     -1003 /* user-supplied fcn failure */
#define SUNADAPTCONTROLLER_OPERATION_FAIL    -1004 /* catchall failure code       */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_ADAPTCONTROLLER_H */
