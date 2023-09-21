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
 *    NONE    - empty controller (does nothing)
 *    H       - controls a single-rate step size
 *    HQ      - controls a single-rate step size and method order
 *    MRI_H   - controls two multirate step sizes
 *    MRI_TOL - controls slow and fast relative tolerances
 * ----------------------------------------------------------------- */

typedef enum
{
  SUN_ADAPTCONTROLLER_NONE,
  SUN_ADAPTCONTROLLER_H,
  SUN_ADAPTCONTROLLER_HQ,
  SUN_ADAPTCONTROLLER_MRI_H,
  SUN_ADAPTCONTROLLER_MRI_TOL
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
  int (*estimatestep)(SUNAdaptController C, realtype h,
                      realtype dsm, realtype* hnew);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_HQ type. */
  int (*estimatestepandorder)(SUNAdaptController C, realtype h, int q,
                              realtype dsm, realtype* hnew, int *qnew);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_MRI_H type. */
  int (*estimatemristeps)(SUNAdaptController C, realtype H, realtype h,
                          realtype DSM, realtype dsm,
                          realtype* Hnew, realtype *hnew);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_MRI_TOL type. */
  int (*estimatesteptol)(SUNAdaptController C, realtype H, realtype tolfac,
                         realtype DSM, realtype dsm, realtype *Hnew,
                         realtype* tolfacnew);

  /* OPTIONAL for all SUNAdaptController implementations. */
  int (*destroy)(SUNAdaptController C);
  int (*reset)(SUNAdaptController C);
  int (*setdefaults)(SUNAdaptController C);
  int (*write)(SUNAdaptController C, FILE* fptr);
  int (*setmethodorder)(SUNAdaptController C, int q);
  int (*setembeddingorder)(SUNAdaptController C, int p);
  int (*seterrorbias)(SUNAdaptController C, realtype bias);
  int (*update)(SUNAdaptController C, realtype h, realtype dsm);
  int (*updatemrih)(SUNAdaptController C, realtype H, realtype h,
                    realtype DSM, realtype dsm);
  int (*updatemritol)(SUNAdaptController C, realtype H, realtype tolfac,
                      realtype DSM, realtype dsm);
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
int SUNAdaptController_EstimateStep(SUNAdaptController C, realtype h,
                                    realtype dsm, realtype* hnew);

/* Combined step size + order controller function.  This is called
   following a time step with size 'h' and order 'q' that has local
   error factor 'dsm'.  The controller should estimate 'hnew' and
   'qnew' so that the ensuing step will have 'dsm' value JUST BELOW 1
   with minimal computational effort. */
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStepAndOrder(SUNAdaptController C, realtype h,
                                            int q, realtype dsm,
                                            realtype* hnew, int *qnew);

/* Combined slow/fast multirate step size controller function.  This
   is called following a slow multirate time step with sizes 'H' and
   'h' (slow and fast, resp.), and error factors 'DSM' and 'dsm'
   (slow and fast, resp.). The controller should estimate slow and
   fast steps 'Hnew' and 'hnew', resp., so that the ensuing step will
   have 'DSM' and 'dsm' values JUST BELOW 1 with minimal computational
   effort. */
SUNDIALS_EXPORT
int SUNAdaptController_EstimateMRISteps(SUNAdaptController C, realtype H,
                                        realtype h, realtype DSM, realtype dsm,
                                        realtype* Hnew, realtype *hnew);

/* Combined slow step/fast tolerance multirate controller function.
   This is called following a slow multirate time step with size 'H'
   and fast/slow relative tolerance ratio 'tolfac', and error factors
   'DSM' and 'dsm' (slow and fast, resp.).  The controller should
   estimate slow stepsize 'Hnew' and updated relative tolerance ratio
   'tolfacnew', so that the ensuing step will have 'DSM' and 'dsm'
   values JUST BELOW 1 with minimal computational effort. */
SUNDIALS_EXPORT
int SUNAdaptController_EstimateStepTol(SUNAdaptController C, realtype H,
                                       realtype tolfac, realtype DSM,
                                       realtype dsm, realtype *Hnew,
                                       realtype* tolfacnew);

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

/* Function to set the asymptotic order of accuracy for the method. */
SUNDIALS_EXPORT
int SUNAdaptController_SetMethodOrder(SUNAdaptController C, int q);

/* Function to set the asymptotic order of accuracy for the embedding. */
SUNDIALS_EXPORT
int SUNAdaptController_SetEmbeddingOrder(SUNAdaptController C, int p);

/* Function to set an error bias factor to use for scaling the local error
   'dsm' factors above. */
SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias(SUNAdaptController C, realtype bias);

/* Function to notify the controller of a successful time step with size
   h and local error factor dsm, indicating that the step size or local
   error factor can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNAdaptController_Update(SUNAdaptController C, realtype h, realtype dsm);

/* Function to notify the controller of a successful multirate time step
   with sizes H and h, and local error factors DSM and dsm, indicating that
   the step sizes or local error factors can be saved for subsequent
   controller functions. */
SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRIH(SUNAdaptController C, realtype H, realtype h,
                                  realtype DSM, realtype dsm);

/* Function to notify the controller of a successful multirate time step
   with size H and fast tolerance factor tolfac, and local error factors
   DSM and dsm, indicating that the step size, tolerance factor, or local
   error factors can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRITol(SUNAdaptController C, realtype H, realtype tolfac,
                                    realtype DSM, realtype dsm);

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
