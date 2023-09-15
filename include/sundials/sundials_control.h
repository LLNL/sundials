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

#ifndef _SUNDIALS_CONTROLLER_H
#define _SUNDIALS_CONTROLLER_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * SUNControl types (currently, only "H" is implemented; others
 * are planned):
 *    NONE    - empty controller (does nothing)
 *    H       - controls a single-rate step size
 *    HQ      - controls a single-rate step size and method order
 *    MRI_H   - controls two multirate step sizes
 *    MRI_TOL - controls slow and fast relative tolerances
 * ----------------------------------------------------------------- */

typedef enum
{
  SUNDIALS_CONTROL_NONE,
  SUNDIALS_CONTROL_H,
  SUNDIALS_CONTROL_HQ,
  SUNDIALS_CONTROL_MRI_H,
  SUNDIALS_CONTROL_MRI_TOL
} SUNControl_Type;

/* -----------------------------------------------------------------
 * Generic definition of SUNControl
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNControl_Ops object */
typedef _SUNDIALS_STRUCT_ generic_SUNControl_Ops_* SUNControl_Ops;

/* Forward reference for pointer to SUNControl object */
typedef _SUNDIALS_STRUCT_ generic_SUNControl_* SUNControl;

/* Structure containing function pointers to controller operations  */
struct generic_SUNControl_Ops_
{
  /* REQUIRED of all controller implementations. */
  SUNControl_Type (*gettype)(SUNControl C);
  int (*destroy)(SUNControl C);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_H type. */
  int (*estimatestep)(SUNControl C, realtype h,
                      realtype dsm, realtype* hnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_HQ type. */
  int (*estimatestepandorder)(SUNControl C, realtype h, int q,
                              realtype dsm, realtype* hnew, int *qnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_MRI_H type. */
  int (*estimatemristeps)(SUNControl C, realtype H, realtype h,
                          realtype DSM, realtype dsm,
                          realtype* Hnew, realtype *hnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_MRI_TOL type. */
  int (*estimatesteptol)(SUNControl C, realtype H, realtype tolfac,
                         realtype DSM, realtype dsm, realtype *Hnew,
                         realtype* tolfacnew);

  /* OPTIONAL for all SUNControl implementations. */
  int (*reset)(SUNControl C);
  int (*setdefaults)(SUNControl C);
  int (*write)(SUNControl C, FILE* fptr);
  int (*setmethodorder)(SUNControl C, int q);
  int (*setembeddingorder)(SUNControl C, int p);
  int (*seterrorbias)(SUNControl C, realtype bias);
  int (*update)(SUNControl C, realtype h, realtype dsm);
  int (*updatemrih)(SUNControl C, realtype H, realtype h,
                    realtype DSM, realtype dsm);
  int (*updatemritol)(SUNControl C, realtype H, realtype tolfac,
                      realtype DSM, realtype dsm);
  int (*space)(SUNControl C, long int *lenrw, long int *leniw);
#ifdef __cplusplus
  generic_SUNControl_Ops_() = default;
#endif

};

/* A SUNControl is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct generic_SUNControl_
{
  void* content;
  SUNControl_Ops ops;
  SUNContext sunctx;
#ifdef __cplusplus
  generic_SUNControl_() = default;
#endif
};

/* -----------------------------------------------------------------
 * Functions exported by SUNControl module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNHeuristics data structure. */
SUNDIALS_EXPORT
SUNControl SUNControl_NewEmpty(SUNContext sunctx);

/* Function to report the type of a SUNHeuristics object. */
SUNDIALS_EXPORT
SUNControl_Type SUNControl_GetType(SUNControl C);

/* Function to deallocate a SUNHeuristics object.

   Any return value other than SUNCONTROL_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNControl_Destroy(SUNControl C);

/* Main step size controller function.  This is called following
   a time step with size 'h' and local error factor 'dsm', and the
   controller should estimate 'hnew' so that the ensuing step
   will have 'dsm' value JUST BELOW 1.

   Any return value other than SUNCONTROL_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNControl_EstimateStep(SUNControl C, realtype h, realtype dsm,
                            realtype* hnew);

/* Combined step size + order controller function.  This is called
   following a time step with size 'h' and order 'q' that has local
   error factor 'dsm'.  The controller should estimate 'hnew' and
   'qnew' so that the ensuing step will have 'dsm' value JUST BELOW 1
   with minimal computational effort. */
SUNDIALS_EXPORT
int SUNControl_EstimateStepAndOrder(SUNControl C, realtype h, int q,
                                    realtype dsm, realtype* hnew,
                                    int *qnew);

/* Combined slow/fast multirate step size controller function.  This
   is called following a slow multirate time step with sizes 'H' and
   'h' (slow and fast, resp.), and error factors 'DSM' and 'dsm'
   (slow and fast, resp.). The controller should estimate slow and
   fast steps 'Hnew' and 'hnew', resp., so that the ensuing step will
   have 'DSM' and 'dsm' values JUST BELOW 1 with minimal computational
   effort. */
SUNDIALS_EXPORT
int SUNControl_EstimateMRISteps(SUNControl C, realtype H, realtype h,
                                realtype DSM, realtype dsm,
                                realtype* Hnew, realtype *hnew);

/* Combined slow step/fast tolerance multirate controller function.
   This is called following a slow multirate time step with size 'H'
   and fast/slow relative tolerance ratio 'tolfac', and error factors
   'DSM' and 'dsm' (slow and fast, resp.).  The controller should
   estimate slow stepsize 'Hnew' and updated relative tolerance ratio
   'tolfacnew', so that the ensuing step will have 'DSM' and 'dsm'
   values JUST BELOW 1 with minimal computational effort. */
SUNDIALS_EXPORT
int SUNControl_EstimateStepTol(SUNControl C, realtype H,
                               realtype tolfac, realtype DSM,
                               realtype dsm, realtype *Hnew,
                               realtype* tolfacnew);

/* Function to reset the controller to its initial state, e.g., if
   it stores a small number of previous dsm or step size values. */
SUNDIALS_EXPORT
int SUNControl_Reset(SUNControl C);

/* Function to set the controller parameters to their default values. */
SUNDIALS_EXPORT
int SUNControl_SetDefaults(SUNControl C);

/* Function to write all controller parameters to the indicated file
   pointer. */
SUNDIALS_EXPORT
int SUNControl_Write(SUNControl C, FILE* fptr);

/* Function to set the asymptotic order of accuracy for the method. */
SUNDIALS_EXPORT
int SUNControl_SetMethodOrder(SUNControl C, int q);

/* Function to set the asymptotic order of accuracy for the embedding. */
SUNDIALS_EXPORT
int SUNControl_SetEmbeddingOrder(SUNControl C, int p);

/* Function to set an error bias factor to use for scaling the local error
   'dsm' factors above. */
SUNDIALS_EXPORT
int SUNControl_SetErrorBias(SUNControl C, realtype bias);

/* Function to notify the controller of a successful time step with size
   h and local error factor dsm, indicating that the step size or local
   error factor can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNControl_Update(SUNControl C, realtype h, realtype dsm);

/* Function to notify the controller of a successful multirate time step
   with sizes H and h, and local error factors DSM and dsm, indicating that
   the step sizes or local error factors can be saved for subsequent
   controller functions. */
SUNDIALS_EXPORT
int SUNControl_UpdateMRIH(SUNControl C, realtype H, realtype h,
                          realtype DSM, realtype dsm);

/* Function to notify the controller of a successful multirate time step
   with size H and fast tolerance factor tolfac, and local error factors
   DSM and dsm, indicating that the step size, tolerance factor, or local
   error factors can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNControl_UpdateMRITol(SUNControl C, realtype H, realtype tolfac,
                            realtype DSM, realtype dsm);

/* Function to return the memory requirements of the controller object. */
SUNDIALS_EXPORT
int SUNControl_Space(SUNControl C, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNControl error codes
 * ----------------------------------------------------------------- */

#define SUNCONTROL_SUCCESS           0     /* function successfull        */
#define SUNCONTROL_ILL_INPUT         -1001 /* illegal function input      */
#define SUNCONTROL_MEM_FAIL          -1002 /* failed memory access/alloc  */
#define SUNCONTROL_USER_FCN_FAIL     -1003 /* user-supplied fcn failure */
#define SUNCONTROL_OPERATION_FAIL    -1004 /* catchall failure code       */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_CONTROLLER_H */
