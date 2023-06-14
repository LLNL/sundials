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
 *    H       - controls a single-rate step size
 *    HQ      - controls a single-rate step size and method order
 *    MRI_H   - controls two multirate step sizes
 *    MRI_TOL - controls slow and fast relative tolerances
 * ----------------------------------------------------------------- */

typedef enum
{
  SUNDIALS_CONTROL_H,
  SUNDIALS_CONTROL_HQ,
  SUNDIALS_CONTROL_MRI_H,
  SUNDIALS_CONTROL_MRI_TOL
} SUNControl_ID;

/* -----------------------------------------------------------------
 * Generic definition of SUNControl
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNControl_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNControl_Ops* SUNControl_Ops;

/* Forward reference for pointer to SUNControl object */
typedef _SUNDIALS_STRUCT_ _generic_SUNControl* SUNControl;

/* Structure containing function pointers to controller operations  */
struct _generic_SUNControl_Ops
{
  /* REQUIRED of all controller implementations. */
  SUNControl_ID (*getid)(SUNControl C);
  void (*destroy)(SUNControl C);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_H type. */
  int (*estimatestep)(SUNControl C, realtype h,
                      realtype dsm, realtype* hnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_HQ type. */
  int (*estimatestepandorder)(SUNControl C, realtype h, int q,
                              realtype dsm, realtype* hnew, int *qnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_MRI_H type. */
  int (*estimatemristeps)(SUNControl C, realtype H, realtype h,
                          realtype DSM, realtype* Hnew, realtype *hnew);

  /* REQUIRED for controllers of SUNDIALS_CONTROL_MRI_TOL type. */
  int (*estimatesteptol)(SUNControl C, realtype H, realtype tolfac,
                         realtype DSM, realtype *Hnew,
                         realtype* tolfacnew);

  /* OPTIONAL for all SUNControl implementations. */
  int (*reset)(SUNControl C);
  int (*setdefaults)(SUNControl C);
  int (*write)(SUNControl C, FILE* fptr);
  int (*setmethodorder)(SUNControl C, int q);
  int (*setembeddingorder)(SUNControl C, int p);
  int (*setsafetyfactor)(SUNControl C, realtype safety);
  int (*seterrorbias)(SUNControl C, realtype bias);
  int (*update)(SUNControl C, realtype h, realtype dsm);
  int (*space)(SUNControl C, long int *lenrw, long int *leniw);
#ifdef __cplusplus
  _generic_SUNControl_Ops() = default;
#endif

};

/* A SUNControl is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct _generic_SUNControl
{
  void* content;
  SUNControl_Ops ops;
  SUNContext sunctx;
#ifdef __cplusplus
  _generic_SUNControl() = default;
#endif
};

/* -----------------------------------------------------------------
 * Functions exported by SUNControl module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNHeuristics data structure. */
SUNDIALS_EXPORT
SUNControl SUNControlNewEmpty(SUNContext sunctx);

/* Function to free an empty SUNHeuristics data structure. */
SUNDIALS_EXPORT
void SUNControlFreeEmpty(SUNControl C);

/* Function to report the ID of a SUNHeuristics object. */
SUNDIALS_EXPORT
SUNControl_ID SUNControlGetID(SUNControl C);

/* Function to deallocate a SUNHeuristics object. */
SUNDIALS_EXPORT
void SUNControlDestroy(SUNControl C);

/* Main step size controller function.  This is called following
   a time step with size 'h' and local error factor 'dsm', and the
   controller should estimate 'hnew' so that the ensuing step
   will have 'dsm' value JUST BELOW 1.

   Any return value other than SUNCONTROL_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNControlEstimateStep(SUNControl C, realtype h, realtype dsm,
                           realtype* hnew);

/* Combined step size + order controller function.  This is called
   following a time step with size 'h' and order 'q' that has local
   error factor 'dsm'.  The controller should estimate 'hnew' and
   'qnew' so that the ensuing step will have 'dsm' value JUST BELOW 1
   with minimal computational effort. */
SUNDIALS_EXPORT
int SUNControlEstimateStepAndOrder(SUNControl C, realtype h, int q,
                                   realtype dsm, realtype* hnew,
                                   int *qnew);

/* Combined slow/fast multirate step size controller function.  This
   is called following a slow multirate time step with sizes 'H' and
   'h' (slow and fast, resp.), and where the slow step has local error
   factor 'dsm' (the MRI controller should request the fast error
   factor directly from the fast integrator).  The controller should
   estimate slow and fast steps 'Hnew' and 'hnew', resp., so that the
   ensuing step will have 'dsm' value JUST BELOW 1 with minimal
   computational effort. */
SUNDIALS_EXPORT
int SUNControlEstimateMRISteps(SUNControl C, realtype H, realtype h,
                               realtype dsm, realtype* Hnew,
                               realtype *hnew);

/* Combined slow step/fast tolerance multirate controller function.
   This is called following a slow multirate time step with size 'H'
   and fast/slow relative tolerance ratio 'tolfac', and where the
   slow step has local error factor 'dsm' (the MRI controller should
   request the fast error factor directly from the fast integrator).
   The controller should estimate slow stepsize 'Hnew' and updated
   relative tolerance ratio 'tolfacnew', so that the ensuing step
   will have 'dsm' value JUST BELOW 1 with minimal computational
   effort. */
SUNDIALS_EXPORT
int SUNControlEstimateStepTol(SUNControl C, realtype H,
                              realtype tolfac, realtype dsm,
                              realtype *Hnew, realtype* tolfacnew);

/* Function to reset the controller to its initial state, e.g., if
   it stores a small number of previous dsm or step size values. */
SUNDIALS_EXPORT
int SUNControlReset(SUNControl C);

/* Function to set the controller parameters to their default values. */
SUNDIALS_EXPORT
int SUNControlSetDefaults(SUNControl C);

/* Function to write all controller parameters to the indicated file
   pointer. */
SUNDIALS_EXPORT
int SUNControlWrite(SUNControl C, FILE* fptr);

/* Function to set the asymptotic order of accuracy for the method. */
SUNDIALS_EXPORT
int SUNControlSetMethodOrder(SUNControl C, int q);

/* Function to set the asymptotic order of accuracy for the embedding. */
SUNDIALS_EXPORT
int SUNControlSetEmbeddingOrder(SUNControl C, int p);

/* Function to set a step size safety factor. */
SUNDIALS_EXPORT
int SUNControlSetSafetyFactor(SUNControl C, realtype safety);

/* Function to set an error bias factor to use for scaling the local error
   'dsm' factors above. */
SUNDIALS_EXPORT
int SUNControlSetErrorBias(SUNControl C, realtype bias);

/* Function to notify the controller of a successful a time step with size
   h and local error factor dsm, indicating the the step size or local error
   factor can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
int SUNControlUpdate(SUNControl C, realtype h, realtype dsm);

/* Function to return the memory requirements of the controller object. */
SUNDIALS_EXPORT
int SUNControlSpace(SUNControl C, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNControl error codes
 * ----------------------------------------------------------------- */

#define SUNCONTROL_SUCCESS           0     /* function successfull        */
#define SUNCONTROL_ILL_INPUT         -1001 /* illegal function input      */
#define SUNCONTROL_MEM_FAIL          -1002 /* failed memory access/alloc  */
#define SUNCONTROL_USER_FCN_FAIL     -1103 /* user-supplied fcn failure */
#define SUNCONTROL_OPERATION_FAIL    -1004 /* catchall failure code       */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_CONTROLLER_H */
