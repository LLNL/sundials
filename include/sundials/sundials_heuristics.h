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
 * SUNDIALS heuristic adaptivity class. These objects constrain step
 * sizes based on solver heuristics.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_HEURISTICS_H
#define _SUNDIALS_HEURISTICS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include "sundials/sundials_types.h"
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * SUNHeuristics types (currently, only "H" is implemented; others
 * are planned):
 *    STD     - standard heuristic stepsize constraints
 *    NULL    - placeholder object that applies no constraints
 * ----------------------------------------------------------------- */

typedef enum
{
  SUNDIALS_HEURISTICS_STD,
  SUNDIALS_HEURISTICS_NULL
} SUNHeuristics_ID;

/* -----------------------------------------------------------------
 * User-supplied function types
 * ----------------------------------------------------------------- */

typedef int (*SUNExpStabFn)(realtype *hstab, void *user_data);

/* -----------------------------------------------------------------
 * Generic definition of SUNHeuristics
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNHeuristics_Ops object */
typedef _SUNDIALS_STRUCT_ generic_SUNHeuristics_Ops_* SUNHeuristics_Ops;

/* Forward reference for pointer to SUNHeuristics object */
typedef _SUNDIALS_STRUCT_ generic_SUNHeuristics_* SUNHeuristics;

/* Structure containing function pointers to heuristics operations  */
struct generic_SUNHeuristics_Ops_
{
  SUNHeuristics_ID (*getid)(SUNHeuristics H);
  int (*destroy)(SUNHeuristics H);
  int (*constrainstep)(SUNHeuristics H, realtype hcur,
                       realtype hnew, realtype *hconstr);
  int (*etestfail)(SUNHeuristics H, realtype hcur,
                   realtype hnew, int nef, realtype *hconstr);
  int (*convfail)(SUNHeuristics H, realtype hcur, realtype *hconstr);
  int (*boundreduction)(SUNHeuristics H, realtype hcur,
                        realtype hnew, realtype *hconstr);
  int (*boundfirststep)(SUNHeuristics H, realtype h0, realtype *h0constr);
  int (*reset)(SUNHeuristics H);
  int (*update)(SUNHeuristics H);
  int (*setdefaults)(SUNHeuristics H);
  int (*write)(SUNHeuristics H, FILE* fptr);
  int (*setmaxstep)(SUNHeuristics H, realtype hmax);
  int (*setminstep)(SUNHeuristics H, realtype hmin);
  int (*setexpstabfn)(SUNHeuristics H, SUNExpStabFn EStab, void* estab_data);
  int (*setcflfraction)(SUNHeuristics H, realtype cfl_frac);
  int (*setsafetyfactor)(SUNHeuristics H, realtype safety);
  int (*setmaxgrowth)(SUNHeuristics H, realtype mx_growth);
  int (*setminreduction)(SUNHeuristics H, realtype eta_min);
  int (*setfixedstepbounds)(SUNHeuristics H, realtype lb, realtype ub);
  int (*setmaxfirstgrowth)(SUNHeuristics H, realtype etamx1);
  int (*setmaxefailgrowth)(SUNHeuristics H, realtype etamxf);
  int (*setsmallnumefails)(SUNHeuristics H, int small_nef);
  int (*setmaxcfailgrowth)(SUNHeuristics H, realtype etacf);
  int (*getnumexpsteps)(SUNHeuristics H, long int* expsteps);
  int (*getnumaccsteps)(SUNHeuristics H, long int* accsteps);
  int (*space)(SUNHeuristics H, long int *lenrw, long int *leniw);
#ifdef __cplusplus
  generic_SUNHeuristics_Ops_() = default;
#endif

};

/* A SUNHeuristics is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct generic_SUNHeuristics_
{
  void* content;
  SUNHeuristics_Ops ops;
  SUNContext sunctx;
#ifdef __cplusplus
  generic_SUNHeuristics_() = default;
#endif
};

/* -----------------------------------------------------------------
 * Functions exported by SUNHeuristics module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNHeuristics data structure. */
SUNDIALS_EXPORT
SUNHeuristics SUNHeuristics_NewEmpty(SUNContext sunctx);

/* Function to report the ID of a SUNHeuristics object. */
SUNDIALS_EXPORT
SUNHeuristics_ID SUNHeuristics_GetID(SUNHeuristics H);

/* Function to deallocate a SUNHeuristics object.

   Any return value other than SUNHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNHeuristics_Destroy(SUNHeuristics H);

/* Main constraint-application function.  This will attempt to
   change the step hcur to hnew, applying any heuristic bounds
   on the step size adjustments.

   Any return value other than SUNHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNHeuristics_ConstrainStep(SUNHeuristics H, realtype hcur,
                                realtype hnew, realtype *hconstr);

/* Function to apply constraints following a step with unacceptable
   temporal error.  Here, 'hnew' has already been predicted by a
   temporal error controller, and 'nef' is the integrator-provided
   counter of how many temporal error test failures have occurred
   on this time step.

   Any return value other than SUNHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNHeuristics_ETestFail(SUNHeuristics H, realtype hcur,
                            realtype hnew,  int nef, realtype *hconstr);

/* This ensures that a step size reduction is within user-prescribed
   bounds.  If a reduction is requested but no reduction from hcur
   is possible, this returns SUNHEURISTICS_CANNOT_DECREASE, otherwise
   it returns SUNHEURISTICS_SUCCESS. */
SUNDIALS_EXPORT
int SUNHeuristics_BoundReduction(SUNHeuristics H, realtype hcur,
                                 realtype hnew, realtype *hconstr);

/* This bounds the initial step by user-provided min/max step values. */
SUNDIALS_EXPORT
int SUNHeuristics_BoundFirstStep(SUNHeuristics H, realtype h0,
                                 realtype *h0constr);

/* Function to apply constraints following a step with an algebraic
   solver convergence failure. */
SUNDIALS_EXPORT
int SUNHeuristics_ConvFail(SUNHeuristics H, realtype hcur,
                           realtype *hconstr);

/* Function to reset the heuristics to its initial state: zeros
   the "accuracy" and "stability" counters (if present), and resets
   any stepsize growth factor bounds. */
SUNDIALS_EXPORT
int SUNHeuristics_Reset(SUNHeuristics H);

/* Function to notify the heuristics that a time step has
   succeeded, indicating e.g. that the stepsize growh factor should
   should be set to its "default" state. */
SUNDIALS_EXPORT
int SUNHeuristics_Update(SUNHeuristics H);

/* Function to set the heuristics parameters to their default values. */
SUNDIALS_EXPORT
int SUNHeuristics_SetDefaults(SUNHeuristics H);

/* Function to write all heuristics parameters to the indicated
   file pointer. */
SUNDIALS_EXPORT
int SUNHeuristics_Write(SUNHeuristics H, FILE* fptr);

/* Function to set the maximum absolute step size allowed
   (hmax <=0 implies infinite). */
SUNDIALS_EXPORT
int SUNHeuristics_SetMaxStep(SUNHeuristics H, realtype hmax);

/* Function to set the minimum absolute step size allowed
   (hmin <= 0 implies no minimum). */
SUNDIALS_EXPORT
int SUNHeuristics_SetMinStep(SUNHeuristics H, realtype hmin);

/* Function to provide a user-supplied function for the
   maximum stable step size (EStab == NULL disables). */
SUNDIALS_EXPORT
int SUNHeuristics_SetExpStabFn(SUNHeuristics H, SUNExpStabFn EStab,
                               void* estab_data);

/* Function to set a CFL safety factor that should be applied to
   the outputs of EStab, above.  Allowable values must be within
   the open interval (0,1).  A non-positive input implies a reset
   to the default value. */
SUNDIALS_EXPORT
int SUNHeuristics_SetCFLFraction(SUNHeuristics H, realtype cfl_frac);

/* Function to set a step size safety factor that should be used
   to constrain the error-controller-recommended step size. */
SUNDIALS_EXPORT
int SUNHeuristics_SetSafetyFactor(SUNHeuristics C, realtype safety);

/* Function to set maximum stepsize growth factor for general steps.
   Note: the first step uses a separate maximum growth factor.
   Allowable values must be >1.  Any illegal value implies a reset
   to the default.*/
SUNDIALS_EXPORT
int SUNHeuristics_SetMaxGrowth(SUNHeuristics H, realtype mx_growth);

/* Function to specify the maximum stepsize growth factor for the
   first internal time step.  Legal values are greater than 1.0.
   Illegal values imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNHeuristics_SetMaxFirstGrowth(SUNHeuristics H, realtype etamx1);

/* Function to specify the step size growth interval within which
   the step size will remain unchanged.  Allowable values must
   enclose the value 1.0.  Any illegal interval implies a reset to
   the default. */
SUNDIALS_EXPORT
int SUNHeuristics_SetFixedStepBounds(SUNHeuristics H, realtype lb,
                                     realtype ub);

/* Function to set a lower bound for the stepsize adjustment factor
   following a temporal error failure. */
SUNDIALS_EXPORT
int SUNHeuristics_SetMinReduction(SUNHeuristics H, realtype eta_min);

/* Function to set an upper bound for the stepsize adjustment factor
   following a temporal error failure (once at least small_nef
   failures have occurred). Legal values are in the interval (0,1].
   Illegal values imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNHeuristics_SetMaxEFailGrowth(SUNHeuristics H, realtype etamxf);

/* Function to specify the step adaptivity constant small_nef (see
   MaxEFailGrowth above).  Legal values are > 0.  Illegal values
   imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNHeuristics_SetSmallNumEFails(SUNHeuristics H, int small_nef);

/* Function to specify an upper bound for the stepsize adjustment
   factor following an algebraic solver convergence failure. Legal
   values are in the interval (0,1].  Illegal values imply a reset
   to the default value. */
SUNDIALS_EXPORT
int SUNHeuristics_SetMaxCFailGrowth(SUNHeuristics H, realtype etacf);

/* Function to return the current number of stability-limited steps. */
SUNDIALS_EXPORT
int SUNHeuristics_GetNumExpSteps(SUNHeuristics H, long int* expsteps);

/* Function to return the current number of accuracy-limited steps. */
SUNDIALS_EXPORT
int SUNHeuristics_GetNumAccSteps(SUNHeuristics H, long int* accsteps);

/* Function to return the memory requirements of the heuristics object. */
SUNDIALS_EXPORT
int SUNHeuristics_Space(SUNHeuristics H, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNHeuristics error codes
 * ----------------------------------------------------------------- */

#define SUNHEURISTICS_SUCCESS         0     /* function successfull      */
#define SUNHEURISTICS_ILL_INPUT       -1101 /* illegal function input    */
#define SUNHEURISTICS_CANNOT_DECREASE -1102 /* no step decrease possible */
#define SUNHEURISTICS_USER_FCN_FAIL   -1103 /* user-supplied fcn failure */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_HEURISTICS_H */
