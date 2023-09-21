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
 * SUNDIALS time step heuristics class. These objects constrain step
 * sizes based on solver heuristics.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_TIMESTEPHEURISTICS_H
#define _SUNDIALS_TIMESTEPHEURISTICS_H

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
 * SUNTimestepHeuristics types (currently, only "H" is implemented; others
 * are planned):
 *    STD     - standard heuristic stepsize constraints
 *    NULL    - placeholder object that applies no constraints
 * ----------------------------------------------------------------- */

typedef enum
{
  SUN_TIMESTEPHEURISTICS_STD,
  SUN_TIMESTEPHEURISTICS_NULL
} SUNTimestepHeuristics_ID;

/* -----------------------------------------------------------------
 * User-supplied function types
 * ----------------------------------------------------------------- */

typedef int (*SUNExpStabFn)(realtype *hstab, void *user_data);

/* -----------------------------------------------------------------
 * Generic definition of SUNTimestepHeuristics
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNTimestepHeuristics_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNTimestepHeuristics_Ops* SUNTimestepHeuristics_Ops;

/* Forward reference for pointer to SUNTimestepHeuristics object */
typedef _SUNDIALS_STRUCT_ _generic_SUNTimestepHeuristics* SUNTimestepHeuristics;

/* Structure containing function pointers to heuristics operations  */
struct _generic_SUNTimestepHeuristics_Ops
{
  SUNTimestepHeuristics_ID (*getid)(SUNTimestepHeuristics H);
  int (*destroy)(SUNTimestepHeuristics H);
  int (*constrainstep)(SUNTimestepHeuristics H, realtype hcur,
                       realtype hnew, realtype *hconstr);
  int (*etestfail)(SUNTimestepHeuristics H, realtype hcur,
                   realtype hnew, int nef, realtype *hconstr);
  int (*convfail)(SUNTimestepHeuristics H, realtype hcur, realtype *hconstr);
  int (*boundreduction)(SUNTimestepHeuristics H, realtype hcur,
                        realtype hnew, realtype *hconstr);
  int (*boundfirststep)(SUNTimestepHeuristics H, realtype h0, realtype *h0constr);
  int (*reset)(SUNTimestepHeuristics H);
  int (*update)(SUNTimestepHeuristics H);
  int (*setdefaults)(SUNTimestepHeuristics H);
  int (*write)(SUNTimestepHeuristics H, FILE* fptr);
  int (*setmaxstep)(SUNTimestepHeuristics H, realtype hmax);
  int (*setminstep)(SUNTimestepHeuristics H, realtype hmin);
  int (*setexpstabfn)(SUNTimestepHeuristics H, SUNExpStabFn EStab, void* estab_data);
  int (*setcflfraction)(SUNTimestepHeuristics H, realtype cfl_frac);
  int (*setsafetyfactor)(SUNTimestepHeuristics H, realtype safety);
  int (*setmaxgrowth)(SUNTimestepHeuristics H, realtype mx_growth);
  int (*setminreduction)(SUNTimestepHeuristics H, realtype eta_min);
  int (*setfixedstepbounds)(SUNTimestepHeuristics H, realtype lb, realtype ub);
  int (*setmaxfirstgrowth)(SUNTimestepHeuristics H, realtype etamx1);
  int (*setmaxefailgrowth)(SUNTimestepHeuristics H, realtype etamxf);
  int (*setsmallnumefails)(SUNTimestepHeuristics H, int small_nef);
  int (*setmaxcfailgrowth)(SUNTimestepHeuristics H, realtype etacf);
  int (*getnumexpsteps)(SUNTimestepHeuristics H, long int* expsteps);
  int (*getnumaccsteps)(SUNTimestepHeuristics H, long int* accsteps);
  int (*space)(SUNTimestepHeuristics H, long int *lenrw, long int *leniw);
};

/* A SUNTimestepHeuristics is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct _generic_SUNTimestepHeuristics
{
  void* content;
  SUNTimestepHeuristics_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNTimestepHeuristics module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNTimestepHeuristics data structure. */
SUNDIALS_EXPORT
SUNTimestepHeuristics SUNTimestepHeuristics_NewEmpty(SUNContext sunctx);

/* Function to report the ID of a SUNTimestepHeuristics object. */
SUNDIALS_EXPORT
SUNTimestepHeuristics_ID SUNTimestepHeuristics_GetID(SUNTimestepHeuristics H);

/* Function to deallocate a SUNTimestepHeuristics object.

   Any return value other than SUNTIMESTEPHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Destroy(SUNTimestepHeuristics H);

/* Main constraint-application function.  This will attempt to
   change the step hcur to hnew, applying any heuristic bounds
   on the step size adjustments.

   Any return value other than SUNTIMESTEPHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ConstrainStep(SUNTimestepHeuristics H, realtype hcur,
                                realtype hnew, realtype *hconstr);

/* Function to apply constraints following a step with unacceptable
   temporal error.  Here, 'hnew' has already been predicted by a
   temporal error controller, and 'nef' is the integrator-provided
   counter of how many temporal error test failures have occurred
   on this time step.

   Any return value other than SUNTIMESTEPHEURISTICS_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ETestFail(SUNTimestepHeuristics H, realtype hcur,
                            realtype hnew,  int nef, realtype *hconstr);

/* This ensures that a step size reduction is within user-prescribed
   bounds.  If a reduction is requested but no reduction from hcur
   is possible, this returns SUNTIMESTEPHEURISTICS_CANNOT_DECREASE, otherwise
   it returns SUNTIMESTEPHEURISTICS_SUCCESS. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_BoundReduction(SUNTimestepHeuristics H, realtype hcur,
                                 realtype hnew, realtype *hconstr);

/* This bounds the initial step by user-provided min/max step values. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_BoundFirstStep(SUNTimestepHeuristics H, realtype h0,
                                 realtype *h0constr);

/* Function to apply constraints following a step with an algebraic
   solver convergence failure. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ConvFail(SUNTimestepHeuristics H, realtype hcur,
                           realtype *hconstr);

/* Function to reset the heuristics to its initial state: zeros
   the "accuracy" and "stability" counters (if present), and resets
   any stepsize growth factor bounds. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Reset(SUNTimestepHeuristics H);

/* Function to notify the heuristics that a time step has
   succeeded, indicating e.g. that the stepsize growh factor should
   should be set to its "default" state. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Update(SUNTimestepHeuristics H);

/* Function to set the heuristics parameters to their default values. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetDefaults(SUNTimestepHeuristics H);

/* Function to write all heuristics parameters to the indicated
   file pointer. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Write(SUNTimestepHeuristics H, FILE* fptr);

/* Function to set the maximum absolute step size allowed
   (hmax <=0 implies infinite). */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxStep(SUNTimestepHeuristics H, realtype hmax);

/* Function to set the minimum absolute step size allowed
   (hmin <= 0 implies no minimum). */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMinStep(SUNTimestepHeuristics H, realtype hmin);

/* Function to provide a user-supplied function for the
   maximum stable step size (EStab == NULL disables). */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetExpStabFn(SUNTimestepHeuristics H, SUNExpStabFn EStab,
                               void* estab_data);

/* Function to set a CFL safety factor that should be applied to
   the outputs of EStab, above.  Allowable values must be within
   the open interval (0,1).  A non-positive input implies a reset
   to the default value. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetCFLFraction(SUNTimestepHeuristics H, realtype cfl_frac);

/* Function to set a step size safety factor that should be used
   to constrain the error-controller-recommended step size. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetSafetyFactor(SUNTimestepHeuristics C, realtype safety);

/* Function to set maximum stepsize growth factor for general steps.
   Note: the first step uses a separate maximum growth factor.
   Allowable values must be >1.  Any illegal value implies a reset
   to the default.*/
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxGrowth(SUNTimestepHeuristics H, realtype mx_growth);

/* Function to specify the maximum stepsize growth factor for the
   first internal time step.  Legal values are greater than 1.0.
   Illegal values imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxFirstGrowth(SUNTimestepHeuristics H, realtype etamx1);

/* Function to specify the step size growth interval within which
   the step size will remain unchanged.  Allowable values must
   enclose the value 1.0.  Any illegal interval implies a reset to
   the default. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetFixedStepBounds(SUNTimestepHeuristics H, realtype lb,
                                     realtype ub);

/* Function to set a lower bound for the stepsize adjustment factor
   following a temporal error failure. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMinReduction(SUNTimestepHeuristics H, realtype eta_min);

/* Function to set an upper bound for the stepsize adjustment factor
   following a temporal error failure (once at least small_nef
   failures have occurred). Legal values are in the interval (0,1].
   Illegal values imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxEFailGrowth(SUNTimestepHeuristics H, realtype etamxf);

/* Function to specify the step adaptivity constant small_nef (see
   MaxEFailGrowth above).  Legal values are > 0.  Illegal values
   imply a reset to the default value. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetSmallNumEFails(SUNTimestepHeuristics H, int small_nef);

/* Function to specify an upper bound for the stepsize adjustment
   factor following an algebraic solver convergence failure. Legal
   values are in the interval (0,1].  Illegal values imply a reset
   to the default value. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxCFailGrowth(SUNTimestepHeuristics H, realtype etacf);

/* Function to return the current number of stability-limited steps. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_GetNumExpSteps(SUNTimestepHeuristics H, long int* expsteps);

/* Function to return the current number of accuracy-limited steps. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_GetNumAccSteps(SUNTimestepHeuristics H, long int* accsteps);

/* Function to return the memory requirements of the heuristics object. */
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Space(SUNTimestepHeuristics H, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNTimestepHeuristics error codes
 * ----------------------------------------------------------------- */

#define SUNTIMESTEPHEURISTICS_SUCCESS         0     /* function successfull      */
#define SUNTIMESTEPHEURISTICS_ILL_INPUT       -1101 /* illegal function input    */
#define SUNTIMESTEPHEURISTICS_CANNOT_DECREASE -1102 /* no step decrease possible */
#define SUNTIMESTEPHEURISTICS_USER_FCN_FAIL   -1103 /* user-supplied fcn failure */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_TIMESTEPHEURISTICS_H */
