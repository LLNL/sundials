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

typedef int (*SUNExpStabFn)(N_Vector y, realtype t,
                            realtype *hstab, void *user_data);

/* -----------------------------------------------------------------
 * Generic definition of SUNHeuristics
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNHeuristics_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNHeuristics_Ops* SUNHeuristics_Ops;

/* Forward reference for pointer to SUNHeuristics object */
typedef _SUNDIALS_STRUCT_ _generic_SUNHeuristics* SUNHeuristics;

/* Structure containing function pointers to heuristics operations  */
struct _generic_SUNHeuristics_Ops
{
  /* REQUIRED of all heuristics implementations. */
  SUNHeuristics_ID (*getid)(SUNHeuristics H);
  void (*destroy)(SUNHeuristics H);

  /* OPTIONAL for all SUNHeuristics implementations. */
  int (*constrainstep)(SUNHeuristics H, realtype hcur,
                       realtype hnew, booleantype efail,
                       booleantype cfail, realtype *hconstr);
  int (*reset)(SUNHeuristics H);
  int (*setdefaults)(SUNHeuristics H);
  int (*write)(SUNHeuristics H, FILE* fptr);
  int (*setexpstabfn)(SUNHeuristics H, SUNExpStabFn EStab, void* estab_data);
  int (*setcflfraction)(SUNHeuristics H, realtype cfl_frac);
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
  _generic_SUNHeuristics_Ops() = default;
#endif

};

/* A SUNHeuristics is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct _generic_SUNHeuristics
{
  void* content;
  SUNHeuristics_Ops ops;
  SUNContext sunctx;
#ifdef __cplusplus
  _generic_SUNHeuristics() = default;
#endif
};

/* -----------------------------------------------------------------
 * Functions exported by SUNHeuristics module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT SUNHeuristics SUNHeuristicsNewEmpty(SUNContext sunctx);
SUNDIALS_EXPORT void SUNHeuristicsFreeEmpty(SUNHeuristics H);
SUNDIALS_EXPORT SUNHeuristics_ID SUNHeuristicsGetID(SUNHeuristics H);
SUNDIALS_EXPORT int SUNHeuristicsConstrainStep(
                       SUNHeuristics H, realtype hcur,
                       realtype hnew, booleantype efail,
                       booleantype cfail, realtype *hconstr);
SUNDIALS_EXPORT int SUNHeuristicsReset(SUNHeuristics H);
SUNDIALS_EXPORT int SUNHeuristicsSetDefaults(SUNHeuristics H);
SUNDIALS_EXPORT int SUNHeuristicsWrite(SUNHeuristics H, FILE* fptr);

SUNDIALS_EXPORT int SUNHeuristicsSetExpStabFn(
                       SUNHeuristics H, SUNExpStabFn EStab, void* estab_data);
SUNDIALS_EXPORT int SUNHeuristicsSetCFLFraction(
                       SUNHeuristics H, realtype cfl_frac);
SUNDIALS_EXPORT int SUNHeuristicsSetMaxGrowth(
                       SUNHeuristics H, realtype mx_growth);
SUNDIALS_EXPORT int SUNHeuristicsSetMinReduction(
                       SUNHeuristics H, realtype eta_min);
SUNDIALS_EXPORT int SUNHeuristicsSetFixedStepBounds(
                       SUNHeuristics H, realtype lb, realtype ub);
SUNDIALS_EXPORT int SUNHeuristicsSetMaxFirstGrowth(
                       SUNHeuristics H, realtype etamx1);
SUNDIALS_EXPORT int SUNHeuristicsSetMaxEFailGrowth(
                       SUNHeuristics H, realtype etamxf);
SUNDIALS_EXPORT int SUNHeuristicsSetSmallNumEFails(
                       SUNHeuristics H, int small_nef);
SUNDIALS_EXPORT int SUNHeuristicsSetMaxCFailGrowth(
                       SUNHeuristics H, realtype etacf);
SUNDIALS_EXPORT int SUNHeuristicsGetNumExpSteps(
                       SUNHeuristics H, long int* expsteps);
SUNDIALS_EXPORT int SUNHeuristicsGetNumAccSteps(
                       SUNHeuristics H, long int* accsteps);
SUNDIALS_EXPORT int SUNHeuristicsSpace(
                       SUNHeuristics H, long int *lenrw, long int *leniw);


/* -----------------------------------------------------------------
 * SUNHeuristics error codes
 * ----------------------------------------------------------------- */

#define SUNHEURISTICS_SUCCESS           0     /* function successfull        */
#define SUNHEURISTICS_ILL_INPUT         -1101 /* illegal function input      */
#define SUNHEURISTICS_MEM_FAIL          -1102 /* failed memory access/alloc  */
#define SUNHEURISTICS_OPERATION_FAIL    -1103 /* catchall failure code       */

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_HEURISTICS_H */
