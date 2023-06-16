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
 * This is the header file for the SUNHeuristics_Default module.
 * -----------------------------------------------------------------*/

#ifndef _SUNHEURISTICS_DEFAULT_H
#define _SUNHEURISTICS_DEFAULT_H

#include <stdio.h>
#include <sundials/sundials_heuristics.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------
 * Default implementation of SUNHeuristics
 * --------------------------------------- */

struct _SUNHeuristicsContent_Default {
  realtype     hmax_inv;   /* inverse of maximum allowable time step     */
  realtype     hmin;       /* minimum allowable time step                */
  realtype     etamax;     /* eta <= etamax                              */
  realtype     etamx1;     /* max step size change on first step         */
  realtype     etamxf;     /* h reduction factor on multiple error fails */
  realtype     etamin;     /* eta >= etamin on error test fail           */
  int          small_nef;  /* bound to determine 'multiple' above        */
  realtype     etacf;      /* h reduction factor on nonlinear conv fail  */
  SUNExpStabFn expstab;    /* step stability function                    */
  void*        estab_data; /* user pointer passed to expstab             */
  realtype     cfl;        /* cfl safety factor                          */
  realtype     safety;     /* step safety factor                         */
  realtype     growth;     /* maximum step growth safety factor          */
  realtype     lbound;     /* eta lower bound to leave h unchanged       */
  realtype     ubound;     /* eta upper bound to leave h unchanged       */
  long int     nst_acc;    /* num accuracy-limited internal steps        */
  long int     nst_exp;    /* num stability-limited internal steps       */
};

typedef struct _SUNHeuristicsContent_Default *SUNHeuristicsContent_Default;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNHeuristics SUNHeuristicsDefault(SUNContext sunctx);
SUNDIALS_EXPORT
SUNHeuristics_ID SUNHeuristicsGetID_Default(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsConstrainStep_Default(SUNHeuristics H, realtype hcur,
                                       realtype hnew, realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsETestFail_Default(SUNHeuristics H, realtype hcur,
                                   realtype hnew, int nef, realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsConvFail_Default(SUNHeuristics H, realtype hcur,
                                  realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsBoundReduction_Default(SUNHeuristics H, realtype hcur,
                                        realtype hnew, realtype* hconstr);
SUNDIALS_EXPORT
int SUNHeuristicsReset_Default(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsUpdate_Default(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsSetDefaults_Default(SUNHeuristics H);
SUNDIALS_EXPORT
int SUNHeuristicsWrite_Default(SUNHeuristics H, FILE* fptr);
SUNDIALS_EXPORT
int SUNHeuristicsSetMaxStep_Default(SUNHeuristics H, realtype hmax);
SUNDIALS_EXPORT
int SUNHeuristicsSetMinStep_Default(SUNHeuristics H, realtype hmin);
SUNDIALS_EXPORT
int SUNHeuristicsSetExpStabFn_Default(SUNHeuristics H, SUNExpStabFn EStab,
                                      void* estab_data);
SUNDIALS_EXPORT
int SUNHeuristicsSetCFLFraction_Default(SUNHeuristics H, realtype cfl_frac);
SUNDIALS_EXPORT
int SUNHeuristicsSetSafetyFactor_Default(SUNHeuristics H, realtype safety);
SUNDIALS_EXPORT
int SUNHeuristicsSetMaxGrowth_Default(SUNHeuristics H, realtype mx_growth);
SUNDIALS_EXPORT
int SUNHeuristicsSetMinReduction_Default(SUNHeuristics H, realtype eta_min);
SUNDIALS_EXPORT
int SUNHeuristicsSetFixedStepBounds_Default(SUNHeuristics H, realtype lb, realtype ub);
SUNDIALS_EXPORT
int SUNHeuristicsSetMaxFirstGrowth_Default(SUNHeuristics H, realtype etamx1);
SUNDIALS_EXPORT
int SUNHeuristicsSetMaxEFailGrowth_Default(SUNHeuristics H, realtype etamxf);
SUNDIALS_EXPORT
int SUNHeuristicsSetSmallNumEFails_Default(SUNHeuristics H, int small_nef);
SUNDIALS_EXPORT
int SUNHeuristicsSetMaxCFailGrowth_Default(SUNHeuristics H, realtype etacf);
SUNDIALS_EXPORT
int SUNHeuristicsGetNumExpSteps_Default(SUNHeuristics H, long int* expsteps);
SUNDIALS_EXPORT
int SUNHeuristicsGetNumAccSteps_Default(SUNHeuristics H, long int* accsteps);
SUNDIALS_EXPORT
int SUNHeuristicsSpace_Default(SUNHeuristics H, long int *lenrw, long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNHEURISTICS_DEFAULT_H */
