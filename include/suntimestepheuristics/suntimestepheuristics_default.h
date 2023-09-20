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
 * This is the header file for the SUNTimestepHeuristics_Default
 * module.
 * -----------------------------------------------------------------*/

#ifndef _SUNTIMESTEPHEURISTICS_DEFAULT_H
#define _SUNTIMESTEPHEURISTICS_DEFAULT_H

#include <stdio.h>
#include <sundials/sundials_timestepheuristics.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------
 * Default implementation of SUNTimestepHeuristics
 * ----------------------------------------------- */

struct SUNTimestepHeuristicsContent_Default_ {
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

typedef struct SUNTimestepHeuristicsContent_Default_ *SUNTimestepHeuristicsContent_Default;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNTimestepHeuristics SUNTimestepHeuristics_Default(SUNContext sunctx);
SUNDIALS_EXPORT
SUNTimestepHeuristics_ID SUNTimestepHeuristics_GetID_Default(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ConstrainStep_Default(SUNTimestepHeuristics H,
                                                realtype hcur,
                                                realtype hnew, realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ETestFail_Default(SUNTimestepHeuristics H, realtype hcur,
                                            realtype hnew, int nef,
                                            realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_ConvFail_Default(SUNTimestepHeuristics H, realtype hcur,
                                           realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_BoundReduction_Default(SUNTimestepHeuristics H,
                                                 realtype hcur,
                                                 realtype hnew, realtype* hconstr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_BoundFirstStep_Default(SUNTimestepHeuristics H, realtype h0,
                                                 realtype* h0constr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Reset_Default(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Update_Default(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetDefaults_Default(SUNTimestepHeuristics H);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Write_Default(SUNTimestepHeuristics H, FILE* fptr);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxStep_Default(SUNTimestepHeuristics H, realtype hmax);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMinStep_Default(SUNTimestepHeuristics H, realtype hmin);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetExpStabFn_Default(SUNTimestepHeuristics H,
                                               SUNExpStabFn EStab, void* estab_data);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetCFLFraction_Default(SUNTimestepHeuristics H,
                                                 realtype cfl_frac);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetSafetyFactor_Default(SUNTimestepHeuristics H,
                                                  realtype safety);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxGrowth_Default(SUNTimestepHeuristics H,
                                               realtype mx_growth);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMinReduction_Default(SUNTimestepHeuristics H,
                                                  realtype eta_min);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetFixedStepBounds_Default(SUNTimestepHeuristics H,
                                                     realtype lb, realtype ub);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxFirstGrowth_Default(SUNTimestepHeuristics H,
                                                    realtype etamx1);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxEFailGrowth_Default(SUNTimestepHeuristics H,
                                                    realtype etamxf);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetSmallNumEFails_Default(SUNTimestepHeuristics H,
                                                    int small_nef);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_SetMaxCFailGrowth_Default(SUNTimestepHeuristics H,
                                                    realtype etacf);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_GetNumExpSteps_Default(SUNTimestepHeuristics H,
                                                 long int* expsteps);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_GetNumAccSteps_Default(SUNTimestepHeuristics H,
                                                 long int* accsteps);
SUNDIALS_EXPORT
int SUNTimestepHeuristics_Space_Default(SUNTimestepHeuristics H, long int *lenrw,
                                        long int *leniw);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNTIMESTEPHEURISTICS_DEFAULT_H */
