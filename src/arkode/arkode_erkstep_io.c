/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKODE ERKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_erkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <suncontrol/suncontrol_pi.h>
#include <sunheuristics/sunheuristics_default.h>


/*===============================================================
  ERKStep Optional input functions (wrappers for generic ARKODE
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int ERKStepSetDenseOrder(void *arkode_mem, int dord) {
  return(ERKStepSetInterpolantDegree(arkode_mem, dord)); }
int ERKStepSetInterpolantDegree(void *arkode_mem, int degree) {
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return(arkSetInterpolantDegree(arkode_mem, degree)); }
int ERKStepSetInterpolantType(void *arkode_mem, int itype) {
  return(arkSetInterpolantType(arkode_mem, itype)); }
int ERKStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data) {
  return(arkSetErrHandlerFn(arkode_mem, ehfun, eh_data)); }
int ERKStepSetErrFile(void *arkode_mem, FILE *errfp) {
  return(arkSetErrFile(arkode_mem, errfp)); }
int ERKStepSetUserData(void *arkode_mem, void *user_data) {
  return(arkSetUserData(arkode_mem, user_data)); }
int ERKStepSetDiagnostics(void *arkode_mem, FILE *diagfp) {
  return(arkSetDiagnostics(arkode_mem, diagfp)); }
int ERKStepSetMaxNumSteps(void *arkode_mem, long int mxsteps) {
  return(arkSetMaxNumSteps(arkode_mem, mxsteps)); }
int ERKStepSetMaxHnilWarns(void *arkode_mem, int mxhnil) {
  return(arkSetMaxHnilWarns(arkode_mem, mxhnil)); }
int ERKStepSetInitStep(void *arkode_mem, realtype hin) {
  return(arkSetInitStep(arkode_mem, hin)); }
int ERKStepSetStopTime(void *arkode_mem, realtype tstop) {
  return(arkSetStopTime(arkode_mem, tstop)); }
int ERKStepClearStopTime(void *arkode_mem) {
  return(arkClearStopTime(arkode_mem)); }
int ERKStepSetRootDirection(void *arkode_mem, int *rootdir) {
  return(arkSetRootDirection(arkode_mem, rootdir)); }
int ERKStepSetNoInactiveRootWarn(void *arkode_mem) {
  return(arkSetNoInactiveRootWarn(arkode_mem)); }
int ERKStepSetConstraints(void *arkode_mem, N_Vector constraints) {
  return(arkSetConstraints(arkode_mem, constraints)); }
int ERKStepSetMaxNumConstrFails(void *arkode_mem, int maxfails) {
  return(arkSetMaxNumConstrFails(arkode_mem, maxfails)); }
int ERKStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessFn ProcessStep) {
  return(arkSetPostprocessStepFn(arkode_mem, ProcessStep)); }
int ERKStepSetPostprocessStageFn(void *arkode_mem,
                                ARKPostProcessFn ProcessStage) {
  return(arkSetPostprocessStageFn(arkode_mem, ProcessStage)); }
int ERKStepSetMaxErrTestFails(void *arkode_mem, int maxnef) {
  return(arkSetMaxErrTestFails(arkode_mem, maxnef)); }
int ERKStepSetFixedStep(void *arkode_mem, realtype hfixed) {
  return(arkSetFixedStep(arkode_mem, hfixed)); }
int ERKStepSetController(void *arkode_mem, SUNControl C) {
  return(arkSetController(arkode_mem, C)); }
int ERKStepSetHeuristics(void *arkode_mem, SUNHeuristics H) {
  return(arkSetHeuristics(arkode_mem, H)); }


/*===============================================================
  ERKStep Optional output functions (wrappers for generic ARKODE
  utility routines). All are documented in arkode_io.c.
  ===============================================================*/

int ERKStepGetNumStepAttempts(void *arkode_mem, long int *nstep_attempts) {
  return(arkGetNumStepAttempts(arkode_mem, nstep_attempts)); }
int ERKStepGetNumSteps(void *arkode_mem, long int *nsteps) {
  return(arkGetNumSteps(arkode_mem, nsteps)); }
int ERKStepGetActualInitStep(void *arkode_mem, realtype *hinused) {
  return(arkGetActualInitStep(arkode_mem, hinused)); }
int ERKStepGetLastStep(void *arkode_mem, realtype *hlast) {
  return(arkGetLastStep(arkode_mem, hlast)); }
int ERKStepGetCurrentStep(void *arkode_mem, realtype *hcur) {
  return(arkGetCurrentStep(arkode_mem, hcur)); }
int ERKStepGetCurrentTime(void *arkode_mem, realtype *tcur) {
  return(arkGetCurrentTime(arkode_mem, tcur)); }
int ERKStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact) {
  return(arkGetTolScaleFactor(arkode_mem, tolsfact)); }
int ERKStepGetErrWeights(void *arkode_mem, N_Vector eweight) {
  return(arkGetErrWeights(arkode_mem, eweight)); }
int ERKStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw) {
  return(arkGetWorkSpace(arkode_mem, lenrw, leniw)); }
int ERKStepGetNumGEvals(void *arkode_mem, long int *ngevals) {
  return(arkGetNumGEvals(arkode_mem, ngevals)); }
int ERKStepGetRootInfo(void *arkode_mem, int *rootsfound) {
  return(arkGetRootInfo(arkode_mem, rootsfound)); }
int ERKStepGetStepStats(void *arkode_mem, long int *nsteps,
                        realtype *hinused, realtype *hlast,
                        realtype *hcur, realtype *tcur) {
  return(arkGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur)); }
int ERKStepGetNumConstrFails(void *arkode_mem, long int *nconstrfails) {
  return(arkGetNumConstrFails(arkode_mem, nconstrfails)); }
int ERKStepGetNumErrTestFails(void *arkode_mem, long int *netfails) {
  return(arkGetNumErrTestFails(arkode_mem, netfails)); }
int ERKStepGetUserData(void *arkode_mem, void** user_data) {
  return(arkGetUserData(arkode_mem, user_data)); }
char *ERKStepGetReturnFlagName(long int flag) {
  return(arkGetReturnFlagName(flag)); }


/*===============================================================
  DEPRECATED ERKStep optional input/output functions
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetMinStep: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMinStep(void *arkode_mem, realtype hmin)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMinStep(ark_mem->hconstraints, hmin);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetMaxStep: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMaxStep(void *arkode_mem, realtype hmax)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMaxStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMaxStep(ark_mem->hconstraints, hmax);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetCFLFraction: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetCFLFraction", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetCFLFraction(ark_mem->hconstraints, cfl_frac);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetSafetyFactor: user should create/attach a
  SUNControl object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetSafetyFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNControlSetSafetyFactor(ark_mem->hcontroller, safety);
  if (retval != SUNCONTROL_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetErrorBias: user should create/attach a
  SUNControl object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetErrorBias", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNControlSetErrorBias(ark_mem->hcontroller, bias);
  if (retval != SUNCONTROL_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetMaxGrowth: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMaxGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMaxGrowth(ark_mem->hconstraints, mx_growth);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetMinReduction: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMinReduction(void *arkode_mem, realtype eta_min)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMinReduction", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMinReduction(ark_mem->hconstraints, eta_min);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetFixedStepBounds: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetFixedStepBounds", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetFixedStepBounds(ark_mem->hconstraints, lb, ub);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetAdaptivityMethod: user should create/attach a
  specific SUNControl object.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault,
                               int pq, realtype adapt_params[3]) {
  return(arkSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params)); }

/*---------------------------------------------------------------
  ERKStepSetAdaptivityFn: user should create/attach a custom
  SUNControl object.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data) {
  return(arkSetAdaptivityFn(arkode_mem, hfun, h_data)); }

/*---------------------------------------------------------------
  ERKStepSetMaxFirstGrowth: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMaxFirstGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMaxFirstGrowth(ark_mem->hconstraints, etamx1);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetMaxEFailGrowth: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetMaxEFailGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetMaxEFailGrowth(ark_mem->hconstraints, etamxf);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetSmallNumEFails: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepSetSmallNumEFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  retval = SUNHeuristicsSetSmallNumEFails(ark_mem->hconstraints, small_nef);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_ILL_INPUT); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetStabilityFn: user should create/attach a
  SUNHeuristics object, and set this directly therein.
  ---------------------------------------------------------------*/
int ERKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data) {
  return(arkSetStabilityFn(arkode_mem, EStab, estab_data)); }




/*===============================================================
  ERKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetDefaults:

  Resets all ERKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int ERKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;
  long int lenrw, leniw;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set default ARKODE infrastructure parameters */
  retval = arkSetDefaults(arkode_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSetDefaults",
                    "Error setting ARKODE infrastructure defaults");
    return(retval);
  }

  /* Remove current SUNHeuristics object, and replace with "Default" */
  retval = SUNHeuristicsSpace(ark_mem->hconstraints, &lenrw, &leniw);
  if (retval == SUNHEURISTICS_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNHeuristicsDestroy(ark_mem->hconstraints);
  ark_mem->hconstraints = NULL;
  ark_mem->hconstraints = SUNHeuristicsDefault(ark_mem->sunctx);
  if (ark_mem->hconstraints == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep",
                    "ERKStepSetDefaults",
                    "SUNHeuristicsDefault allocation failure");
    return(ARK_MEM_FAIL);
  }
  retval = SUNHeuristicsSpace(ark_mem->hconstraints, &lenrw, &leniw);
  if (retval == SUNHEURISTICS_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }

  /* Remove current SUNControl object, and replace with "PI" */
  retval = SUNControlSpace(ark_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNControlDestroy(ark_mem->hcontroller);
  ark_mem->hcontroller = NULL;
  ark_mem->hcontroller = SUNControlPI(ark_mem->sunctx);
  if (ark_mem->hcontroller == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep",
                    "ERKStepSetDefaults",
                    "SUNControlPI allocation failure");
    return(ARK_MEM_FAIL);
  }

  /* Set default values for integrator optional inputs
     (overwrite some adaptivity params for ERKStep use) */
  step_mem->q = Q_DEFAULT;                     /* method order */
  step_mem->p = 0;                             /* embedding order */
  step_mem->stages = 0;                        /* no stages */
  step_mem->B = NULL;                          /* no Butcher table */
  (void) SUNHeuristicsSetMaxEFailGrowth(ark_mem->hconstraints, RCONST(0.3));
  (void) SUNControlSetSafetyFactor(ark_mem->hcontroller, RCONST(0.99));
  (void) SUNControlSetErrorBias(ark_mem->hcontroller, RCONST(1.2));
  (void) SUNHeuristicsSetMaxGrowth(ark_mem->hconstraints, RCONST(25.0));
  (void) SUNControlPI_SetParams(ark_mem->hcontroller, SUNFALSE,
                                RCONST(0.8), RCONST(0.31));
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int ERKStepSetOrder(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetOrder",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    step_mem->q = Q_DEFAULT;
  } else {
    step_mem->q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  step_mem->stages = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.

  If d==NULL, then the method is automatically flagged as a
  fixed-step method; a user MUST also call either
  ERKStepSetFixedStep or ERKStepSetInitStep to set the desired
  time step size.
  ---------------------------------------------------------------*/
int ERKStepSetTable(void *arkode_mem, ARKodeButcherTable B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for legal inputs */
  if (B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSetTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /* set the relevant parameters */
  step_mem->stages = B->stages;
  step_mem->q = B->q;
  step_mem->p = B->p;

  /* copy the table into step memory */
  step_mem->B = ARKodeButcherTable_Copy(B);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSetTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetTableNum:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the integer flag passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int ERKStepSetTableNum(void *arkode_mem, ARKODE_ERKTableID etable)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that argument specifies an explicit table */
  if (etable<ARKODE_MIN_ERK_NUM || etable>ARKODE_MAX_ERK_NUM) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSetTableNum",
                    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /* fill in table based on argument */
  step_mem->B = ARKodeButcherTable_LoadERK(etable);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSetTableNum",
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  step_mem->stages = step_mem->B->stages;
  step_mem->q = step_mem->B->q;
  step_mem->p = step_mem->B->p;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetTableName:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the string passed to ARKodeButcherTable_LoadERKByNmae()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int ERKStepSetTableName(void *arkode_mem, const char *etable)
{
  return ERKStepSetTableNum(arkode_mem,
                            arkButcherTableERKNameToID(etable));
}

/*===============================================================
  ERKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int ERKStepGetNumRhsEvals(void *arkode_mem, long int *fevals)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get values from step_mem */
  *fevals = step_mem->nfe;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetCurrentButcherTable:

  Sets pointers to the Butcher table currently in use.
  ---------------------------------------------------------------*/
int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                  ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetCurrentButcherTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get tables from step_mem */
  *B = step_mem->B;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetEstLocalErrors: (updated to the correct vector, but
  need to verify that it is unchanged between filling the
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ERKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetEstLocalErrors",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ERKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                               long int *accsteps, long int *attempts,
                               long int *fevals, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetTimestepperStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set integrator outputs */
  *attempts = ark_mem->nst_attempts;
  *fevals   = step_mem->nfe;
  *netfails = ark_mem->netf;

  /* get expsteps and accsteps from heuristics object */
  retval = SUNHeuristicsGetNumExpSteps(ark_mem->hconstraints, expsteps);
  if (retval != ARK_SUCCESS) { return(ARK_HEURISTICS_ERR); }
  retval = SUNHeuristicsGetNumAccSteps(ark_mem->hconstraints, accsteps);
  if (retval != ARK_SUCCESS) { return(ARK_HEURISTICS_ERR); }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumExpSteps:

  Returns the current number of stability-limited steps.

  This is a convenience routine, that merely calls
  SUNHeuristicsGetNumExpSteps for the result.
  ---------------------------------------------------------------*/
int ERKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepGetNumExpSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (SUNHeuristicsGetNumExpSteps(ark_mem->hconstraints, nsteps)
      != SUNHEURISTICS_SUCCESS) { return(ARK_HEURISTICS_ERR); }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps.

  This is a convenience routine, that merely calls
  SUNHeuristicsGetNumAccSteps for the result.
  ---------------------------------------------------------------*/
int ERKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ERKStep",
                    "ERKStepGetNumAccSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (SUNHeuristicsGetNumAccSteps(ark_mem->hconstraints, nsteps)
      != SUNHEURISTICS_SUCCESS) { return(ARK_HEURISTICS_ERR); }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepPrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int ERKStepPrintAllStats(void *arkode_mem, FILE *outfile, SUNOutputFormat fmt)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepPrintAllStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  retval = arkPrintAllStats(arkode_mem, outfile, fmt);
  if (retval != ARK_SUCCESS) return(retval);

  switch(fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "RHS fn evals                 = %ld\n", step_mem->nfe);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, "\n");
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ERKStepPrintAllStats",
                    "Invalid formatting option.");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ERKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKODE infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepWriteParameters",
                    "Error writing ARKODE infrastructure parameters");
    return(retval);
  }

  /* print integrator parameters to file */
  fprintf(fp, "ERKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",step_mem->q);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepWriteButcher", "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* print Butcher table to file */
  fprintf(fp, "\nERKStep Butcher table (stages = %i):\n", step_mem->stages);
  ARKodeButcherTable_Write(step_mem->B, fp);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
