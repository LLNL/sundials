/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for ARKODE's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "arkode/arkode_sprkstep.h"
#include "arkode_impl.h"
#include "arkode_sprkstep_impl.h"
#include "arkode_interp_impl.h"

/*===============================================================
  SPRKStep Exported functions -- Required
  ===============================================================*/

void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, realtype t0, N_Vector y0, SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  booleantype nvectorOK;
  int retval;

  /* Check that f1 and f2 are supplied */
  if (f1 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  if (f2 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  if (!sunctx) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_NULL_SUNCTX);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = sprkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeSPRKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeSPRKStepMem) malloc(sizeof(struct ARKodeSPRKStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::SPRKStep",
                    "SPRKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeSPRKStepMemRec));

  /* Allocate vectors in stepper mem */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata))) {
    SPRKStepFree((void**) &ark_mem);
    return(NULL); 
  }

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol   = NULL;
  ark_mem->step_attachmasssol  = NULL;
  ark_mem->step_disablelsetup  = NULL;
  ark_mem->step_disablemsetup  = NULL;
  ark_mem->step_getlinmem      = NULL;
  ark_mem->step_getmassmem     = NULL;
  ark_mem->step_getimplicitrhs = NULL;
  ark_mem->step_mmult          = NULL;
  ark_mem->step_getgammas      = NULL;
  ark_mem->step_init           = sprkStep_Init;
  ark_mem->step_fullrhs        = sprkStep_FullRHS;
  ark_mem->step                = sprkStep_TakeStep;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for SPRKStep optional inputs */
  retval = SPRKStepSetDefaults((void *)ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::SPRKStep",
                    "SPRKStepCreate",
                    "Error setting default solver options");
    SPRKStepFree((void**) &ark_mem);  return(NULL);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Update the ARKODE workspace requirements */
  // ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  // ark_mem->lrw += 10;

  /* Initialize the counters */
  step_mem->nf1 = 0;
  step_mem->nf2 = 0;
  step_mem->istage = 0;

  // /* Initialize external polynomial forcing data */
  // step_mem->expforcing = SUNFALSE;
  // step_mem->impforcing = SUNFALSE;
  // step_mem->forcing    = NULL;
  // step_mem->nforcing   = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::SPRKStep", "SPRKStepCreate",
                    "Unable to initialize main ARKODE infrastructure");
    SPRKStepFree((void**) &ark_mem);  return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
  SPRKStepResize:

  This routine resizes the memory within the SPRKStep module.
  It first resizes the main ARKODE infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int SPRKStepResize(void *arkode_mem, N_Vector y0, realtype hscale,
                  realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepResize",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* resize ARKODE infrastructure memory */
  retval = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::SPRKStep", "SPRKStepResize",
                    "Unable to resize main ARKODE infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  SPRKStepReInit:

  This routine re-initializes the SPRKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE::SPRKStep",
                    "SPRKStepReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check that f1 and f2 are supplied */
  if (f1 == NULL || f2 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::SPRKStep",
                    "SPRKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::SPRKStep", "SPRKStepReInit",
                    "Unable to reinitialize main ARKODE infrastructure");
    return(retval);
  }

  /* Initialize the counters */
  step_mem->nf1 = 0;
  step_mem->nf2 = 0;
  step_mem->istage = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  SPRKStepReset:

  This routine resets the SPRKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int SPRKStepReset(void* arkode_mem, realtype tR, N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepReset",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::SPRKStep", "SPRKStepReset",
                    "Unable to initialize main ARKODE infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  SPRKStepSStolerances, SPRKStepSVtolerances, SPRKStepWFtolerances,
  SPRKStepResStolerance, SPRKStepResVtolerance, SPRKStepResFtolerance:

  These routines set integration tolerances (wrappers for general
  ARKODE utility routines)
  ---------------------------------------------------------------*/

// int SPRKStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
// {
//   return(arkSStolerances((ARKodeMem) arkode_mem, reltol, abstol));
// }

// int SPRKStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
// {
//   return(arkSVtolerances((ARKodeMem) arkode_mem, reltol, abstol));
// }

// int SPRKStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
// {
//   return(arkWFtolerances((ARKodeMem) arkode_mem, efun));
// }

// int SPRKStepResStolerance(void *arkode_mem, realtype rabstol)
// {
//   return(arkResStolerance((ARKodeMem) arkode_mem, rabstol));
// }

// int SPRKStepResVtolerance(void *arkode_mem, N_Vector rabstol)
// {
//   return(arkResVtolerance((ARKodeMem) arkode_mem, rabstol));
// }

// int SPRKStepResFtolerance(void *arkode_mem, ARKRwtFn rfun)
// {
//   return(arkResFtolerance((ARKodeMem) arkode_mem, rfun));
// }


/*---------------------------------------------------------------
  SPRKStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKODE utility routine)
  ---------------------------------------------------------------*/
int SPRKStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  int retval;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkEvolve((ARKodeMem) arkode_mem, tout, yout, tret, itask);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}


/*---------------------------------------------------------------
  SPRKStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKODE utility routine)
  ---------------------------------------------------------------*/
int SPRKStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  int retval;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkGetDky((ARKodeMem) arkode_mem, t, k, dky);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}

/*---------------------------------------------------------------
  SPRKStepFree frees all SPRKStep memory, and then calls an ARKODE
  utility routine to free the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void SPRKStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL SPRKStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {
    step_mem = (ARKodeSPRKStepMem) ark_mem->step_mem;
    
    if (step_mem->sdata != NULL) {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }
    
    ARKodeSPRKMem_Free(step_mem->method);

    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }

  /* free memory for overall ARKODE infrastructure */
  arkFree(arkode_mem);
}


/*===============================================================
  SPRKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKODE
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
  sprkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  For all initialization types, this routine sets the relevant
  TakeStep routine based on the current problem configuration.

  With initialization type FIRST_INIT this routine:

  With initialization type FIRST_INIT or RESIZE_INIT, this routine:

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
int sprkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  int j, retval;
  booleantype reset_efun;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return(ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT) {

    if (!step_mem->method) {
      switch (step_mem->q) {
        case 1:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_1);
          break;
        case 2:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_2);
          break;
        case 3:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_3);
          break;
        case 4:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_4);
          break;
        case 5:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_5);
          break;
        case 6:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_6);
          break;
        case 7:
        case 8:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_8);
          break;
        case 9:
        case 10:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_10);
        default:
          step_mem->method = ARKodeSPRKMem_Load(SPRKSTEP_DEFAULT_4);
          break;
      }
    }
   
  }

  return(ARK_SUCCESS);
}

/* Utility to call f1 and increment the counter */
int sprkStep_f1(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur, N_Vector f1, void* user_data)
{
  int retval = step_mem->f1(tcur, ycur, f1, user_data);
  step_mem->nf1++;
  return retval;
}

/* Utility to call f2 and increment the counter */
int sprkStep_f2(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur, N_Vector f2, void* user_data)
{
  int retval = step_mem->f2(tcur, ycur, f2, user_data);
  step_mem->nf2++;
  return retval;
}

/*---------------------------------------------------------------
  sprkStep_FullRHS:

  Rewriting the problem
    My' = fe(t,y) + fi(t,y)
  in the form
    y' = M^{-1}*[ fe(t,y) + fi(t,y) ],
  this routine computes the full right-hand side vector,
    f = M^{-1}*[ fe(t,y) + fi(t,y) ]

  This will be called in one of three 'modes':
    ARK_FULLRHS_START -> called at the beginning of a simulation
                         or after post processing at step
    ARK_FULLRHS_END   -> called at the end of a successful step
    ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If it is called in ARK_FULLRHS_START mode, we store the vectors
  fe(t,y) and fi(t,y) in Fe[0] and Fi[0] for possible reuse in the
  first stage of the subsequent time step.

  If it is called in ARK_FULLRHS_END mode and the ARK method
  coefficients support it, we may just copy vectors Fe[stages] and
  Fi[stages] to fill f instead of calling fe() and fi().

  ARK_FULLRHS_OTHER mode is only called for dense output in-between
  steps, or when estimating the initial time step size, so we strive to
  store the intermediate parts so that they do not interfere
  with the other two modes.
  ---------------------------------------------------------------*/
int sprkStep_FullRHS(void* arkode_mem, realtype t, N_Vector y, N_Vector f,
                     int mode)
{
  // ARKodeMem ark_mem;
  // ARKodeSPRKStepMem step_mem;
  // int nvec, retval;
  // booleantype recomputeRHS;
  // realtype* cvals;
  // N_Vector* Xvecs;

  // /* access ARKodeSPRKStepMem structure */
  // retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_FullRHS",
  //                                &ark_mem, &step_mem);
  // if (retval != ARK_SUCCESS)  return(retval);

  // /* local shortcuts for use with fused vector operations */
  // cvals = step_mem->cvals;
  // Xvecs = step_mem->Xvecs;

  // /* setup mass-matrix if required (use output f as a temporary) */
  // if ((step_mem->mass_type == MASS_TIMEDEP) && (step_mem->msetup != NULL)) {
  //   retval = step_mem->msetup((void *) ark_mem, t, f,
  //                             ark_mem->tempv2, ark_mem->tempv3);
  //   if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
  // }

  // /* perform RHS functions contingent on 'mode' argument */
  // switch(mode) {

  // /* ARK_FULLRHS_START: called at the beginning of a simulation
  //    Store the vectors fe(t,y) and fi(t,y) in Fe[0] and Fi[0] for
  //    possible reuse in the first stage of the subsequent time step */
  // case ARK_FULLRHS_START:

  //   /* call fe if the problem has an explicit component */
  //   if (step_mem->explicit) {
  //     retval = sprkStep_Fe(step_mem, t, y, step_mem->Fe[0], ark_mem->user_data);
  //     if (retval != 0) {
  //       arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                       "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //       return(ARK_RHSFUNC_FAIL);
  //     }
  //     /* apply external polynomial forcing */
  //     if (step_mem->expforcing) {
  //       cvals[0] = ONE;
  //       Xvecs[0] = step_mem->Fe[0];
  //       nvec     = 1;
  //       sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //       N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fe[0]);
  //     }
  //   }

  //   /* call fi if the problem has an implicit component */
  //   if (step_mem->implicit) {
  //     retval = sprkStep_Fi(step_mem, t, y, step_mem->Fi[0], ark_mem->user_data);
  //     if (retval != 0) {
  //       arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                       "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //       return(ARK_RHSFUNC_FAIL);
  //     }
  //     /* apply external polynomial forcing */
  //     if (step_mem->impforcing) {
  //       cvals[0] = ONE;
  //       Xvecs[0] = step_mem->Fi[0];
  //       nvec     = 1;
  //       sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //       N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fi[0]);
  //     }
  //   }

  //   /* combine RHS vector(s) into output */
  //   if (step_mem->explicit && step_mem->implicit) { /* ImEx */
  //     N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
  //   } else if (step_mem->implicit) {                   /* implicit */
  //     N_VScale(ONE, step_mem->Fi[0], f);
  //   } else {                                           /* explicit */
  //     N_VScale(ONE, step_mem->Fe[0], f);
  //   }

  //   break;


  // /* ARK_FULLRHS_END: called at the end of a successful step
  //    If the ARK method coefficients support it, we just copy the last stage RHS
  //    vectors to fill f instead of calling fe() and fi().
  //    Copy the results to Fe[0] and Fi[0] if the ARK coefficients support it. */
  // case ARK_FULLRHS_END:

  //   /* determine if explicit/implicit RHS functions need to be recomputed */
  //   recomputeRHS = SUNFALSE;
  //   if ( step_mem->explicit && (SUNRabs(step_mem->Be->c[step_mem->stages-1]-ONE)>TINY) )
  //     recomputeRHS = SUNTRUE;
  //   if ( step_mem->implicit && (SUNRabs(step_mem->Bi->c[step_mem->stages-1]-ONE)>TINY) )
  //     recomputeRHS = SUNTRUE;

  //   /* base RHS calls on recomputeRHS argument */
  //   if (recomputeRHS) {

  //     /* call fe if the problem has an explicit component */
  //     if (step_mem->explicit) {
  //       retval = sprkStep_Fe(step_mem, t, y, step_mem->Fe[0], ark_mem->user_data);
  //       if (retval != 0) {
  //         arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                         "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //         return(ARK_RHSFUNC_FAIL);
  //       }
  //       /* apply external polynomial forcing */
  //       if (step_mem->expforcing) {
  //         cvals[0] = ONE;
  //         Xvecs[0] = step_mem->Fe[0];
  //         nvec     = 1;
  //         sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //         N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fe[0]);
  //       }
  //     }

  //     /* call fi if the problem has an implicit component */
  //     if (step_mem->implicit) {
  //       retval = sprkStep_Fi(step_mem, t, y, step_mem->Fi[0], ark_mem->user_data);
  //       if (retval != 0) {
  //         arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                         "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //         return(ARK_RHSFUNC_FAIL);
  //       }
  //       /* apply external polynomial forcing */
  //       if (step_mem->impforcing) {
  //         cvals[0] = ONE;
  //         Xvecs[0] = step_mem->Fi[0];
  //         nvec     = 1;
  //         sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //         N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fi[0]);
  //       }
  //     }
  //   } else {
  //     if (step_mem->explicit)
  //       N_VScale(ONE, step_mem->Fe[step_mem->stages-1], step_mem->Fe[0]);
  //     if (step_mem->implicit)
  //       N_VScale(ONE, step_mem->Fi[step_mem->stages-1], step_mem->Fi[0]);
  //   }

  //   /* combine RHS vector(s) into output */
  //   if (step_mem->explicit && step_mem->implicit) { /* ImEx */
  //     N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
  //   } else if (step_mem->implicit) {                   /* implicit */
  //     N_VScale(ONE, step_mem->Fi[0], f);
  //   } else {                                           /* explicit */
  //     N_VScale(ONE, step_mem->Fe[0], f);
  //   }

  //   break;

  // /* ARK_FULLRHS_OTHER: called for dense output in-between steps or for
  //    estimation of the initial time step size, store the intermediate
  //    calculations in such a way as to not interfere with the other two modes */
  // case ARK_FULLRHS_OTHER:

  //   /* call fe if the problem has an explicit component (store in ark_tempv2) */
  //   if (step_mem->explicit) {
  //     retval = sprkStep_Fe(step_mem, t, y, ark_mem->tempv2, ark_mem->user_data);
  //     if (retval != 0) {
  //       arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                       "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //       return(ARK_RHSFUNC_FAIL);
  //     }
  //     /* apply external polynomial forcing */
  //     if (step_mem->expforcing) {
  //       cvals[0] = ONE;
  //       Xvecs[0] = ark_mem->tempv2;
  //       nvec     = 1;
  //       sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //       N_VLinearCombination(nvec, cvals, Xvecs, ark_mem->tempv2);
  //     }
  //   }

  //   /* call fi if the problem has an implicit component (store in sdata) */
  //   if (step_mem->implicit) {
  //     retval = sprkStep_Fi(step_mem, t, y, step_mem->sdata, ark_mem->user_data);
  //     if (retval != 0) {
  //       arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                       "sprkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
  //       return(ARK_RHSFUNC_FAIL);
  //     }
  //     /* apply external polynomial forcing */
  //     if (step_mem->impforcing) {
  //       cvals[0] = ONE;
  //       Xvecs[0] = step_mem->sdata;
  //       nvec     = 1;
  //       sprkStep_ApplyForcing(step_mem, t, ONE, &nvec);
  //       N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
  //     }
  //   }

  //   /* combine RHS vector(s) into output */
  //   if (step_mem->explicit && step_mem->implicit) { /* ImEx */
  //     N_VLinearSum(ONE, step_mem->sdata, ONE, ark_mem->tempv2, f);
  //   } else if (step_mem->implicit) {                   /* implicit */
  //     N_VScale(ONE, step_mem->sdata, f);
  //   } else {                                           /* explicit */
  //     N_VScale(ONE, ark_mem->tempv2, f);
  //   }

  //   break;

  // default:
  //   /* return with RHS failure if unknown mode is passed */
  //   arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::SPRKStep",
  //                   "sprkStep_FullRHS", "Unknown full RHS mode");
  //   return(ARK_RHSFUNC_FAIL);
  // }

  // /* if M != I, then update f = M^{-1}*f */
  // if (step_mem->mass_type != MASS_IDENTITY) {
  //   retval = step_mem->msolve((void *) ark_mem, f,
  //                             step_mem->nlscoef/ark_mem->h);
  //   if (retval != ARK_SUCCESS) {
  //     arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE::SPRKStep",
  //                     "sprkStep_FullRHS", "Mass matrix solver failure");
  //     return(ARK_MASSSOLVE_FAIL);
  //   }
  // }

  return(ARK_SUCCESS);
}

int sprkStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  ARKodeSPRKMem method_mem;
  N_Vector prev_stage;
  int retval, is;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_TakeStep",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  method_mem = step_mem->method;

  prev_stage = ark_mem->yn;
  for (is = 0; is < method_mem->stages; is++) {
    retval = sprkStep_SPRKStage(ark_mem, step_mem, prev_stage,
                                method_mem->b[is], method_mem->a[is], ark_mem->ycur);
    if (retval != ARK_SUCCESS) return(retval);
    prev_stage = ark_mem->ycur;
  }

  *nflagPtr = 0;
  *dsmPtr = 0;

  return ARK_SUCCESS;
}

/* Increment SPRK algorithm with compensated summation */
int sprkStep_TakeStep_Compensated(void* arkode_mem, realtype *dsmPtr, int *nflagPtr)
{
  ARKodeMem ark_mem;
  ARKodeSPRKStepMem step_mem;
  ARKodeSPRKMem method;
  int retval, is;
  N_Vector delta_Yi, yn_plus_delta_Yi;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_TakeStep_SPRK",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  method = step_mem->method;

  /* other shortcuts */
  delta_Yi = ark_mem->tempv1;
  yn_plus_delta_Yi = ark_mem->tempv2;

  /* [ \Delta Q_0 ] = [ 0 ]
     [ \Delta P_0 ] = [ 0 ] */
  N_VConst(ZERO, delta_Yi);

  /* loop over internal stages to the step */
  for (is=0; is < method->stages; is++) {
    /* store current stage index */
    step_mem->istage = is;

    /* set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + method->b[is]*ark_mem->h;

    /* [ q_n ] + [ \Delta Q_i ]
       [     ] + [            ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    /* evaluate Fi with previous stage increment */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to set other outputs to zero */
    retval = sprkStep_f1(step_mem, ark_mem->tcur, yn_plus_delta_Yi,
                         step_mem->sdata, ark_mem->user_data);
    if (retval != 0) return(ARK_RHSFUNC_FAIL);

    /* update the implicit stage
       [            ] = [                ] + [       ]
       [ \Delta P_i ] = [ \Delta P_{i-1} ] + [ sdata ] */
    N_VLinearSum(ONE, delta_Yi,  ark_mem->h * method->b[is], step_mem->sdata, delta_Yi);

    /* [     ] + [            ]
       [ p_n ] + [ \Delta P_i ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    /* evaluate Fe with the current p_n + \Delta P_i */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to set other outputs to zero */
    retval = sprkStep_f2(step_mem, ark_mem->tn + method->a[is]*ark_mem->h,
                         yn_plus_delta_Yi, step_mem->sdata, ark_mem->user_data);
    if (retval != 0) return(ARK_RHSFUNC_FAIL);

    /* update the explicit stage
       [ \Delta Q_i ] = [ \Delta Q_{i-1} ] + [ sdata ]
       [            ] = [                ] + [       ] */
    N_VLinearSum(ONE, delta_Yi, ark_mem->h * method->a[is], step_mem->sdata, delta_Yi);
  }

  /*
    Now we compute the step solution via compensated summation.
     [ q_{n+1} ] = [ q_n ] + [ \Delta Q_i ]
     [ p_{n+1} ] = [ p_n ] + [ \Delta P_i ] */
  N_VLinearSum(ONE, delta_Yi, -ONE, ark_mem->yerr, delta_Yi);
  N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->yn, ark_mem->tempv3);
  N_VLinearSum(ONE, ark_mem->tempv3, -ONE, delta_Yi, ark_mem->yerr);

  *nflagPtr = 0;
  *dsmPtr = 0;

  return 0;
}

/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
  sprkStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int sprkStep_AccessStepMem(void* arkode_mem, const char *fname,
                           ARKodeMem *ark_mem, ARKodeSPRKStepMem *step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::SPRKStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKODE::SPRKStep",
                    fname, MSG_SPRKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSPRKStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  sprkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype sprkStep_CheckNVector(N_Vector tmpl)
{
  if ( (tmpl->ops->nvclone     == NULL) ||
       (tmpl->ops->nvdestroy   == NULL) ||
       (tmpl->ops->nvlinearsum == NULL) ||
       (tmpl->ops->nvconst     == NULL) ||
       (tmpl->ops->nvscale     == NULL) ||
       (tmpl->ops->nvwrmsnorm  == NULL) )
    return(SUNFALSE);
  return(SUNTRUE);
}

int sprkStep_SPRKStage(ARKodeMem ark_mem, ARKodeSPRKStepMem step_mem, N_Vector prev_stage,
                       sunrealtype bi, sunrealtype ai, N_Vector stage_result)
{
  int retval = 0;

  /* set current stage time(s) */
  ark_mem->tcur = ark_mem->tn + bi*ark_mem->h;

  /* evaluate f_1 at the previous stage value */
  N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to set other outputs to zero */
  retval = sprkStep_f2(step_mem, ark_mem->tcur, prev_stage, step_mem->sdata, ark_mem->user_data);
  if (retval != 0) return ARK_RHSFUNC_FAIL;

  /* update ycur with the q stage */
  N_VLinearSum(ONE, prev_stage, ark_mem->h*bi, step_mem->sdata, stage_result);

  /* evaluate f_2 with the stage value for q */
  N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to set other outputs to zero */
  retval = sprkStep_f1(step_mem, ark_mem->tn + ai*ark_mem->h, stage_result, step_mem->sdata, ark_mem->user_data);
  if (retval != 0) return ARK_RHSFUNC_FAIL;

  /* update ycur with the stage value for p */
  N_VLinearSum(ONE, stage_result, ark_mem->h*ai, step_mem->sdata, stage_result);

  /* keep track of the stage number */
  step_mem->istage++;

  return ARK_SUCCESS;
}
