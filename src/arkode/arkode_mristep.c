/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for ARKode's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

#define NO_DEBUG_OUTPUT
/* #define DEBUG_OUTPUT */
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)



/*===============================================================
  MRIStep Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  Create MRIStep integrator memory struct
  ---------------------------------------------------------------*/
void* MRIStepCreate(ARKRhsFn fs, realtype t0, N_Vector y0,
                    MRISTEP_ID inner_step_id, void* inner_mem)
{
  void *arkode_mem; /* ARKode MRIStep memory */
  int   retval;     /* return value          */

  /* Create the ARKode MRIStep memory */
  arkode_mem = mriStep_Create(fs, t0, y0);
  if (arkode_mem == NULL) return(NULL);

  /* Set the inner integrator */
  switch (inner_step_id) {
  case MRISTEP_ARKSTEP:
    retval = mriStep_AttachARK(arkode_mem, inner_mem);
    break;
  default:
    arkProcessError((ARKodeMem) arkode_mem, ARK_ILL_INPUT,
                    "ARKode::MRIStep", "MRIStepCreate",
                    "Invalid inner integrator option");
    MRIStepFree(&arkode_mem);
    return(NULL);
  }

  /* check if inner integrator was attached successfully */
  if (retval != ARK_SUCCESS) return(NULL);

  /* return ARKode memory */
  return(arkode_mem);
}


/*---------------------------------------------------------------
  MRIStepResize:

  This routine resizes the memory within the MRIStep module.
  It first resizes the main ARKode infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int MRIStepResize(void *arkode_mem, N_Vector y0, realtype t0,
                  ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval, i, flag;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepResize",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKode infrastructure memory (use hscale = 1.0) */
  flag = arkResize(ark_mem, y0, RCONST(1.0), t0, resize, resize_data);
  if (flag != ARK_SUCCESS) {
    arkProcessError(ark_mem, flag, "ARKode::MRIStep", "MRIStepResize",
                    "Unable to resize main ARKode infrastructure");
    return(flag);
  }

  /* Resize the inner forcing vector */
  if (step_mem->inner_forcing != NULL) {
    for (i = 0; i < step_mem->inner_num_forcing; i++) {
      retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                            liw_diff, y0, &(step_mem->inner_forcing[i]));
      if (retval != ARK_SUCCESS) return(retval);
    }
  }

  /* Resize the RHS vectors */
  for (i=0; i<step_mem->stages; i++) {
    retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                          liw_diff, y0, &step_mem->F[i]);
    if (retval != ARK_SUCCESS)  return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepReInit:

  This routine re-initializes the MRIStep module to solve a new
  problem of the same size as was previously solved.

  NOTE: the inner stepper needs to be reinitialized before
  calling this function.
  ---------------------------------------------------------------*/
int MRIStepReInit(void* arkode_mem, ARKRhsFn fs, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check that fs is supplied */
  if (fs == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  retval = arkReInit(arkode_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(retval);
  }

  /* Copy the input parameters into ARKode state */
  step_mem->fs = fs;

  /* Initialize all the counters */
  step_mem->nfs = 0;

  /* Reattach the inner integrator */
  switch (step_mem->inner_stepper_id) {
  case MRISTEP_ARKSTEP:
    retval = mriStep_AttachARK(arkode_mem, step_mem->inner_mem);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT,
                    "ARKode::MRIStep", "MRIStepReInit",
                    "Invalid inner integrator option");
    return(ARK_ILL_INPUT);
  }

  /* check if inner integrator was attached successfully */
  if (retval != ARK_SUCCESS) return(ARK_INNERSTEP_ATTACH_ERR);

  /* Reattach user data */
  if (ark_mem->user_data != NULL) {
    retval = MRIStepSetUserData(arkode_mem, ark_mem->user_data);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepReInit",
                      "Unable to re-attach user_data");
      return(retval);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepRootInit:

  Initialize (attach) a rootfinding problem to the stepper
  (wrappers for general ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


/*---------------------------------------------------------------
  MRIStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkEvolve(ark_mem, tout, yout, tret, itask));
}


/*---------------------------------------------------------------
  MRIStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetDky(ark_mem, t, k, dky));
}


/*---------------------------------------------------------------
  MRIStepFree frees all MRIStep memory, and then calls an ARKode
  utility routine to free the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void MRIStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL MRIStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

    /* free the Butcher table */
    if (step_mem->B != NULL) {
      ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->B);
      step_mem->B = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }

    /* free the inner forcing vector */
    if (step_mem->inner_forcing != NULL) {
      for (j = 0; j < step_mem->inner_num_forcing; j++) {
        arkFreeVec(ark_mem, &(step_mem->inner_forcing[j]));
        step_mem->inner_forcing[j] = NULL;
      }
    }

    /* free the RHS vectors */
    if (step_mem->F != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->F[j]);
      free(step_mem->F);
      step_mem->F = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL) {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= (step_mem->stages + 1);
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (step_mem->stages + 1);
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }

  /* free memory for overall ARKode infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
  MRIStepPrintMem:

  This routine outputs the memory from the MRIStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void MRIStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepPrintMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* output data from main ARKode infrastructure */
  fprintf(outfile,"MRIStep Slow Stepper Mem:\n");
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"MRIStep: q = %i\n", step_mem->q);
  fprintf(outfile,"MRIStep: p = %i\n", step_mem->p);
  fprintf(outfile,"MRIStep: stages = %i\n", step_mem->stages);

  /* output long integer quantities */
  fprintf(outfile,"MRIStep: nfs = %li\n", step_mem->nfs);

  /* output realtype quantities */
  fprintf(outfile,"MRIStep: Butcher table:\n");
  ARKodeButcherTable_Write(step_mem->B, outfile);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */
  for (i=0; i<step_mem->stages; i++) {
    fprintf(outfile,"MRIStep: F[%i]:\n", i);
    N_VPrint_Serial(step_mem->F[i]);
  }
#endif

  return;
}



/*===============================================================
  MRIStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Create outer memory structure
  ---------------------------------------------------------------*/

void* mriStep_Create(ARKRhsFn fs, realtype t0, N_Vector y0)
{
  ARKodeMem        ark_mem;         /* outer ARKode memory  */
  ARKodeMRIStepMem step_mem;        /* outer stepper memory */
  booleantype      nvectorOK;
  int              retval;

  /* Check that fs is supplied */
  if (fs == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = mriStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate();
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeMRIStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeMRIStepMem) malloc(sizeof(struct ARKodeMRIStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeMRIStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init    = mriStep_Init;
  ark_mem->step_fullrhs = mriStep_FullRHS;
  ark_mem->step         = mriStep_TakeStep;
  ark_mem->step_mem     = (void*) step_mem;

  /* Set default values for MRIStep optional inputs */
  retval = MRIStepSetDefaults((void *) ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep",
                    "MRIStepCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Allocate the general MRI stepper vectors using y0 as a template */
  /* NOTE: F, inner_forcing, cvals and Xvecs will be allocated later on
     (based on the MRI method) */

  /* Copy the slow RHS function into stepper memory */
  step_mem->fs = fs;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 11;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 1;

  /* Initialize all the counters */
  step_mem->nfs = 0;

  /* Initialize pre and post inner evolve functions */
  step_mem->pre_inner_evolve  = NULL;
  step_mem->post_inner_evolve = NULL;

  /* Initialize main ARKode infrastructure (allocates vectors) */
  retval = arkInit(ark_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
  Attach ARKStep inner integrator
  ---------------------------------------------------------------*/

int mriStep_AttachARK(void* arkode_mem, void* inner_mem)
{
  ARKodeMem         ark_mem;         /* outer ARKode memory  */
  ARKodeMRIStepMem  step_mem;        /* outer stepper memory */
  ARKodeMem         inner_ark_mem;   /* inner ARKode memory  */
  ARKodeARKStepMem  inner_step_mem;  /* inner stepper memory */
  int               retval;          /* return value         */

  /* Access MRIStep memory */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_AttachARK",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { MRIStepFree(&arkode_mem); return(-1); }

  /* Access ARKStep memory */
  retval = arkStep_AccessStepMem(inner_mem, "mriStep_AttachARK",
                                 &inner_ark_mem, &inner_step_mem);
  if (retval != ARK_SUCCESS) { MRIStepFree(&arkode_mem); return(-1); }

  /* Set RHS function wrappers */
  if (inner_step_mem->explicit) {
    step_mem->ff1      = inner_step_mem->fe;
    inner_step_mem->fe = mriStep_InnerRhsFnForcing;
    if (inner_step_mem->implicit) {
      step_mem->ff2      = inner_step_mem->fi;
      inner_step_mem->fi = mriStep_InnerRhsFn;
    }
  } else {
    step_mem->ff1      = inner_step_mem->fi;
    inner_step_mem->fi = mriStep_InnerRhsFnForcing;
  }

  /* attach the inner stepper and initialize the inner stepper return flag */
  step_mem->inner_mem        = inner_mem;
  step_mem->inner_stepper_id = MRISTEP_ARKSTEP;
  step_mem->inner_retval     = ARK_SUCCESS;
  step_mem->inner_evolve     = mriStep_EvolveInnerARK;

  /* attach outer stepper mem to inner stepper as user data */
  inner_ark_mem->user_data = arkode_mem;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup (init_type == 0) or arkPostResizeSetup
  (init_type == 1).

  With init_type == 0, this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK
    stages, method order, or solver options

  With init_type == 1, this routine does nothing.
  ---------------------------------------------------------------*/
int mriStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval, j;

  /* immediately return if init_type == 1 */
  if (init_type == 1)  return(ARK_SUCCESS);

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* assume fixed outer step size */
  if (!ark_mem->fixedstep) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                    "Adaptive outer time stepping is not currently supported");
    return(ARK_ILL_INPUT);
  }

  /* Create Butcher table (if not already set) */
  retval = mriStep_SetButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                    "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher table is OK */
  retval = mriStep_CheckButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_Init", "Error in Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* Allocate MRI RHS vector memory, update storage requirements */
  /*   Allocate F[0] ... F[stages-1] if needed */
  if (step_mem->F == NULL)
    step_mem->F = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
  for (j=0; j<step_mem->stages; j++) {
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->F[j])))
      return(ARK_MEM_FAIL);
  }
  ark_mem->liw += step_mem->stages;  /* pointers */

  /* Allocate MRI forcing vectors if needed */
  step_mem->inner_num_forcing = 1; /* only support 1 forcing vector at this time */
  if (step_mem->inner_forcing == NULL) {
    step_mem->inner_forcing = (N_Vector *) calloc(step_mem->inner_num_forcing,
                                                  sizeof(N_Vector));
    for (j = 0; j < step_mem->inner_num_forcing; j++) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->inner_forcing[j])))
        return(ARK_MEM_FAIL);
    }
  }

  /* Allocate reusable arrays for fused vector interface */
  if (step_mem->cvals == NULL) {
    step_mem->cvals = (realtype *) calloc(step_mem->stages+1, sizeof(realtype));
    if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
    ark_mem->lrw += (step_mem->stages + 1);
  }
  if (step_mem->Xvecs == NULL) {
    step_mem->Xvecs = (N_Vector *) calloc(step_mem->stages+1, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
    ark_mem->liw += (step_mem->stages + 1);   /* pointers */
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS functions,
  f(t,y) = fs(t,y) + ff(t,y).

  This will be called in one of three 'modes':
    0 -> called at the beginning of a simulation
    1 -> called at the end of a successful step
    2 -> called elsewhere (e.g. for dense output)

  If it is called in mode 0, we store the vectors f(t,y) in F[0]
  for possible reuse in the first stage of the subsequent time step.

  If it is called in mode 1, we reevauate f(t,y). At this time no
  checks are made to see if the method coefficient support copying
  vectors F[stages] to fill f instead of calling f().

  Mode 2 is only called for dense output in-between steps, so we
  strive to store the intermediate parts so that they do not
  interfere with the other two modes.
  ---------------------------------------------------------------*/
int mriStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vector fs(t,y) in F[0] for possible reuse
     in the first stage of the subsequent time step */
  case 0:

    /* call fs */
    retval = step_mem->fs(t, y, step_mem->F[0], ark_mem->user_data);
    step_mem->nfs++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* call ff */
    retval = mriStep_InnerFullRhs(arkode_mem, t, y, f);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    N_VLinearSum(ONE, step_mem->F[0], ONE, f, f);

    break;


  /* Mode 1: called at the end of a successful step
     This always recomputes the full RHS (i.e., this is the
     same as case 0). */
  case 1:

    /* call fs */
    retval = step_mem->fs(t, y, step_mem->F[0], ark_mem->user_data);
    step_mem->nfs++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* call ff */
    retval = mriStep_InnerFullRhs(arkode_mem, t, y, f);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    N_VLinearSum(ONE, step_mem->F[0], ONE, f, f);

    break;

  /*  Mode 2: called for dense output in-between steps
      store the intermediate calculations in such a way as to not
      interfere with the other two modes */
  default:

    /* call fs */
    retval = step_mem->fs(t, y, ark_mem->tempv2, ark_mem->user_data);
    step_mem->nfs++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* call ff */
    retval = mriStep_InnerFullRhs(arkode_mem, t, y, f);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    N_VLinearSum(ONE, ark_mem->tempv2, ONE, f, f);

    break;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_TakeStep:

  This routine serves the primary purpose of the MRIStep module:
  it performs a single successful MRI step (if possible).
  Multiple attempts may be taken in this process -- once a step
  completes with successful (non)linear solves at each stage and
  passes the error estimate, the routine returns successfully.
  If it cannot do so, it returns with an appropriate error flag.
  ---------------------------------------------------------------*/
int mriStep_TakeStep(void* arkode_mem)
{
  ARKodeMem         ark_mem;   /* outer ARKode memory   */
  ARKodeMRIStepMem  step_mem;  /* outer stepper memory  */

  realtype dsm, rcdiff, t0;
  int retval, is, eflag, js, nvec;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access the MRIStep mem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_TakeStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  dsm   = ZERO;
  eflag = ARK_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {

#ifdef DEBUG_OUTPUT
    printf("stage 0 RHS:\n");
    N_VPrint_Serial(step_mem->F[0]);
#endif

    /* Loop over internal stages to the step; since the method is explicit
       the first stage RHS is just the slow RHS from the start of the step */
    for (is=1; is<step_mem->stages; is++) {

      /* Set current stage time */
      ark_mem->tcur = ark_mem->tn + step_mem->B->c[is]*ark_mem->h;

#ifdef DEBUG_OUTPUT
      printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
             ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report)
        fprintf(ark_mem->diagfp, "MRIStep  step  %li  %"RSYM"  %i  %"RSYM"\n",
                ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

      /* compute forcing vector of inner steps (assumes c[is] != c[is-1]) */
      rcdiff = ONE / (step_mem->B->c[is] - step_mem->B->c[is-1]);
      nvec = 0;
      for (js=0; js<is; js++) {
        cvals[js] = rcdiff * (step_mem->B->A[is][js] - step_mem->B->A[is-1][js]);
        Xvecs[js] = step_mem->F[js];
        nvec++;
      }

      retval = N_VLinearCombination(nvec, cvals, Xvecs,
                                    step_mem->inner_forcing[0]);
      if (retval != 0) return(ARK_VECTOROP_ERR);

      /* initial time for inner integrator */
      t0 = ark_mem->tn + step_mem->B->c[is-1]*ark_mem->h;

      /* pre inner evolve function */
      if (step_mem->pre_inner_evolve) {
        retval = step_mem->pre_inner_evolve(t0, step_mem->inner_forcing,
                                            step_mem->inner_num_forcing,
                                            ark_mem->user_data);
        if (retval != 0) return(ARK_OUTERTOINNER_FAIL);
      }

      /* advance inner method in time */
      step_mem->inner_retval =
        step_mem->inner_evolve(step_mem->inner_mem, t0, ark_mem->ycur,
                               ark_mem->tcur);
      if (step_mem->inner_retval < 0) return(ARK_INNERSTEP_FAIL);

      /* post inner evolve function */
      if (step_mem->post_inner_evolve) {
        retval = step_mem->post_inner_evolve(ark_mem->tcur, ark_mem->ycur,
                                             ark_mem->user_data);
        if (retval != 0) return(ARK_INNERTOOUTER_FAIL);
      }

      /* compute updated slow RHS */
      retval = step_mem->fs(ark_mem->tcur, ark_mem->ycur,
                            step_mem->F[is], ark_mem->user_data);
      step_mem->nfs++;
      if (retval < 0)  return(ARK_RHSFUNC_FAIL);
      if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef DEBUG_OUTPUT
      printf("RHS:\n");
      N_VPrint_Serial(step_mem->F[is]);
#endif

    } /* loop over stages */

    /* Compute time step solution */

    /* set current time for solution */
    ark_mem->tcur = ark_mem->tn + ark_mem->h;

    /* compute forcing vector of inner steps (assumes c[stages-1] != 1) */
    rcdiff = ONE / (ONE - step_mem->B->c[step_mem->stages-1]);
    nvec = 0;
    for (js=0; js<step_mem->stages; js++) {
      cvals[js] = rcdiff * (step_mem->B->b[js] - step_mem->B->A[step_mem->stages-1][js]);
      Xvecs[js] = step_mem->F[js];
      nvec++;
    }

    retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->inner_forcing[0]);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* initial time for inner integrator */
    t0 = ark_mem->tn + step_mem->B->c[step_mem->stages-1]*ark_mem->h;

    /* pre inner evolve function */
    if (step_mem->pre_inner_evolve) {
      retval = step_mem->pre_inner_evolve(t0, step_mem->inner_forcing,
                                          step_mem->inner_num_forcing,
                                          ark_mem->user_data);
      if (retval != 0) return(ARK_OUTERTOINNER_FAIL);
    }

    /* advance inner method in time */
    step_mem->inner_retval =
      step_mem->inner_evolve(step_mem->inner_mem, t0, ark_mem->ycur,
                             ark_mem->tcur);
    if (step_mem->inner_retval < 0) return(ARK_INNERSTEP_FAIL);

    /* post inner evolve function */
    if (step_mem->post_inner_evolve) {
      retval = step_mem->post_inner_evolve(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0) return(ARK_INNERTOOUTER_FAIL);
    }

#ifdef DEBUG_OUTPUT
    printf("error estimate = %"RSYM"\n", dsm);
    printf("updated solution:\n");
    N_VPrint_Serial(ark_mem->ycur);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "MRIStep  etest  %li  %"RSYM"  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, dsm);

#ifdef DEBUG_OUTPUT
    printf("error test flag = %i\n", eflag);
#endif

    /* Restart step attempt (recompute all stages) if error test fails recoverably */
    if (eflag == TRY_AGAIN)  continue;

    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS)  return(eflag);

    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  } /* loop over step attempts */


  /* The step has completed successfully, clean up and
     consider change of step size */
  retval = mriStep_PrepareNextStep(ark_mem, dsm);
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    fname, MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeMRIStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}



/*---------------------------------------------------------------
  mriStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype mriStep_CheckNVector(N_Vector tmpl)
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


/*---------------------------------------------------------------
  mriStep_SetButcherTable

  This routine determines the MRI method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
int mriStep_SetButcherTable(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int itable;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_SetButcherTable", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* if table has already been specified, just return */
  if (step_mem->B != NULL) return(ARK_SUCCESS);

  /* initialize table number to an illegal value */
  itable = -1;

  /* select method based on order */
  switch (step_mem->q) {
  case(3):
    itable = DEFAULT_MRI_TABLE_3;
    break;
  default:    /* no available method, set default */
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_SetButcherTable",
                    "No explicit MRI method at requested order, using q=3.");
    itable = DEFAULT_MRI_TABLE_3;
    break;
  }

  /* load the Butcher table */
  if (itable > -1) step_mem->B = ARKodeButcherTable_LoadERK(itable);

  /* set [redundant] stored values for stage numbers and method orders */
  if (step_mem->B != NULL) {
    step_mem->stages = step_mem->B->stages;
    step_mem->q = step_mem->B->q;
    step_mem->p = step_mem->B->p;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_CheckButcherTable

  This routine runs through the outer Butcher table to ensure
  that it meets all necessary requirements, including:
    strictly lower-triangular (ERK)
    method order q > 0 (all)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
  ---------------------------------------------------------------*/
int mriStep_CheckButcherTable(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeMRIStepMem step_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that MRI table is strictly lower triangular */
  okay = SUNTRUE;
  for (i=0; i<step_mem->stages; i++)
    for (j=i; j<step_mem->stages; j++)
      if (SUNRabs(step_mem->B->A[i][j]) > tol)
        okay = SUNFALSE;
  if (!okay) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "As Butcher table is implicit!");
    return(ARK_ILL_INPUT);
  }

  /* check that stages times are unique and ordered */
  okay = SUNTRUE;
  for (i=1; i<step_mem->stages; i++) {
    if (SUNRabs(step_mem->B->c[i] - step_mem->B->c[i-1]) < tol)
      okay = SUNFALSE;
    else if (step_mem->B->c[i] - step_mem->B->c[i-1] < ZERO)
      okay = SUNFALSE;
  }
  if (!okay) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "Stage times must be unique and ordered.");
    return(ARK_ILL_INPUT);
  }

  /* check that the last stage is not at or beyond the final time */
  okay = SUNTRUE;
  if (SUNRabs(ONE - step_mem->B->c[step_mem->stages-1]) < tol)
    okay = SUNFALSE;
  else if (ONE - step_mem->B->c[step_mem->stages-1] < ZERO)
    okay = SUNFALSE;

  if (!okay) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_CheckButcherTable",
                    "Final stage time must be less than 1.");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_PrepareNextStep

  This routine handles MRI-specific updates following a successful
  step: copying the MRI result to the current solution vector,
  updating the error/step history arrays, and setting the
  prospective step size, hprime, for the next step.  Along with
  hprime, it sets the ratio eta=hprime/h.  It also updates other
  state variables related to a change of step size.
  ---------------------------------------------------------------*/
int mriStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  /* If fixed time-stepping requested, defer
     step size changes until next step */
  if (ark_mem->fixedstep) {
    ark_mem->hprime = ark_mem->h;
    ark_mem->eta = ONE;
    return(ARK_SUCCESS);
  }

  /* Set hprime value for next step size */
  ark_mem->hprime = ark_mem->h * ark_mem->eta;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Inner stepper evolve functions
  ---------------------------------------------------------------*/

int mriStep_EvolveInnerARK(void* inner_mem, realtype t0,
                           N_Vector y0, realtype tout)
{
  ARKodeMem         inner_ark_mem;   /* inner ARKode memory        */
  ARKodeARKStepMem  inner_step_mem;  /* inner stepper memory       */
  realtype          hf, hi;          /* time step sizes            */
  realtype          tspan;           /* integration time span      */
  realtype          tret;            /* return time                */
  int               retval;          /* return value               */

  /* initialize hf to avoid compiler warning */
  hf = ZERO;
  
  /* access the inner ARKStep mem structure */
  retval = arkStep_AccessStepMem(inner_mem,
                                 "mriStep_EvolveInnerARK",
                                 &inner_ark_mem, &inner_step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* if fixed time stepping adjust the step size */
  /* >>>>>>> SHOULD BE ABLE TO REMOVE THIS IF TSTOP IS UPDATED <<<<<<< */
  if (inner_ark_mem->fixedstep) {

    /* save a copy of the inner step size */
    hf = inner_ark_mem->hin;

    /* compute inner step size */
    tspan = tout - t0;
    hi    = tspan / SUNRceil(tspan / hf);

    /* set inner time step size */
    retval = ARKStepSetFixedStep(inner_mem, hi);
    if (retval != ARK_SUCCESS) return(retval);
  }

  /* set stop time */
  retval = ARKStepSetStopTime(inner_mem, tout);
  if (retval != ARK_SUCCESS) return(retval);

  /* evolve inner ODE */
  retval = ARKStepEvolve(inner_mem, tout, y0, &tret, ARK_NORMAL);
  if (retval < 0) return(retval);

  /* if fixed time stepping restore fast step size */
  if (inner_ark_mem->fixedstep) {
    retval = ARKStepSetFixedStep(inner_mem, hf);
    if (retval != ARK_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Wrappers for inner stepper RHS functions
  ---------------------------------------------------------------*/


/* Inner RHS with added forcing */
int mriStep_InnerRhsFnForcing(realtype t, N_Vector y, N_Vector ydot,
                              void *user_data)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access outer integrator memory */
  ark_mem  = (ARKodeMem) user_data;
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* call user function */
  retval = step_mem->ff1(t, y, ydot, ark_mem->user_data);
  if (retval != 0) return(retval);

  /* add forcing from the outer integrator */
  N_VLinearSum(ONE, step_mem->inner_forcing[0], ONE, ydot, ydot);

  /* successfully computed RHS */
  return(0);
}


/* Inner RHS without forcing */
int mriStep_InnerRhsFn(realtype t, N_Vector y, N_Vector ydot,
                       void *user_data)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* access outer integrator memory */
  ark_mem  = (ARKodeMem) user_data;
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* call user function */
  return(step_mem->ff2(t, y, ydot, ark_mem->user_data));
}


/* Inner full RHS without forcing */
int mriStep_InnerFullRhs(void *arkode_mem, realtype t, N_Vector y,
                         N_Vector ydot)
{
  ARKodeMem         ark_mem;         /* outer ARKode memory   */
  ARKodeMRIStepMem  step_mem;        /* outer stepper memory  */
  ARKodeMem         inner_ark_mem;   /* inner ARKode memory   */
  ARKodeARKStepMem  inner_step_mem;  /* inner stepper memory  */
  int               retval;          /* return value          */

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Access ARKStep memory */
  retval = arkStep_AccessStepMem(step_mem->inner_mem, "mriStep_AttachARK",
                                 &inner_ark_mem, &inner_step_mem);
  if (retval != ARK_SUCCESS) { MRIStepFree(&arkode_mem); return(-1); }

  /* call user RHS function */
  if (inner_step_mem->explicit && inner_step_mem->implicit) {

    /* explicit */
    retval = step_mem->ff1(t, y, ark_mem->tempv2, ark_mem->user_data);
    inner_step_mem->nfe++;
    if (retval != 0) return(retval);

    /* implicit */
    retval = step_mem->ff2(t, y, ydot, ark_mem->user_data);
    inner_step_mem->nfi++;
    if (retval != 0) return(retval);

    /* sum fe and fi */
    N_VLinearSum(ONE, ark_mem->tempv2, ONE, ydot, ydot);

  } else if (inner_step_mem->explicit) {

    /* explicit */
    retval = step_mem->ff1(t, y, ydot, ark_mem->user_data);
    inner_step_mem->nfe++;
    if (retval != 0) return(retval);

  } else {

    /* implicit */
    retval = step_mem->ff1(t, y, ydot, ark_mem->user_data);
    inner_step_mem->nfi++;
    if (retval != 0) return(retval);

  }

  /* successfully computed RHS */
  return(ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
