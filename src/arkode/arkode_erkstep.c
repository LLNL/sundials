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
 * This is the implementation file for ARKODE's ERK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode/arkode_butcher.h"
#include "arkode_impl.h"
#include "arkode_erkstep_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>


/*===============================================================
  ERKStep Exported functions -- Required
  ===============================================================*/

void* ERKStepCreate(ARKRhsFn f, realtype t0, N_Vector y0, SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  booleantype nvectorOK;
  int retval;

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  if (!sunctx) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_NULL_SUNCTX);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = erkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeERKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeERKStepMem) malloc(sizeof(struct ARKodeERKStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep",
                    "ERKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeERKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init    = erkStep_Init;
  ark_mem->step_fullrhs = erkStep_FullRHS;
  ark_mem->step         = erkStep_TakeStep;
  ark_mem->step_mem     = (void*) step_mem;

  /* Set default values for ERKStep optional inputs */
  retval = ERKStepSetDefaults((void *) ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::ERKStep", "ERKStepCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Allocate the general ERK stepper vectors using y0 as a template */
  /* NOTE: F, cvals and Xvecs will be allocated later on
     (based on the number of ERK stages) */

  /* Copy the input parameters into ARKODE state */
  step_mem->f = f;

  /* Update the ARKODE workspace requirements -- UPDATE */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* Initialize all the counters */
  step_mem->nfe = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::ERKStep", "ERKStepCreate",
                    "Unable to initialize main ARKODE infrastructure");
    return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
  ERKStepResize:

  This routine resizes the memory within the ERKStep module.
  It first resizes the main ARKODE infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int ERKStepResize(void *arkode_mem, N_Vector y0, realtype hscale,
                  realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepReSize",
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

  /* resize ARKODE infrastructure memory */
  retval = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::ERKStep", "ERKStepResize",
                    "Unable to resize main ARKODE infrastructure");
    return(retval);
  }

  /* Resize the RHS vectors */
  for (i=0; i<step_mem->stages; i++) {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                      liw_diff, y0, &step_mem->F[i])) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep", "ERKStepResize",
                      "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepReInit:

  This routine re-initializes the ERKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE::ERKStep",
                    "ERKStepReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f = f;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::ERKStep", "ERKStepReInit",
                    "Unable to initialize main ARKODE infrastructure");
    return(retval);
  }

  /* Initialize all the counters */
  step_mem->nfe = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepReset:

  This routine resets the ERKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int ERKStepReset(void* arkode_mem, realtype tR, N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepReset",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKODE::ERKStep", "ERKStepReset",
                    "Unable to initialize main ARKODE infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSStolerances, ERKStepSVtolerances, ERKStepWFtolerances:

  These routines set integration tolerances (wrappers for general
  ARKODE utility routines)
  ---------------------------------------------------------------*/
int ERKStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}


int ERKStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}


int ERKStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}


int ERKStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


int ERKStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkEvolve(ark_mem, tout, yout, tret, itask);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}


int ERKStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkGetDky(ark_mem, t, k, dky);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}


/*---------------------------------------------------------------
  ERKStepFree frees all ERKStep memory, and then calls an ARKODE
  utility routine to free the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void ERKStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL ERKStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

    /* free the Butcher table */
    if (step_mem->B != NULL) {
      ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->B);
      step_mem->B = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
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

  /* free memory for overall ARKODE infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
  ERKStepPrintMem:

  This routine outputs the memory from the ERKStep structure and
  the main ARKODE infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void ERKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepPrintMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* output data from main ARKODE infrastructure */
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"ERKStep: q = %i\n", step_mem->q);
  fprintf(outfile,"ERKStep: p = %i\n", step_mem->p);
  fprintf(outfile,"ERKStep: stages = %i\n", step_mem->stages);

  /* output long integer quantities */
  fprintf(outfile,"ERKStep: nfe = %li\n", step_mem->nfe);

  /* output realtype quantities */
  fprintf(outfile,"ERKStep: Butcher table:\n");
  ARKodeButcherTable_Write(step_mem->B, outfile);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  for (i=0; i<step_mem->stages; i++) {
    fprintf(outfile,"ERKStep: F[%i]:\n", i);
    N_VPrintFile(step_mem->F[i], outfile);
  }
#endif
}



/*===============================================================
  ERKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKODE
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  erkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK
    stages, method order, or solver options
  - sets the call_fullrhs flag

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int erkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval, j;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
    return(ARK_SUCCESS);

  /* enforce use of arkEwtSmallReal if using a fixed step size
     and an internal error weight function */
  if ( ark_mem->fixedstep && !ark_mem->user_efun ) {
    ark_mem->user_efun = SUNFALSE;
    ark_mem->efun      = arkEwtSetSmallReal;
    ark_mem->e_data    = ark_mem;
  }

  /* Create Butcher table (if not already set) */
  retval = erkStep_SetButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", "erkStep_Init",
                    "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher table are OK */
  retval = erkStep_CheckButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "erkStep_Init", "Error in Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Retrieve/store method and embedding orders now that table is finalized */
  step_mem->q = ark_mem->hadapt_mem->q = step_mem->B->q;
  step_mem->p = ark_mem->hadapt_mem->p = step_mem->B->p;

  /* Ensure that if adaptivity is enabled, then method includes embedding coefficients */
  if (!ark_mem->fixedstep && (step_mem->p == 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", "erkStep_Init",
                    "Adaptive timestepping cannot be performed without embedding coefficients");
    return(ARK_ILL_INPUT);
  }

  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate F[0] ... F[stages-1] if needed */
  if (step_mem->F == NULL)
    step_mem->F = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
  for (j=0; j<step_mem->stages; j++) {
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->F[j])))
      return(ARK_MEM_FAIL);
  }
  ark_mem->liw += step_mem->stages;  /* pointers */

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

  /* Limit max interpolant degree (negative input only overwrites the current
     interpolant degree if it is greater than abs(input). */
  if (ark_mem->interp != NULL)
  {
    if (step_mem->q > 1)
    {
      /* Limit max degree to at most one less than the method global order */
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q - 1));
    }
    else
    {
      /* Allow for linear interpolant with first order methods to ensure
         solution values are returned at the time interval end points */
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q));
    }

    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", "erkStep_Init",
                      "Unable to update interpolation polynomial degree");
      return (ARK_ILL_INPUT);
    }
  }

  /* Signal to shared arkode module that full RHS evaluations are required */
  ark_mem->call_fullrhs = SUNTRUE;

  return(ARK_SUCCESS);
}


/*------------------------------------------------------------------------------
  erkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function, f(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If this function is called in ARK_FULLRHS_START or ARK_FULLRHS_END mode and
  evaluating the RHS functions is necessary, we store the vector f(t,y) in Fe[0]
  for reuse in the first stage of the subsequent time step.

  In ARK_FULLRHS_END mode we check if the method is "stiffly accurate" and, if
  appropriate, copy the vector F[stages - 1] to F[0] for reuse in the first
  stage of the subsequent time step.

  ARK_FULLRHS_OTHER mode is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to store the
  intermediate parts so that they do not interfere with the other two modes.
  ----------------------------------------------------------------------------*/
int erkStep_FullRHS(void* arkode_mem, realtype t, N_Vector y, N_Vector f,
                    int mode)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  booleantype recomputeRHS;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  case ARK_FULLRHS_START:

    /* compute the RHS */
    if (!(ark_mem->fn_is_current))
    {
      retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                        "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* copy RHS vector into output */
    N_VScale(ONE, step_mem->F[0], f);

    break;

  case ARK_FULLRHS_END:

    /* determine if RHS function needs to be recomputed */
    if (!(ark_mem->fn_is_current))
    {
      recomputeRHS = !ARKodeButcherTable_IsStifflyAccurate(step_mem->B);

      /* First Same As Last methods are not FSAL when relaxation is enabled */
      if (ark_mem->relax_enabled) { recomputeRHS = SUNTRUE; }

      /* base RHS calls on recomputeRHS argument */
      if (recomputeRHS)
      {
        /* call f */
        retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
        step_mem->nfe++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                        "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
          return(ARK_RHSFUNC_FAIL);
        }
      }
      else
      {
        N_VScale(ONE, step_mem->F[step_mem->stages-1], step_mem->F[0]);
      }
    }

    /* copy RHS vector into output */
    N_VScale(ONE, step_mem->F[0], f);

    break;

  case ARK_FULLRHS_OTHER:

    /* call f */
    retval = step_mem->f(t, y, f, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                      "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    break;

  default:
    /* return with RHS failure if unknown mode is passed */
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                    "erkStep_FullRHS", "Unknown full RHS mode");
    return(ARK_RHSFUNC_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_TakeStep:

  This routine serves the primary purpose of the ERKStep module:
  it performs a single ERK step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int erkStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr)
{
  int retval, is, js, nvec, mode;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success */
  *nflagPtr = ARK_SUCCESS;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_TakeStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::erkStep_TakeStep", "start-stage",
                     "step = %li, stage = 0, h = %"RSYM", tcur = %"RSYM,
                     ark_mem->nst, ark_mem->h, ark_mem->tcur);
#endif

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::erkStep_TakeStep", "stage",
                     "z[0] =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::erkStep_TakeStep", "stage RHS",
                     "F[0] =", "");
  N_VPrintFile(step_mem->F[0], ARK_LOGGER->debug_fp);
#endif

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current))
  {
    mode = (ark_mem->initsetup) ? ARK_FULLRHS_START : ARK_FULLRHS_END;
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                   ark_mem->fn, mode);
    if (retval) { return ARK_RHSFUNC_FAIL; }
    ark_mem->fn_is_current = SUNTRUE;
  }

  /* Loop over internal stages to the step; since the method is explicit
     the first stage RHS is just the full RHS from the start of the step */
  for (is=1; is<step_mem->stages; is++) {

    /* Set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + step_mem->B->c[is]*ark_mem->h;

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::erkStep_TakeStep", "start-stage",
                       "step = %li, stage = %i, h = %"RSYM", tcur = %"RSYM,
                       ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

    /* Set ycur to current stage solution */
    nvec = 0;
    for (js=0; js<is; js++) {
      cvals[nvec] = ark_mem->h * step_mem->B->A[is][js];
      Xvecs[nvec] = step_mem->F[js];
      nvec += 1;
    }
    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec += 1;

    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->ProcessStage != NULL) {
      retval = ark_mem->ProcessStage(ark_mem->tcur,
                                     ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) return(ARK_POSTPROCESS_STAGE_FAIL);
    }

    /* compute updated RHS */
    retval = step_mem->f(ark_mem->tcur, ark_mem->ycur,
                         step_mem->F[is], ark_mem->user_data);
    step_mem->nfe++;
    if (retval < 0)  return(ARK_RHSFUNC_FAIL);
    if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::erkStep_TakeStep", "stage RHS",
                       "F[%i] =", is);
    N_VPrintFile(step_mem->F[is], ARK_LOGGER->debug_fp);
#endif

  } /* loop over stages */

  /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
  retval = erkStep_ComputeSolutions(ark_mem, dsmPtr);
  if (retval < 0)  return(retval);

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::erkStep_TakeStep", "updated solution",
                     "ycur =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::erkStep_TakeStep", "error-test",
                     "step = %li, h = %"RSYM", dsm = %"RSYM,
                     ark_mem->nst, ark_mem->h, *dsmPtr);
#endif

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  erkStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int erkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeERKStepMem *step_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    fname, MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeERKStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype erkStep_CheckNVector(N_Vector tmpl)
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
  erkStep_SetButcherTable

  This routine determines the ERK method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
int erkStep_SetButcherTable(ARKodeMem ark_mem)
{
  int etable;
  ARKodeERKStepMem step_mem;
  sunindextype Bliw, Blrw;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_SetButcherTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* if table has already been specified, just return */
  if (step_mem->B != NULL)
    return(ARK_SUCCESS);

  /* initialize table number to illegal values */
  etable = -1;

  /* select method based on order */
  switch (step_mem->q) {
  case(2):
    etable = ERKSTEP_DEFAULT_2;
    break;
  case(3):
    etable = ERKSTEP_DEFAULT_3;
    break;
  case(4):
    etable = ERKSTEP_DEFAULT_4;
    break;
  case(5):
    etable = ERKSTEP_DEFAULT_5;
    break;
  case(6):
    etable = ERKSTEP_DEFAULT_6;
    break;
  case(7):
    etable = ERKSTEP_DEFAULT_7;
    break;
  case(8):
    etable = ERKSTEP_DEFAULT_8;
    break;
  case(9):
    etable = ERKSTEP_DEFAULT_9;
    break;
  default:    /* no available method, set default */
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "erkStep_SetButcherTable",
                    "No explicit method at requested order, using q=9.");
    etable = ERKSTEP_DEFAULT_9;
    break;
  }

  if (etable > -1)
    step_mem->B = ARKodeButcherTable_LoadERK(etable);

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* set [redundant] stored values for stage numbers and method orders */
  if (step_mem->B != NULL) {
    step_mem->stages = step_mem->B->stages;
    step_mem->q = step_mem->B->q;
    step_mem->p = step_mem->B->p;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_CheckButcherTable

  This routine runs through the explicit Butcher table to ensure
  that it meets all necessary requirements, including:
    strictly lower-triangular (ERK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int erkStep_CheckButcherTable(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeERKStepMem step_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable",
                    "stages < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable",
                    "method order < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable",
                    "embedding order < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that embedding exists */
  if ((step_mem->p > 0) && (!ark_mem->fixedstep)) {
    if (step_mem->B->d == NULL) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                      "erkStep_CheckButcherTable",
                      "no embedding!");
      return(ARK_INVALID_TABLE);
    }
  }

  /* check that ERK table is strictly lower triangular */
  okay = SUNTRUE;
  for (i=0; i<step_mem->stages; i++)
    for (j=i; j<step_mem->stages; j++)
      if (SUNRabs(step_mem->B->A[i][j]) > tol)
        okay = SUNFALSE;
  if (!okay) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable",
                    "Ae Butcher table is implicit!");
    return(ARK_INVALID_TABLE);
  }

  /* check if all b values are positive for relaxation */
  if (ark_mem->relax_enabled)
  {
    if (step_mem->q < 2)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                      "erkStep_CheckButcherTables",
                      "The Butcher table must be at least second order!");
      return ARK_INVALID_TABLE;
    }

    for (i = 0; i < step_mem->stages; i++)
    {
      if (step_mem->B->b[i] < ZERO)
      {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKODE::ERKStep",
                        "erkStep_CheckButcherTable",
                        "The Butcher table has a negative b value!");
        return ARK_INVALID_TABLE;
      }
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_tempv1, in case the
  calling routine wishes to examine the error locations.

  Note: at this point in the step, the vector ark_tempv1 may be
  used as a temporary vector.
  ---------------------------------------------------------------*/
int erkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_ComputeSolutions", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;


  /* Compute time step solution */
  /*   set arrays for fused vector operation */
  nvec = 0;
  for (j=0; j<step_mem->stages; j++) {
    cvals[nvec] = ark_mem->h * step_mem->B->b[j];
    Xvecs[nvec] = step_mem->F[j];
    nvec += 1;
  }
  cvals[nvec] = ONE;
  Xvecs[nvec] = ark_mem->yn;
  nvec += 1;

  /*   call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep) {

    /* set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<step_mem->stages; j++) {
      cvals[nvec] = ark_mem->h * (step_mem->B->b[j] - step_mem->B->d[j]);
      Xvecs[nvec] = step_mem->F[j];
      nvec += 1;
    }

    /* call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* fill error norm */
    *dsmPtr = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return(ARK_SUCCESS);
}


/* -----------------------------------------------------------------------------
 * erkStep_RelaxDeltaE
 *
 * Computes the change in the relaxation functions for use in relaxation methods
 * delta_e = h * sum_i b_i * <rjac(z_i), f_i>
 * ---------------------------------------------------------------------------*/

int erkStep_RelaxDeltaE(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                        long int* num_relax_jac_evals, sunrealtype* delta_e_out)
{
  int i, j, nvec, retval;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;
  N_Vector z_stage = ark_mem->tempv2;
  N_Vector J_relax = ark_mem->tempv3;

  /* Access the stepper memory structure */
  if (!(ark_mem->step_mem))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_RelaxDeltaE", MSG_ERKSTEP_NO_MEM);
    return ARK_MEM_NULL;
  }
  step_mem = (ARKodeERKStepMem)(ark_mem->step_mem);

  /* Initialize output */
  *delta_e_out = ZERO;

  /* Set arrays for fused vector operation */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  for (i = 0; i < step_mem->stages; i++)
  {
    /* Construct stages z[i] = y_n + h * sum_j Ae[i,j] Fe[j] + Ai[i,j] Fi[j] */
    nvec = 0;

    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec++;

    for (j = 0; j < i; j++)
    {
      cvals[nvec] = ark_mem->h * step_mem->B->A[i][j];
      Xvecs[nvec] = step_mem->F[j];
      nvec++;
    }

    retval = N_VLinearCombination(nvec, cvals, Xvecs, z_stage);
    if (retval) return ARK_VECTOROP_ERR;

    /* Evaluate the Jacobian at z_i */
    retval = relax_jac_fn(z_stage, J_relax, ark_mem->user_data);
    (*num_relax_jac_evals)++;
    if (retval < 0) { return ARK_RELAX_JAC_FAIL; }
    if (retval > 0) { return ARK_RELAX_JAC_RECV; }

    /* Update estimates */
    if (J_relax->ops->nvdotprodlocal && J_relax->ops->nvdotprodmultiallreduce)
    {
      *delta_e_out += step_mem->B->b[i] * N_VDotProdLocal(J_relax,
                                                          step_mem->F[i]);
    }
    else
    {
      *delta_e_out += step_mem->B->b[i] * N_VDotProd(J_relax,
                                                     step_mem->F[i]);
    }
  }

  if (J_relax->ops->nvdotprodlocal && J_relax->ops->nvdotprodmultiallreduce)
  {
    retval = N_VDotProdMultiAllReduce(1, J_relax, delta_e_out);
    if (retval) { return ARK_VECTOROP_ERR; }
  }

  *delta_e_out *= ark_mem->h;

  return ARK_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * erkStep_GetOrder
 *
 * Returns the method order
 * ---------------------------------------------------------------------------*/

int erkStep_GetOrder(ARKodeMem ark_mem)
{
  ARKodeERKStepMem step_mem = (ARKodeERKStepMem)(ark_mem->step_mem);
  return step_mem->q;
}

/*===============================================================
  EOF
  ===============================================================*/
