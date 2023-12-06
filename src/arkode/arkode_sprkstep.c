/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the implementation file for ARKODE's SPRK time stepper
 * module.
 *--------------------------------------------------------------*/

#include "arkode/arkode_sprkstep.h"

#include <arkode/arkode.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_sprkstep_impl.h"

/*===============================================================
  SPRKStep Exported functions -- Required
  ===============================================================*/

void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0, N_Vector y0,
                     SUNContext sunctx)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  sunbooleantype nvectorOK   = 0;
  int retval                 = 0;

  /* Check that f1 and f2 are supplied */
  if (!f1)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  if (!f2)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  /* Check for legal input parameters */
  if (!y0)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (NULL);
  }

  if (!sunctx)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return (NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = sprkStep_CheckNVector(y0);
  if (!nvectorOK)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_NVECTOR);
    return (NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (NULL);
  }

  /* Allocate ARKodeSPRKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeSPRKStepMem)malloc(sizeof(struct ARKodeSPRKStepMemRec));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    return (NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeSPRKStepMemRec));

  /* Allocate vectors in stepper mem */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata)))
  {
    SPRKStepFree((void**)&ark_mem);
    return (NULL);
  }

  if (ark_mem->use_compensated_sums)
  {
    if (!arkAllocVec(ark_mem, y0, &(step_mem->yerr)))
    {
      SPRKStepFree((void**)&ark_mem);
      return (NULL);
    }
  }
  else { step_mem->yerr = NULL; }
  ark_mem->step_init    = sprkStep_Init;
  ark_mem->step_fullrhs = sprkStep_FullRHS;
  ark_mem->step         = sprkStep_TakeStep;
  ark_mem->step_mem     = (void*)step_mem;

  /* Set default values for SPRKStep optional inputs */
  retval = SPRKStepSetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    SPRKStepFree((void**)&ark_mem);
    return (NULL);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Initialize the counters */
  step_mem->nf1    = 0;
  step_mem->nf2    = 0;
  step_mem->istage = 0;

  /* Zero yerr for compensated summation */
  if (ark_mem->use_compensated_sums) { N_VConst(ZERO, step_mem->yerr); }

  /* SPRKStep uses Lagrange interpolation by default, since Hermite is
     less compatible with these methods. */
  arkSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    SPRKStepFree((void**)&ark_mem);
    return (NULL);
  }

  return ((void*)ark_mem);
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
int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,
                   N_Vector y0)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepReInit", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that f1 and f2 are supplied */
  if (!f1 || !f2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (!y0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to reinitialize main ARKODE infrastructure");
    return (retval);
  }

  /* Initialize the counters */
  step_mem->nf1    = 0;
  step_mem->nf2    = 0;
  step_mem->istage = 0;

  /* Zero yerr for compensated summation */
  N_VConst(ZERO, step_mem->yerr);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepReset:

  This routine resets the SPRKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStepReset", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return (retval);
  }

  N_VConst(SUN_RCONST(0.0), step_mem->yerr);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  SPRKStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKODE utility routine)
  ---------------------------------------------------------------*/
int SPRKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                   sunrealtype* tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem = NULL;
  int retval        = 0;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkEvolve(ark_mem, tout, yout, tret, itask);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return (retval);
}

/*---------------------------------------------------------------
  SPRKStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKODE utility routine)
  ---------------------------------------------------------------*/
int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem = NULL;
  int retval        = 0;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkGetDky(ark_mem, t, k, dky);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return (retval);
}

/*---------------------------------------------------------------
  SPRKStepFree frees all SPRKStep memory, and then calls an ARKODE
  utility routine to free the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void SPRKStepFree(void** arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL) { return; }

  /* conditional frees on non-NULL SPRKStep module */
  ark_mem = (ARKodeMem)(*arkode_mem);
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeSPRKStepMem)ark_mem->step_mem;

    if (step_mem->sdata != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }

    if (step_mem->yerr != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->yerr);
      step_mem->yerr = NULL;
    }

    ARKodeSPRKTable_Free(step_mem->method);

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

  With initialization type FIRST_INIT or RESIZE_INIT, this routine
  this routines loads the default method of the selected order
  if necessary.

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
int sprkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_Init", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return (ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    if (!step_mem->method)
    {
      switch (step_mem->q)
      {
      case 1:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_1);
        break;
      case 2:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_2);
        break;
      case 3:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_3);
        break;
      case 4:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4);
        break;
      case 5:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_5);
        break;
      case 6:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_6);
        break;
      case 7:
      case 8:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_8);
        break;
      case 9:
      case 10:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_10);
        break;
      default:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4);
        break;
      }
    }
  }

  /* Limit max interpolant degree (negative input only overwrites the current
     interpolant degree if it is greater than abs(input). */
  if (ark_mem->interp != NULL)
  {
    if (step_mem->method->q > 1)
    {
      /* Limit max degree to at most one less than the method global order */
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp,
                                  -(step_mem->method->q - 1));
    }
    else
    {
      /* Allow for linear interpolant with first order methods to ensure
         solution values are returned at the time interval end points */
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp,
                                  -(step_mem->method->q));
    }

    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Unable to update interpolation polynomial degree");
      return (ARK_ILL_INPUT);
    }
  }

  return (ARK_SUCCESS);
}

int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem = NULL;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkRootInit(ark_mem, nrtfn, g));
}

/* Utility to call f1 and increment the counter */
SUNDIALS_INLINE
int sprkStep_f1(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur,
                N_Vector f1, void* user_data)
{
  int retval = step_mem->f1(tcur, ycur, f1, user_data);
  step_mem->nf1++;
  return retval;
}

/* Utility to call f2 and increment the counter */
SUNDIALS_INLINE
int sprkStep_f2(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur,
                N_Vector f2, void* user_data)
{
  int retval = step_mem->f2(tcur, ycur, f2, user_data);
  step_mem->nf2++;
  return retval;
}

/*------------------------------------------------------------------------------
  sprkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS,
  f1(t,y) + f2(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  Since RHS values are not stored in SPRKStep we evaluate the RHS functions for
  all modes.
  ----------------------------------------------------------------------------*/
int sprkStep_FullRHS(void* arkode_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode)
{
  int retval                 = 0;
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "SPRKStep_FullRHS", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:
  case ARK_FULLRHS_END:
  case ARK_FULLRHS_OTHER:

    /* Since f1 and f2 do not have overlapping outputs and so the f vector is
       passed to both RHS functions. */

    retval = sprkStep_f1(step_mem, t, y, f, ark_mem->user_data);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }

    retval = sprkStep_f2(step_mem, t, y, f, ark_mem->user_data);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }

    break;

  default:
    /* return with RHS failure if unknown mode is passed */
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    "Unknown full RHS mode");
    return (ARK_RHSFUNC_FAIL);
  }

  return (ARK_SUCCESS);
}

/* Standard formulation of SPRK.
   This requires only 2 vectors in principle, but we use three
   since we persist the stage data. Only the stage data vector
   belongs to SPRKStep, the other two are reused from the ARKODE core. */
int sprkStep_TakeStep(void* arkode_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  N_Vector prev_stage        = NULL;
  N_Vector curr_stage        = NULL;
  sunrealtype ci             = SUN_RCONST(0.0);
  sunrealtype chati          = SUN_RCONST(0.0);
  int is                     = 0;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_TakeStep", &ark_mem,
                                  &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  prev_stage = ark_mem->yn;
  curr_stage = ark_mem->ycur;
  for (is = 0; is < step_mem->method->stages; is++)
  {
    /* load/compute coefficients */
    sunrealtype ai    = step_mem->method->a[is];
    sunrealtype ahati = step_mem->method->ahat[is];

    ci += ai;
    chati += ahati;

    /* store current stage index */
    step_mem->istage = is;

    /* evaluate p' with the previous velocity */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                        set other outputs to zero */
    retval = sprkStep_f1(step_mem, ark_mem->tn + chati * ark_mem->h, prev_stage,
                         step_mem->sdata, ark_mem->user_data);
    if (retval != 0) { return ARK_RHSFUNC_FAIL; }

    /* position update */
    N_VLinearSum(ONE, prev_stage, ark_mem->h * ahati, step_mem->sdata,
                 curr_stage);

    /* set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + chati * ark_mem->h;

    /* evaluate q' with the current positions */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                        set other outputs to zero */
    retval = sprkStep_f2(step_mem, ark_mem->tn + ci * ark_mem->h, curr_stage,
                         step_mem->sdata, ark_mem->user_data);
    if (retval != 0) { return ARK_RHSFUNC_FAIL; }

    /* velocity update */
    N_VLinearSum(ONE, curr_stage, ark_mem->h * ai, step_mem->sdata, curr_stage);

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->ProcessStage != NULL)
    {
      retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) { return (ARK_POSTPROCESS_STAGE_FAIL); }
    }

    /* keep track of the stage number */
    step_mem->istage++;

    prev_stage = curr_stage;
  }

  *nflagPtr = 0;
  *dsmPtr   = 0;

  return ARK_SUCCESS;
}

/* Increment SPRK algorithm with compensated summation.
   This algorithm requires 6 vectors, but 5 of them are reused
   from the ARKODE core. */
int sprkStep_TakeStep_Compensated(void* arkode_mem, sunrealtype* dsmPtr,
                                  int* nflagPtr)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  N_Vector delta_Yi          = NULL;
  N_Vector yn_plus_delta_Yi  = NULL;
  N_Vector diff              = NULL;
  sunrealtype ci             = SUN_RCONST(0.0);
  sunrealtype chati          = SUN_RCONST(0.0);
  int is                     = 0;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(arkode_mem, "sprkStep_TakeStep_SPRK",
                                  &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Vector shortcuts */
  delta_Yi         = ark_mem->tempv1;
  yn_plus_delta_Yi = ark_mem->tempv2;
  diff             = ark_mem->tempv3;

  /* [ \Delta P_0 ] = [ 0 ]
     [ \Delta Q_0 ] = [ 0 ] */
  N_VConst(ZERO, delta_Yi);

  /* loop over internal stages to the step */
  for (is = 0; is < step_mem->method->stages; is++)
  {
    /* load/compute coefficients */
    sunrealtype ai    = step_mem->method->a[is];
    sunrealtype ahati = step_mem->method->ahat[is];

    ci += ai;
    chati += ahati;

    /* store current stage index */
    step_mem->istage = is;

    /* [     ] + [            ]
       [ q_n ] + [ \Delta Q_i ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    /* Evaluate p' with the previous velocity */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                        set other outputs to zero */
    retval = sprkStep_f1(step_mem, ark_mem->tn + chati * ark_mem->h,
                         yn_plus_delta_Yi, step_mem->sdata, ark_mem->user_data);
    if (retval != 0) { return (ARK_RHSFUNC_FAIL); }

    /* Incremental position update:
       [ \Delta P_i ] = [ \Delta P_{i-1} ] + [ sdata ]
       [            ] = [                ] + [       ] */
    N_VLinearSum(ONE, delta_Yi, ark_mem->h * ahati, step_mem->sdata, delta_Yi);

    /* [ p_n ] + [ \Delta P_i ]
       [     ] + [            ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    /* set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + chati * ark_mem->h;

    /* Evaluate q' with the current positions */
    N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                        set other outputs to zero */
    retval = sprkStep_f2(step_mem, ark_mem->tn + ci * ark_mem->h,
                         yn_plus_delta_Yi, step_mem->sdata, ark_mem->user_data);
    if (retval != 0) { return (ARK_RHSFUNC_FAIL); }

    /* Incremental velocity update:
       [            ] = [                ] + [       ]
       [ \Delta Q_i ] = [ \Delta Q_{i-1} ] + [ sdata ] */
    N_VLinearSum(ONE, delta_Yi, ark_mem->h * ai, step_mem->sdata, delta_Yi);

    /* if user-supplied stage postprocessing function, we error out since it
     * wont work with the increment form */
    if (ark_mem->ProcessStage != NULL)
    {
      arkProcessError(ark_mem, ARK_POSTPROCESS_STAGE_FAIL, __LINE__, __func__,
                      __FILE__,
                      "Compensated summation is not compatible with stage "
                      "PostProcessing!\n");
      return (ARK_POSTPROCESS_STAGE_FAIL);
    }
  }

  /*
    Now we compute the step solution via compensated summation.
     [ p_{n+1} ] = [ p_n ] + [ \Delta P_i ]
     [ q_{n+1} ] = [ q_n ] + [ \Delta Q_i ] */
  N_VLinearSum(ONE, delta_Yi, -ONE, step_mem->yerr, delta_Yi);
  N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->yn, diff);
  N_VLinearSum(ONE, diff, -ONE, delta_Yi, step_mem->yerr);

  *nflagPtr = 0;
  *dsmPtr   = SUN_RCONST(0.0);

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
int sprkStep_AccessStepMem(void* arkode_mem, const char* fname,
                           ARKodeMem* ark_mem, ARKodeSPRKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_SPRKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSPRKStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  sprkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype sprkStep_CheckNVector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvscale == NULL) || (tmpl->ops->nvwrmsnorm == NULL))
  {
    return (SUNFALSE);
  }
  return (SUNTRUE);
}
