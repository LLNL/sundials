/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for ARKODE's forcing method
 *--------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_forcingstep.h>
#include <sundials/sundials_nvector.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include "arkode_forcingstep_impl.h"

/*---------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem.
  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
static int forcingStep_AccessStepMem(const ARKodeMem ark_mem,
                                       const char* const fname,
                                       ARKodeForcingStepMem* const step_mem)
{
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeForcingStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
static sunbooleantype forcingStep_CheckNVector(const N_Vector y)
{
  return y->ops->nvlinearsum != NULL;
}

/*---------------------------------------------------------------
  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets the forcing for stepper2

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
static int forcingStep_Init(const ARKodeMem ark_mem, const int init_type)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* assume fixed outer step size */
  if (!ark_mem->fixedstep)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Adaptive outer time stepping is not currently supported");
    return ARK_ILL_INPUT;
  }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return ARK_SUCCESS;
  }

  SUNStepper_SetForcing(step_mem->stepper2, 1, ark_mem->yn);
  ark_mem->interp_degree = 1;

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  In ForcingStep, we accumulate the RHS functions in ARK_FULLRHS_OTHER mode.
  Generally, inner steppers will not have the correct yn when this function is
  called and will not be able to reuse a function evaluation since their state
  resets at the next ARKodeEvolve call.
  ----------------------------------------------------------------------------*/
static int forcingStep_FullRHS(const ARKodeMem ark_mem, const sunrealtype t,
                                 const N_Vector y, const N_Vector f,
                                 SUNDIALS_MAYBE_UNUSED const int mode)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = step_mem->stepper1->ops->fullrhs(step_mem->stepper1, t, y, f,
                                               ARK_FULLRHS_OTHER);
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return (ARK_RHSFUNC_FAIL);
  }

  retval = step_mem->stepper2->ops->fullrhs(step_mem->stepper2, t, y,
                                            ark_mem->tempv1,
                                            ARK_FULLRHS_OTHER);
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return (ARK_RHSFUNC_FAIL);
  }
  N_VLinearSum(SUN_RCONST(1.0), f, SUN_RCONST(1.0), ark_mem->tempv1, f);

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine performs a single step of the forcing method.
  ---------------------------------------------------------------*/
static int forcingStep_TakeStep(const ARKodeMem ark_mem,
                                  sunrealtype* const dsmPtr, int* const nflagPtr)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
  *dsmPtr   = ZERO;        /* No error estimate */

  const SUNStepper s1 = step_mem->stepper1;
  const sunrealtype tout = ark_mem->tn + ark_mem->h;
  sunrealtype tret = 0;
  int stop_reason = 0;

  /* Evolve partition 1 on its own */
  SUNErrCode err = SUNStepper_Reset(s1, ark_mem->tn, ark_mem->yn);
  if (err != SUN_SUCCESS) { return err; }
  err = SUNStepper_SetStopTime(s1, tout);
  if (err != SUN_SUCCESS) { return err; }
  err = SUNStepper_Evolve(s1, ark_mem->tn, tout, ark_mem->ycur, &tret, &stop_reason);
  if (err != SUN_SUCCESS) { return err; }
  if (stop_reason < 0) { return stop_reason; }
  step_mem->n_stepper_evolves++;

  /* Write tendency into stepper 2 forcing */
  const SUNStepper s2 = step_mem->stepper2;
  const sunrealtype hinv = SUN_RCONST(1.0) / ark_mem->h;
  N_VLinearSum(hinv, ark_mem->ycur, -hinv, ark_mem->yn, s2->forcing[0]);

  /* No need to reset stepper 2 since its trajectory is continuous and in sync
   * with the outer forcing method */
  err = SUNStepper_SetStopTime(s2, tout);
  if (err != SUN_SUCCESS) { return err; }

  /* Eolve partition 2 with the forcing */
  err = SUNStepper_Evolve(s2, ark_mem->tn, tout, ark_mem->ycur, &tret, &stop_reason);
  if (err != SUN_SUCCESS) { return err; }
  if (stop_reason < 0) { return stop_reason; }
  step_mem->n_stepper_evolves++;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Prints integrator statistics
  ---------------------------------------------------------------*/
static int forcingStep_PrintAllStats(const ARKodeMem ark_mem,
                                       FILE* const outfile,
                                       const SUNOutputFormat fmt)
{
  // TODO(SBR): update when https://github.com/LLNL/sundials/pull/517 merged
  ARKodeForcingStepMem step_mem = NULL;
  const int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Stepper evolves              = %ld\n",
            step_mem->n_stepper_evolves);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, "Stepper evolves,%ld\n", step_mem->n_stepper_evolves);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Frees all ForcingStep memory.
  ---------------------------------------------------------------*/
static void forcingStep_Free(const ARKodeMem ark_mem)
{
  free(ark_mem->step_mem);
  ark_mem->step_mem = NULL;
}

/*---------------------------------------------------------------
  This routine outputs the memory from the ForcingStep
  structure to a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
static void forcingStep_PrintMem(const ARKodeMem ark_mem, FILE* const outfile)
{
  ARKodeForcingStepMem step_mem = NULL;
  const int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output long integer quantities */
  fprintf(outfile, "ForcingStep: n_stepper_evolves = %li\n",
          step_mem->n_stepper_evolves);
}

/*---------------------------------------------------------------
  Resets all ForcingStep optional inputs to their default
  values. Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
static int forcingStep_SetDefaults(const ARKodeMem ark_mem)
{
  ark_mem->interp_type = ARK_INTERP_LAGRANGE;

  return ARK_SUCCESS;
}

static int forcingStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* Reset stepper 2 to be consistent. No need to reset stepper 1 since that
   * happens at the start of a step anyway. */
  retval = SUNStepper_Reset(step_mem->stepper2, tR, yR);
  if (retval != ARK_SUCCESS) { return ARK_INNERSTEP_FAIL; }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Creates the ForcingStep integrator
  ---------------------------------------------------------------*/
void* ForcingStepCreate(SUNStepper stepper1,
                        SUNStepper stepper2,
                        sunrealtype t0, N_Vector y0,
                        SUNContext sunctx)
{
  if (stepper1 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper1 = NULL illegal.");
    return NULL;
  }

  if (stepper2 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper2 = NULL illegal.");
    return NULL;
  }

  if (y0 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return NULL;
  }

  if (sunctx == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return NULL;
  }

  /* Test if all required vector operations are implemented */
  if (!forcingStep_CheckNVector(y0))
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_NVECTOR);
    return (NULL);
  }

  /* Create ark_mem structure and set default values */
  ARKodeMem ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (NULL);
  }

  const ARKodeForcingStepMem step_mem =
    (ARKodeForcingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->stepper1 = stepper1;
  step_mem->stepper2 = stepper2;
  step_mem->n_stepper_evolves = 0;

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init            = forcingStep_Init;
  ark_mem->step_fullrhs         = forcingStep_FullRHS;
  ark_mem->step                 = forcingStep_TakeStep;
  ark_mem->step_printallstats   = forcingStep_PrintAllStats;
  ark_mem->step_free            = forcingStep_Free;
  ark_mem->step_printmem        = forcingStep_PrintMem;
  ark_mem->step_setdefaults     = forcingStep_SetDefaults;
  ark_mem->step_reset           = forcingStep_Reset;
  ark_mem->step_mem             = (void*)step_mem;

  /* Set default values for ARKStep optional inputs */
  int retval = forcingStep_SetDefaults(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  return ark_mem;
}
