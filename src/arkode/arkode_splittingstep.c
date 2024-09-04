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
 * This is the implementation file for ARKODE's operator
 * splitting module
 *--------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <arkode/execution_policy/arkode_execution_policy_serial.h>
#include <sundials/sundials_nvector.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include "arkode_splittingstep_impl.h"

/*---------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem.
  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
static int splittingStep_AccessStepMem(const ARKodeMem ark_mem,
                                       const char* const fname,
                                       ARKodeSplittingStepMem* const step_mem)
{
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
static int splittingStep_AccessARKODEStepMem(void* const arkode_mem,
                                             const char* const fname,
                                             ARKodeMem* const ark_mem,
                                             ARKodeSplittingStepMem* const step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  return splittingStep_AccessStepMem(*ark_mem, __func__, step_mem);
}

/*---------------------------------------------------------------
  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
static sunbooleantype splittingStep_CheckNVector(const N_Vector y)
{
  // TODO(SBR): check all ops are correct
  return y->ops->nvclone != NULL && y->ops->nvdestroy != NULL &&
         y->ops->nvlinearsum != NULL && y->ops->nvscale != NULL;
}

/*---------------------------------------------------------------
  This routine determines the splitting coefficients to use,
  based on the desired accuracy.
  ---------------------------------------------------------------*/
static int splittingStep_SetCoefficients(const ARKodeMem ark_mem,
                                         const ARKodeSplittingStepMem step_mem)
{
  if (step_mem->coefficients != NULL) { return ARK_SUCCESS; }

  if (step_mem->order <= 1)
  {
    /* Lie-Trotter is the default (order < 1) */
    step_mem->coefficients =
      SplittingStepCoefficients_LieTrotter(step_mem->partitions);
  }
  else if (step_mem->order == 3)
  {
    step_mem->coefficients =
      SplittingStepCoefficients_ThirdOrderSuzuki(step_mem->partitions);
  }
  else if (step_mem->order % 2 == 0)
  {
    /* Triple jump only works for even order */
    step_mem->coefficients =
      SplittingStepCoefficients_TripleJump(step_mem->partitions, step_mem->order);
  }
  else
  {
    /* Bump the order up to be even but with an error */
    const int new_order = step_mem->order + 1;
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "No splitting method at requested order, using q=%i.",
                    new_order);
    step_mem->coefficients =
      SplittingStepCoefficients_TripleJump(step_mem->partitions, new_order);
  }

  if (step_mem->coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to allocate splitting coefficients");
    return ARK_MEM_FAIL;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the splitting coefficients to be used
  - sets/checks the execution policy to be used

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
static int splittingStep_Init(const ARKodeMem ark_mem, const int init_type)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return ARK_SUCCESS; }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    retval = splittingStep_SetCoefficients(ark_mem, step_mem);
    if (retval != ARK_SUCCESS) { return retval; }

    if (step_mem->policy == NULL)
    {
      step_mem->policy = ARKodeSplittingExecutionPolicy_New_Serial();
      if (step_mem->policy == NULL)
      {
        if (step_mem->coefficients == NULL)
        {
          arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                          "Failed to allocate execution policy");
          return ARK_MEM_FAIL;
        }
      }
      step_mem->own_policy = SUNTRUE;
    }
  }

  ark_mem->interp_degree =
    SUNMAX(1, SUNMIN(step_mem->order - 1, ark_mem->interp_degree));

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y) + ... + f^P(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  In SplittingStep, we accumulate the RHS functions in ARK_FULLRHS_OTHER mode.
  Generally, inner steppers will not have the correct yn when this function is
  called and will not be able to reuse a function evaluation since their state
  resets at the next ARKodeEvolve call.
  ----------------------------------------------------------------------------*/
static int splittingStep_FullRHS(const ARKodeMem ark_mem, const sunrealtype t,
                                 const N_Vector y, const N_Vector f,
                                 SUNDIALS_MAYBE_UNUSED const int mode)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = step_mem->steppers[0]->ops->fullrhs(step_mem->steppers[0], t, y, f,
                                               ARK_FULLRHS_OTHER);
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return (ARK_RHSFUNC_FAIL);
  }

  for (int i = 1; i < step_mem->partitions; i++)
  {
    retval = step_mem->steppers[i]->ops->fullrhs(step_mem->steppers[i], t, y,
                                                 ark_mem->tempv1,
                                                 ARK_FULLRHS_OTHER);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }
    N_VLinearSum(SUN_RCONST(1.0), f, SUN_RCONST(1.0), ark_mem->tempv1, f);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine performs a sequential operator splitting method
  ---------------------------------------------------------------*/
static int splittingStep_SequentialMethod(const int i, const N_Vector y,
                                          void* const user_data)
{
  ARKodeMem ark_mem               = (ARKodeMem)user_data;
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  const SplittingStepCoefficients coefficients = step_mem->coefficients;

  for (int j = 0; j < coefficients->stages; j++)
  {
    for (int k = 0; k < coefficients->partitions; k++)
    {
      const sunrealtype t_start = ark_mem->tn +
                                  coefficients->beta[i][j][k] * ark_mem->h;
      const sunrealtype t_end = ark_mem->tn +
                                coefficients->beta[i][j + 1][k] * ark_mem->h;

      if (t_start == t_end) { continue; }

      const SUNStepper stepper = step_mem->steppers[k];
      SUNErrCode err           = stepper->ops->reset(stepper, t_start, y);
      if (err != SUN_SUCCESS) { return err; }

      err = stepper->ops->setstoptime(stepper, t_end);
      if (err != SUN_SUCCESS) { return err; }

      sunrealtype tret = 0;
      int stop_reason  = 0;
      err = SUNStepper_Evolve(stepper, t_start, t_end, y, &tret, &stop_reason);
      step_mem->n_stepper_evolves++;

      if (err != SUN_SUCCESS) { return err; }
      if (stop_reason < 0) { return stop_reason; }
    }
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine performs a single step of the splitting method.
  ---------------------------------------------------------------*/
static int splittingStep_TakeStep(const ARKodeMem ark_mem,
                                  sunrealtype* const dsmPtr, int* const nflagPtr)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
  *dsmPtr   = ZERO;        /* No error estimate */

  const SplittingStepCoefficients coefficients = step_mem->coefficients;

  return step_mem->policy->execute(step_mem->policy,
                                   splittingStep_SequentialMethod, ark_mem->yn,
                                   ark_mem->ycur, ark_mem->tempv1,
                                   coefficients->alpha,
                                   coefficients->sequential_methods, ark_mem);
}

/*---------------------------------------------------------------
  Prints integrator statistics
  ---------------------------------------------------------------*/
static int splittingStep_PrintAllStats(const ARKodeMem ark_mem,
                                       FILE* const outfile,
                                       const SUNOutputFormat fmt)
{
  // TODO(SBR): update when https://github.com/LLNL/sundials/pull/517 merged
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
static int splittingStep_WriteParameters(const ARKodeMem ark_mem, FILE* const fp)
{
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  fprintf(fp, "SplittingStep time step module parameters:\n  Method order %i\n\n",
          step_mem->order);

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  This routine resizes the memory within the SplittingStep module.
  ---------------------------------------------------------------*/
static int splittingStep_Resize(SUNDIALS_MAYBE_UNUSED const ARKodeMem ark_mem,
                                SUNDIALS_MAYBE_UNUSED const N_Vector ynew,
                                SUNDIALS_MAYBE_UNUSED const sunrealtype hscale,
                                SUNDIALS_MAYBE_UNUSED const sunrealtype t0,
                                SUNDIALS_MAYBE_UNUSED const ARKVecResizeFn resize,
                                SUNDIALS_MAYBE_UNUSED void* const resize_data)
{
  /* Nothing to do since the step_mem has no vectors. Users are responsible for
   * resizing the SUNSteppers. */
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Frees all SplittingStep memory.
  ---------------------------------------------------------------*/
static void splittingStep_Free(const ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  if (step_mem == NULL) { return; }

  free(step_mem->steppers);
  if (step_mem->own_policy)
  {
    ARKodeSplittingExecutionPolicy_Free(&step_mem->policy);
  }

  SplittingStepCoefficients_Free(step_mem->coefficients);
  free(step_mem);
  ark_mem->step_mem = NULL;
}

/*---------------------------------------------------------------
  This routine outputs the memory from the SplittingStep
  structure to a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
static void splittingStep_PrintMem(const ARKodeMem ark_mem, FILE* const outfile)
{
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output integer quantities */
  fprintf(outfile, "SplittingStep: partitions = %i\n", step_mem->partitions);
  fprintf(outfile, "SplittingStep: order = %i\n", step_mem->order);

  /* output long integer quantities */
  fprintf(outfile, "SplittingStep: n_stepper_evolves = %li\n",
          step_mem->n_stepper_evolves);

  /* output sunrealtype quantities */
  fprintf(outfile, "SplittingStep: Coefficients:\n");
  SplittingStepCoefficients_Write(step_mem->coefficients, outfile);
}

/*---------------------------------------------------------------
  Specifies the method order
  ---------------------------------------------------------------*/
static int splittingStep_SetOrder(const ARKodeMem ark_mem, const int order)
{
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user-provided value, or default, depending on argument */
  step_mem->order = order <= 0 ? 1 : order;

  SplittingStepCoefficients_Free(step_mem->coefficients);
  step_mem->coefficients = NULL;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Resets all SplittingStep optional inputs to their default
  values. Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
static int splittingStep_SetDefaults(const ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = splittingStep_SetOrder(ark_mem, 0);
  if (retval != ARK_SUCCESS) { return retval; }

  if (step_mem->own_policy) { free(step_mem->policy); }
  step_mem->own_policy = SUNFALSE;

  ark_mem->interp_type = ARK_INTERP_LAGRANGE;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Creates the SplittingStep integrator
  ---------------------------------------------------------------*/
void* SplittingStepCreate(SUNStepper* const steppers, const int partitions,
                          const sunrealtype t0, const N_Vector y0,
                          const SUNContext sunctx)
{
  if (steppers == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "steppers = NULL illegal.");
    return NULL;
  }

  for (int i = 0; i < partitions; i++)
  {
    if (steppers[i] == NULL)
    {
      arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "steppers[i] = NULL illegal.", i);
      return NULL;
    }
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
  if (!splittingStep_CheckNVector(y0))
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

  const ARKodeSplittingStepMem step_mem =
    (ARKodeSplittingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->steppers = malloc(partitions * sizeof(*steppers));
  if (step_mem->steppers == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  memcpy(step_mem->steppers, steppers, partitions * sizeof(*steppers));

  step_mem->coefficients      = NULL;
  step_mem->policy            = NULL;
  step_mem->n_stepper_evolves = 0;
  step_mem->partitions        = partitions;
  step_mem->order             = 0;
  step_mem->own_policy        = SUNFALSE;

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init            = splittingStep_Init;
  ark_mem->step_fullrhs         = splittingStep_FullRHS;
  ark_mem->step                 = splittingStep_TakeStep;
  ark_mem->step_printallstats   = splittingStep_PrintAllStats;
  ark_mem->step_writeparameters = splittingStep_WriteParameters;
  ark_mem->step_resize          = splittingStep_Resize;
  ark_mem->step_free            = splittingStep_Free;
  ark_mem->step_printmem        = splittingStep_PrintMem;
  ark_mem->step_setdefaults     = splittingStep_SetDefaults;
  ark_mem->step_setorder        = splittingStep_SetOrder;
  ark_mem->step_mem             = (void*)step_mem;

  /* Set default values for ARKStep optional inputs */
  int retval = splittingStep_SetDefaults(ark_mem);
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

/*---------------------------------------------------------------
  Sets the SplittingStep coefficients.
  ---------------------------------------------------------------*/
int SplittingStep_SetCoefficients(void* arkode_mem,
                                  SplittingStepCoefficients coefficients)
{
  ARKodeMem ark_mem               = NULL;
  ARKodeSplittingStepMem step_mem = NULL;
  const int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__,
                                                       &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Splitting coefficients must be non-NULL");
    return ARK_ILL_INPUT;
  }

  if (step_mem->partitions != coefficients->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The splitting method has %i partitions but the coefficients have %i.",
                    step_mem->partitions, coefficients->partitions);
    return ARK_ILL_INPUT;
  }

  SplittingStepCoefficients_Free(step_mem->coefficients);
  step_mem->coefficients = SplittingStepCoefficients_Copy(coefficients);
  if (step_mem->coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to copy splitting coefficients");
    return (ARK_MEM_NULL);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Sets the exection policy.
  ---------------------------------------------------------------*/
int SplittingStep_SetExecutionPolicy(void* arkode_mem,
                                     ARKodeSplittingExecutionPolicy policy)
{
  ARKodeMem ark_mem               = NULL;
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (policy == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The execution policy must be non-NULL");
    return ARK_ILL_INPUT;
  }

  if (step_mem->own_policy)
  {
    ARKodeSplittingExecutionPolicy_Free(&step_mem->policy);
  }
  step_mem->policy     = policy;
  step_mem->own_policy = SUNFALSE;

  return ARK_SUCCESS;
}
