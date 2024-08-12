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
 * TODO
 *--------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <sundials/sundials_nvector.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include "arkode_splittingstep_impl.h"

#define DEFAULT_ORDER 2

// TODO: check the 2nd argument of arkProcessError is correct
// TODO: workspace size

static int splittingStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                                       ARKodeSplittingStepMem* step_mem)
{
  /* access ARKodeSPRKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

static int splittingStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                             ARKodeMem* ark_mem,
                                             ARKodeSplittingStepMem* step_mem)
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

static sunbooleantype splittingStep_CheckNVector(const N_Vector y)
{
  // TODO: check all ops are correct
  return y->ops->nvclone != NULL && y->ops->nvdestroy != NULL &&
         y->ops->nvlinearsum != NULL && y->ops->nvscale != NULL;
}

static int splittingStep_Init(ARKodeMem ark_mem, const int init_type)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return ARK_SUCCESS; }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    if (step_mem->coefficients == NULL)
    {
      if (step_mem->order <= 1)
      {
        step_mem->coefficients =
          ARKodeSplittingCoefficients_LieTrotter(step_mem->partitions);
        // TODO: check if non-NULL?
      }
      else
      {
        // TODO consider extrapolation for odd orders
        const int next_even_order = step_mem->order + (step_mem->order % 2);
        step_mem->coefficients =
          ARKodeSplittingCoefficients_TripleJump(step_mem->partitions,
                                                 next_even_order);
        // TODO: check if non-NULL?
      }
    }

    if (step_mem->policy == NULL)
    {
      step_mem->policy     = ARKodeSplittingExecutionPolicy_Serial();
      step_mem->own_policy = SUNTRUE;
    }
  }

  ark_mem->interp_degree =
    SUNMAX(1, SUNMIN(step_mem->order - 1, ark_mem->interp_degree));

  return ARK_SUCCESS;
}

static int splittingStep_FullRHS(ARKodeMem ark_mem, const sunrealtype t,
                                 const N_Vector y, const N_Vector f,
                                 const int mode)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  // retval = innerStepper_FullRhs(step_mem->steppers[0], t, y, f,
  //                               ARK_FULLRHS_OTHER);
  if (retval != 0)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return (ARK_RHSFUNC_FAIL);
  }

  for (int i = 1; i < step_mem->partitions; i++)
  {
    // retval = innerStepper_FullRhs(step_mem->steppers[i], t, y, ark_mem->tempv1,
    //                               ARK_FULLRHS_OTHER);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }
    N_VLinearSum(ONE, f, ONE, ark_mem->tempv1, f);
  }

  return ARK_SUCCESS;
}

static int splittingStep_Stage(const int i, const N_Vector y, void* user_data)
{
  // Note: this should not write to user_data to avoid race conditions when executed in parallel

  ARKodeMem ark_mem               = (ARKodeMem)user_data;
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  ARKodeSplittingCoefficients coefficients = step_mem->coefficients;

  for (int j = 0; j < coefficients->stages; j++)
  {
    for (int k = 0; k < coefficients->partitions; k++)
    {
      const sunrealtype tcur = ark_mem->tn +
                               coefficients->beta[i][j][k] * ark_mem->h;
      const sunrealtype tnext = ark_mem->tn + coefficients->beta[i][j + 1][k] *
                                ark_mem->h;
      if (tcur != tnext)
      {
        SUNStepper stepper = step_mem->steppers[k];
        if (SUNTRUE)
        { // TODO: condition for FSAL. Consider special cases where the
          // coefficients have changed, it's the first step, or reset is called
          // TODO: what about backward integration?
          stepper->ops->reset(stepper, tcur, y);
        }
        retval = stepper->ops->evolve(stepper, tcur, tnext, y);

        if (retval != ARK_SUCCESS) { return retval; }
      }
    }
  }

  return ARK_SUCCESS;
}

static int splittingStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                  int* nflagPtr)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  ARKodeSplittingCoefficients coefficients = step_mem->coefficients;

  return step_mem->policy->execute(step_mem->policy, splittingStep_Stage,
                                     ark_mem->yn, ark_mem->ycur,
                                     ark_mem->tempv1, coefficients->alpha,
                                     coefficients->sequential_methods, ark_mem);
}

int splittingStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile,
                                SUNOutputFormat fmt)
{
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Stepper evolves = %ld\n", step_mem->n_stepper_evolves);
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

int splittingStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  fprintf(fp, "SplittingStep time step module parameters:\n  Method order %i\n\n",
          step_mem->order);

  return ARK_SUCCESS;
}

static int splittingStep_Resize(ARKodeMem ark_mem, N_Vector ynew,
                                sunrealtype hscale, sunrealtype t0,
                                ARKVecResizeFn resize, void* resize_data)
{
  // TODO: update lrw, liw?
  return ARK_SUCCESS;
}

static void splittingStep_Free(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  if (step_mem == NULL) { return; }

  // TODO keep track of workspace size
  if (step_mem->own_policy)
  {
    ARKodeSplittingExecutionPolicyFree(&step_mem->policy);
  }

  ARKodeSplittingCoefficients_Free(step_mem->coefficients);
  free(step_mem);
  ark_mem->step_mem = NULL;
}

void splittingStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeSplittingStepMem step_mem;
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
  ARKodeSplittingCoefficients_Write(step_mem->coefficients, outfile);

  //TODO: consider printing policy
}

int splittingStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem;
  const int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  step_mem->order        = DEFAULT_ORDER;
  step_mem->coefficients = NULL;
  step_mem->policy       = NULL;
  return ARK_SUCCESS;
}

static int splittingStep_SetOrder(ARKodeMem ark_mem, const int order)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user-provided value, or default, depending on argument */
  if (order <= 0) { step_mem->order = 1; }
  else { step_mem->order = order; }

  if (step_mem->coefficients)
  {
    ARKodeSplittingCoefficients_Free(step_mem->coefficients);
    step_mem->coefficients = NULL;
  }

  return ARK_SUCCESS;
}

void* SplittingStepCreate(SUNStepper* steppers, const int partitions,
                          const sunrealtype t0, N_Vector y0, SUNContext sunctx)
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

  ARKodeSplittingStepMem step_mem =
    (ARKodeSplittingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->steppers          = steppers;
  step_mem->coefficients      = NULL;
  step_mem->policy            = NULL; // TODO: when to initialize?
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
  ark_mem->interp_type          = ARK_INTERP_LAGRANGE;

  /* Set default values for ARKStep optional inputs */
  int retval = splittingStep_SetDefaults(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  // TODO update liw, lrw

  /* Initialize main ARKODE infrastructure */
  // TODO will this create interpolation then have it immediately replaced?
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

int SplittingStep_SetCoefficients(void* arkode_mem,
                                  ARKodeSplittingCoefficients coefficients)
{
  ARKodeMem ark_mem;
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (step_mem->partitions != coefficients->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The splitting method has %i partitions but the coefficients have %i.",
                    step_mem->partitions, coefficients->partitions);
    return ARK_ILL_INPUT;
  }

  ARKodeSplittingCoefficients_Free(step_mem->coefficients);
  step_mem->coefficients = ARKodeSplittingCoefficients_Copy(coefficients);
  // TODO adjust workspace size

  return ARK_SUCCESS;
}

int SplittingStep_SetExecutionPolicy(void* arkode_mem,
                                     ARKodeSplittingExecutionPolicy policy)
{
  ARKodeMem ark_mem;
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  // TODO Error if policy is NULL?

  if (step_mem->own_policy)
  {
    ARKodeSplittingExecutionPolicyFree(&step_mem->policy);
  }
  step_mem->policy = policy;

  // TODO: adjust workspace size

  return ARK_SUCCESS;
}