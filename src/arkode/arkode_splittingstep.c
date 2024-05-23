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
#include "arkode_splittingstep_impl.h"

#define DEFAULT_ORDER 2

// TODO: which function go in arkode_splittingstep_io.c? It seems like static functions in this file makes more sense
// TODO: check the 2nd argument of arkProcessError is correct
// TODO: workspace size
// TODO: use Lagrange interpolation by default to avoid full RHS calls?

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
  return (ARK_SUCCESS);
}

static sunbooleantype splittingStep_CheckNVector(N_Vector y) // TODO: check all ops are correct
{
  return y->ops->nvclone == NULL || y->ops->nvdestroy == NULL ||
         y->ops->nvlinearsum == NULL;
}

static int splittingStep_Init(ARKodeMem ark_mem, int init_type)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return (ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT && step_mem->coefficients == NULL)
  {
    if (step_mem->order <= 1)
    {
      step_mem->coefficients = ARKodeSplittingCoefficients_LieTrotter(step_mem->partitions);
    }
    else
    {
      // TODO consider extrapolation for odd orders
      const int next_even_order = (step_mem->order + 1) / 2;
      step_mem->coefficients = ARKodeSplittingCoefficients_TripleJump(step_mem->partitions,
                                                          next_even_order);
    }
  }

  // TODO: set Lagrange interpolation and order
  // TODO: initialize policy if NULL

  return (ARK_SUCCESS);
}

static int splittingStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y,
                                 N_Vector f, int mode)
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

  return 0;
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
  ark_mem->tcur                = ark_mem->tn;

  // for (int j = 0; j < coefficients->stages; j++)
  // {
  //   // const sunrealtype tout = ark_mem->tcur +
  //   //                          coefficients->
  //   for (int k = 0; k < coefficients->partitions; k++)
  //   {
  //     if (coefficients->beta[i][j][k] != 0)
  //     {
  //       // SUNStepper stepper = step_mem->steppers[k];
  //       // retval             = stepper->evolve(stepper, ark_mem->tcur, tout, y);
  //     }
  //   }
  //   // ark_mem->tcur = tout;
  // }

  return 0;
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

  return (ARK_SUCCESS);
}

int splittingStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  fprintf(fp, "SplittingStep time step module parameters:\n  Method order %i\n\n",
          step_mem->order);

  return (ARK_SUCCESS);
}

static int splittingStep_Resize(ARKodeMem ark_mem, N_Vector ynew,
                                sunrealtype hscale, sunrealtype t0,
                                ARKVecResizeFn resize, void* resize_data)
{
  // TODO: update lrw, liw?
  return (ARK_SUCCESS);
}

static void splittingStep_Free(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  if (step_mem == NULL) { return; }

  // TODO keep track of workspace size
  if (step_mem->own_policy)
  {
    // TODO figure out name
    // SplittingPolicyFree(step_mem->policy);
  }

  ARKodeSplittingCoefficients_Free(step_mem->coefficients);
}

void splittingStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeSplittingStepMem step_mem;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
}

int splittingStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem;

  /* access ARKodeLSRKStepMem structure */
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set default values for integrator optional inputs
     (overwrite some adaptivity params for LSRKStep use) */
  step_mem->order  = DEFAULT_ORDER; /* method order */
  step_mem->coefficients = NULL;          /* no Butcher table */
  return (ARK_SUCCESS);
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

  return (ARK_SUCCESS);
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
  step_mem->coefficients            = NULL;
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

int SplittingStep_SetCoefficients(void *arkode_mem, ARKodeSplittingCoefficients coefficients) {
  int retval;
  ARKodeMem ark_mem;
  ARKodeSplittingStepMem step_mem;
  // TODO set memory

  if (step_mem->partitions != coefficients->partitions) {
    // TODO error
  }

  step_mem->coefficients = ARKodeSplittingCoefficients_Copy(coefficients);
  // TODO adjust workspace size

  return ARK_SUCCESS;
}

int SplittingStep_SetExecutionPolicy(void* arkode_mem, ARKodeSplittingExecutionPolicy policy) {
  int retval;
  ARKodeMem ark_mem;
  ARKodeSplittingStepMem step_mem;
  // TODO set memory

  // TODO: make copy?
  step_mem->policy = policy;

  // TODO: adjust workspace size

  return ARK_SUCCESS;
}