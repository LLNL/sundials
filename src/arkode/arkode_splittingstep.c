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

static int splittingStep_Init(void* arkode_mem, int init_type) { return 0; }

static int splittingStep_FullRHS(void* arkode_mem, sunrealtype t, N_Vector y,
                                 N_Vector f, int mode)
{
  return 0;
}

static int splittignStep_TakeStep(void* arkode_mem, sunrealtype* dsm, int* nflag) {
  return 0;
}

void* ARKStepCreateSplitting(SUNStepper* steppers, const int partitions,
                             const sunrealtype t0, N_Vector y0, SUNContext sunctx)
{
  if (steppers == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "TODO"); // TODO: make message macro
    return NULL;
  }

  for (int i = 0; i < partitions; i++)
  {
    if (steppers[i] == NULL)
    {
      arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "TODO"); // TODO: make message macro
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
  if (!arkStep_CheckNVector(y0))
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

  ARKodeSplittingStepMem step_mem = (ARKodeSplittingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKStepFree((void**)&ark_mem);
    return NULL;
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
  ark_mem->step_init           = splittingStep_Init;
  ark_mem->step_fullrhs        = splittingStep_FullRHS;
  ark_mem->step                = splittignStep_TakeStep;
  ark_mem->step_mem            = (void*)step_mem;

  /* Set default values for ARKStep optional inputs */
  retval = ARKStepSetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }

  /* Initialize main ARKODE infrastructure */
  retval = splittingStep_Init(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }

  return ((void*)ark_mem);
}