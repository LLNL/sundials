/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for MRIStep's multirate adaptivity controller
 * layer.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

/* SUNAdaptController_MRIStep heuristic constant: */
/*   minimum estimated temporal error for inner solver */
#define INNER_MIN_DSM SUNRsqrt(SUN_UNIT_ROUNDOFF)

/*--------------------------------------------
  MRIStep SUNAdaptController wrapper functions
  --------------------------------------------*/

SUNAdaptController SUNAdaptController_MRIStep(void* arkode_mem,
                                              SUNAdaptController CMRI)
{
  SUNAdaptController C;
  mriStepControlContent content;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* Return with failure if input controller is NULL or has
     unsupported type */
  if (CMRI == NULL) { return (NULL); }
  if (SUNAdaptController_GetType(CMRI) != SUN_ADAPTCONTROLLER_MRI_TOL)
  {
    return (NULL);
  }

  /* Return with failure if stepper is inaccessible */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return (NULL);

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(CMRI->sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_MRIStep;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_MRIStep;
  C->ops->reset        = SUNAdaptController_Reset_MRIStep;
  C->ops->write        = SUNAdaptController_Write_MRIStep;
  C->ops->updateh      = SUNAdaptController_UpdateH_MRIStep;
  C->ops->space        = SUNAdaptController_Space_MRIStep;

  /* Create content */
  content = NULL;
  content = (mriStepControlContent)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach ARKODE memory, MRI stepper memory and MRI controller objects */
  content->ark_mem  = ark_mem;
  content->step_mem = step_mem;
  content->C        = CMRI;

  /* Attach content and return */
  C->content = content;
  return (C);
}

int SUNAdaptController_EstimateStep_MRIStep(SUNAdaptController C, sunrealtype H,
                                            int P, sunrealtype DSM,
                                            sunrealtype* Hnew)
{
  /* Shortcuts to ARKODE and MRIStep memory */
  int retval;
  ARKodeMem ark_mem         = MRICONTROL_A(C);
  ARKodeMRIStepMem step_mem = MRICONTROL_S(C);
  if ((ark_mem == NULL) || (step_mem == NULL)) { return SUN_ERR_MEM_FAIL; }

  /* Enforce bound on inner_dsm */
  step_mem->inner_dsm = SUNMAX(step_mem->inner_dsm, INNER_MIN_DSM);

  /* Estimate slow stepsize from MRI controller */
  retval = SUNAdaptController_EstimateStepTol(MRICONTROL_C(C), H,
                                              step_mem->inner_control,
                                              P, DSM, step_mem->inner_dsm, Hnew,
                                              &(step_mem->inner_control_new));
  return retval;
}

int SUNAdaptController_UpdateH_MRIStep(SUNAdaptController C, sunrealtype H,
                                       sunrealtype DSM)
{
  /* Shortcuts to ARKODE and MRIStep memory */
  ARKodeMem ark_mem         = MRICONTROL_A(C);
  ARKodeMRIStepMem step_mem = MRICONTROL_S(C);
  if ((ark_mem == NULL) || (step_mem == NULL)) { return SUN_ERR_MEM_FAIL; }

  /* Update MRI controller */
  int retval;
  retval = SUNAdaptController_UpdateMRITol(MRICONTROL_C(C), H,
                                           step_mem->inner_control, DSM,
                                           step_mem->inner_dsm);
  if (retval != SUN_SUCCESS) { return (retval); }

  /* Update inner controller parameter to most-recent prediction */
  step_mem->inner_control = step_mem->inner_control_new;

  /* return with success*/
  return SUN_SUCCESS;
}

SUNAdaptController_Type SUNAdaptController_GetType_MRIStep(SUNAdaptController C)
{
  return SUNAdaptController_GetType(MRICONTROL_C(C));
}

SUNErrCode SUNAdaptController_Reset_MRIStep(SUNAdaptController C)
{
  return SUNAdaptController_Reset(MRICONTROL_C(C));
}

SUNErrCode SUNAdaptController_Write_MRIStep(SUNAdaptController C, FILE* fptr)
{
  return SUNAdaptController_Write(MRICONTROL_C(C), fptr);
}

SUNErrCode SUNAdaptController_Space_MRIStep(SUNAdaptController C,
                                            long int* lenrw, long int* leniw)
{
  return SUNAdaptController_Space(MRICONTROL_C(C), lenrw, leniw);
}

/*===============================================================
  EOF
  ===============================================================*/
