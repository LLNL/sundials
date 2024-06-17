/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This is the implementation file for the optional input and
 * output functions for the ARKODE LSRKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_lsrkstep_impl.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepSetTableName:

  Specifies to use a pre-existing LSRK table for the problem.
  ---------------------------------------------------------------*/
int LSRKStepSetTableName(void* arkode_mem, const char* etable)
{
  /* Fill this in */
  return ARK_SUCCESS;
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepGetNumRhsEvals:

  Returns the current number of calls to f
  ---------------------------------------------------------------*/
int LSRKStepGetNumRhsEvals(void* arkode_mem, long int* fevals)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *fevals = step_mem->nfe;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int LSRKStepGetTimestepperStats(void* arkode_mem, long int* expsteps,
                               long int* accsteps, long int* attempts,
                               long int* fevals, long int* netfails)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  // /* set expsteps and accsteps from adaptivity structure */
  // *expsteps = ark_mem->hadapt_mem->nst_exp;
  // *accsteps = ark_mem->hadapt_mem->nst_acc;

  // /* set remaining outputs */
  // *attempts = ark_mem->nst_attempts;
  // *fevals   = step_mem->nfe;
  // *netfails = ark_mem->netf;

  printf("\nLSRKStepGetTimestepperStats is not ready yet!\n");

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_SetDefaults:

  Resets all LSRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int lsrkStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeLSRKStepMem step_mem;
  int retval;
  long int lenrw, leniw;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Remove current SUNAdaptController object, and replace with "PI" */
  retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                    &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  if (ark_mem->hadapt_mem->owncontroller)
  {
    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }

  printf("\nCheck hcontroller in %s at line: %d\n\n", __func__, __LINE__);
  // ark_mem->hadapt_mem->hcontroller = NULL;
  // ark_mem->hadapt_mem->hcontroller = SUNAdaptController_PI(ark_mem->sunctx);
  // if (ark_mem->hadapt_mem->hcontroller == NULL)
  // {
  //   arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
  //                   "SUNAdaptControllerPI allocation failure");
  //   return (ARK_MEM_FAIL);
  // }
  // ark_mem->hadapt_mem->owncontroller = SUNTRUE;
  // retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
  //                                   &leniw);
  // if (retval == SUN_SUCCESS)
  // {
  //   ark_mem->liw += leniw;
  //   ark_mem->lrw += lenrw;
  // }


  /* Set default values for integrator optional inputs
     (overwrite some adaptivity params for LSRKStep use) */
  step_mem->reqstages = 0;                       /* no stages */
  step_mem->absh = 1.0/ark_mem->hmax_inv;        /* set absh = hmax */ 
  step_mem->err = 0;
  step_mem->errold = 0;

  /* Counters and stats*/
  step_mem->nfe = 0;
  step_mem->sprnfe = 0;
  step_mem->stagemax = 0;
  step_mem->stagemaxlimit = SUNMAX(2, round(SUNRsqrt(ark_mem->reltol/(10.0*ark_mem->uround))));
  step_mem->nreject = 0;
  step_mem->nstsig = 0;

  /* Spectral radius info */  
  step_mem->sprad = 1;
  step_mem->sprmax = 0;
  step_mem->sprmin = 0;
  step_mem->sprsfty = 1.1;
  step_mem->sprupdatepar = 25;

  /* Flags */
  step_mem->newspr = SUNTRUE;
  step_mem->jacatt = SUNFALSE;

  ark_mem->hadapt_mem->etamxf = SUN_RCONST(0.3); /* max change on error-failed step */
  ark_mem->hadapt_mem->safety = SUN_RCONST(0.99); /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->growth = SUN_RCONST(25.0); /* step adaptivity growth factor */
  
  printf("\nCheck SUNAdaptController_SetErrorBias in %s at line: %d\n\n", __func__, __LINE__);
  // (void)SUNAdaptController_SetErrorBias(ark_mem->hadapt_mem->hcontroller,
  //                                       SUN_RCONST(1.2));
  // (void)SUNAdaptController_SetParams_PI(ark_mem->hadapt_mem->hcontroller,
  //                                       SUN_RCONST(0.8), -SUN_RCONST(0.31));
  
  printf("\nlsrkStep_SetDefaults is not ready yet!\n");

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int lsrkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "RHS fn evals                 = %ld\n", step_mem->nfe);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, "\n");
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

  printf("\nlsrkStep_PrintAllStats is not ready yet!\n");

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int lsrkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* print integrator parameters to file */
  fprintf(fp, "LSRKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", NULL);
  fprintf(fp, "\n");

  printf("\nlsrkStep_WriteParameters is not ready yet!\n");

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
