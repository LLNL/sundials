/*---------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include "arkode_lsrkstep_impl.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepSetMethod sets method
    ARKODE_LSRK_RKC_2
    ARKODE_LSRK_RKL_2
    ARKODE_LSRK_RKG_2
    ARKODE_LSRK_SSP_S_2
    ARKODE_LSRK_SSP_S_3
    ARKODE_LSRK_SSP_10_4
  ---------------------------------------------------------------*/
int LSRKStepSetMethod(void* arkode_mem, ARKODE_LSRKMethodType method)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  switch (method)
  {
  case ARKODE_LSRK_RKC_2:
    ark_mem->step          = lsrkStep_TakeStepRKC;
    step_mem->LSRKmethod   = ARKODE_LSRK_RKC_2;
    step_mem->nfusedopvecs = 5;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    break;
  case ARKODE_LSRK_RKL_2:
    ark_mem->step          = lsrkStep_TakeStepRKL;
    step_mem->LSRKmethod   = ARKODE_LSRK_RKL_2;
    step_mem->nfusedopvecs = 5;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    break;
  // case ARKODE_LSRK_RKG_2:
  //   ark_mem->step = lsrkStep_TakeStepRKG;
  //   step_mem->LSRKmethod = ARKODE_LSRK_RKG_2;
  //   step_mem->nfusedopvecs = 5;
  //   step_mem->q = ark_mem->hadapt_mem->q = 2;
  //   step_mem->p = ark_mem->hadapt_mem->p = 2;
  //   break;
  case ARKODE_LSRK_SSP_S_2:
    ark_mem->step               = lsrkStep_TakeStepSSPs2;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSP_S_2;
    step_mem->is_SSP            = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->req_stages        = 10;
    step_mem->nfusedopvecs      = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 1;
    break;
  case ARKODE_LSRK_SSP_S_3:
    ark_mem->step               = lsrkStep_TakeStepSSPs3;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSP_S_3;
    step_mem->is_SSP            = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->req_stages        = 9;
    step_mem->nfusedopvecs      = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 3;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    break;
  case ARKODE_LSRK_SSP_10_4:
    ark_mem->step               = lsrkStep_TakeStepSSP104;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSP_10_4;
    step_mem->is_SSP            = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->req_stages        = 10;
    step_mem->nfusedopvecs      = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 4;
    step_mem->p = ark_mem->hadapt_mem->p = 3;
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return (ARK_ILL_INPUT);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFn specifies the dom_eig function.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set the dom_eig routine pointer, and update relevant flags */
  if (dom_eig != NULL)
  {
    step_mem->is_ext_dom_eig = SUNTRUE;
    step_mem->extDomEig      = dom_eig;

    return (ARK_SUCCESS);
  }
  else
  {
    step_mem->is_ext_dom_eig = SUNFALSE;
    step_mem->extDomEig      = NULL;

    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Internal dom_eig is not supported yet!");
    return (ARK_ILL_INPUT);
  }
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFrequency sets dom_eig computation frequency -
  Dominated Eigenvalue is recomputed after "nsteps" successful steps.

    nsteps = 0 refers to Constant Jacobian
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigFrequency(void* arkode_mem, int nsteps)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (nsteps < 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "nsteps must be greater than or equal to 0");
    return (ARK_ILL_INPUT);
  }

  if (nsteps == 0)
  {
    step_mem->const_Jac    = SUNTRUE;
    step_mem->dom_eig_freq = 1;
  }
  else
  {
    step_mem->dom_eig_freq = nsteps;
    step_mem->const_Jac    = SUNFALSE;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetMaxNumStages sets the maximum number of stages allowed.
  ---------------------------------------------------------------*/
int LSRKStepSetMaxNumStages(void* arkode_mem, int stage_max_limit)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (stage_max_limit < 2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stage_max_limit must be greater than or equal to 2");
    return (ARK_ILL_INPUT);
  }

  step_mem->stage_max_limit = stage_max_limit;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigSafetyFactor sets the safety factor for the DomEigs.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigSafetyFactor(void* arkode_mem, sunrealtype dom_eig_sfty)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (dom_eig_sfty < SUN_RCONST(1.0))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "dom_eig_sfty must be greater than or equal to 1");
    return (ARK_ILL_INPUT);
  }

  step_mem->dom_eig_sfty = dom_eig_sfty;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetSSPStageNum sets the number of stages in the following
  SSP methods:

      ARKODE_LSRK_SSP_S_2  -- num_of_stages must be greater than or equal to 2
      ARKODE_LSRK_SSP_S_3  -- num_of_stages must be a perfect square greater than or equal to 9
      ARKODE_LSRK_SSP_10_4 -- num_of_stages must be equal to 10 - no need to call!

  This set routine must be called after calling LSRKStepSetMethod with an SSP method
  ---------------------------------------------------------------*/
int LSRKStepSetSSPStageNum(void* arkode_mem, int num_of_stages)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  switch (step_mem->LSRKmethod)
  {
  case ARKODE_LSRK_SSP_S_2:
    if (num_of_stages < 2)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "num_of_stages must be greater than or equal to 2");
      return (ARK_ILL_INPUT);
    }
    break;

  case ARKODE_LSRK_SSP_S_3:
    int root = SUNRsqrt(num_of_stages);
    if ((num_of_stages < 9) || ((root * root) != num_of_stages))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "num_of_stages must be a perfect square greater than or equal to 9");
      return (ARK_ILL_INPUT);
    }
    break;

  case ARKODE_LSRK_SSP_10_4:
    if (num_of_stages != 10)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "SSP10_4 method has a prefixed num_of_stages = 10");
      return (ARK_ILL_INPUT);
    }
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Call LSRKStepSetMethod to declare SSP method type first!");
    return (ARK_ILL_INPUT);
    break;
  }
  step_mem->req_stages = num_of_stages;

  return (ARK_SUCCESS);
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepGetNumRhsEvals:

  Returns the current number of calls to f
  ---------------------------------------------------------------*/
int LSRKStepGetNumRhsEvals(void* arkode_mem, long int* fe_evals,
                           long int* fi_evals)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *fe_evals = step_mem->nfe;
  *fi_evals = step_mem->nfi;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepGetNumDomEigUpdates:

  Returns the number of dominant eigenvalue updates
  ---------------------------------------------------------------*/
int LSRKStepGetNumDomEigUpdates(void* arkode_mem, long int* num_dom_eig_updates)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *num_dom_eig_updates = step_mem->num_dom_eig_updates;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepGetMaxNumStages:

  Returns the max number of stages taken
  ---------------------------------------------------------------*/
int LSRKStepGetMaxNumStages(void* arkode_mem, int* stage_max)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *stage_max = step_mem->stage_max;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepGetAverageStageNum:

  Returns the average number of stages taken
  ---------------------------------------------------------------*/
int LSRKStepGetAverageStageNum(void* arkode_mem, sunrealtype* averstage)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *averstage = ((sunrealtype)step_mem->nfe) /
               ((sunrealtype)ark_mem->nst_attempts);

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
  if (ark_mem->hadapt_mem->owncontroller)
  {
    retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                      &leniw);
    if (retval == SUN_SUCCESS)
    {
      ark_mem->liw -= leniw;
      ark_mem->lrw -= lenrw;
    }
    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = SUNAdaptController_PI(ark_mem->sunctx);
  if (ark_mem->hadapt_mem->hcontroller == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "SUNAdaptControllerPI allocation failure");
    return (ARK_MEM_FAIL);
  }
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;
  retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                    &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }

  /* Set default values for integrator optional inputs
     (overwrite some adaptivity params for LSRKStep use) */
  step_mem->req_stages = 0; /* no stages */

  /* Spectral info */
  step_mem->lambdaR      = 0;
  step_mem->lambdaI      = 0;
  step_mem->sprad        = 0;
  step_mem->spr_max      = 0;
  step_mem->spr_min      = 0;
  step_mem->dom_eig_sfty = 1.01;
  step_mem->dom_eig_freq = 25;

  /* Flags */
  step_mem->new_dom_eig = SUNTRUE;
  step_mem->const_Jac   = SUNFALSE;
  step_mem->jacatt      = SUNFALSE;
  step_mem->is_SSP      = SUNFALSE;

  ark_mem->hadapt_mem->adjust = 0;               /* set default adjustment */
  ark_mem->hadapt_mem->etamxf = SUN_RCONST(0.3); /* max change on error-failed step */
  ark_mem->hadapt_mem->safety = SUN_RCONST(0.99); /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->growth = SUN_RCONST(25.0); /* step adaptivity growth factor */

  (void)SUNAdaptController_SetErrorBias(ark_mem->hadapt_mem->hcontroller,
                                        SUN_RCONST(1.2));
  (void)SUNAdaptController_SetParams_PI(ark_mem->hadapt_mem->hcontroller,
                                        SUN_RCONST(0.8), -SUN_RCONST(0.31));
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
int lsrkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele)
{
  int retval;
  ARKodeLSRKStepMem step_mem;
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* return an error if local truncation error is not computed */
  if (ark_mem->fixedstep) { return (ARK_STEPPER_UNSUPPORTED); }

  /* otherwise, copy local truncation error vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_PrintAllStats:

  Prints integrator statistics for STS methods
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
    fprintf(outfile, "RHS fn evals                  = %ld\n", step_mem->nfe);
    fprintf(outfile, "RHS fn evals for Dom. Eigs.   = %ld\n",
            step_mem->dom_eig_nfe);
    fprintf(outfile, "Number of dom_eig update calls = %ld\n",
            step_mem->num_dom_eig_updates);
    fprintf(outfile, "Max. num. of stages taken     = %d\n", step_mem->stage_max);
    fprintf(outfile, "Max. num. of stages allowed   = %d\n",
            step_mem->stage_max_limit);

    fprintf(outfile, "Max. spectral radius         = %.2f\n", step_mem->spr_max);
    fprintf(outfile, "Min. spectral radius         = %.2f\n", step_mem->spr_min);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, ",RHS fn evals for Dom. Eigs.,%ld", step_mem->dom_eig_nfe);
    fprintf(outfile, ",Number of dom_eig update calls,%ld",
            step_mem->num_dom_eig_updates);
    fprintf(outfile, ",Max. num. of stages taken,%d", step_mem->stage_max);
    fprintf(outfile, ",Max. num. of stages allowed,%d",
            step_mem->stage_max_limit);

    fprintf(outfile, ",Max. spectral radius,%.2f", step_mem->spr_max);
    fprintf(outfile, ",Min. spectral radius,%.2f", step_mem->spr_min);
    fprintf(outfile, "\n");
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkSSPStep_PrintAllStats:

  Prints integrator statistics for SSP methods
  ---------------------------------------------------------------*/
int lsrkSSPStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile,
                              SUNOutputFormat fmt)
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
    fprintf(outfile, "Current stages taken         = %d\n", step_mem->req_stages);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, ",Current stages taken,%d", step_mem->req_stages);

    fprintf(outfile, "\n");
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

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
  switch (step_mem->LSRKmethod)
  {
  case ARKODE_LSRK_RKC_2:
    fprintf(fp, "LSRKStep RKC time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 2);
    break;
  case ARKODE_LSRK_RKL_2:
    fprintf(fp, "LSRKStep RKL time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 2);
    break;
  // case ARKODE_LSRK_RKG_2:
  //   fprintf(fp, "LSRKStep RKG time step module parameters:\n");
  //   fprintf(fp, "  Method order %i\n", 2);
  //   break;
  case ARKODE_LSRK_SSP_S_2:
    fprintf(fp, "LSRKStep SSP(s,2) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 2);
    break;
  case ARKODE_LSRK_SSP_S_3:
    fprintf(fp, "LSRKStep SSP(s,3) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 3);
    break;
  case ARKODE_LSRK_SSP_10_4:
    fprintf(fp, "LSRKStep SSP(10,4) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 4);
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return (ARK_ILL_INPUT);
  }
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
