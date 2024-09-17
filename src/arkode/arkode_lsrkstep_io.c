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
    ARKODE_LSRK_SSPs_2
    ARKODE_LSRK_SSPs_3
    ARKODE_LSRK_SSP10_4
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
  case ARKODE_LSRK_SSPs_2:
    ark_mem->step               = lsrkStep_TakeStepSSPs2;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSPs_2;
    step_mem->isSSP             = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->reqstages         = 10;
    step_mem->nfusedopvecs      = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 1;
    break;
  case ARKODE_LSRK_SSPs_3:
    ark_mem->step               = lsrkStep_TakeStepSSPs3;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSPs_3;
    step_mem->isSSP             = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->reqstages         = 9;
    step_mem->nfusedopvecs      = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 3;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    break;
  case ARKODE_LSRK_SSP10_4:
    ark_mem->step               = lsrkStep_TakeStepSSP104;
    step_mem->LSRKmethod        = ARKODE_LSRK_SSP10_4;
    step_mem->isSSP             = SUNTRUE;
    ark_mem->step_printallstats = lsrkSSPStep_PrintAllStats;
    step_mem->reqstages         = 10;
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
  LSRKStepSetDomEigFn specifies the DomEig function.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn DomEig)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set the DomEig routine pointer, and update relevant flags */
  if (DomEig != NULL)
  {
    step_mem->isextDomEig = SUNTRUE;
    step_mem->extDomEig   = DomEig;

    return (ARK_SUCCESS);
  }
  else
  {
    step_mem->isextDomEig = SUNFALSE;
    step_mem->extDomEig   = NULL;

    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Internal DomEig is not supported yet!");
    return (ARK_ILL_INPUT);
  }
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFrequency sets DomEig computation frequency -
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
    step_mem->constJac   = SUNTRUE;
    step_mem->domeigfreq = 1;
  }
  else { step_mem->domeigfreq = nsteps; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetMaxNumStages sets the maximum number of stages allowed.
  ---------------------------------------------------------------*/
int LSRKStepSetMaxNumStages(void* arkode_mem, int stagemaxlimit)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (stagemaxlimit < 2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stagemaxlimit must be greater than or equal to 2");
    return (ARK_ILL_INPUT);
  }

  step_mem->stagemaxlimit = stagemaxlimit;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigSafetyFactor sets the safety factor for the DomEigs.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigSafetyFactor(void* arkode_mem, sunrealtype domeigsfty)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (domeigsfty < SUN_RCONST(1.0))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "domeigsfty must be greater than or equal to 1");
    return (ARK_ILL_INPUT);
  }

  step_mem->domeigsfty = domeigsfty;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepSetSSPStageNum sets the number of stages in the following
  SSP methods:

      ARKODE_LSRK_SSPs_2  -- numofstages must be greater than or equal to 2
      ARKODE_LSRK_SSPs_3  -- numofstages must be a full-square greater than or equal to 9
      ARKODE_LSRK_SSP10_4 -- numofstages must be equal to 10 - no need to call!

  This set routine must be called after calling LSRKStepSetMethod with an SSP method
  ---------------------------------------------------------------*/
int LSRKStepSetSSPStageNum(void* arkode_mem, int numofstages)
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
  case ARKODE_LSRK_SSPs_2:
    if (numofstages < 2)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "numofstages must be greater than or equal to 2");
      return (ARK_ILL_INPUT);
    }
    step_mem->reqstages = numofstages;
    break;

  case ARKODE_LSRK_SSPs_3:
    if (numofstages < 9 ||
        (ceil(SUNRsqrt(numofstages)) != floor(SUNRsqrt(numofstages))))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "numofstages must be a full-square greater than or equal to 9");
      return (ARK_ILL_INPUT);
    }
    else { step_mem->reqstages = numofstages; }

    break;

  case ARKODE_LSRK_SSP10_4:
    if (numofstages != 10)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "SSP10_4 method has a prefixed numofstages = 10");
      return (ARK_ILL_INPUT);
    }
    else { step_mem->reqstages = numofstages; }

    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Call LSRKStepSetMethod to declare SSP method type first!");
    return (ARK_ILL_INPUT);
    break;
  }

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
int LSRKStepGetNumDomEigUpdates(void* arkode_mem, long int* ndomeigupdates)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *ndomeigupdates = step_mem->ndomeigupdates;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  LSRKStepGetMaxNumStages:

  Returns the max number of stages taken
  ---------------------------------------------------------------*/
int LSRKStepGetMaxNumStages(void* arkode_mem, int* stagemax)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get values from step_mem */
  *stagemax = step_mem->stagemax;

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
  ark_mem->hadapt_mem->hcontroller = NULL;
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
  step_mem->reqstages = 0; /* no stages */

  /* Counters and stats*/
  step_mem->nfe            = 0;
  step_mem->domeignfe      = 0;
  step_mem->stagemax       = 0;
  step_mem->ndomeigupdates = 0;
  step_mem->stagemaxlimit =
    (int)round(SUNRsqrt(ark_mem->reltol / (10.0 * ark_mem->uround)));
  step_mem->stagemaxlimit =
    (step_mem->stagemaxlimit > 2) ? step_mem->stagemaxlimit : 2;
  step_mem->nstsig = 0;

  /* Spectral info */
  step_mem->lambdaR    = 0;
  step_mem->lambdaI    = 0;
  step_mem->sprad      = 0;
  step_mem->sprmax     = 0;
  step_mem->sprmin     = 0;
  step_mem->domeigsfty = 1.01;
  step_mem->domeigfreq = 25;

  /* Flags */
  step_mem->newdomeig = SUNTRUE;
  step_mem->constJac  = SUNFALSE;
  step_mem->jacatt    = SUNFALSE;
  step_mem->isSSP     = SUNFALSE;

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
            step_mem->domeignfe);
    fprintf(outfile, "Number of DomEig update calls = %ld\n",
            step_mem->ndomeigupdates);
    fprintf(outfile, "Max. num. of stages taken     = %d\n", step_mem->stagemax);
    fprintf(outfile, "Max. num. of stages allowed   = %d\n",
            step_mem->stagemaxlimit);

    fprintf(outfile, "Max. spectral radius         = %.2f\n", step_mem->sprmax);
    fprintf(outfile, "Min. spectral radius         = %.2f\n", step_mem->sprmin);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, ",RHS fn evals for Dom. Eigs.,%ld", step_mem->domeignfe);
    fprintf(outfile, ",Number of DomEig update calls,%ld",
            step_mem->ndomeigupdates);
    fprintf(outfile, ",Max. num. of stages taken,%d", step_mem->stagemax);
    fprintf(outfile, ",Max. num. of stages allowed,%d", step_mem->stagemaxlimit);

    fprintf(outfile, ",Max. spectral radius,%.2f", step_mem->sprmax);
    fprintf(outfile, ",Min. spectral radius,%.2f", step_mem->sprmin);
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
    fprintf(outfile, "Current stages taken         = %d\n", step_mem->reqstages);
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, ",RHS fn evals,%ld", step_mem->nfe);
    fprintf(outfile, ",Current stages taken,%d", step_mem->reqstages);

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
    fprintf(fp, "\n");
    break;
  case ARKODE_LSRK_RKL_2:
    fprintf(fp, "LSRKStep RKL time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 2);
    fprintf(fp, "\n");
    break;
  // case ARKODE_LSRK_RKG_2:
  //   fprintf(fp, "LSRKStep RKG time step module parameters:\n");
  //   fprintf(fp, "  Method order %i\n", 2);
  //   fprintf(fp, "\n");
  //   break;
  case ARKODE_LSRK_SSPs_2:
    fprintf(fp, "LSRKStep SSP(s,2) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 2);
    fprintf(fp, "\n");
    break;
  case ARKODE_LSRK_SSPs_3:
    fprintf(fp, "LSRKStep SSP(s,3) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 3);
    fprintf(fp, "\n");
    break;
  case ARKODE_LSRK_SSP10_4:
    fprintf(fp, "LSRKStep SSP(10,4) time step module parameters:\n");
    fprintf(fp, "  Method order %i\n", 4);
    fprintf(fp, "\n");
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return (ARK_ILL_INPUT);
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
