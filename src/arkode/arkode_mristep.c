/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 *                Rujeko Chinomona @ SMU
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
 * This is the implementation file for ARKODE's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_mristep_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0, N_Vector y0,
                    MRIStepInnerStepper stepper, SUNContext sunctx)
{
  ARKodeMem ark_mem;         /* outer ARKODE memory   */
  ARKodeMRIStepMem step_mem; /* outer stepper memory  */
  SUNNonlinearSolver NLS;    /* default nonlin solver */
  sunbooleantype nvectorOK;
  int retval;

  /* Check that at least one of fse, fsi is supplied and is to be used*/
  if (fse == NULL && fsi == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (NULL);
  }

  /* Check that stepper is supplied */
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The inner stepper memory is NULL");
    return (NULL);
  }

  /* Check that context is supplied */
  if (!sunctx)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return (NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = mriStep_CheckNVector(y0);
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

  /* Allocate ARKodeMRIStepMem structure, and initialize to zero */
  step_mem = (ARKodeMRIStepMem)calloc(1, sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol              = mriStep_AttachLinsol;
  ark_mem->step_disablelsetup             = mriStep_DisableLSetup;
  ark_mem->step_getlinmem                 = mriStep_GetLmem;
  ark_mem->step_getimplicitrhs            = mriStep_GetImplicitRHS;
  ark_mem->step_getgammas                 = mriStep_GetGammas;
  ark_mem->step_init                      = mriStep_Init;
  ark_mem->step_fullrhs                   = mriStep_FullRHS;
  ark_mem->step                           = mriStep_TakeStepMRIGARK;
  ark_mem->step_setuserdata               = mriStep_SetUserData;
  ark_mem->step_printallstats             = mriStep_PrintAllStats;
  ark_mem->step_writeparameters           = mriStep_WriteParameters;
  ark_mem->step_resize                    = mriStep_Resize;
  ark_mem->step_reset                     = mriStep_Reset;
  ark_mem->step_free                      = mriStep_Free;
  ark_mem->step_printmem                  = mriStep_PrintMem;
  ark_mem->step_setdefaults               = mriStep_SetDefaults;
  ark_mem->step_computestate              = mriStep_ComputeState;
  ark_mem->step_setorder                  = mriStep_SetOrder;
  ark_mem->step_setnonlinearsolver        = mriStep_SetNonlinearSolver;
  ark_mem->step_setlinear                 = mriStep_SetLinear;
  ark_mem->step_setnonlinear              = mriStep_SetNonlinear;
  ark_mem->step_setnlsrhsfn               = mriStep_SetNlsRhsFn;
  ark_mem->step_setdeduceimplicitrhs      = mriStep_SetDeduceImplicitRhs;
  ark_mem->step_setnonlincrdown           = mriStep_SetNonlinCRDown;
  ark_mem->step_setnonlinrdiv             = mriStep_SetNonlinRDiv;
  ark_mem->step_setdeltagammamax          = mriStep_SetDeltaGammaMax;
  ark_mem->step_setlsetupfrequency        = mriStep_SetLSetupFrequency;
  ark_mem->step_setpredictormethod        = mriStep_SetPredictorMethod;
  ark_mem->step_setmaxnonliniters         = mriStep_SetMaxNonlinIters;
  ark_mem->step_setnonlinconvcoef         = mriStep_SetNonlinConvCoef;
  ark_mem->step_setstagepredictfn         = mriStep_SetStagePredictFn;
  ark_mem->step_getnumlinsolvsetups       = mriStep_GetNumLinSolvSetups;
  ark_mem->step_getcurrentgamma           = mriStep_GetCurrentGamma;
  ark_mem->step_getestlocalerrors         = mriStep_GetEstLocalErrors;
  ark_mem->step_getnonlinearsystemdata    = mriStep_GetNonlinearSystemData;
  ark_mem->step_getnumnonlinsolviters     = mriStep_GetNumNonlinSolvIters;
  ark_mem->step_getnumnonlinsolvconvfails = mriStep_GetNumNonlinSolvConvFails;
  ark_mem->step_getnonlinsolvstats        = mriStep_GetNonlinSolvStats;
  ark_mem->step_supports_adaptive         = SUNTRUE;
  ark_mem->step_supports_implicit         = SUNTRUE;
  ark_mem->step_supports_relaxation       = SUNTRUE;
  ark_mem->step_mem                       = (void*)step_mem;

  /* Set default values for optional inputs */
  retval = mriStep_SetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Allocate the general MRI stepper vectors using y0 as a template */
  /* NOTE: Fse, Fsi, inner_forcing, cvals, Xvecs, sdata, zpred and zcor will
     be allocated later on (based on the MRI method) */

  /* Copy the slow RHS functions into stepper memory */
  step_mem->fse = fse;
  step_mem->fsi = fsi;

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit_rhs = (fse == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit_rhs = (fsi == NULL) ? SUNFALSE : SUNTRUE;

  /* Update the ARKODE workspace requirements */
  ark_mem->liw += 49; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
  ark_mem->lrw += 14;

  /* Create a default Newton NLS object (just in case; will be deleted if
     the user attaches a nonlinear solver) */
  step_mem->NLS    = NULL;
  step_mem->ownNLS = SUNFALSE;

  if (step_mem->implicit_rhs)
  {
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (!NLS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error creating default Newton solver");
      ARKodeFree((void**)&ark_mem);
      return (NULL);
    }
    retval = ARKodeSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error attaching default Newton solver");
      ARKodeFree((void**)&ark_mem);
      return (NULL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit  = NULL;
  step_mem->lsetup = NULL;
  step_mem->lsolve = NULL;
  step_mem->lfree  = NULL;
  step_mem->lmem   = NULL;

  /* Initialize initial error norm  */
  step_mem->eRNrm = ONE;

  /* Initialize all the counters */
  step_mem->nfse      = 0;
  step_mem->nfsi      = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;
  step_mem->nls_fails = 0;

  /* Initialize fused op work space */
  step_mem->cvals = NULL;
  step_mem->Xvecs = NULL;

  /* Initialize adaptivity parameters */
  step_mem->inner_control     = ONE;
  step_mem->inner_dsm         = ONE;
  step_mem->inner_control_new = ONE;

  /* Initialize pre and post inner evolve functions */
  step_mem->pre_inner_evolve  = NULL;
  step_mem->post_inner_evolve = NULL;

  /* Initialize main ARKODE infrastructure (allocates vectors) */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Attach the inner stepper memory */
  step_mem->stepper = stepper;

  /* Check for required stepper functions */
  retval = mriStepInnerStepper_HasRequiredOps(step_mem->stepper);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "A required inner stepper function is NULL");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* return ARKODE memory */
  return ((void*)ark_mem);
}

/*---------------------------------------------------------------
  MRIStepReInit:

  This routine re-initializes the MRIStep module to solve a new
  problem of the same size as was previously solved (all counter
  values are set to 0).

  NOTE: the inner stepper needs to be reinitialized before
  calling this function.
  ---------------------------------------------------------------*/
int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0,
                  N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  SUNNonlinearSolver NLS;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that at least one of fse, fsi is supplied and is to be used */
  if (fse == NULL && fsi == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit_rhs = (fse == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit_rhs = (fsi == NULL) ? SUNFALSE : SUNTRUE;

  /* Create a default Newton NLS object (just in case; will be deleted if
     the user attaches a nonlinear solver) */
  if (step_mem->implicit_rhs && !(step_mem->NLS))
  {
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (!NLS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error creating default Newton solver");
      ARKodeFree((void**)&ark_mem);
      return (ARK_MEM_FAIL);
    }
    retval = ARKodeSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error attaching default Newton solver");
      ARKodeFree((void**)&ark_mem);
      return (ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* ReInitialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to reinitialize main ARKODE infrastructure");
    return (retval);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->fse = fse;
  step_mem->fsi = fsi;

  /* Initialize all the counters */
  step_mem->nfse      = 0;
  step_mem->nfsi      = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;

  return (ARK_SUCCESS);
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_Resize:

  This routine resizes the memory within the MRIStep module.
  ---------------------------------------------------------------*/
int mriStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                   sunrealtype t0, ARKVecResizeFn resize, void* resize_data)
{
  ARKodeMRIStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Determine change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Resize Fse */
  if (step_mem->Fse)
  {
    if (!arkResizeVecArray(resize, resize_data, step_mem->nstages_allocated, y0,
                           &(step_mem->Fse), lrw_diff, &(ark_mem->lrw),
                           liw_diff, &(ark_mem->liw)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
    if (step_mem->unify_Fs) { step_mem->Fsi = step_mem->Fse; }
  }

  /* Resize Fsi */
  if (step_mem->Fsi && !step_mem->unify_Fs)
  {
    if (!arkResizeVecArray(resize, resize_data, step_mem->nstages_allocated, y0,
                           &(step_mem->Fsi), lrw_diff, &(ark_mem->lrw),
                           liw_diff, &(ark_mem->liw)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
  }

  /* Resize the nonlinear solver interface vectors (if applicable) */
  if (step_mem->sdata != NULL)
  {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                      &step_mem->sdata))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
  }
  if (step_mem->zpred != NULL)
  {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                      &step_mem->zpred))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
  }
  if (step_mem->zcor != NULL)
  {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                      &step_mem->zcor))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
  }

  /* If a NLS object was previously used, destroy and recreate default Newton
     NLS object (can be replaced by user-defined object if desired) */
  if ((step_mem->NLS != NULL) && (step_mem->ownNLS))
  {
    /* destroy existing NLS object */
    retval = SUNNonlinSolFree(step_mem->NLS);
    if (retval != ARK_SUCCESS) { return (retval); }
    step_mem->NLS    = NULL;
    step_mem->ownNLS = SUNFALSE;

    /* create new Newton NLS object */
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (NLS == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error creating default Newton solver");
      return (ARK_MEM_FAIL);
    }

    /* attach new Newton NLS object */
    retval = ARKodeSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error attaching default Newton solver");
      return (ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Resize the inner stepper vectors */
  retval = mriStepInnerStepper_Resize(step_mem->stepper, resize, resize_data,
                                      lrw_diff, liw_diff, y0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  /* reset nonlinear solver counters */
  if (step_mem->NLS != NULL) { step_mem->nsetups = 0; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_Reset:

  This routine resets the MRIStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).  It is called after the main ARKODE
  infrastructure is reset.
  ---------------------------------------------------------------*/
int mriStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Reset the inner integrator with this same state */
  retval = mriStepInnerStepper_Reset(step_mem->stepper, tR, yR);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to reset the inner stepper");
    return (ARK_INNERSTEP_FAIL);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_ComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
int mriStep_ComputeState(ARKodeMem ark_mem, N_Vector zcor, N_Vector z)
{
  int retval;
  ARKodeMRIStepMem step_mem;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_Free frees all MRIStep memory.
  ---------------------------------------------------------------*/
void mriStep_Free(ARKodeMem ark_mem)
{
  sunindextype Cliw, Clrw;
  ARKodeMRIStepMem step_mem;

  /* nothing to do if ark_mem is already NULL */
  if (ark_mem == NULL) { return; }

  /* conditional frees on non-NULL MRIStep module */
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;

    /* free the coupling structure and derived quantities */
    if (step_mem->MRIC != NULL)
    {
      MRIStepCoupling_Space(step_mem->MRIC, &Cliw, &Clrw);
      MRIStepCoupling_Free(step_mem->MRIC);
      step_mem->MRIC = NULL;
      ark_mem->liw -= Cliw;
      ark_mem->lrw -= Clrw;
      if (step_mem->stagetypes)
      {
        free(step_mem->stagetypes);
        step_mem->stagetypes = NULL;
        ark_mem->liw -= (step_mem->stages + 1);
      }
      if (step_mem->stage_map)
      {
        free(step_mem->stage_map);
        step_mem->stage_map = NULL;
        ark_mem->liw -= step_mem->stages;
      }
      if (step_mem->Ae_row)
      {
        free(step_mem->Ae_row);
        step_mem->Ae_row = NULL;
        ark_mem->lrw -= step_mem->stages;
      }
      if (step_mem->Ai_row)
      {
        free(step_mem->Ai_row);
        step_mem->Ai_row = NULL;
        ark_mem->lrw -= step_mem->stages;
      }
    }

    /* free the nonlinear solver memory (if applicable) */
    if ((step_mem->NLS != NULL) && (step_mem->ownNLS))
    {
      SUNNonlinSolFree(step_mem->NLS);
      step_mem->ownNLS = SUNFALSE;
    }
    step_mem->NLS = NULL;

    /* free the linear solver memory */
    if (step_mem->lfree != NULL)
    {
      step_mem->lfree((void*)ark_mem);
      step_mem->lmem = NULL;
    }

    /* free the sdata, zpred and zcor vectors */
    if (step_mem->sdata != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }
    if (step_mem->zpred != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->zpred);
      step_mem->zpred = NULL;
    }
    if (step_mem->zcor != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->zcor);
      step_mem->zcor = NULL;
    }

    /* free the RHS vectors */
    if (step_mem->Fse)
    {
      arkFreeVecArray(step_mem->nstages_allocated, &(step_mem->Fse),
                      ark_mem->lrw1, &(ark_mem->lrw), ark_mem->liw1,
                      &(ark_mem->liw));
      if (step_mem->unify_Fs) { step_mem->Fsi = NULL; }
    }

    if (step_mem->Fsi)
    {
      arkFreeVecArray(step_mem->nstages_allocated, &(step_mem->Fsi),
                      ark_mem->lrw1, &(ark_mem->lrw), ark_mem->liw1,
                      &(ark_mem->liw));
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL)
    {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= (step_mem->nfusedopvecs);
    }
    if (step_mem->Xvecs != NULL)
    {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (step_mem->nfusedopvecs);
    }
    step_mem->nfusedopvecs = 0;

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }
}

/*---------------------------------------------------------------
  mriStep_PrintMem:

  This routine outputs the memory from the MRIStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
void mriStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeMRIStepMem step_mem;
  int i, retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output integer quantities */
  fprintf(outfile, "MRIStep: q = %i\n", step_mem->q);
  fprintf(outfile, "MRIStep: p = %i\n", step_mem->p);
  fprintf(outfile, "MRIStep: istage = %i\n", step_mem->istage);
  fprintf(outfile, "MRIStep: stages = %i\n", step_mem->stages);
  fprintf(outfile, "MRIStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile, "MRIStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile, "MRIStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile, "MRIStep: convfail = %i\n", step_mem->convfail);
  fprintf(outfile, "MRIStep: stagetypes =");
  for (i = 0; i <= step_mem->stages; i++)
  {
    fprintf(outfile, " %i", step_mem->stagetypes[i]);
  }

  /* output long integer quantities */
  fprintf(outfile, "MRIStep: nfse = %li\n", step_mem->nfse);
  fprintf(outfile, "MRIStep: nfsi = %li\n", step_mem->nfsi);
  fprintf(outfile, "MRIStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile, "MRIStep: nstlp = %li\n", step_mem->nstlp);
  fprintf(outfile, "MRIStep: nls_iters = %li\n", step_mem->nls_iters);

  /* output boolean quantities */
  fprintf(outfile, "MRIStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile, "MRIStep: user_linear_timedep = %i\n",
          step_mem->linear_timedep);
  fprintf(outfile, "MRIStep: user_explicit = %i\n", step_mem->explicit_rhs);
  fprintf(outfile, "MRIStep: user_implicit = %i\n", step_mem->implicit_rhs);
  fprintf(outfile, "MRIStep: jcur = %i\n", step_mem->jcur);
  fprintf(outfile, "MRIStep: ownNLS = %i\n", step_mem->ownNLS);

  /* output sunrealtype quantities */
  fprintf(outfile, "MRIStep: Coupling structure:\n");
  MRIStepCoupling_Write(step_mem->MRIC, outfile);

  fprintf(outfile, "MRIStep: gamma = %" RSYM "\n", step_mem->gamma);
  fprintf(outfile, "MRIStep: gammap = %" RSYM "\n", step_mem->gammap);
  fprintf(outfile, "MRIStep: gamrat = %" RSYM "\n", step_mem->gamrat);
  fprintf(outfile, "MRIStep: crate = %" RSYM "\n", step_mem->crate);
  fprintf(outfile, "MRIStep: delp = %" RSYM "\n", step_mem->delp);
  fprintf(outfile, "MRIStep: eRNrm = %" RSYM "\n", step_mem->eRNrm);
  fprintf(outfile, "MRIStep: nlscoef = %" RSYM "\n", step_mem->nlscoef);
  fprintf(outfile, "MRIStep: crdown = %" RSYM "\n", step_mem->crdown);
  fprintf(outfile, "MRIStep: rdiv = %" RSYM "\n", step_mem->rdiv);
  fprintf(outfile, "MRIStep: dgmax = %" RSYM "\n", step_mem->dgmax);
  fprintf(outfile, "MRIStep: Ae_row =");
  for (i = 0; i < step_mem->nstages_active; i++)
  {
    fprintf(outfile, " %" RSYM, step_mem->Ae_row[i]);
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "MRIStep: Ai_row =");
  for (i = 0; i < step_mem->nstages_active; i++)
  {
    fprintf(outfile, " %" RSYM, step_mem->Ai_row[i]);
  }
  fprintf(outfile, "\n");

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "MRIStep: sdata:\n");
  N_VPrintFile(step_mem->sdata, outfile);
  fprintf(outfile, "MRIStep: zpred:\n");
  N_VPrintFile(step_mem->zpred, outfile);
  fprintf(outfile, "MRIStep: zcor:\n");
  N_VPrintFile(step_mem->zcor, outfile);
  if (step_mem->Fse)
    for (i = 0; i < step_mem->nstages_active; i++)
    {
      fprintf(outfile, "MRIStep: Fse[%i]:\n", i);
      N_VPrintFile(step_mem->Fse[i], outfile);
    }
  if (step_mem->Fsi && !step_mem->unify_Fs)
    for (i = 0; i < step_mem->nstages_active; i++)
    {
      fprintf(outfile, "MRIStep: Fsi[%i]:\n", i);
      N_VPrintFile(step_mem->Fsi[i], outfile);
    }
#endif

  /* print the inner stepper memory */
  mriStepInnerStepper_PrintMem(step_mem->stepper, outfile);
  return;
}

/*---------------------------------------------------------------
  mriStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  MRIStep module.
  ---------------------------------------------------------------*/
int mriStep_AttachLinsol(ARKodeMem ark_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup, ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type, void* lmem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* free any existing system solver */
  if (step_mem->lfree != NULL) { step_mem->lfree(ark_mem); }

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit  = linit;
  step_mem->lsetup = lsetup;
  step_mem->lsolve = lsolve;
  step_mem->lfree  = lfree;
  step_mem->lmem   = lmem;

  /* Reset all linear solver counters */
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  MRIStep module.
  ---------------------------------------------------------------*/
void mriStep_DisableLSetup(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* nullify the lsetup function pointer */
  step_mem->lsetup = NULL;
}

/*---------------------------------------------------------------
  mriStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.
  ---------------------------------------------------------------*/
void* mriStep_GetLmem(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure, and return lmem */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (NULL); }
  return (step_mem->lmem);
}

/*---------------------------------------------------------------
  mriStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
ARKRhsFn mriStep_GetImplicitRHS(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure, and return fi */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (NULL); }
  if (step_mem->implicit_rhs) { return (step_mem->fsi); }
  else { return (NULL); }
}

/*---------------------------------------------------------------
  mriStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
int mriStep_GetGammas(ARKodeMem ark_mem, sunrealtype* gamma, sunrealtype* gamrat,
                      sunbooleantype** jcur, sunbooleantype* dgamma_fail)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set outputs */
  *gamma       = step_mem->gamma;
  *gamrat      = step_mem->gamrat;
  *jcur        = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the coefficient tables to be used
  - allocates any internal memory that depend on the MRI method
    structure or solver options

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int mriStep_Init(ARKodeMem ark_mem, int init_type)
{
  ARKodeMRIStepMem step_mem;
  int retval, j;
  sunbooleantype reset_efun;
  SUNAdaptController_Type adapt_type;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return (ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    /* enforce use of arkEwtSmallReal if using a fixed step size for
       an explicit method, an internal error weight function, and not performing
       accumulated temporal error estimation */
    reset_efun = SUNTRUE;
    if (step_mem->implicit_rhs) { reset_efun = SUNFALSE; }
    if (!ark_mem->fixedstep) { reset_efun = SUNFALSE; }
    if (ark_mem->user_efun) { reset_efun = SUNFALSE; }
    if (ark_mem->AccumErrorType >= 0) { reset_efun = SUNFALSE; }
    if (reset_efun)
    {
      ark_mem->user_efun = SUNFALSE;
      ark_mem->efun      = arkEwtSetSmallReal;
      ark_mem->e_data    = ark_mem;
    }

    /* Create coupling structure (if not already set) */
    retval = mriStep_SetCoupling(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Could not create coupling table");
      return (ARK_ILL_INPUT);
    }

    /* Check that coupling structure is OK */
    retval = mriStep_CheckCoupling(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in coupling table");
      return (ARK_ILL_INPUT);
    }

    /* Attach correct TakeStep routine for this coupling table */
    switch (step_mem->MRIC->type)
    {
    case MRISTEP_EXPLICIT:
    case MRISTEP_IMPLICIT:
    case MRISTEP_IMEX: ark_mem->step = mriStep_TakeStepMRIGARK; break;
    case MRISTEP_MERK: ark_mem->step = mriStep_TakeStepMERK; break;
    case MRISTEP_MRISR: ark_mem->step = mriStep_TakeStepMRISR; break;
    }

    /* Retrieve/store method and embedding orders now that tables are finalized */
    step_mem->stages = step_mem->MRIC->stages;
    step_mem->q = ark_mem->hadapt_mem->q = step_mem->MRIC->q;
    step_mem->p = ark_mem->hadapt_mem->p = step_mem->MRIC->p;

    /* Ensure that if adaptivity or error accumulation is enabled, then
       method includes embedding coefficients */
    if ((!ark_mem->fixedstep || (ark_mem->AccumErrorType >= 0)) &&
        (step_mem->p == 0))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Temporal error estimation cannot be performed without embedding coefficients");
      return (ARK_ILL_INPUT);
    }

    /* allocate/fill derived quantities from MRIC structure */

    /* stage map */
    if (step_mem->stage_map)
    {
      free(step_mem->stage_map);
      ark_mem->liw -= step_mem->stages;
    }
    step_mem->stage_map = (int*)calloc(step_mem->MRIC->stages,
                                       sizeof(*step_mem->stage_map));
    if (step_mem->stage_map == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->liw += step_mem->MRIC->stages;
    retval = mriStepCoupling_GetStageMap(step_mem->MRIC, step_mem->stage_map,
                                         &(step_mem->nstages_active));
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in coupling table");
      return (ARK_ILL_INPUT);
    }

    /* stage types */
    if (step_mem->stagetypes)
    {
      free(step_mem->stagetypes);
      ark_mem->liw -= step_mem->stages;
    }
    step_mem->stagetypes = (int*)calloc(step_mem->MRIC->stages + 1,
                                        sizeof(*step_mem->stagetypes));
    if (step_mem->stagetypes == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->liw += (step_mem->MRIC->stages + 1);
    step_mem->stagetypes[0] =
      MRISTAGE_ERK_NOFAST; /* 1st stage is always ERK_NOFAST */
    for (j = 1; j <= step_mem->MRIC->stages; j++)
    {
      step_mem->stagetypes[j] = mriStepCoupling_GetStageType(step_mem->MRIC, j);
    }

    /* explicit RK coefficient row */
    if (step_mem->Ae_row)
    {
      free(step_mem->Ae_row);
      ark_mem->lrw -= step_mem->stages;
    }
    step_mem->Ae_row = (sunrealtype*)calloc(step_mem->MRIC->stages,
                                            sizeof(*step_mem->Ae_row));
    if (step_mem->Ae_row == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->lrw += step_mem->MRIC->stages;

    /* implicit RK coefficient row */
    if (step_mem->Ai_row)
    {
      free(step_mem->Ai_row);
      ark_mem->lrw -= step_mem->stages;
    }
    step_mem->Ai_row = (sunrealtype*)calloc(step_mem->MRIC->stages,
                                            sizeof(*step_mem->Ai_row));
    if (step_mem->Ai_row == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->lrw += step_mem->MRIC->stages;

    /* Allocate reusable arrays for fused vector interface */
    if (step_mem->cvals)
    {
      free(step_mem->cvals);
      ark_mem->lrw -= step_mem->nfusedopvecs;
    }
    if (step_mem->Xvecs)
    {
      free(step_mem->Xvecs);
      ark_mem->liw -= step_mem->nfusedopvecs;
    }
    step_mem->nfusedopvecs = 2 * step_mem->MRIC->stages + 2;
    step_mem->cvals        = (sunrealtype*)calloc(step_mem->nfusedopvecs,
                                                  sizeof(*step_mem->cvals));
    if (step_mem->cvals == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->lrw += step_mem->nfusedopvecs;

    step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs,
                                        sizeof(*step_mem->Xvecs));
    if (step_mem->Xvecs == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
    ark_mem->liw += step_mem->nfusedopvecs;

    /* Retrieve/store method and embedding orders now that tables are finalized */
    step_mem->stages = step_mem->MRIC->stages;
    step_mem->q      = step_mem->MRIC->q;
    step_mem->p      = step_mem->MRIC->p;

    /* If an MRISR method is applied to a non-ImEx problem, we "unify"
       the Fse and Fsi vectors to point at the same memory */
    step_mem->unify_Fs = SUNFALSE;
    if ((step_mem->MRIC->type == MRISTEP_MRISR) &&
        ((step_mem->explicit_rhs && !step_mem->implicit_rhs) ||
         (!step_mem->explicit_rhs && step_mem->implicit_rhs)))
    {
      step_mem->unify_Fs = SUNTRUE;
    }

    /* Allocate MRI RHS vector memory, update storage requirements */
    /*   Allocate Fse[0] ... Fse[nstages_active - 1] and           */
    /*   Fsi[0] ... Fsi[nstages_active - 1] if needed              */
    if (step_mem->nstages_allocated < step_mem->nstages_active)
    {
      if (step_mem->nstages_allocated)
      {
        if (step_mem->explicit_rhs)
        {
          arkFreeVecArray(step_mem->nstages_allocated, &(step_mem->Fse),
                          ark_mem->lrw1, &(ark_mem->lrw), ark_mem->liw1,
                          &(ark_mem->liw));
          if (step_mem->unify_Fs) { step_mem->Fsi = NULL; }
        }
        if (step_mem->implicit_rhs)
        {
          arkFreeVecArray(step_mem->nstages_allocated, &(step_mem->Fsi),
                          ark_mem->lrw1, &(ark_mem->lrw), ark_mem->liw1,
                          &(ark_mem->liw));
          if (step_mem->unify_Fs) { step_mem->Fse = NULL; }
        }
      }
      if (step_mem->explicit_rhs && !step_mem->unify_Fs)
      {
        if (!arkAllocVecArray(step_mem->nstages_active, ark_mem->ewt,
                              &(step_mem->Fse), ark_mem->lrw1, &(ark_mem->lrw),
                              ark_mem->liw1, &(ark_mem->liw)))
        {
          return (ARK_MEM_FAIL);
        }
      }
      if (step_mem->implicit_rhs && !step_mem->unify_Fs)
      {
        if (!arkAllocVecArray(step_mem->nstages_active, ark_mem->ewt,
                              &(step_mem->Fsi), ark_mem->lrw1, &(ark_mem->lrw),
                              ark_mem->liw1, &(ark_mem->liw)))
        {
          return (ARK_MEM_FAIL);
        }
      }
      if (step_mem->unify_Fs)
      {
        if (!arkAllocVecArray(step_mem->nstages_active, ark_mem->ewt,
                              &(step_mem->Fse), ark_mem->lrw1, &(ark_mem->lrw),
                              ark_mem->liw1, &(ark_mem->liw)))
        {
          return (ARK_MEM_FAIL);
        }
        step_mem->Fsi = step_mem->Fse;
      }

      step_mem->nstages_allocated = step_mem->nstages_active;
    }

    /* if any slow stage is implicit, allocate sdata, zpred, zcor vectors;
       if all stages explicit, free default NLS object, and detach all
       linear solver routines.  Note: step_mem->implicit_rhs will only equal
       SUNTRUE if an implicit table has been user-provided. */
    if (step_mem->implicit_rhs)
    {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->sdata)))
      {
        return (ARK_MEM_FAIL);
      }
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->zpred)))
      {
        return (ARK_MEM_FAIL);
      }
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->zcor)))
      {
        return (ARK_MEM_FAIL);
      }
    }
    else
    {
      if ((step_mem->NLS != NULL) && (step_mem->ownNLS))
      {
        SUNNonlinSolFree(step_mem->NLS);
        step_mem->NLS    = NULL;
        step_mem->ownNLS = SUNFALSE;
      }
      step_mem->linit  = NULL;
      step_mem->lsetup = NULL;
      step_mem->lsolve = NULL;
      step_mem->lfree  = NULL;
      step_mem->lmem   = NULL;
    }

    /* Free existing reusable arrays if they are the wrong size */
    if (step_mem->nfusedopvecs != 2 * step_mem->stages + 2)
    {
      if (step_mem->cvals) { free(step_mem->cvals); }
      if (step_mem->Xvecs) { free(step_mem->Xvecs); }
      step_mem->cvals = NULL;
      step_mem->Xvecs = NULL;
    }

    /* Allocate reusable arrays for fused vector interface */
    step_mem->nfusedopvecs = 2 * step_mem->stages + 2;
    if (step_mem->cvals == NULL)
    {
      step_mem->cvals = (sunrealtype*)calloc(step_mem->nfusedopvecs,
                                             sizeof(sunrealtype));
      if (step_mem->cvals == NULL) { return (ARK_MEM_FAIL); }
      ark_mem->lrw += (step_mem->nfusedopvecs);
    }
    if (step_mem->Xvecs == NULL)
    {
      step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs,
                                          sizeof(N_Vector));
      if (step_mem->Xvecs == NULL) { return (ARK_MEM_FAIL); }
      ark_mem->liw += (step_mem->nfusedopvecs); /* pointers */
    }

    /* Allocate inner stepper data */
    retval = mriStepInnerStepper_AllocVecs(step_mem->stepper,
                                           step_mem->MRIC->nmat, ark_mem->ewt);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error allocating inner stepper memory");
      return (ARK_MEM_FAIL);
    }

    /* Limit max interpolant degree (negative input only overwrites the current
       interpolant degree if it is greater than abs(input). */
    if (ark_mem->interp != NULL)
    {
      if (step_mem->q > 1)
      {
        /* Limit max degree to at most one less than the method global order */
        retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q - 1));
      }
      else
      {
        /* Allow for linear interpolant with first order methods to ensure
           solution values are returned at the time interval end points */
        retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q));
      }

      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        "Unable to update interpolation polynomial degree");
        return (ARK_ILL_INPUT);
      }
    }
  }

  /* Call linit (if it exists) */
  if (step_mem->linit)
  {
    retval = step_mem->linit(ark_mem);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_LINIT_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_LINIT_FAIL);
      return (ARK_LINIT_FAIL);
    }
  }

  /* Initialize the nonlinear solver object (if it exists) */
  if (step_mem->NLS)
  {
    retval = mriStep_NlsInit(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to initialize SUNNonlinearSolver object");
      return (ARK_NLS_INIT_FAIL);
    }
  }

  /* get timestep adaptivity type, and return an error if an
     incompatible type is detected */
  adapt_type = SUNAdaptController_GetType(ark_mem->hadapt_mem->hcontroller);
  if ((adapt_type != SUN_ADAPTCONTROLLER_MRI_H) &&
      (adapt_type != SUN_ADAPTCONTROLLER_MRI_TOL) &&
      (adapt_type != SUN_ADAPTCONTROLLER_H))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "SUNAdaptController type is unsupported by MRIStep");
    return (ARK_ILL_INPUT);
  }

  /*** Perform timestep adaptivity checks and initial setup ***/

  if (ark_mem->fixedstep)
  {
    /* Non-adaptive controller: user must have supplied initial step
       size, and indicated fixed time stepping */
    if ((ark_mem->hin == ZERO) || (ark_mem->fixedstep == SUNFALSE))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Timestep adaptivity disabled, but missing user-defined fixed stepsize");
      return (ARK_ILL_INPUT);
    }
  }
  else
  {
    /* Controller provides adaptivity (at least at the slow time scale):
       - verify that the MRI method includes an embedding, and
       - estimate initial slow step size (store in ark_mem->hin) */
    if (step_mem->MRIC->p == 0)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Timestep adaptivity enabled, but non-embedded MRI table specified");
      return (ARK_ILL_INPUT);
    }
    if (ark_mem->hin == ZERO)
    {
      /*   tempv1 = fslow(t0, y0) */
      if (mriStep_SlowRHS(ark_mem, ark_mem->tcur, ark_mem->yn, ark_mem->tempv1,
                          ARK_FULLRHS_START) != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        "error calling slow RHS function(s)");
        return (ARK_RHSFUNC_FAIL);
      }
      retval = mriStep_Hin(ark_mem, ark_mem->tcur, ark_mem->tout, ark_mem->yn,
                           ark_mem->tempv1, ark_mem->ycur, ark_mem->tempv2,
                           ark_mem->tempv3, mriStep_SlowRHS, &(ark_mem->hin));
      if (retval != ARK_SUCCESS)
      {
        retval = arkHandleFailure(ark_mem, retval);
        return (retval);
      }
    }
  }

  /* Perform additional setup for (H,tol) controller */
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL)
  {
    /* Verify that adaptivity type is supported by inner stepper */
    if (!mriStepInnerStepper_SupportsRTolAdaptivity(step_mem->stepper))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "MRI-TOL SUNAdaptController provided, but unsupported by inner stepper");
      return (ARK_ILL_INPUT);
    }

    /* initialize fast stepper to use the same relative tolerance as MRIStep */
    step_mem->inner_control = ONE;
  }

  /* Perform additional setup for (H,h) controller */
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H)
  {
    /* verify that adaptivity type is supported by inner stepper */
    if (!mriStepInnerStepper_SupportsStepAdaptivity(step_mem->stepper))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "MRI-H SUNAdaptController provided, but unsupported by inner stepper");
      return (ARK_ILL_INPUT);
    }

    /* If user has left inner_hfactor unset, reset to our default to indicate
       that we **do not trust** fixed-stepsize fast error estimates. */
    if (step_mem->inner_hfactor < ZERO)
    {
      step_mem->inner_hfactor = INNER_HFACTOR;
    }

    /* initialize fast stepper fixed step size (store in step_mem->inner_control) */
    if (step_mem->stepper->ops->fullrhs)
    {
      /* tempv1 = ffast(t0, y0) */
      if (mriStep_FastRHS(ark_mem, ark_mem->tcur, ark_mem->yn, ark_mem->tempv1,
                          ARK_FULLRHS_START) != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        "error calling fast RHS function(s)");
        return (ARK_RHSFUNC_FAIL);
      }

      /* compute initial fast step size */
      retval = arkHin(ark_mem, ark_mem->tcur, ark_mem->tout, ark_mem->yn,
                      ark_mem->tempv1, ark_mem->ycur, ark_mem->tempv2,
                      ark_mem->tempv3, mriStep_FastRHS,
                      &(step_mem->inner_control));
      if (retval != ARK_SUCCESS)
      {
        retval = arkHandleFailure(ark_mem, retval);
        return (retval);
      }
    }
    else
    {
      /* set step_mem->inner_control to equal H/100 */
      step_mem->inner_control = SUN_RCONST(0.01) * ark_mem->hin;
    }
    /* Pass fixed stepsize to inner stepper */
    retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                              step_mem->inner_control);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "error setting fast fixed step");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  else
  {
    /* If user has left inner_hfactor unset, reset to zero to indicate
       either that we **trust** tolerance-based adaptive fast error
       estimates, or that the method does not use fast error estimates. */
    if (step_mem->inner_hfactor < ZERO) { step_mem->inner_hfactor = ZERO; }
  }

  /* /\* Signal to shared arkode module that fullrhs is required after each step *\/ */
  /* /\* TO-DO: verify whether this is actually required *\/ */
  /* ark_mem->call_fullrhs = SUNTRUE; */
  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  mriStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y)  + ff(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If this function is called in ARK_FULLRHS_START or ARK_FULLRHS_END mode and
  evaluating the RHS functions is necessary, we store the vectors fse(t,y) and
  fsi(t,y) in Fse[0] and Fsi[0] for possible reuse in the first stage of the
  subsequent time step.

  ARK_FULLRHS_OTHER mode is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to store the
  intermediate parts so that they do not interfere with the other two modes.

  Presently ff(t,y) is always called with ARK_FULLRHS_OTHER mode.
  ----------------------------------------------------------------------------*/
int mriStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode)
{
  ARKodeMRIStepMem step_mem;
  int i, sa_stage, retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* ensure that inner stepper provides fullrhs function */
  if (!(step_mem->stepper->ops->fullrhs))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_MISSING_FULLRHS);
    return ARK_ILL_INPUT;
  }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:

    /* compute the full RHS */
    if (!(ark_mem->fn_is_current))
    {
      /* compute the explicit component */
      if (step_mem->explicit_rhs)
      {
        retval = step_mem->fse(t, y, step_mem->Fse[0], ark_mem->user_data);
        step_mem->nfse++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
          return (ARK_RHSFUNC_FAIL);
        }
      }

      /* compute the implicit component */
      if (step_mem->implicit_rhs)
      {
        retval = step_mem->fsi(t, y, step_mem->Fsi[0], ark_mem->user_data);
        step_mem->nfsi++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
          return (ARK_RHSFUNC_FAIL);
        }
      }

      /* compute the fast component (force new RHS computation) */
      retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                           ARK_FULLRHS_OTHER);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs)
    {
      /* ImEx */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    }
    else if (step_mem->implicit_rhs)
    {
      /* implicit */
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    }
    else
    {
      /* explicit */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
    }

    break;

  case ARK_FULLRHS_END:

    /* compute the full RHS */
    if (!(ark_mem->fn_is_current))
    {
      /* if the method had a stiffly-accurate internal stage, use the
         already-computed RHS vectors */
      sa_stage = -1;
      for (i = 1; i < step_mem->stages; i++)
      {
        if (step_mem->stagetypes[i] == MRISTAGE_STIFF_ACC)
        {
          sa_stage = step_mem->stage_map[i - 1];
        }
      }
      if (sa_stage > -1)
      {
        /* copy the explicit component */
        if (step_mem->explicit_rhs)
        {
          N_VScale(ONE, step_mem->Fse[sa_stage], step_mem->Fse[0]);
        }

        /* copy the implicit component */
        if (step_mem->implicit_rhs)
        {
          N_VScale(ONE, step_mem->Fsi[sa_stage], step_mem->Fsi[0]);
        }
      }
      else
      {
        /* compute the explicit component */
        if (step_mem->explicit_rhs)
        {
          retval = step_mem->fse(t, y, step_mem->Fse[0], ark_mem->user_data);
          step_mem->nfse++;
          if (retval != 0)
          {
            arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                            __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
            return (ARK_RHSFUNC_FAIL);
          }
        }

        /* compute the implicit component */
        if (step_mem->implicit_rhs)
        {
          retval = step_mem->fsi(t, y, step_mem->Fsi[0], ark_mem->user_data);
          step_mem->nfsi++;
          if (retval != 0)
          {
            arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                            __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
            return (ARK_RHSFUNC_FAIL);
          }
        }
      }

      /* compute the fast component (force new RHS computation) */
      retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                           ARK_FULLRHS_OTHER);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs)
    {
      /* ImEx */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    }
    else if (step_mem->implicit_rhs)
    {
      /* implicit */
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    }
    else
    {
      /* explicit */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
    }

    break;

  case ARK_FULLRHS_OTHER:

    /* compute the explicit component and store in ark_tempv2 */
    if (step_mem->explicit_rhs)
    {
      retval = step_mem->fse(t, y, ark_mem->tempv2, ark_mem->user_data);
      step_mem->nfse++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* compute the implicit component and store in sdata */
    if (step_mem->implicit_rhs)
    {
      retval = step_mem->fsi(t, y, step_mem->sdata, ark_mem->user_data);
      step_mem->nfsi++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* compute the fast component (force new RHS computation) */
    retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                         ARK_FULLRHS_OTHER);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs)
    {
      /* ImEx */
      N_VLinearSum(ONE, ark_mem->tempv2, ONE, f, f);
      N_VLinearSum(ONE, step_mem->sdata, ONE, f, f);
    }
    else if (step_mem->implicit_rhs)
    {
      /* implicit */
      N_VLinearSum(ONE, step_mem->sdata, ONE, f, f);
    }
    else
    {
      /* explicit */
      N_VLinearSum(ONE, ark_mem->tempv2, ONE, f, f);
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

/*---------------------------------------------------------------
  mriStep_TakeStepMRIGARK:

  This routine serves the primary purpose of the MRIStep module:
  it performs a single MRI step (with embedding, if possible).

  The vector ark_mem->yn holds the previous time-step solution
  on input, and the vector ark_mem->ycur should hold the result
  of this step on output.

  If timestep adaptivity is enabled, this routine also computes
  the error estimate y-ytilde, where ytilde is the
  embedded solution, and the norm weights come from ark_ewt.
  This esimate is stored in ark_mem->tempv1, in case the calling
  routine wishes to examine the error locations.

  The output variable dsmPtr should contain a scalar-valued
  estimate of the temporal error from this step, ||y-ytilde||_WRMS
  if timestep adaptivity is enabled; otherwise it should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure
  Since the fast-scale evolution could be considered a different
  type of "algebraic solver", we similarly report any fast-scale
  evolution error as a recoverable nflagPtr value.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int mriStep_TakeStepMRIGARK(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeMRIStepMem step_mem;          /* outer stepper memory       */
  int is;                             /* current stage index        */
  int retval;                         /* reusable return flag       */
  N_Vector tmp;                       /* N_Vector pointer           */
  SUNAdaptController_Type adapt_type; /* timestep adaptivity type   */
  sunrealtype t0, tf;                 /* start/end of each stage    */
  sunbooleantype calc_fslow;

  /* access the MRIStep mem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* initialize algebraic solver convergence flag to success;
     error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* initial time and state for step */
  N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
  t0 = ark_mem->tn;

  /* if MRI adaptivity is enabled: reset fast accumulated error,
     and send appropriate control parameter to the fast integrator */
  adapt_type = SUNAdaptController_GetType(ark_mem->hadapt_mem->hcontroller);
  if ((adapt_type == SUN_ADAPTCONTROLLER_MRI_H) ||
      (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL))
  {
    step_mem->inner_dsm = ZERO;
    retval              = mriStepInnerStepper_ResetError(step_mem->stepper);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper error estimate");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H)
  {
    retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                              step_mem->inner_control);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set a fixed stepsize in the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL)
  {
    retval =
      mriStepInnerStepper_SetRTol(step_mem->stepper,
                                  step_mem->inner_control * ark_mem->reltol);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set the inner stepper tolerance");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* if any temporal adaptivity is enabled: reset the inner integrator
     to the beginning of this step (in case of recomputation) */
  if (!ark_mem->fixedstep)
  {
    retval = mriStepInnerStepper_Reset(step_mem->stepper, t0, ark_mem->yn);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* call nonlinear solver setup if it exists */
  if (step_mem->NLS)
  {
    if ((step_mem->NLS)->ops->setup)
    {
      N_VConst(ZERO, ark_mem->tempv3); /* set guess to 0 */
      retval = SUNNonlinSolSetup(step_mem->NLS, ark_mem->tempv3, ark_mem);
      if (retval < 0) { return (ARK_NLS_SETUP_FAIL); }
      if (retval > 0) { return (ARK_NLS_SETUP_RECVR); }
    }
  }

  /* Evaluate the slow RHS functions if needed. NOTE: We do not use the full RHS
     function here (unlike ERKStep and ARKStep) since it does not need to check
     for FSAL or SA methods and thus avoids potentially unnecessary evaluations
     of the inner (fast) RHS function */

  if (!(ark_mem->fn_is_current))
  {
    /* compute the explicit component */
    if (step_mem->explicit_rhs)
    {
      retval = step_mem->fse(t0, ark_mem->yn, step_mem->Fse[0],
                             ark_mem->user_data);
      step_mem->nfse++;
      if (retval) { return ARK_RHSFUNC_FAIL; }
    }

    /* compute the implicit component */
    if (step_mem->implicit_rhs)
    {
      retval = step_mem->fsi(t0, ark_mem->yn, step_mem->Fsi[0],
                             ark_mem->user_data);
      step_mem->nfsi++;
      if (retval) { return ARK_RHSFUNC_FAIL; }
    }

    ark_mem->fn_is_current = SUNTRUE;
  }

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep step %li,  stage 0,  h = %" RSYM ",  t_n = %" RSYM "\n",
         ark_mem->nst, ark_mem->h, t0);
#endif

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRIGARK", "slow stage",
                     "z[0] =", "");
  N_VPrintFile(ark_mem->yn, ARK_LOGGER->debug_fp);

  if (step_mem->explicit_rhs)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "slow explicit RHS",
                       "Fse[0] =", "");
    N_VPrintFile(step_mem->Fse[0], ARK_LOGGER->debug_fp);
  }
  if (step_mem->implicit_rhs)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "slow implicit RHS",
                       "Fsi[0] =", "");
    N_VPrintFile(step_mem->Fsi[0], ARK_LOGGER->debug_fp);
  }
#endif

  /* The first stage is the previous time-step solution, so its RHS
     is the [already-computed] slow RHS from the start of the step */

  /* Loop over remaining stages */
  for (is = 1; is < step_mem->stages; is++)
  {
    /* Set relevant stage times (including desired stage time for implicit solves) */
    t0 = ark_mem->tn + step_mem->MRIC->c[is - 1] * ark_mem->h;
    tf = ark_mem->tcur = ark_mem->tn + step_mem->MRIC->c[is] * ark_mem->h;

    /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "start-stage",
                       "step = %li, stage = %i, stage type = %d, h = %" RSYM
                       ", tcur = %" RSYM,
                       ark_mem->nst, is, step_mem->stagetypes[is], ark_mem->h,
                       ark_mem->tcur);
#endif

    /* Determine current stage type, and call corresponding routine; the
       vector ark_mem->ycur stores the previous stage solution on input, and
       should store the result of this stage solution on output. */
    switch (step_mem->stagetypes[is])
    {
    case (MRISTAGE_ERK_FAST):
      retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, is, t0, tf);
      if (retval != ARK_SUCCESS)
      {
        *nflagPtr = CONV_FAIL;
        break;
      }
      retval = mriStep_StageERKFast(ark_mem, step_mem, is, t0, tf, ark_mem->ycur,
                                    ark_mem->tempv2, SUNFALSE, SUNTRUE);
      if (retval != ARK_SUCCESS) { *nflagPtr = CONV_FAIL; }
      break;
    case (MRISTAGE_ERK_NOFAST):
      retval = mriStep_StageERKNoFast(ark_mem, step_mem, is);
      break;
    case (MRISTAGE_DIRK_NOFAST):
      retval = mriStep_StageDIRKNoFast(ark_mem, step_mem, is, nflagPtr);
      break;
    case (MRISTAGE_DIRK_FAST):
      retval = mriStep_StageDIRKFast(ark_mem, step_mem, is, nflagPtr);
      break;
    case (MRISTAGE_STIFF_ACC): retval = ARK_SUCCESS; break;
    }
    if (retval != ARK_SUCCESS) { return (retval); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "slow stage",
                       "z[%i] =", is);
    N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

    /* apply user-supplied stage postprocessing function (if supplied) */
    if ((ark_mem->ProcessStage != NULL) &&
        (step_mem->stagetypes[is] != MRISTAGE_STIFF_ACC))
    {
      retval = ark_mem->ProcessStage(tf, ark_mem->ycur, ark_mem->user_data);
      if (retval != 0) { return (ARK_POSTPROCESS_STAGE_FAIL); }
    }

    /* conditionally reset the inner integrator with the modified stage solution */
    if (step_mem->stagetypes[is] != MRISTAGE_STIFF_ACC)
    {
      if ((step_mem->stagetypes[is] != MRISTAGE_ERK_FAST) ||
          (ark_mem->ProcessStage != NULL))
      {
        retval = mriStepInnerStepper_Reset(step_mem->stepper, tf, ark_mem->ycur);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to reset the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }
      }
    }

    /* Compute updated slow RHS, except:
       1. at last stage which is the new solution, since that RHS occurs in arkCompleteStep.
       2. if the stage is excluded from stage_map
       3. if the next stage has "STIFF_ACC" type, and temporal estimation is disabled */
    calc_fslow = SUNTRUE;
    if (is == step_mem->stages - 1) { calc_fslow = SUNFALSE; }
    if (step_mem->stage_map[is] == -1) { calc_fslow = SUNFALSE; }
    if (is < step_mem->stages - 1)
    {
      if (ark_mem->fixedstep && (ark_mem->AccumErrorType < 0) &&
          (step_mem->stagetypes[is + 1] == MRISTAGE_STIFF_ACC))
      {
        calc_fslow = SUNFALSE;
      }
    }
    if (calc_fslow)
    {
      /* store explicit slow rhs */
      if (step_mem->explicit_rhs)
      {
        retval = step_mem->fse(tf, ark_mem->ycur,
                               step_mem->Fse[step_mem->stage_map[is]],
                               ark_mem->user_data);
        step_mem->nfse++;
        if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
        if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRIGARK",
                           "slow explicit RHS", "Fse[%i] =", is);
        N_VPrintFile(step_mem->Fse[step_mem->stage_map[is]],
                     ARK_LOGGER->debug_fp);
#endif
      }

      /* store implicit slow rhs  */
      if (step_mem->implicit_rhs)
      {
        if (!step_mem->deduce_rhs ||
            (step_mem->stagetypes[is] != MRISTAGE_DIRK_NOFAST))
        {
          retval = step_mem->fsi(tf, ark_mem->ycur,
                                 step_mem->Fsi[step_mem->stage_map[is]],
                                 ark_mem->user_data);
          step_mem->nfsi++;
        }
        else
        {
          N_VLinearSum(ONE / step_mem->gamma, step_mem->zcor,
                       -ONE / step_mem->gamma, step_mem->sdata,
                       step_mem->Fsi[step_mem->stage_map[is]]);
        }

        if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
        if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRIGARK",
                           "slow implicit RHS", "Fsi[%i] =", is);
        N_VPrintFile(step_mem->Fsi[step_mem->stage_map[is]],
                     ARK_LOGGER->debug_fp);
#endif
      }

      /* if this is the next-to-last internal stage, and if temporal error estimation
         is enabled, save ycur into ark_mem->tempv4 (may be needed for embedding) */
      if ((is == step_mem->stages - 2) &&
          (!ark_mem->fixedstep || (ark_mem->AccumErrorType >= 0)))
      {
        N_VScale(ONE, ark_mem->ycur, ark_mem->tempv4);
      }

    } /* compute slow RHS */
  }   /* loop over stages */

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRIGARK", "updated solution",
                     "ycur =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

  /* if temporal error estimation is enabled: compute slow time scale error and store in dsmPtr */
  if (!ark_mem->fixedstep || (ark_mem->AccumErrorType >= 0))
  {
    /* Temporarily swap ark_mem->ycur and ark_mem->tempv4 pointers, so that
       during this embedding "stage":
         - ark_mem->ycur starts out as what it was for the last internal stage,
           and ends up as the embedded solution.
         - ark_mem->tempv4 archives the current time step solution vector. */
    tmp             = ark_mem->ycur;
    ark_mem->ycur   = ark_mem->tempv4;
    ark_mem->tempv4 = tmp;

    /* Reset ark_mem->tcur as the time value corresponding with the end of the step */
    /* Set relevant stage times (including desired stage time for implicit solves) */
    t0 = ark_mem->tn + step_mem->MRIC->c[step_mem->stages - 2] * ark_mem->h;
    tf = ark_mem->tcur = ark_mem->tn +
                         step_mem->MRIC->c[step_mem->stages - 1] * ark_mem->h;

    /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "embedding-stage",
                       "step = %li, stage = %i, stage type = %d, h = %" RSYM
                       ", tcur = %" RSYM,
                       ark_mem->nst, step_mem->stages,
                       step_mem->stagetypes[step_mem->stages], ark_mem->h,
                       ark_mem->tcur);
#endif

    /* Determine embedding stage type, and call corresponding routine; the
       vector ark_mem->ycur stores the previous stage solution on input, and
       should store the result of this stage solution on output. */
    switch (step_mem->stagetypes[step_mem->stages])
    {
    case (MRISTAGE_ERK_FAST):
      retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, step_mem->stages,
                                           t0, tf);
      if (retval != ARK_SUCCESS)
      {
        *nflagPtr = CONV_FAIL;
        break;
      }
      retval = mriStep_StageERKFast(ark_mem, step_mem, step_mem->stages, t0, tf,
                                    ark_mem->ycur, ark_mem->tempv2, SUNTRUE,
                                    SUNFALSE);
      if (retval != ARK_SUCCESS) { *nflagPtr = CONV_FAIL; }
      break;
    case (MRISTAGE_ERK_NOFAST):
      retval = mriStep_StageERKNoFast(ark_mem, step_mem, step_mem->stages);
      break;
    case (MRISTAGE_DIRK_NOFAST):
      retval = mriStep_StageDIRKNoFast(ark_mem, step_mem, step_mem->stages,
                                       nflagPtr);
      break;
    case (MRISTAGE_DIRK_FAST):
      retval = mriStep_StageDIRKFast(ark_mem, step_mem, step_mem->stages,
                                     nflagPtr);
      break;
    }
    if (retval != ARK_SUCCESS) { return (retval); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRIGARK", "embedded solution",
                       "ytilde =", "");
    N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

    /* Compute temporal error estimate via difference between step
       solution and embedding, store in ark_mem->tempv1, and take norm. */
    N_VLinearSum(ONE, ark_mem->tempv4, -ONE, ark_mem->ycur, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);

    /* Swap back ark_mem->ycur with ark_mem->tempv4. */
    tmp             = ark_mem->ycur;
    ark_mem->ycur   = ark_mem->tempv4;
    ark_mem->tempv4 = tmp;
  }

  /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRIGARK", "end-step",
                     "step = %li, t = %" RSYM ", h = %" RSYM ", dsm = %" RSYM,
                     ark_mem->nst, ark_mem->tn, ark_mem->h, *dsmPtr);
#endif

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_TakeStepMRISR:

  This routine performs a single MRISR step.

  The vector ark_mem->yn holds the previous time-step solution
  on input, and the vector ark_mem->ycur should hold the result
  of this step on output.

  If timestep adaptivity is enabled, this routine also computes
  the error estimate y-ytilde, where ytilde is the
  embedded solution, and the norm weights come from ark_ewt.
  This esimate is stored in ark_mem->tempv1, in case the calling
  routine wishes to examine the error locations.

  The output variable dsmPtr should contain a scalar-valued
  estimate of the temporal error from this step, ||y-ytilde||_WRMS
  if timestep adaptivity is enabled; otherwise it should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure
  Since the fast-scale evolution could be considered a different
  type of "algebraic solver", we similarly report any fast-scale
  evolution error as a recoverable nflagPtr value.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int mriStep_TakeStepMRISR(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeMRIStepMem step_mem;          /* outer stepper memory       */
  int stage, j;                       /* stage indices              */
  int retval;                         /* reusable return flag       */
  N_Vector ytilde;                    /* embedded solution          */
  N_Vector ytemp;                     /* temporary vector           */
  SUNAdaptController_Type adapt_type; /* timestep adaptivity type   */
  sunbooleantype embedding;           /* flag indicating embedding  */
  sunbooleantype solution;            /*   or solution stages       */
  sunbooleantype impl_corr;           /* is slow correct. implicit? */
  sunbooleantype store_imprhs;        /* temporary storage          */
  sunbooleantype store_exprhs;
  sunrealtype cstage; /* current stage abscissa     */
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* access the MRIStep mem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* initialize algebraic solver convergence flag to success;
     error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* set N_Vector shortcuts */
  ytilde = ark_mem->tempv4;
  ytemp  = ark_mem->tempv2;

  /* if MRI adaptivity is enabled: reset fast accumulated error,
     and send appropriate control parameter to the fast integrator */
  adapt_type = SUNAdaptController_GetType(ark_mem->hadapt_mem->hcontroller);
  if ((adapt_type == SUN_ADAPTCONTROLLER_MRI_H) ||
      (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL))
  {
    step_mem->inner_dsm = ZERO;
    retval              = mriStepInnerStepper_ResetError(step_mem->stepper);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper error estimate");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H)
  {
    retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                              step_mem->inner_control);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set a fixed stepsize in the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL)
  {
    retval =
      mriStepInnerStepper_SetRTol(step_mem->stepper,
                                  step_mem->inner_control * ark_mem->reltol);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set the inner stepper tolerance");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* if any temporal adaptivity is enabled: reset the inner integrator
     to the beginning of this step (in case of recomputation) */
  if (!ark_mem->fixedstep)
  {
    retval = mriStepInnerStepper_Reset(step_mem->stepper, ark_mem->tn,
                                       ark_mem->yn);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* Evaluate the slow RHS functions if needed. NOTE: We do not use the full RHS
     function here (unlike ERKStep and ARKStep) since it does not need to check
     for FSAL or SA methods and thus avoids potentially unnecessary evaluations
     of the inner (fast) RHS function */
  if (!(ark_mem->fn_is_current))
  {
    /* compute the explicit component */
    if (step_mem->explicit_rhs)
    {
      retval = step_mem->fse(ark_mem->tn, ark_mem->yn, step_mem->Fse[0],
                             ark_mem->user_data);
      step_mem->nfse++;
      if (retval) { return ARK_RHSFUNC_FAIL; }
    }

    /* compute the implicit component */
    if (step_mem->implicit_rhs)
    {
      retval = step_mem->fsi(ark_mem->tn, ark_mem->yn, step_mem->Fsi[0],
                             ark_mem->user_data);
      step_mem->nfsi++;
      if (retval) { return ARK_RHSFUNC_FAIL; }
    }

    /* combine both RHS into Fse for ImEx problems */
    if (step_mem->implicit_rhs && step_mem->explicit_rhs)
    {
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, step_mem->Fsi[0],
                   step_mem->Fse[0]);
    }

    ark_mem->fn_is_current = SUNTRUE;
  }

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep step %li,  stage 0,  h = %" RSYM ",  t_n = %" RSYM "\n",
         ark_mem->nst, ark_mem->h, ark_mem->tn);
#endif

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRISR", "slow stage", "z[0] =", "");
  N_VPrintFile(ark_mem->yn, ARK_LOGGER->debug_fp);

  if (step_mem->explicit_rhs)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRISR", "slow explicit RHS",
                       "Fse[0] =", "");
    N_VPrintFile(step_mem->Fse[0], ARK_LOGGER->debug_fp);
  }

  if (step_mem->implicit_rhs)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRISR", "slow implicit RHS",
                       "Fsi[0] =", "");
    N_VPrintFile(step_mem->Fsi[0], ARK_LOGGER->debug_fp);
  }
#endif

  /* The first stage is the previous time-step solution, so its RHS
     is the [already-computed] slow RHS from the start of the step */

  /* Loop over stages */
  for (stage = 1; stage <= step_mem->stages; stage++)
  {
    /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRISR", "start-stage",
                       "step = %li, stage = %i, h = %" RSYM, ark_mem->nst,
                       stage, ark_mem->h);
#endif

    /* Determine if this is an "embedding" or "solution" stage */
    solution  = (stage == step_mem->stages - 1);
    embedding = (stage == step_mem->stages);

    /* Skip the embedding if we're using fixed time-stepping and
       temporal error estimation is disabled */
    if (ark_mem->fixedstep && embedding && (ark_mem->AccumErrorType < 0))
    {
      break;
    }

    /* Set initial condition for this stage */
    N_VScale(ONE, ark_mem->yn, ark_mem->ycur);

    /* Set current stage abscissa */
    cstage = (embedding) ? ONE : step_mem->MRIC->c[stage];

    /* Compute forcing function for inner solver: temporarily
       force explicit_rhs and disable implicit_rhs flag since MRISR
       methods ignore the G coefficients in the forcing function. */
    store_imprhs           = step_mem->implicit_rhs;
    store_exprhs           = step_mem->explicit_rhs;
    step_mem->implicit_rhs = SUNFALSE;
    step_mem->explicit_rhs = SUNTRUE;
    retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, stage, ark_mem->tn,
                                         ark_mem->tn + cstage * ark_mem->h);
    if (retval != ARK_SUCCESS)
    {
      *nflagPtr = CONV_FAIL;
      break;
    }
    step_mem->implicit_rhs = store_imprhs;
    step_mem->explicit_rhs = store_exprhs;

    /* Evolve fast IVP for this stage:
         force reset due to "stage-restart" structure
         get inner dsm on all non-embedding stages */
    retval = mriStep_StageERKFast(ark_mem, step_mem, stage, ark_mem->tn,
                                  ark_mem->tn + cstage * ark_mem->h,
                                  ark_mem->ycur, ytemp, SUNTRUE, !embedding);
    if (retval != ARK_SUCCESS) { *nflagPtr = CONV_FAIL; }

    /* set current stage time for implicit correction, postprocessing
       and RHS calls */
    ark_mem->tcur = ark_mem->tn + cstage * ark_mem->h;

    /* perform MRISR slow/implicit correction */
    impl_corr = SUNFALSE;
    if (step_mem->implicit_rhs)
    {
      /* determine whether implicit RHS correction will require an implicit solve */
      impl_corr = SUNRabs(step_mem->MRIC->G[0][stage][stage]) > tol;

      /* perform implicit solve for correction */
      if (impl_corr)
      {
        /* store current stage index (for an "embedded" stage, subtract 1) */
        step_mem->istage = (stage == step_mem->stages) ? stage - 1 : stage;

        /* Call predictor for current stage solution (result placed in zpred) */
        retval = mriStep_Predict(ark_mem, step_mem->istage, step_mem->zpred);
        if (retval != ARK_SUCCESS) { return (retval); }

        /* If a user-supplied predictor routine is provided, call that here
           Note that mriStep_Predict is *still* called, so this user-supplied
           routine can just "clean up" the built-in prediction, if desired. */
        if (step_mem->stage_predict)
        {
          retval = step_mem->stage_predict(ark_mem->tcur, step_mem->zpred,
                                           ark_mem->user_data);
          if (retval < 0) { return (ARK_USER_PREDICT_FAIL); }
          if (retval > 0) { return (TRY_AGAIN); }
        }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRISR", "predictor",
                           "zpred =", "");
        N_VPrintFile(step_mem->zpred, ARK_LOGGER->debug_fp);
#endif

        /* fill sdata with explicit contributions to correction */
        step_mem->cvals[0] = ONE;
        step_mem->Xvecs[0] = ark_mem->ycur;
        step_mem->cvals[1] = -ONE;
        step_mem->Xvecs[1] = step_mem->zpred;
        for (j = 0; j < stage; j++)
        {
          step_mem->cvals[j + 2] = ark_mem->h * step_mem->MRIC->G[0][stage][j];
          step_mem->Xvecs[j + 2] = step_mem->Fsi[j];
        }
        retval = N_VLinearCombination(stage + 2, step_mem->cvals,
                                      step_mem->Xvecs, step_mem->sdata);
        if (retval != 0) { return (ARK_VECTOROP_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRISR", "rhs data",
                           "sdata =", "");
        N_VPrintFile(step_mem->sdata, ARK_LOGGER->debug_fp);
#endif

        /* Update gamma for implicit solver */
        step_mem->gamma = ark_mem->h * step_mem->MRIC->G[0][stage][stage];
        if (ark_mem->firststage) { step_mem->gammap = step_mem->gamma; }
        step_mem->gamrat =
          (ark_mem->firststage) ? ONE : step_mem->gamma / step_mem->gammap;

        /* perform implicit solve (result is stored in ark_mem->ycur); return
           with positive value on anything but success */
        *nflagPtr = mriStep_Nls(ark_mem, *nflagPtr);
        if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }
      }
      /* perform explicit update for correction */
      else
      {
        step_mem->cvals[0] = ONE;
        step_mem->Xvecs[0] = ark_mem->ycur;
        for (j = 0; j < stage; j++)
        {
          step_mem->cvals[j + 1] = ark_mem->h * step_mem->MRIC->G[0][stage][j];
          step_mem->Xvecs[j + 1] = step_mem->Fsi[j];
        }
        retval = N_VLinearCombination(stage + 1, step_mem->cvals,
                                      step_mem->Xvecs, ark_mem->ycur);
        if (retval != 0) { return (ARK_VECTOROP_ERR); }
      }
    }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMRISR", "slow stage",
                       "z[%i] =", stage);
    N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

    /* apply user-supplied stage postprocessing function (if supplied),
       and reset the inner integrator with the modified stage solution */
    if (ark_mem->ProcessStage != NULL)
    {
      retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) { return (ARK_POSTPROCESS_STAGE_FAIL); }
      retval = mriStepInnerStepper_Reset(step_mem->stepper, ark_mem->tcur,
                                         ark_mem->ycur);
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                        __FILE__, "Unable to reset the inner stepper");
        return (ARK_INNERSTEP_FAIL);
      }
    }

    /* Compute updated slow RHS (except for final solution or embedding) */
    if ((!solution) && (!embedding))
    {
      /* store explicit slow rhs */
      if (step_mem->explicit_rhs)
      {
        retval = step_mem->fse(ark_mem->tcur, ark_mem->ycur,
                               step_mem->Fse[stage], ark_mem->user_data);
        step_mem->nfse++;
        if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
        if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRISR", "slow explicit RHS",
                           "Fse[%i] =", stage);
        N_VPrintFile(step_mem->Fse[stage], ARK_LOGGER->debug_fp);
#endif
      }

      /* store implicit slow rhs */
      if (step_mem->implicit_rhs)
      {
        if (!step_mem->deduce_rhs || !impl_corr)
        {
          retval = step_mem->fsi(ark_mem->tcur, ark_mem->ycur,
                                 step_mem->Fsi[stage], ark_mem->user_data);
          step_mem->nfsi++;
        }
        else
        {
          N_VLinearSum(ONE / step_mem->gamma, step_mem->zcor,
                       -ONE / step_mem->gamma, step_mem->sdata,
                       step_mem->Fsi[stage]);
        }

        if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
        if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMRISR", "slow implicit RHS",
                           "Fsi[%i] =", stage);
        N_VPrintFile(step_mem->Fsi[stage], ARK_LOGGER->debug_fp);
#endif
      }

      /* combine both RHS into Fse for ImEx problems */
      if (step_mem->implicit_rhs && step_mem->explicit_rhs)
      {
        N_VLinearSum(ONE, step_mem->Fse[stage], ONE, step_mem->Fsi[stage],
                     step_mem->Fse[stage]);
      }
    }

    /* If this is the solution stage, archive for error estimation */
    if (solution) { N_VScale(ONE, ark_mem->ycur, ytilde); }

  } /* loop over stages */

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRISR", "updated solution",
                     "ycur =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

  /* if temporal error estimation is enabled: compute estimate via difference between
     step solution and embedding, store in ark_mem->tempv1, store norm in dsmPtr, and
     copy solution back to ycur */
  if (!ark_mem->fixedstep || (ark_mem->AccumErrorType >= 0))
  {
    N_VLinearSum(ONE, ytilde, -ONE, ark_mem->ycur, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
    N_VScale(ONE, ytilde, ark_mem->ycur);
  }

  /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMRISR", "end-step",
                     "step = %li, t = %" RSYM ", h = %" RSYM ", dsm = %" RSYM,
                     ark_mem->nst, ark_mem->tn, ark_mem->h, *dsmPtr);
#endif

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_TakeStepMERK:

  This routine performs a single MERK step.

  The vector ark_mem->yn holds the previous time-step solution
  on input, and the vector ark_mem->ycur should hold the result
  of this step on output.

  If timestep adaptivity is enabled, this routine also computes
  the error estimate y-ytilde, where ytilde is the
  embedded solution, and the norm weights come from ark_ewt.
  This esimate is stored in ark_mem->tempv1, in case the calling
  routine wishes to examine the error locations.

  The output variable dsmPtr should contain a scalar-valued
  estimate of the temporal error from this step, ||y-ytilde||_WRMS
  if timestep adaptivity is enabled; otherwise it should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure
  Since the fast-scale evolution could be considered a different
  type of "algebraic solver", we similarly report any fast-scale
  evolution error as a recoverable nflagPtr value.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int mriStep_TakeStepMERK(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeMRIStepMem step_mem;          /* outer stepper memory       */
  int ig;                             /* current stage group index  */
  int is;                             /* stage index in group       */
  int stage, nextstage;               /* current/next stages        */
  int retval;                         /* reusable return flag       */
  N_Vector ytilde;                    /* embedded solution          */
  N_Vector ytemp;                     /* temporary vector           */
  SUNAdaptController_Type adapt_type; /* timestep adaptivity type   */
  sunrealtype t0, tf;                 /* start/end of each stage    */
  sunbooleantype embedding;           /* flag indicating embedding  */
  sunbooleantype solution;            /*   or solution stages       */
  sunrealtype cstage;                 /* current stage abscissa     */

  /* access the MRIStep mem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* initialize algebraic solver convergence flag to success;
     error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* set N_Vector shortcuts */
  ytilde = ark_mem->tempv4;
  ytemp  = ark_mem->tempv2;

  /* initial time for step */
  t0 = ark_mem->tn;

  /* if MRI adaptivity is enabled: reset fast accumulated error,
     and send appropriate control parameter to the fast integrator */
  adapt_type = SUNAdaptController_GetType(ark_mem->hadapt_mem->hcontroller);
  if ((adapt_type == SUN_ADAPTCONTROLLER_MRI_H) ||
      (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL))
  {
    step_mem->inner_dsm = ZERO;
    retval              = mriStepInnerStepper_ResetError(step_mem->stepper);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper error estimate");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H)
  {
    retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                              step_mem->inner_control);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set a fixed stepsize in the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }
  if (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL)
  {
    retval =
      mriStepInnerStepper_SetRTol(step_mem->stepper,
                                  step_mem->inner_control * ark_mem->reltol);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to set the inner stepper tolerance");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* if any temporal adaptivity is enabled: reset the inner integrator
     to the beginning of this step (in case of recomputation) */
  if (!ark_mem->fixedstep)
  {
    retval = mriStepInnerStepper_Reset(step_mem->stepper, t0, ark_mem->yn);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* Evaluate the slow RHS if needed. NOTE: We do not use the full RHS here
     (unlike ERKStep and ARKStep) since it does not need to check for FSAL
     or SA methods and thus avoids potentially unnecessary evaluations of the
     inner (fast) RHS function */
  if (!(ark_mem->fn_is_current))
  {
    retval = step_mem->fse(t0, ark_mem->yn, step_mem->Fse[0], ark_mem->user_data);
    step_mem->nfse++;
    if (retval) { return ARK_RHSFUNC_FAIL; }
    ark_mem->fn_is_current = SUNTRUE;
  }

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep step %li,  stage 0,  h = %" RSYM ",  t_n = %" RSYM "\n",
         ark_mem->nst, ark_mem->h, t0);
#endif

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMERK", "slow stage", "z[0] =", "");
  N_VPrintFile(ark_mem->yn, ARK_LOGGER->debug_fp);

  if (step_mem->explicit_rhs)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMERK", "slow explicit RHS",
                       "Fse[0] =", "");
    N_VPrintFile(step_mem->Fse[0], ARK_LOGGER->debug_fp);
  }
#endif

  /* The first stage is the previous time-step solution, so its RHS
     is the [already-computed] slow RHS from the start of the step */

  /* Loop over stage groups */
  for (ig = 0; ig < step_mem->MRIC->ngroup; ig++)
  {
    /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_TakeStepMERK", "start-stage",
                       "step = %li, stage group = %i, h = %" RSYM, ark_mem->nst,
                       ig, ark_mem->h);
#endif

    /* Set up fast RHS for this stage group */
    retval = mriStep_ComputeInnerForcing(ark_mem, step_mem,
                                         step_mem->MRIC->group[ig][0],
                                         ark_mem->tn, ark_mem->tn + ark_mem->h);
    if (retval != ARK_SUCCESS) { return (retval); }

    /* Set initial condition for this stage group */
    N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
    t0 = ark_mem->tn;

    /* Evolve fast IVP over each subinterval in stage group */
    for (is = 0; is < step_mem->stages; is++)
    {
      /* Get stage index from group; skip to the next group if
         we've reached the end of this one */
      stage = step_mem->MRIC->group[ig][is];
      if (stage < 0) { break; }
      nextstage = -1;
      if (stage < step_mem->stages)
      {
        nextstage = step_mem->MRIC->group[ig][is + 1];
      }

      /* Determine if this is an "embedding" or "solution" stage */
      embedding = solution = SUNFALSE;
      if (ig == step_mem->MRIC->ngroup - 2)
      {
        if ((stage >= 0) && (nextstage < 0)) { embedding = SUNTRUE; }
      }
      if (ig == step_mem->MRIC->ngroup - 1)
      {
        if ((stage >= 0) && (nextstage < 0)) { solution = SUNTRUE; }
      }

      /* Skip the embedding if we're using fixed time-stepping and
         temporal error estimation is disabled */
      if (ark_mem->fixedstep && embedding && (ark_mem->AccumErrorType < 0))
      {
        break;
      }

      /* Set current stage abscissa */
      cstage = (stage >= step_mem->stages) ? ONE : step_mem->MRIC->c[stage];

      /* Set desired output time for subinterval */
      tf = ark_mem->tn + cstage * ark_mem->h;

      /* Evolve fast IVP for this stage:
           force reset on first stage in group
           get inner dsm on all non-embedding stages */
      retval = mriStep_StageERKFast(ark_mem, step_mem, stage, t0, tf,
                                    ark_mem->ycur, ytemp, is == 0, !embedding);
      if (retval != ARK_SUCCESS) { *nflagPtr = CONV_FAIL; }

      /* Update "initial time" for next stage in group */
      t0 = tf;

      /* set current stage time for postprocessing and RHS calls */
      ark_mem->tcur = tf;

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::mriStep_TakeStepMERK", "slow stage",
                         "z[%i] =", stage);
      N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

      /* apply user-supplied stage postprocessing function (if supplied),
         and reset the inner integrator with the modified stage solution */
      if (ark_mem->ProcessStage != NULL)
      {
        retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                       ark_mem->user_data);
        if (retval != 0) { return (ARK_POSTPROCESS_STAGE_FAIL); }
        retval = mriStepInnerStepper_Reset(step_mem->stepper, ark_mem->tcur,
                                           ark_mem->ycur);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to reset the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }
      }

      /* Compute updated slow RHS (except for final solution or embedding) */
      if ((!solution) && (!embedding))
      {
        /* store explicit slow rhs */
        retval = step_mem->fse(ark_mem->tcur, ark_mem->ycur,
                               step_mem->Fse[stage], ark_mem->user_data);
        step_mem->nfse++;
        if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
        if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::mriStep_TakeStepMERK", "slow explicit RHS",
                           "Fse[%i] =", is);
        N_VPrintFile(step_mem->Fse[stage], ARK_LOGGER->debug_fp);
#endif
      }

      /* If this is the embedding stage, archive solution for error estimation */
      if (embedding) { N_VScale(ONE, ark_mem->ycur, ytilde); }

    } /* loop over stages */

  } /* loop over stage groups */

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMERK", "updated solution",
                     "ycur =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

  /* if temporal error estimation is enabled: compute estimate via difference between
     step solution and embedding, store in ark_mem->tempv1, and store norm in dsmPtr */
  if (!ark_mem->fixedstep || (ark_mem->AccumErrorType >= 0))
  {
    N_VLinearSum(ONE, ytilde, -ONE, ark_mem->ycur, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }

  /* Solver diagnostics reporting */
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_TakeStepMERK", "end-step",
                     "step = %li, t = %" RSYM ", h = %" RSYM ", dsm = %" RSYM,
                     ark_mem->nst, ark_mem->tn, ark_mem->h, *dsmPtr);
#endif

  return (ARK_SUCCESS);
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_AccessARKODEStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int mriStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeMRIStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  /* access ARKodeMRIStepMem structure */
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeMRIStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_AccessStepMem:

  Shortcut routine to unpack step_mem structure from ark_mem.
  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int mriStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeMRIStepMem* step_mem)
{
  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype mriStep_CheckNVector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvscale == NULL) || (tmpl->ops->nvwrmsnorm == NULL))
  {
    return (SUNFALSE);
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  mriStep_SetCoupling

  This routine determines the MRI method to use, based on the
  desired accuracy and fixed/adaptive time stepping choice.
  ---------------------------------------------------------------*/
int mriStep_SetCoupling(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  sunindextype Cliw, Clrw;
  int q_actual;
  ARKODE_MRITableID table_id;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;
  q_actual = step_mem->q;

  /* if coupling has already been specified, just return */
  if (step_mem->MRIC != NULL) { return (ARK_SUCCESS); }

  if (q_actual < 1 || q_actual > 4)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "No MRI method at requested order, using q=3.");
    q_actual = 3;
  }

  /* select method based on order and type */

  /**** ImEx methods ****/
  if (step_mem->implicit_rhs && step_mem->explicit_rhs)
  {
    if (!ark_mem->fixedstep)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "No adaptive ImEx MRI methods are yet available.");
      return (ARK_INVALID_TABLE);
    }
    switch (q_actual)
    {
    case 1: table_id = MRISTEP_DEFAULT_IMEX_SD_1; break;
    case 2: table_id = MRISTEP_DEFAULT_IMEX_SD_2; break;
    case 3: table_id = MRISTEP_DEFAULT_IMEX_SD_3; break;
    case 4: table_id = MRISTEP_DEFAULT_IMEX_SD_4; break;
    }
    /**** implicit methods ****/
  }
  else if (step_mem->implicit_rhs)
  {
    switch (q_actual)
    {
    case 1: table_id = MRISTEP_DEFAULT_IMPL_SD_1; break;
    case 2: table_id = MRISTEP_DEFAULT_IMPL_SD_2; break;
    case 3: table_id = MRISTEP_DEFAULT_IMPL_SD_3; break;
    case 4: table_id = MRISTEP_DEFAULT_IMPL_SD_4; break;
    }

    /**** explicit methods ****/
  }
  else
  {
    switch (q_actual)
    {
    case 1: table_id = MRISTEP_DEFAULT_EXPL_1; break;
    case 2: table_id = MRISTEP_DEFAULT_EXPL_2; break;
    case 3: table_id = MRISTEP_DEFAULT_EXPL_3; break;
    case 4: table_id = MRISTEP_DEFAULT_EXPL_4; break;
    case 5: table_id = MRISTEP_DEFAULT_EXPL_5_AD; break;
    }
  }

  step_mem->MRIC = MRIStepCoupling_LoadTable(table_id);
  if (step_mem->MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "An error occurred in constructing coupling table.");
    return (ARK_INVALID_TABLE);
  }

  /* note coupling structure space requirements */
  MRIStepCoupling_Space(step_mem->MRIC, &Cliw, &Clrw);
  ark_mem->liw += Cliw;
  ark_mem->lrw += Clrw;

  /* set [redundant] stored values for stage numbers and
     method/embedding orders */
  step_mem->stages = step_mem->MRIC->stages;
  step_mem->q      = step_mem->MRIC->q;
  step_mem->p      = step_mem->MRIC->p;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_CheckCoupling

  This routine runs through the MRI coupling structure to ensure
  that it meets all necessary requirements, including:
    sorted abscissae, with c[0] = 0 and c[end] = 1
    lower-triangular (i.e., ERK or DIRK)
    all DIRK stages are solve-decoupled [temporarily]
    method order q > 0 (all)
    stages > 0 (all)

  Returns ARK_SUCCESS if it passes, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int mriStep_CheckCoupling(ARKodeMem ark_mem)
{
  int i, j, k;
  sunbooleantype okay;
  ARKodeMRIStepMem step_mem;
  sunrealtype Gabs, Wabs;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->MRIC->stages < 1)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "stages < 1!");
    return (ARK_INVALID_TABLE);
  }

  /* check that method order q > 0 */
  if (step_mem->MRIC->q < 1)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "method order < 1");
    return (ARK_INVALID_TABLE);
  }

  /* check that embedding order p > 0 (if adaptive) */
  if ((step_mem->MRIC->p < 1) && (!ark_mem->fixedstep))
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "embedding order < 1");
    return (ARK_INVALID_TABLE);
  }

  /* Check that coupling table has compatible type */
  if (step_mem->implicit_rhs && step_mem->explicit_rhs &&
      (step_mem->MRIC->type != MRISTEP_IMEX) &&
      (step_mem->MRIC->type != MRISTEP_MRISR))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid coupling table for an IMEX problem!");
    return (ARK_ILL_INPUT);
  }
  if (step_mem->explicit_rhs && (step_mem->MRIC->type != MRISTEP_EXPLICIT) &&
      (step_mem->MRIC->type != MRISTEP_IMEX) &&
      (step_mem->MRIC->type != MRISTEP_MERK) &&
      (step_mem->MRIC->type != MRISTEP_MRISR))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid coupling table for an explicit problem!");
    return (ARK_ILL_INPUT);
  }
  if (step_mem->implicit_rhs && (step_mem->MRIC->type != MRISTEP_IMPLICIT) &&
      (step_mem->MRIC->type != MRISTEP_IMEX) &&
      (step_mem->MRIC->type != MRISTEP_MRISR))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid coupling table for an implicit problem!");
    return (ARK_ILL_INPUT);
  }

  /* Check that the matrices are defined appropriately */
  if (step_mem->MRIC->type == MRISTEP_IMEX)
  {
    /* ImEx */
    if (!(step_mem->MRIC->W) || !(step_mem->MRIC->G))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Invalid coupling table for an IMEX problem!");
      return (ARK_ILL_INPUT);
    }
  }
  else if (step_mem->MRIC->type == MRISTEP_EXPLICIT)
  {
    /* Explicit */
    if (!(step_mem->MRIC->W) || step_mem->MRIC->G)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Invalid coupling table for an explicit problem!");
      return (ARK_ILL_INPUT);
    }
  }
  else if (step_mem->MRIC->type == MRISTEP_IMPLICIT)
  {
    /* Implicit */
    if (step_mem->MRIC->W || !(step_mem->MRIC->G))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Invalid coupling table for an implicit problem!");
      return (ARK_ILL_INPUT);
    }
  }

  /* Check that W tables are strictly lower triangular */
  if (step_mem->MRIC->W)
  {
    Wabs = SUN_RCONST(0.0);
    for (k = 0; k < step_mem->MRIC->nmat; k++)
    {
      for (i = 0; i < step_mem->MRIC->stages; i++)
      {
        for (j = i; j < step_mem->MRIC->stages; j++)
        {
          Wabs += SUNRabs(step_mem->MRIC->W[k][i][j]);
        }
      }
    }
    if (Wabs > tol)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Coupling can be up to ERK (at most)!");
      return (ARK_INVALID_TABLE);
    }
  }

  /* Check that G tables are lower triangular */
  if (step_mem->MRIC->G)
  {
    Gabs = SUN_RCONST(0.0);
    for (k = 0; k < step_mem->MRIC->nmat; k++)
    {
      for (i = 0; i < step_mem->MRIC->stages; i++)
      {
        for (j = i + 1; j < step_mem->MRIC->stages; j++)
        {
          Gabs += SUNRabs(step_mem->MRIC->G[k][i][j]);
        }
      }
    }
    if (Gabs > tol)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Coupling can be up to DIRK (at most)!");
      return (ARK_INVALID_TABLE);
    }
  }

  /* Check that no stage has MRISTAGE_DIRK_FAST type (for now) */
  okay = SUNTRUE;
  for (i = 0; i < step_mem->MRIC->stages; i++)
  {
    if (mriStepCoupling_GetStageType(step_mem->MRIC, i) == MRISTAGE_DIRK_FAST)
    {
      okay = SUNFALSE;
    }
  }
  if (!okay)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "solve-coupled DIRK stages not currently supported");
    return (ARK_INVALID_TABLE);
  }

  /* check that MRI-GARK stage times are sorted */
  if ((step_mem->MRIC->type == MRISTEP_IMPLICIT) ||
      (step_mem->MRIC->type == MRISTEP_EXPLICIT) ||
      (step_mem->MRIC->type == MRISTEP_IMEX))
  {
    okay = SUNTRUE;
    for (i = 1; i < step_mem->MRIC->stages; i++)
    {
      if ((step_mem->MRIC->c[i] - step_mem->MRIC->c[i - 1]) < -tol)
      {
        okay = SUNFALSE;
      }
    }
    if (!okay)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Stage times must be sorted.");
      return (ARK_INVALID_TABLE);
    }
  }

  /* check that the first stage is just the old step solution */
  Gabs = SUNRabs(step_mem->MRIC->c[0]);
  for (k = 0; k < step_mem->MRIC->nmat; k++)
  {
    for (j = 0; j < step_mem->MRIC->stages; j++)
    {
      if (step_mem->MRIC->W) { Gabs += SUNRabs(step_mem->MRIC->W[k][0][j]); }
      if (step_mem->MRIC->G) { Gabs += SUNRabs(step_mem->MRIC->G[k][0][j]); }
    }
  }
  if (Gabs > tol)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "First stage must equal old solution.");
    return (ARK_INVALID_TABLE);
  }

  /* check that the last stage is at the final time */
  if (SUNRabs(ONE - step_mem->MRIC->c[step_mem->MRIC->stages - 1]) > tol)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "Final stage time must be equal 1.");
    return (ARK_INVALID_TABLE);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_StageERKFast

  This routine performs a single MRI stage, is, with explicit
  slow time scale and fast time scale that requires evolution.

  On input, ycur is the initial condition for the fast IVP at t0.
  On output, ycur is the solution of the fast IVP at tf.
  The vector ytemp is only used if temporal adaptivity is enabled,
  and the fast error is not provided by the fast integrator.

  force_reset indicates whether ycur differs from the result of
  the previous fast evolution, in which case the inner integrator
  needs to be reset.

  get_inner_dsm indicates whether this stage is one that should
  accumulate an inner temporal error estimate.
  ---------------------------------------------------------------*/
int mriStep_StageERKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, int is,
                         sunrealtype t0, sunrealtype tf, N_Vector ycur,
                         N_Vector ytemp, sunbooleantype force_reset,
                         sunbooleantype get_inner_dsm)
{
  int retval;                         /* reusable return flag */
  SUNAdaptController_Type adapt_type; /* timestep adaptivity type */

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep ERK fast stage\n");
#endif

  /* pre inner evolve function (if supplied) */
  if (step_mem->pre_inner_evolve)
  {
    retval = step_mem->pre_inner_evolve(t0, step_mem->stepper->forcing,
                                        step_mem->stepper->nforcing,
                                        ark_mem->user_data);
    if (retval != 0) { return (ARK_OUTERTOINNER_FAIL); }
  }

  /* if the input state has changed since last called, reset the inner integrator */
  if (force_reset)
  {
    retval = mriStepInnerStepper_Reset(step_mem->stepper, t0, ycur);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reset the inner stepper");
      return (ARK_INNERSTEP_FAIL);
    }
  }

  /* Get the adaptivity type (if applicable) */
  adapt_type = (get_inner_dsm)
                 ? SUNAdaptController_GetType(ark_mem->hadapt_mem->hcontroller)
                 : SUN_ADAPTCONTROLLER_NONE;

  /* if we'll need to manually accumulate the fast error estimate, archive the
     current state to provide an initial condition for the ``fast embedding'' */
  if ((get_inner_dsm) && (adapt_type == SUN_ADAPTCONTROLLER_MRI_H) &&
      (step_mem->inner_hfactor != ZERO))
  {
    N_VScale(ONE, ycur, ytemp);
  }

  /* advance inner method in time */
  retval = mriStepInnerStepper_Evolve(step_mem->stepper, t0, tf, ycur);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                    "Failure when evolving the inner stepper");
    return (ARK_INNERSTEP_FAIL);
  }

  /* for normal stages with MRI adaptivity enabled, get an
     estimate for the fast time scale error */
  if (get_inner_dsm)
  {
    /* if the fast integrator uses adaptive steps, retrieve the error estimate */
    if (adapt_type == SUN_ADAPTCONTROLLER_MRI_TOL)
    {
      retval = mriStepInnerStepper_GetError(step_mem->stepper,
                                            &(step_mem->inner_dsm));
      if (retval != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                        __FILE__, "Unable to get accumulated error from the inner stepper");
        return (ARK_INNERSTEP_FAIL);
      }
    }

    /* the fast integrator uses fixed steps */
    if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H)
    {
      /* if we trust inner integrator accumulated error estimate, call it here */
      if (step_mem->inner_hfactor == ZERO)
      {
        retval = mriStepInnerStepper_GetError(step_mem->stepper,
                                              &(step_mem->inner_dsm));
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to get accumulated error from the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }
      }
      else /* manually accumulate fast error estimate */
      {
        /* reset fast integrator for time interval */
        retval = mriStepInnerStepper_Reset(step_mem->stepper, t0, ytemp);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to reset the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }

        /* update fixed-step size for fast integrator */
        retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                                  step_mem->inner_control *
                                                    step_mem->inner_hfactor);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to set a fixed stepsize in the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }

        /* evolve fast integrator using modified step size */
        retval = mriStepInnerStepper_Evolve(step_mem->stepper, t0, tf, ytemp);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Failure when evolving the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }

        /* compute solution difference */
        N_VLinearSum(ONE, ytemp, -ONE, ycur, ytemp);

        /* accumulate fast error estimate */
        step_mem->inner_dsm =
          SUNMAX(step_mem->inner_dsm,
                 ark_mem->reltol * N_VWrmsNorm(ytemp, ark_mem->ewt));

        /* reset fast integrator to the main evolution result */
        retval = mriStepInnerStepper_Reset(step_mem->stepper, tf, ycur);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to reset the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }
        retval = mriStepInnerStepper_SetFixedStep(step_mem->stepper,
                                                  step_mem->inner_control);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__,
                          __FILE__, "Unable to set a fixed stepsize in the inner stepper");
          return (ARK_INNERSTEP_FAIL);
        }
      } /* if (step_mem->inner_hfactor == ZERO) */
    }   /* if (adapt_type == SUN_ADAPTCONTROLLER_MRI_H) */
  }     /* if (get_inner_dsm) */

  /* post inner evolve function (if supplied) */
  if (step_mem->post_inner_evolve)
  {
    retval = step_mem->post_inner_evolve(tf, ycur, ark_mem->user_data);
    if (retval != 0) { return (ARK_INNERTOOUTER_FAIL); }
  }

  /* return with success */
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_StageERKNoFast

  This routine performs a single MRI stage with explicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
int mriStep_StageERKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, int is)
{
  int retval, j, nvec;

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep ERK stage\n");
#endif

  /* determine effective ERK coefficients (store in Ae_row and Ai_row) */
  retval = mriStep_RKCoeffs(step_mem->MRIC, is, step_mem->stage_map,
                            step_mem->Ae_row, step_mem->Ai_row);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* call fused vector operation to perform ERK update -- bound on
     j needs "SUNMIN" to handle the case of an "embedding" stage */
  step_mem->cvals[0] = ONE;
  step_mem->Xvecs[0] = ark_mem->ycur;
  nvec               = 1;
  for (j = 0; j < SUNMIN(is, step_mem->stages); j++)
  {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1)
    {
      step_mem->cvals[nvec] = ark_mem->h *
                              step_mem->Ae_row[step_mem->stage_map[j]];
      step_mem->Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[j]];
      nvec += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1)
    {
      step_mem->cvals[nvec] = ark_mem->h *
                              step_mem->Ai_row[step_mem->stage_map[j]];
      step_mem->Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[j]];
      nvec += 1;
    }
  }
  /* Is there a case where we have an explicit update with Fsi? */

  retval = N_VLinearCombination(nvec, step_mem->cvals, step_mem->Xvecs,
                                ark_mem->ycur);
  if (retval != 0) { return (ARK_VECTOROP_ERR); }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_StageDIRKFast

  This routine performs a single stage of a "solve coupled"
  MRI method, i.e. a stage that is DIRK on the slow time scale
  and involves evolution of the fast time scale, in a
  fully-coupled fashion.
  ---------------------------------------------------------------*/
int mriStep_StageDIRKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, int is,
                          int* nflagPtr)
{
#ifdef SUNDIALS_DEBUG
  printf("    MRIStep DIRK fast stage\n");
#endif

  /* this is not currently implemented */
  arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                  "This routine is not yet implemented.");
  return (ARK_INVALID_TABLE);
}

/*---------------------------------------------------------------
  mriStep_StageDIRKNoFast

  This routine performs a single MRI stage with implicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
int mriStep_StageDIRKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                            int is, int* nflagPtr)
{
  int retval;

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep DIRK stage\n");
#endif

  /* store current stage index (for an "embedded" stage, subtract 1) */
  step_mem->istage = (is == step_mem->stages) ? is - 1 : is;

  /* Call predictor for current stage solution (result placed in zpred) */
  retval = mriStep_Predict(ark_mem, step_mem->istage, step_mem->zpred);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* If a user-supplied predictor routine is provided, call that here
     Note that mriStep_Predict is *still* called, so this user-supplied
     routine can just "clean up" the built-in prediction, if desired. */
  if (step_mem->stage_predict)
  {
    retval = step_mem->stage_predict(ark_mem->tcur, step_mem->zpred,
                                     ark_mem->user_data);
    if (retval < 0) { return (ARK_USER_PREDICT_FAIL); }
    if (retval > 0) { return (TRY_AGAIN); }
  }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_StageDIRKNoFast", "predictor",
                     "zpred =", "");
  N_VPrintFile(step_mem->zpred, ARK_LOGGER->debug_fp);
#endif

  /* determine effective DIRK coefficients (store in cvals) */
  retval = mriStep_RKCoeffs(step_mem->MRIC, is, step_mem->stage_map,
                            step_mem->Ae_row, step_mem->Ai_row);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set up data for evaluation of DIRK stage residual (data stored in sdata) */
  retval = mriStep_StageSetup(ark_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                     "ARKODE::mriStep_StageDIRKNoFast", "rhs data",
                     "sdata =", "");
  N_VPrintFile(step_mem->sdata, ARK_LOGGER->debug_fp);
#endif

  /* perform implicit solve (result is stored in ark_mem->ycur); return
     with positive value on anything but success */
  *nflagPtr = mriStep_Nls(ark_mem, *nflagPtr);
  if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_ComputeInnerForcing

  Constructs the 'coefficient' vectors for the forcing polynomial
  for a 'fast' outer MRI-GARK stage i:

  p_i(theta) = sum_{k=0}^{n-1} forcing[k] * theta^k

  where theta = (t - t0) / (tf-t0) is the mapped 'time' for
  each 'fast' MRIStep evolution, with:
  * t0 -- the start of this outer MRIStep stage
  * tf-t0, the temporal width of this MRIStep stage
  * n -- shorthand for MRIC->nmat

  Defining cdiff = (tf-t0)/h, explicit and solve-decoupled
  implicit or IMEX MRI-based methods define this forcing polynomial
  for each outer stage i > 0:

  p_i(theta) = w_i,0(theta) * fse_0 + ... + w_i,{i-1}(theta) * fse_{i-1}
             + g_i,0(theta) * fsi_0 + ... + g_i,{i-1}(theta) * fsi_{i-1}

  where

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = 1/cdiff * MRIC->W[k][i][j]

  and

  g_i,j(theta) = g_0,i,j + g_1,i,j * theta + ... + g_n,i,j * theta^{n-1},
  g_k,i,j = 1/cdiff * MRIC->G[k][i][j]

  Converting to the appropriate form, we have

  p_i(theta) = ( w_0,i,0 * fse_0 + ... + w_0,i,{i-1} * fse_{i-1} +
                 g_0,i,0 * fsi_0 + ... + g_0,i,{i-1} * fsi_{i-1} ) * theta^0
             + ( w_1,i,0 * fse_0 + ... + w_1,i,{i-1} * fse_{i-1} +
                 g_1,i,0 * fsi_0 + ... + g_1,i,{i-1} * fsi_{i-1} ) * theta^1
                                    .
                                    .
                                    .
             + ( w_n,i,0 * fse_0 + ... + w_n,i,{i-1} * fse_{i-1} +
                 g_n,i,0 * fsi_0 + ... + g_n,i,{i-1} * fsi_{i-1} ) * theta^{n-1}

  Thus we define the forcing vectors for k = 0,...,nmat - 1

  forcing[k] = w_k,i,0 * fse_0 + ... + w_k,i,{i-1} * fse_{i-1}
             + g_k,i,0 * fsi_0 + ... + g_k,i,{i-1} * fsi_{i-1}

             = 1 / cdiff *
               ( W[k][i][0] * fse_0 + ... + W[k][i][i-1] * fse_{i-1} +
               ( G[k][i][0] * fsi_0 + ... + G[k][i][i-1] * fsi_{i-1} )

  We may use an identical formula for MERK methods, so long as we set t0=tn,
  tf=tn+h, stage_map[j]=j (identity map), and implicit_rhs=SUNFALSE.
  With this configuration: tf-t0=h, theta = (t-tn)/h, and cdiff=1.  MERK methods
  define the forcing polynomial for each outer stage i > 0 as:

  p_i(theta) = w_i,0(theta) * fse_0 + ... + w_i,{i-1}(theta) * fse_{i-1}

  where

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = MRIC->W[k][i][j]

  which is equivalent to the formula above.

  We may use a similar formula for MRISR methods, so long as we set t0=tn,
  tf=tn+h*ci, stage_map[j]=j (identity map), and implicit_rhs=SUNFALSE.
  With this configuration: tf-t0=ci*h, theta = (t-tn)/(ci*h), and cdiff=1/ci.
  MRISR methods define the forcing polynomial for each outer stage i > 0 as:

  p_i(theta) = w_i,0(theta) * fs_0 + ... + w_i,{i-1}(theta) * fs_{i-1}

  where fs_j = fse_j + fsi_j and

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = 1/ci * MRIC->W[k][i][j]

  which is equivalent to the formula above, so long as the stage RHS vectors
  Fse[j] are repurposed to instead store (fse_j + fsi_j).

  This routine additionally returns a success/failure flag:
     ARK_SUCCESS -- successful evaluation
  ---------------------------------------------------------------*/

int mriStep_ComputeInnerForcing(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                                int stage, sunrealtype t0, sunrealtype tf)
{
  sunrealtype rcdiff;
  int j, k, nmat, nstore, retval;
  sunrealtype* cvals;
  N_Vector* Xvecs;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Set inner forcing time normalization constants */
  step_mem->stepper->tshift = t0;
  step_mem->stepper->tscale = tf - t0;

  /* compute inner forcing vectors (assumes cdiff != 0) */
  nstore = 0;
  for (j = 0; j < SUNMIN(stage, step_mem->stages); j++)
  {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1)
    {
      Xvecs[nstore] = step_mem->Fse[step_mem->stage_map[j]];
      nstore += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1)
    {
      Xvecs[nstore] = step_mem->Fsi[step_mem->stage_map[j]];
      nstore += 1;
    }
  }

  nmat   = step_mem->MRIC->nmat;
  rcdiff = ark_mem->h / (tf - t0);

  for (k = 0; k < nmat; k++)
  {
    nstore = 0;
    for (j = 0; j < SUNMIN(stage, step_mem->stages); j++)
    {
      if (step_mem->stage_map[j] > -1)
      {
        if (step_mem->explicit_rhs && step_mem->implicit_rhs)
        {
          /* ImEx */
          cvals[nstore] = rcdiff * step_mem->MRIC->W[k][stage][j];
          nstore += 1;
          cvals[nstore] = rcdiff * step_mem->MRIC->G[k][stage][j];
          nstore += 1;
        }
        else if (step_mem->explicit_rhs)
        {
          /* explicit only */
          cvals[nstore] = rcdiff * step_mem->MRIC->W[k][stage][j];
          nstore += 1;
        }
        else
        {
          /* implicit only */
          cvals[nstore] = rcdiff * step_mem->MRIC->G[k][stage][j];
          nstore += 1;
        }
      }
    }

    retval = N_VLinearCombination(nstore, cvals, Xvecs,
                                  step_mem->stepper->forcing[k]);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }
  }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  for (k = 0; k < nmat; k++)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::mriStep_ComputeInnerForcing", "forcing",
                       "forcing[%i] =", k);
    N_VPrintFile(step_mem->stepper->forcing[k], ARK_LOGGER->debug_fp);
  }
#endif

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  Compute/return the effective RK coefficients for a "nofast"
  stage.  We may assume that "A" has already been allocated.
  ---------------------------------------------------------------*/

int mriStep_RKCoeffs(MRIStepCoupling MRIC, int is, int* stage_map,
                     sunrealtype* Ae_row, sunrealtype* Ai_row)
{
  int j, k;
  sunrealtype kconst;

  if (is < 1 || is > MRIC->stages || !stage_map || !Ae_row || !Ai_row)
  {
    return ARK_INVALID_TABLE;
  }

  /* initialize RK coefficient array */
  for (j = 0; j < MRIC->stages; j++)
  {
    Ae_row[j] = ZERO;
    Ai_row[j] = ZERO;
  }

  /* compute RK coefficients -- note that bounds on j need
     "SUNMIN" to handle the case of an "embedding" stage */
  for (k = 0; k < MRIC->nmat; k++)
  {
    kconst = ONE / (k + ONE);
    if (MRIC->W)
    {
      for (j = 0; j < SUNMIN(is, MRIC->stages - 1); j++)
      {
        if (stage_map[j] > -1)
        {
          Ae_row[stage_map[j]] += (MRIC->W[k][is][j] * kconst);
        }
      }
    }
    if (MRIC->G)
    {
      for (j = 0; j <= SUNMIN(is, MRIC->stages - 1); j++)
      {
        if (stage_map[j] > -1)
        {
          Ai_row[stage_map[j]] += (MRIC->G[k][is][j] * kconst);
        }
      }
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
int mriStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess)
{
  int i, retval, jstage, nvec;
  sunrealtype tau;
  sunrealtype h;
  ARKodeMRIStepMem step_mem;
  sunrealtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Interpolation structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->initsetup)
  {
    N_VScale(ONE, ark_mem->yn, yguess);
    return (ARK_SUCCESS);
  }

  /* set evaluation time tau as relative shift from previous successful time */
  tau = step_mem->MRIC->c[istage] * ark_mem->h / ark_mem->hold;

  /* use requested predictor formula */
  switch (step_mem->predictor)
  {
  case 1:

    /***** Interpolatory Predictor 1 -- all to max order *****/
    retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT) { return (retval); }
    break;

  case 2:

    /***** Interpolatory Predictor 2 -- decrease order w/ increasing level of extrapolation *****/
    retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT) { return (retval); }
    break;

  case 3:

    /***** Cutoff predictor: max order interpolatory output for stages "close"
           to previous step, first-order predictor for subsequent stages *****/
    retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT) { return (retval); }
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i,
           construct a quadratic Hermite interpolant for prediction; otherwise
           use the trivial predictor.  The actual calculations are performed in
           arkPredict_Bootstrap, but here we need to determine the appropriate
           stage, c_j, to use. *****/

    /* determine if any previous stages in step meet criteria */
    jstage = -1;
    for (i = 0; i < istage; i++)
    {
      jstage = (step_mem->MRIC->c[i] != ZERO) ? i : jstage;
    }

    /* if using the trivial predictor, break */
    if (jstage == -1) { break; }

    /* find the "optimal" previous stage to use */
    for (i = 0; i < istage; i++)
    {
      if ((step_mem->MRIC->c[i] > step_mem->MRIC->c[jstage]) &&
          (step_mem->MRIC->c[i] != ZERO) && step_mem->stage_map[i] > -1)
      {
        jstage = i;
      }
    }

    /* set stage time, stage RHS and interpolation values */
    h    = ark_mem->h * step_mem->MRIC->c[jstage];
    tau  = ark_mem->h * step_mem->MRIC->c[istage];
    nvec = 0;
    if (step_mem->implicit_rhs)
    { /* Implicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[jstage]];
      nvec += 1;
    }
    if (step_mem->explicit_rhs)
    { /* Explicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[jstage]];
      nvec += 1;
    }

    /* call predictor routine */
    retval = arkPredict_Bootstrap(ark_mem, h, tau, nvec, cvals, Xvecs, yguess);
    if (retval != ARK_ILL_INPUT) { return (retval); }
    break;
  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->yn, yguess);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_StageSetup

  This routine sets up the stage data for computing the
  solve-decoupled MRI stage residual, along with the step- and
  method-related factors gamma, gammap and gamrat.

  At the ith stage, we compute the residual vector for
  z=z_i=zp+zc:
    r = z - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zp + zc) - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zc - gamma*F(z)) - data,
  where data = (z_{i-1} - zp + h*sum_{j=0}^{i-1} A(i,j)*F(z_j))
  corresponds to existing information.  This routine computes
  this 'data' vector and stores in step_mem->sdata.

  Note: on input, this row A(i,:) is already stored in rkcoeffs.
  ---------------------------------------------------------------*/
int mriStep_StageSetup(ARKodeMem ark_mem)
{
  /* local data */
  ARKodeMRIStepMem step_mem;
  int retval, i, j, nvec;
  sunrealtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = step_mem->istage;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Update gamma (if the method contains an implicit component) */
  step_mem->gamma = ark_mem->h * step_mem->Ai_row[step_mem->stage_map[i]];

  if (ark_mem->firststage) { step_mem->gammap = step_mem->gamma; }
  step_mem->gamrat = (ark_mem->firststage) ? ONE
                                           : step_mem->gamma / step_mem->gammap;

  /* set cvals and Xvecs for setting stage data */
  cvals[0] = ONE;
  Xvecs[0] = ark_mem->ycur;
  cvals[1] = -ONE;
  Xvecs[1] = step_mem->zpred;
  nvec     = 2;

  for (j = 0; j < i; j++)
  {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1)
    {
      cvals[nvec] = ark_mem->h * step_mem->Ae_row[step_mem->stage_map[j]];
      Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[j]];
      nvec += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1)
    {
      cvals[nvec] = ark_mem->h * step_mem->Ai_row[step_mem->stage_map[j]];
      Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[j]];
      nvec += 1;
    }
  }

  /* call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
  if (retval != 0) { return (ARK_VECTOROP_ERR); }

  /* return with success */
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SlowRHS:

  Wrapper routine to call the user-supplied slow RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y), with API matching
  ARKTimestepFullRHSFn.  This is only used to determine an
  initial slow time-step size to use when one is not specified
  by the user (i.e., mode should correspond with
  ARK_FULLRHS_OTHER.
  ---------------------------------------------------------------*/
int mriStep_SlowRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* call fse if the problem has an explicit component */
  if (step_mem->explicit_rhs)
  {
    retval = step_mem->fse(t, y, step_mem->Fse[0], ark_mem->user_data);
    step_mem->nfse++;
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }
  }

  /* call fsi if the problem has an implicit component */
  if (step_mem->implicit_rhs)
  {
    retval = step_mem->fsi(t, y, step_mem->Fsi[0], ark_mem->user_data);
    step_mem->nfsi++;
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }
  }

  /* combine RHS vectors into output */
  if (step_mem->explicit_rhs && step_mem->implicit_rhs) /* ImEx */
  {
    N_VLinearSum(ONE, step_mem->Fse[0], ONE, step_mem->Fsi[0], f);
  }
  else
  {
    if (step_mem->implicit_rhs) { N_VScale(ONE, step_mem->Fsi[0], f); }
    else { N_VScale(ONE, step_mem->Fse[0], f); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_FastRHS:

  This is just a wrapper to call the fast RHS function,
  f(t,y) = ff(t,y), with API matching ARKTimestepFullRHSFn.  This
  is only used to determine an initial fast time-step size to use
  when one is not specified by the user and H-H MRI time step
  adaptivity is enabled.
  ---------------------------------------------------------------*/
int mriStep_FastRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* call ff */
  retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                       ARK_FULLRHS_OTHER);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return (ARK_RHSFUNC_FAIL);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_Hin

  This routine computes a tentative initial step size h0.  This
  employs the same safeguards as ARKODE's arkHin utility routine,
  but employs a simpler algorithm that estimates the first step
  such that an explicit Euler step (for only the slow RHS
  routine(s)) would be within user-specified tolerances of the
  initial condition.
  ---------------------------------------------------------------*/
int mriStep_Hin(ARKodeMem ark_mem, sunrealtype tcur, sunrealtype tout,
                N_Vector ycur, N_Vector fcur, N_Vector ytmp, N_Vector temp1,
                N_Vector temp2, ARKTimestepFullRHSFn rhs, sunrealtype* h)
{
  int sign;
  sunrealtype tdiff, tdist, tround, fnorm, h0_inv;

  /* If tout is too close to tn, give up */
  if ((tdiff = tout - tcur) == ZERO) { return (ARK_TOO_CLOSE); }
  sign   = (tdiff > ZERO) ? 1 : -1;
  tdist  = SUNRabs(tdiff);
  tround = ark_mem->uround * SUNMAX(SUNRabs(tcur), SUNRabs(tout));
  if (tdist < TWO * tround) { return (ARK_TOO_CLOSE); }

  /* h0 should bound the change due to a forward Euler step, and
     include safeguard against "too-small" ||f(t0,y0)||: */
  fnorm  = N_VWrmsNorm(fcur, ark_mem->ewt) / H0_BIAS;
  h0_inv = SUNMAX(ONE / H0_UBFACTOR / tdist, fnorm);
  *h     = sign / h0_inv;
  return (ARK_SUCCESS);
}

/*===============================================================
  User-callable functions for a custom inner integrator
  ===============================================================*/

int MRIStepInnerStepper_Create(SUNContext sunctx, MRIStepInnerStepper* stepper)
{
  if (!sunctx) { return ARK_ILL_INPUT; }

  *stepper = NULL;
  *stepper = (MRIStepInnerStepper)malloc(sizeof(**stepper));
  if (*stepper == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    return (ARK_MEM_FAIL);
  }
  memset(*stepper, 0, sizeof(**stepper));

  (*stepper)->ops = (MRIStepInnerStepper_Ops)malloc(sizeof(*((*stepper)->ops)));
  if ((*stepper)->ops == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    free(*stepper);
    return (ARK_MEM_FAIL);
  }
  memset((*stepper)->ops, 0, sizeof(*((*stepper)->ops)));

  /* initialize stepper data */
  (*stepper)->last_flag = ARK_SUCCESS;
  (*stepper)->sunctx    = sunctx;

  return (ARK_SUCCESS);
}

int MRIStepInnerStepper_Free(MRIStepInnerStepper* stepper)
{
  if (*stepper == NULL) { return ARK_SUCCESS; }

  /* free the inner forcing and fused op workspace vector */
  mriStepInnerStepper_FreeVecs(*stepper);

  /* free operations structure */
  free((*stepper)->ops);

  /* free inner stepper mem */
  free(*stepper);
  *stepper = NULL;

  return (ARK_SUCCESS);
}

int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper, void* content)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }
  stepper->content = content;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper, void** content)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }
  *content = stepper->content;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper,
                                    MRIStepInnerEvolveFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->evolve = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper,
                                     MRIStepInnerFullRhsFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->fullrhs = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper,
                                   MRIStepInnerResetFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->reset = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetAccumulatedErrorGetFn(MRIStepInnerStepper stepper,
                                                 MRIStepInnerGetAccumulatedError fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->geterror = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetAccumulatedErrorResetFn(
  MRIStepInnerStepper stepper, MRIStepInnerResetAccumulatedError fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->reseterror = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetFixedStepFn(MRIStepInnerStepper stepper,
                                       MRIStepInnerSetFixedStep fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->setfixedstep = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_SetRTolFn(MRIStepInnerStepper stepper,
                                  MRIStepInnerSetRTol fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->setrtol = fn;

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper, sunrealtype t,
                                   N_Vector f)
{
  sunrealtype tau, taui;
  int i;

  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  /* always append the constant forcing term */
  stepper->vals[0] = ONE;
  stepper->vecs[0] = f;

  /* compute normalized time tau and initialize tau^i */
  tau  = (t - stepper->tshift) / (stepper->tscale);
  taui = ONE;

  for (i = 0; i < stepper->nforcing; i++)
  {
    stepper->vals[i + 1] = taui;
    stepper->vecs[i + 1] = stepper->forcing[i];
    taui *= tau;
  }

  N_VLinearCombination(stepper->nforcing + 1, stepper->vals, stepper->vecs, f);

  return ARK_SUCCESS;
}

int MRIStepInnerStepper_GetForcingData(MRIStepInnerStepper stepper,
                                       sunrealtype* tshift, sunrealtype* tscale,
                                       N_Vector** forcing, int* nforcing)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  *tshift   = stepper->tshift;
  *tscale   = stepper->tscale;
  *forcing  = stepper->forcing;
  *nforcing = stepper->nforcing;

  return ARK_SUCCESS;
}

/*===============================================================
  Private inner integrator functions
  ===============================================================*/

/* Check for required operations */
int mriStepInnerStepper_HasRequiredOps(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

  if (stepper->ops->evolve) { return ARK_SUCCESS; }
  else { return ARK_ILL_INPUT; }
}

/* Check whether stepper supports fast/slow stepsize adaptivity */
sunbooleantype mriStepInnerStepper_SupportsStepAdaptivity(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) { return SUNFALSE; }
  if (stepper->ops == NULL) { return SUNFALSE; }

  if (stepper->ops->geterror && stepper->ops->reseterror &&
      stepper->ops->setfixedstep)
  {
    return SUNTRUE;
  }
  else { return SUNFALSE; }
}

/* Check whether stepper supports fast/slow tolerance adaptivity */
sunbooleantype mriStepInnerStepper_SupportsRTolAdaptivity(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) { return SUNFALSE; }
  if (stepper->ops == NULL) { return SUNFALSE; }

  if (stepper->ops->geterror && stepper->ops->reseterror && stepper->ops->setrtol)
  {
    return SUNTRUE;
  }
  else { return SUNFALSE; }
}

/* Evolve the inner (fast) ODE */
int mriStepInnerStepper_Evolve(MRIStepInnerStepper stepper, sunrealtype t0,
                               sunrealtype tout, N_Vector y)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops->evolve == NULL) { return ARK_ILL_INPUT; }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(stepper->sunctx->logger, SUN_LOGLEVEL_INFO,
                     "ARKODE::mriStepInnerStepper_Evolve", "start-inner-evolve",
                     "t0 = %" RSYM ", tout = %" RSYM, t0, tout);
#endif

  stepper->last_flag = stepper->ops->evolve(stepper, t0, tout, y);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(stepper->sunctx->logger, SUN_LOGLEVEL_INFO,
                     "ARKODE::mriStepInnerStepper_Evolve", "end-inner-evolve",
                     "flag = %i", stepper->last_flag);
#endif

  return stepper->last_flag;
}

/* Compute the full RHS for inner (fast) time scale TODO(DJG): This function can
   be made optional when fullrhs is not called unconditionally by the ARKODE
   infrastructure e.g., in arkInitialSetup, arkYddNorm, and arkCompleteStep. */
int mriStepInnerStepper_FullRhs(MRIStepInnerStepper stepper, sunrealtype t,
                                N_Vector y, N_Vector f, int mode)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops->fullrhs == NULL) { return ARK_ILL_INPUT; }

  stepper->last_flag = stepper->ops->fullrhs(stepper, t, y, f, mode);
  return stepper->last_flag;
}

/* Reset the inner (fast) stepper state */
int mriStepInnerStepper_Reset(MRIStepInnerStepper stepper, sunrealtype tR,
                              N_Vector yR)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(stepper->sunctx->logger, SUN_LOGLEVEL_INFO,
                     "ARKODE::mriStepInnerStepper_Reset", "reset-inner-state",
                     "tR = %" RSYM, tR);
#endif

  if (stepper->ops->reset)
  {
    stepper->last_flag = stepper->ops->reset(stepper, tR, yR);
    return stepper->last_flag;
  }
  else
  {
    /* assume stepper uses input state and does not need to be reset */
    return ARK_SUCCESS;
  }
}

/* Gets the inner (fast) stepper accumulated error */
int mriStepInnerStepper_GetError(MRIStepInnerStepper stepper,
                                 sunrealtype* accum_error)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

  if (stepper->ops->geterror)
  {
    stepper->last_flag = stepper->ops->geterror(stepper, accum_error);
    return stepper->last_flag;
  }
  else
  {
    /* assume stepper provides exact solution */
    *accum_error = SUN_RCONST(0.0);
    return ARK_SUCCESS;
  }
}

/* Resets the inner (fast) stepper accumulated error */
int mriStepInnerStepper_ResetError(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

  if (stepper->ops->geterror)
  {
    stepper->last_flag = stepper->ops->reseterror(stepper);
    return stepper->last_flag;
  }
  else
  {
    /* assume stepper provides exact solution and needs no reset */
    return ARK_SUCCESS;
  }
}

/* Sets the inner (fast) stepper fixed step size */
int mriStepInnerStepper_SetFixedStep(MRIStepInnerStepper stepper, sunrealtype h)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

  if (stepper->ops->setfixedstep)
  {
    stepper->last_flag = stepper->ops->setfixedstep(stepper, h);
    return stepper->last_flag;
  }
  else
  {
    /* assume stepper provides exact solution using infinitesimally small step */
    return ARK_SUCCESS;
  }
}

/* Sets the inner (fast) stepper relative tolerance scaling factor */
int mriStepInnerStepper_SetRTol(MRIStepInnerStepper stepper, sunrealtype rtol)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }
  if (stepper->ops == NULL) { return ARK_ILL_INPUT; }

  if (stepper->ops->setrtol)
  {
    stepper->last_flag = stepper->ops->setrtol(stepper, rtol);
    return stepper->last_flag;
  }
  else
  {
    /* assume stepper provides exact solution */
    return ARK_SUCCESS;
  }
}

/* Allocate MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_AllocVecs(MRIStepInnerStepper stepper, int count,
                                  N_Vector tmpl)
{
  sunindextype lrw1, liw1;

  if (stepper == NULL) { return ARK_ILL_INPUT; }

  /* Set space requirements for one N_Vector */
  if (tmpl->ops->nvspace) { N_VSpace(tmpl, &lrw1, &liw1); }
  else
  {
    lrw1 = 0;
    liw1 = 0;
  }
  stepper->lrw1 = lrw1;
  stepper->liw1 = liw1;

  /* Set the number of forcing vectors and allocate vectors */
  stepper->nforcing = count;

  if (stepper->nforcing_allocated < stepper->nforcing)
  {
    if (stepper->nforcing_allocated)
    {
      arkFreeVecArray(stepper->nforcing_allocated, &(stepper->forcing),
                      stepper->lrw1, &(stepper->lrw), stepper->liw1,
                      &(stepper->liw));
    }
    if (!arkAllocVecArray(stepper->nforcing, tmpl, &(stepper->forcing),
                          stepper->lrw1, &(stepper->lrw), stepper->liw1,
                          &(stepper->liw)))
    {
      mriStepInnerStepper_FreeVecs(stepper);
      return (ARK_MEM_FAIL);
    }
    stepper->nforcing_allocated = stepper->nforcing;
  }

  /* Allocate fused operation workspace arrays */
  if (stepper->vecs == NULL)
  {
    stepper->vecs = (N_Vector*)calloc(count + 1, sizeof(*stepper->vecs));
    if (stepper->vecs == NULL)
    {
      mriStepInnerStepper_FreeVecs(stepper);
      return (ARK_MEM_FAIL);
    }
  }

  if (stepper->vals == NULL)
  {
    stepper->vals = (sunrealtype*)calloc(count + 1, sizeof(*stepper->vals));
    if (stepper->vals == NULL)
    {
      mriStepInnerStepper_FreeVecs(stepper);
      return (ARK_MEM_FAIL);
    }
  }

  return (ARK_SUCCESS);
}

/* Resize MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_Resize(MRIStepInnerStepper stepper, ARKVecResizeFn resize,
                               void* resize_data, sunindextype lrw_diff,
                               sunindextype liw_diff, N_Vector tmpl)
{
  int retval;

  if (stepper == NULL) { return ARK_ILL_INPUT; }

  retval = arkResizeVecArray(resize, resize_data, stepper->nforcing_allocated,
                             tmpl, &(stepper->forcing), lrw_diff,
                             &(stepper->lrw), liw_diff, &(stepper->liw));
  if (retval != ARK_SUCCESS) { return (ARK_MEM_FAIL); }

  return (ARK_SUCCESS);
}

/* Free MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_FreeVecs(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) { return ARK_ILL_INPUT; }

  arkFreeVecArray(stepper->nforcing_allocated, &(stepper->forcing),
                  stepper->lrw1, &(stepper->lrw), stepper->liw1, &(stepper->liw));

  if (stepper->vecs != NULL)
  {
    free(stepper->vecs);
    stepper->vecs = NULL;
  }

  if (stepper->vals != NULL)
  {
    free(stepper->vals);
    stepper->vals = NULL;
  }

  return (ARK_SUCCESS);
}

/* Print forcing vectors to output file */
void mriStepInnerStepper_PrintMem(MRIStepInnerStepper stepper, FILE* outfile)
{
#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif
  if (stepper == NULL) { return; }

  /* output data from the inner stepper */
  fprintf(outfile, "MRIStepInnerStepper Mem:\n");
  fprintf(outfile, "MRIStepInnerStepper: inner_nforcing = %i\n",
          stepper->nforcing);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  if (stepper->forcing != NULL)
  {
    for (i = 0; i < stepper->nforcing; i++)
    {
      fprintf(outfile, "MRIStep: inner_forcing[%i]:\n", i);
      N_VPrintFile(stepper->forcing[i], outfile);
    }
  }
#endif

  return;
}

/*===============================================================
  EOF
  ===============================================================*/
