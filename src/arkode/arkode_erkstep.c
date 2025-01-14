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
 * This is the implementation file for ARKODE's ERK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>

#include "arkode/arkode_butcher.h"
#include "arkode_erkstep_impl.h"
#include "arkode_impl.h"
#include "arkode_interp_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

void* ERKStepCreate(ARKRhsFn f, sunrealtype t0, N_Vector y0, SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunbooleantype nvectorOK;
  int retval;

  /* Check that f is supplied */
  if (f == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (NULL);
  }

  if (!sunctx)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return (NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = erkStep_CheckNVector(y0);
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

  /* Allocate ARKodeERKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeERKStepMem)malloc(sizeof(struct ARKodeERKStepMemRec));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeERKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init                = erkStep_Init;
  ark_mem->step_fullrhs             = erkStep_FullRHS;
  ark_mem->step                     = erkStep_TakeStep;
  ark_mem->step_printallstats       = erkStep_PrintAllStats;
  ark_mem->step_writeparameters     = erkStep_WriteParameters;
  ark_mem->step_resize              = erkStep_Resize;
  ark_mem->step_free                = erkStep_Free;
  ark_mem->step_printmem            = erkStep_PrintMem;
  ark_mem->step_setdefaults         = erkStep_SetDefaults;
  ark_mem->step_setrelaxfn          = erkStep_SetRelaxFn;
  ark_mem->step_setorder            = erkStep_SetOrder;
  ark_mem->step_getnumrhsevals      = erkStep_GetNumRhsEvals;
  ark_mem->step_getestlocalerrors   = erkStep_GetEstLocalErrors;
  ark_mem->step_setforcing          = erkStep_SetInnerForcing;
  ark_mem->step_supports_adaptive   = SUNTRUE;
  ark_mem->step_supports_relaxation = SUNTRUE;
  ark_mem->step_mem                 = (void*)step_mem;

  /* Set default values for optional inputs */
  retval = erkStep_SetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Allocate the general ERK stepper vectors using y0 as a template */
  /* NOTE: F, cvals and Xvecs will be allocated later on
     (based on the number of ERK stages) */

  /* Copy the input parameters into ARKODE state */
  step_mem->f = f;

  /* Update the ARKODE workspace requirements -- UPDATE */
  ark_mem->liw += 41; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
  ark_mem->lrw += 10;

  /* Initialize all the counters */
  step_mem->nfe = 0;

  /* Initialize fused op work space */
  step_mem->cvals        = NULL;
  step_mem->Xvecs        = NULL;
  step_mem->nfusedopvecs = 0;

  /* Initialize external polynomial forcing data */
  step_mem->forcing  = NULL;
  step_mem->nforcing = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  return ((void*)ark_mem);
}

/*---------------------------------------------------------------
  ERKStepReInit:

  This routine re-initializes the ERKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int ERKStepReInit(void* arkode_mem, ARKRhsFn f, sunrealtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that f is supplied */
  if (f == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (ARK_ILL_INPUT);
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f = f;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return (retval);
  }

  /* Initialize all the counters */
  step_mem->nfe = 0;

  return (ARK_SUCCESS);
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_Resize:

  This routine resizes the memory within the ERKStep module.
  ---------------------------------------------------------------*/
int erkStep_Resize(ARKodeMem ark_mem, N_Vector y0,
                   SUNDIALS_MAYBE_UNUSED sunrealtype hscale,
                   SUNDIALS_MAYBE_UNUSED sunrealtype t0, ARKVecResizeFn resize,
                   void* resize_data)
{
  ARKodeERKStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Determine change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Resize the RHS vectors */
  for (i = 0; i < step_mem->stages; i++)
  {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                      &step_mem->F[i]))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_Free frees all ERKStep memory.
  ---------------------------------------------------------------*/
void erkStep_Free(ARKodeMem ark_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeERKStepMem step_mem;

  /* nothing to do if ark_mem is already NULL */
  if (ark_mem == NULL) { return; }

  /* conditional frees on non-NULL ERKStep module */
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeERKStepMem)ark_mem->step_mem;

    /* free the Butcher table */
    if (step_mem->B != NULL)
    {
      ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->B);
      step_mem->B = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }

    /* free the RHS vectors */
    if (step_mem->F != NULL)
    {
      for (j = 0; j < step_mem->stages; j++)
      {
        arkFreeVec(ark_mem, &step_mem->F[j]);
      }
      free(step_mem->F);
      step_mem->F = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL)
    {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= step_mem->nfusedopvecs;
    }
    if (step_mem->Xvecs != NULL)
    {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= step_mem->nfusedopvecs;
    }
    step_mem->nfusedopvecs = 0;

    /* free work arrays for MRI forcing */
    if (step_mem->stage_times)
    {
      free(step_mem->stage_times);
      step_mem->stage_times = NULL;
      ark_mem->lrw -= step_mem->stages;
    }

    if (step_mem->stage_coefs)
    {
      free(step_mem->stage_coefs);
      step_mem->stage_coefs = NULL;
      ark_mem->lrw -= step_mem->stages;
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }
}

/*---------------------------------------------------------------
  erkStep_PrintMem:

  This routine outputs the memory from the ERKStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
void erkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeERKStepMem step_mem;
  int retval;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output integer quantities */
  fprintf(outfile, "ERKStep: q = %i\n", step_mem->q);
  fprintf(outfile, "ERKStep: p = %i\n", step_mem->p);
  fprintf(outfile, "ERKStep: stages = %i\n", step_mem->stages);

  /* output long integer quantities */
  fprintf(outfile, "ERKStep: nfe = %li\n", step_mem->nfe);

  /* output sunrealtype quantities */
  fprintf(outfile, "ERKStep: Butcher table:\n");
  ARKodeButcherTable_Write(step_mem->B, outfile);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  for (i = 0; i < step_mem->stages; i++)
  {
    fprintf(outfile, "ERKStep: F[%i]:\n", i);
    N_VPrintFile(step_mem->F[i], outfile);
  }
#endif
}

/*---------------------------------------------------------------
  erkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK
    stages, method order, or solver options
  - sets the call_fullrhs flag

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int erkStep_Init(ARKodeMem ark_mem, SUNDIALS_MAYBE_UNUSED sunrealtype tout,
                 int init_type)
{
  ARKodeERKStepMem step_mem;
  sunbooleantype reset_efun;
  int retval, j;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return (ARK_SUCCESS);
  }

  /* enforce use of arkEwtSmallReal if using a fixed step size,
     an internal error weight function, and not performing accumulated
     temporal error estimation */
  reset_efun = SUNTRUE;
  if (!ark_mem->fixedstep) { reset_efun = SUNFALSE; }
  if (ark_mem->user_efun) { reset_efun = SUNFALSE; }
  if (ark_mem->AccumErrorType != ARK_ACCUMERROR_NONE) { reset_efun = SUNFALSE; }
  if (reset_efun)
  {
    ark_mem->user_efun = SUNFALSE;
    ark_mem->efun      = arkEwtSetSmallReal;
    ark_mem->e_data    = ark_mem;
  }

  /* Create Butcher table (if not already set) */
  retval = erkStep_SetButcherTable(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Could not create Butcher table");
    return (ARK_ILL_INPUT);
  }

  /* Check that Butcher table are OK */
  retval = erkStep_CheckButcherTable(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Error in Butcher table");
    return (ARK_ILL_INPUT);
  }

  /* Retrieve/store method and embedding orders now that table is finalized */
  step_mem->q = ark_mem->hadapt_mem->q = step_mem->B->q;
  step_mem->p = ark_mem->hadapt_mem->p = step_mem->B->p;

  /* Ensure that if adaptivity or error accumulation is enabled, then
       method includes embedding coefficients */
  if ((!ark_mem->fixedstep || (ark_mem->AccumErrorType != ARK_ACCUMERROR_NONE)) &&
      (step_mem->p == 0))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Temporal error estimation cannot be performed without embedding coefficients");
    return (ARK_ILL_INPUT);
  }

  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate F[0] ... F[stages-1] if needed */
  if (step_mem->F == NULL)
  {
    step_mem->F = (N_Vector*)calloc(step_mem->stages, sizeof(N_Vector));
  }
  for (j = 0; j < step_mem->stages; j++)
  {
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->F[j])))
    {
      return (ARK_MEM_FAIL);
    }
  }
  ark_mem->liw += step_mem->stages; /* pointers */

  /* Allocate reusable arrays for fused vector interface */
  step_mem->nfusedopvecs = 2 * step_mem->stages + 2 + step_mem->nforcing;
  if (step_mem->cvals == NULL)
  {
    step_mem->cvals = (sunrealtype*)calloc(step_mem->nfusedopvecs,
                                           sizeof(sunrealtype));
    if (step_mem->cvals == NULL) { return (ARK_MEM_FAIL); }
    ark_mem->lrw += step_mem->nfusedopvecs;
  }
  if (step_mem->Xvecs == NULL)
  {
    step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL) { return (ARK_MEM_FAIL); }
    ark_mem->liw += step_mem->nfusedopvecs; /* pointers */
  }

  /* Allocate workspace for MRI forcing -- need to allocate here as the
     number of stages may not bet set before this point and we assume
     SetInnerForcing has been called before the first step i.e., methods
     start with a fast integration */
  if (step_mem->nforcing > 0)
  {
    if (!(step_mem->stage_times))
    {
      step_mem->stage_times = (sunrealtype*)calloc(step_mem->stages,
                                                   sizeof(sunrealtype));
      ark_mem->lrw += step_mem->stages;
    }

    if (!(step_mem->stage_coefs))
    {
      step_mem->stage_coefs = (sunrealtype*)calloc(step_mem->stages,
                                                   sizeof(sunrealtype));
      ark_mem->lrw += step_mem->stages;
    }
  }

  /* Override the interpolant degree (if needed), used in arkInitialSetup */
  if (step_mem->q > 1 && ark_mem->interp_degree > (step_mem->q - 1))
  {
    /* Limit max degree to at most one less than the method global order */
    ark_mem->interp_degree = step_mem->q - 1;
  }
  else if (step_mem->q == 1 && ark_mem->interp_degree > 1)
  {
    /* Allow for linear interpolant with first order methods to ensure
       solution values are returned at the time interval end points */
    ark_mem->interp_degree = 1;
  }

  /* Signal to shared arkode module that full RHS evaluations are required */
  ark_mem->call_fullrhs = SUNTRUE;

  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  erkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function, f(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called in the following circumstances:
                          (a) at the beginning of a simulation i.e., at
                              (tn, yn) = (t0, y0) or (tR, yR),
                          (b) when transitioning between time steps t_{n-1}
                              \to t_{n} to fill f_{n-1} within the Hermite
                              interpolation module, or
                          (c) by ERKStep at the start of the first internal step.

                          In each case, we may check the fn_is_current flag to
                          know whether the values stored in F[0] are up-to-date,
                          allowing us to copy those values instead of recomputing.
                          If these values are not current, then the RHS should be
                          stored in F[0] for reuse later, before copying the values
                          into the output vector.

     ARK_FULLRHS_END   -> called in the following circumstances:
                          (a) when temporal root-finding is enabled, this will be
                              called in-between steps t_{n-1} \to t_{n} to fill f_{n},
                          (b) when high-order dense output is requested from the
                              Hermite interpolation module in-between steps t_{n-1}
                              \to t_{n} to fill f_{n}, or
                          (c) by ERKStep when starting a time step t_{n} \to t_{n+1}
                              and when using an FSAL method.

                          Again, we may check the fn_is_current flag to know whether
                          ARKODE believes that the values stored in F[0] are
                          up-to-date, and may just be copied.  If the values stored
                          in F[0] are not current, then the only instance where
                          recomputation is not needed is (c), since the values in
                          F[stages - 1] may be copied into F[0].  In all other cases,
                          the RHS should be recomputed and stored in F[0] for reuse
                          later, before copying the values into the output vector.

     ARK_FULLRHS_OTHER -> called in the following circumstances:
                          (a) when estimating the initial time step size,
                          (b) for high-order dense output with the Hermite
                              interpolation module, or
                          (c) by an "outer" stepper when ERKStep is used as an
                              inner solver).

                          All of these instances will occur in-between ERKStep time
                          steps, but the (t,y) input does not correspond to an
                          "official" time step, thus the RHS should always be
                          evaluated, with the values *not* stored in F[0].
  ----------------------------------------------------------------------------*/
int erkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode)
{
  int nvec, retval;
  ARKodeERKStepMem step_mem;
  sunbooleantype recomputeRHS;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  sunrealtype stage_coefs = ONE;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:

    /* compute the RHS if needed */
    if (!(ark_mem->fn_is_current))
    {
      retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* copy RHS into output */
    N_VScale(ONE, step_mem->F[0], f);

    /* apply external polynomial forcing */
    if (step_mem->nforcing > 0)
    {
      cvals[0] = ONE;
      Xvecs[0] = f;
      nvec     = 1;
      erkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
      N_VLinearCombination(nvec, cvals, Xvecs, f);
    }

    break;

  case ARK_FULLRHS_END:

    /* determine if RHS function needs to be recomputed */
    if (!(ark_mem->fn_is_current))
    {
      recomputeRHS = !ARKodeButcherTable_IsStifflyAccurate(step_mem->B);

      /* First Same As Last methods are not FSAL when relaxation is enabled */
      if (ark_mem->relax_enabled) { recomputeRHS = SUNTRUE; }

      /* base RHS call on recomputeRHS argument */
      if (recomputeRHS)
      {
        /* call f */
        retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
        step_mem->nfe++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
          return (ARK_RHSFUNC_FAIL);
        }
      }
      else { N_VScale(ONE, step_mem->F[step_mem->stages - 1], step_mem->F[0]); }

      /* copy RHS vector into output */
      N_VScale(ONE, step_mem->F[0], f);

      /* apply external polynomial forcing */
      if (step_mem->nforcing > 0)
      {
        cvals[0] = ONE;
        Xvecs[0] = f;
        nvec     = 1;
        erkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, f);
      }
    }

    break;

  case ARK_FULLRHS_OTHER:

    /* call f */
    retval = step_mem->f(t, y, f, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }
    /* apply external polynomial forcing */
    if (step_mem->nforcing > 0)
    {
      cvals[0] = ONE;
      Xvecs[0] = f;
      nvec     = 1;
      erkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
      N_VLinearCombination(nvec, cvals, Xvecs, f);
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
  erkStep_TakeStep:

  This routine serves the primary purpose of the ERKStep module:
  it performs a single ERK step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int erkStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval, is, js, nvec, mode;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success */
  *nflagPtr = ARK_SUCCESS;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  SUNLogInfo(ARK_LOGGER, "begin-stage", "stage = 0, tcur = %" RSYM,
             ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current))
  {
    mode   = (ark_mem->initsetup) ? ARK_FULLRHS_START : ARK_FULLRHS_END;
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                   ark_mem->fn, mode);
    if (retval)
    {
      SUNLogInfo(ARK_LOGGER, "end-stage",
                 "status = failed rhs eval, retval = %i", retval);
      return ARK_RHSFUNC_FAIL;
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->F[0], "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stage", "status = success");

  /* Loop over internal stages to the step; since the method is explicit
     the first stage RHS is just the full RHS from the start of the step */
  for (is = 1; is < step_mem->stages; is++)
  {
    /* Set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + step_mem->B->c[is] * ark_mem->h;

    SUNLogInfo(ARK_LOGGER, "begin-stage", "stage = %i, tcur = %" RSYM, is,
               ark_mem->tcur);

    /* Set ycur to current stage solution */
    nvec = 0;
    for (js = 0; js < is; js++)
    {
      cvals[nvec] = ark_mem->h * step_mem->B->A[is][js];
      Xvecs[nvec] = step_mem->F[js];
      nvec += 1;
    }
    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec += 1;

    /* apply external polynomial forcing */
    if (step_mem->nforcing > 0)
    {
      for (js = 0; js < is; js++)
      {
        step_mem->stage_times[js] = ark_mem->tn + step_mem->B->c[js] * ark_mem->h;
        step_mem->stage_coefs[js] = ark_mem->h * step_mem->B->A[is][js];
      }
      erkStep_ApplyForcing(step_mem, step_mem->stage_times,
                           step_mem->stage_coefs, is, &nvec);
    }

    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stage",
                 "status = failed vector op, retval = %i", retval);
      return (ARK_VECTOROP_ERR);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->ProcessStage != NULL)
    {
      retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stage",
                   "status = failed postprocess stage, retval = %i", retval);
        return (ARK_POSTPROCESS_STAGE_FAIL);
      }
    }

    /* compute updated RHS */
    retval = step_mem->f(ark_mem->tcur, ark_mem->ycur, step_mem->F[is],
                         ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->F[is],
                        "F_%i(:) =", is);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stage",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
    if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }

    SUNLogInfo(ARK_LOGGER, "end-stage", "status = success");

  } /* loop over stages */

  SUNLogInfo(ARK_LOGGER, "begin-compute-solution", "");

  /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
  retval = erkStep_ComputeSolutions(ark_mem, dsmPtr);
  if (retval < 0)
  {
    SUNLogInfo(ARK_LOGGER, "end-compute-solution",
               "status = failed compute solution, retval = %i", retval);
    return (retval);
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");
  SUNLogInfo(ARK_LOGGER, "end-compute-solution", "status = success");

  return (ARK_SUCCESS);
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_AccessARKODEStepMem:

  Shortcut routine to unpack both ark_mem and step_mem structures
  from void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int erkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeERKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  /* access ARKodeERKStepMem structure */
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeERKStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_AccessStepMem:

  Shortcut routine to unpack the step_mem structure from
  ark_mem.  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int erkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeERKStepMem* step_mem)
{
  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeERKStepMem)ark_mem->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype erkStep_CheckNVector(N_Vector tmpl)
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
  erkStep_SetButcherTable

  This routine determines the ERK method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
int erkStep_SetButcherTable(ARKodeMem ark_mem)
{
  int etable;
  ARKodeERKStepMem step_mem;
  sunindextype Bliw, Blrw;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem)ark_mem->step_mem;

  /* if table has already been specified, just return */
  if (step_mem->B != NULL) { return (ARK_SUCCESS); }

  /* initialize table number to illegal values */
  etable = -1;

  /* select method based on order */
  switch (step_mem->q)
  {
  case (1): etable = ERKSTEP_DEFAULT_1; break;
  case (2): etable = ERKSTEP_DEFAULT_2; break;
  case (3): etable = ERKSTEP_DEFAULT_3; break;
  case (4): etable = ERKSTEP_DEFAULT_4; break;
  case (5): etable = ERKSTEP_DEFAULT_5; break;
  case (6): etable = ERKSTEP_DEFAULT_6; break;
  case (7): etable = ERKSTEP_DEFAULT_7; break;
  case (8): etable = ERKSTEP_DEFAULT_8; break;
  case (9): etable = ERKSTEP_DEFAULT_9; break;
  default: /* no available method, set default */
    arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                    "No explicit method at requested order, using q=9.");
    etable = ERKSTEP_DEFAULT_9;
    break;
  }

  if (etable > -1) { step_mem->B = ARKodeButcherTable_LoadERK(etable); }

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* set [redundant] stored values for stage numbers and method orders */
  if (step_mem->B != NULL)
  {
    step_mem->stages = step_mem->B->stages;
    step_mem->q      = step_mem->B->q;
    step_mem->p      = step_mem->B->p;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_CheckButcherTable

  This routine runs through the explicit Butcher table to ensure
  that it meets all necessary requirements, including:
    strictly lower-triangular (ERK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int erkStep_CheckButcherTable(ARKodeMem ark_mem)
{
  int i, j;
  sunbooleantype okay;
  ARKodeERKStepMem step_mem;
  sunrealtype tol = SUN_RCONST(1.0e-12);

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem)ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->stages < 1)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "stages < 1!");
    return (ARK_INVALID_TABLE);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "method order < 1!");
    return (ARK_INVALID_TABLE);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep))
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "embedding order < 1!");
    return (ARK_INVALID_TABLE);
  }

  /* check that embedding exists */
  if ((step_mem->p > 0) && (!ark_mem->fixedstep))
  {
    if (step_mem->B->d == NULL)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "no embedding!");
      return (ARK_INVALID_TABLE);
    }
  }

  /* check that ERK table is strictly lower triangular */
  okay = SUNTRUE;
  for (i = 0; i < step_mem->stages; i++)
  {
    for (j = i; j < step_mem->stages; j++)
    {
      if (SUNRabs(step_mem->B->A[i][j]) > tol) { okay = SUNFALSE; }
    }
  }
  if (!okay)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "Ae Butcher table is implicit!");
    return (ARK_INVALID_TABLE);
  }

  /* check if all b values are positive for relaxation */
  if (ark_mem->relax_enabled)
  {
    if (step_mem->q < 2)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "The Butcher table must be at least second order!");
      return ARK_INVALID_TABLE;
    }

    for (i = 0; i < step_mem->stages; i++)
    {
      if (step_mem->B->b[i] < ZERO)
      {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__,
                        __FILE__, "The Butcher table has a negative b value!");
        return ARK_INVALID_TABLE;
      }
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_tempv1, in case the
  calling routine wishes to examine the error locations.

  Note: at this point in the step, the vector ark_tempv1 may be
  used as a temporary vector.
  ---------------------------------------------------------------*/
int erkStep_ComputeSolutions(ARKodeMem ark_mem, sunrealtype* dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem)ark_mem->step_mem;

  /* set N_Vector shortcuts */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;

  /* Compute time step solution */
  /*   set arrays for fused vector operation */
  nvec = 0;
  for (j = 0; j < step_mem->stages; j++)
  {
    cvals[nvec] = ark_mem->h * step_mem->B->b[j];
    Xvecs[nvec] = step_mem->F[j];
    nvec += 1;
  }
  cvals[nvec] = ONE;
  Xvecs[nvec] = ark_mem->yn;
  nvec += 1;

  /* apply external polynomial forcing */
  if (step_mem->nforcing > 0)
  {
    for (j = 0; j < step_mem->stages; j++)
    {
      step_mem->stage_times[j] = ark_mem->tn + step_mem->B->c[j] * ark_mem->h;
      step_mem->stage_coefs[j] = ark_mem->h * step_mem->B->b[j];
    }
    erkStep_ApplyForcing(step_mem, step_mem->stage_times, step_mem->stage_coefs,
                         step_mem->stages, &nvec);
  }

  /*   call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
  if (retval != 0) { return (ARK_VECTOROP_ERR); }

  /* Compute yerr (if step adaptivity or error accumulation enabled) */
  if (!ark_mem->fixedstep || (ark_mem->AccumErrorType != ARK_ACCUMERROR_NONE))
  {
    /* set arrays for fused vector operation */
    nvec = 0;
    for (j = 0; j < step_mem->stages; j++)
    {
      cvals[nvec] = ark_mem->h * (step_mem->B->b[j] - step_mem->B->d[j]);
      Xvecs[nvec] = step_mem->F[j];
      nvec += 1;
    }

    /* apply external polynomial forcing */
    if (step_mem->nforcing > 0)
    {
      for (j = 0; j < step_mem->stages; j++)
      {
        step_mem->stage_times[j] = ark_mem->tn + step_mem->B->c[j] * ark_mem->h;
        step_mem->stage_coefs[j] = ark_mem->h *
                                   (step_mem->B->b[j] - step_mem->B->d[j]);
      }
      erkStep_ApplyForcing(step_mem, step_mem->stage_times,
                           step_mem->stage_coefs, step_mem->stages, &nvec);
    }

    /* call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    /* fill error norm */
    *dsmPtr = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  Internal utility routines for relaxation
  ===============================================================*/

/* -----------------------------------------------------------------------------
 * erkStep_RelaxDeltaE
 *
 * Computes the change in the relaxation functions for use in relaxation methods
 * delta_e = h * sum_i b_i * <rjac(z_i), f_i>
 * ---------------------------------------------------------------------------*/
int erkStep_RelaxDeltaE(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                        long int* num_relax_jac_evals, sunrealtype* delta_e_out)
{
  int i, j, nvec, retval;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;
  N_Vector z_stage = ark_mem->tempv2;
  N_Vector J_relax = ark_mem->tempv3;

  /* Access the stepper memory structure */
  if (!(ark_mem->step_mem))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ERKSTEP_NO_MEM);
    return ARK_MEM_NULL;
  }
  step_mem = (ARKodeERKStepMem)(ark_mem->step_mem);

  /* Initialize output */
  *delta_e_out = ZERO;

  /* Set arrays for fused vector operation */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  for (i = 0; i < step_mem->stages; i++)
  {
    /* Construct stages z[i] = y_n + h * sum_j Ae[i,j] Fe[j] + Ai[i,j] Fi[j] */
    nvec = 0;

    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec++;

    for (j = 0; j < i; j++)
    {
      cvals[nvec] = ark_mem->h * step_mem->B->A[i][j];
      Xvecs[nvec] = step_mem->F[j];
      nvec++;
    }

    retval = N_VLinearCombination(nvec, cvals, Xvecs, z_stage);
    if (retval) { return ARK_VECTOROP_ERR; }

    /* Evaluate the Jacobian at z_i */
    retval = relax_jac_fn(z_stage, J_relax, ark_mem->user_data);
    (*num_relax_jac_evals)++;
    if (retval < 0) { return ARK_RELAX_JAC_FAIL; }
    if (retval > 0) { return ARK_RELAX_JAC_RECV; }

    /* Update estimates */
    if (J_relax->ops->nvdotprodlocal && J_relax->ops->nvdotprodmultiallreduce)
    {
      *delta_e_out += step_mem->B->b[i] *
                      N_VDotProdLocal(J_relax, step_mem->F[i]);
    }
    else
    {
      *delta_e_out += step_mem->B->b[i] * N_VDotProd(J_relax, step_mem->F[i]);
    }
  }

  if (J_relax->ops->nvdotprodlocal && J_relax->ops->nvdotprodmultiallreduce)
  {
    retval = N_VDotProdMultiAllReduce(1, J_relax, delta_e_out);
    if (retval) { return ARK_VECTOROP_ERR; }
  }

  *delta_e_out *= ark_mem->h;

  return ARK_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * erkStep_GetOrder
 *
 * Returns the method order
 * ---------------------------------------------------------------------------*/
int erkStep_GetOrder(ARKodeMem ark_mem)
{
  ARKodeERKStepMem step_mem = (ARKodeERKStepMem)(ark_mem->step_mem);
  return step_mem->q;
}

/*---------------------------------------------------------------
  Utility routines for ERKStep to serve as an MRIStepInnerStepper
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  erkStep_ApplyForcing

  Determines the linear combination coefficients and vectors to apply forcing
  at a given value of the independent variable (t).  This occurs through
  appending coefficients and N_Vector pointers to the underlying cvals and Xvecs
  arrays in the step_mem structure.  The dereferenced input *nvec should indicate
  the next available entry in the cvals/Xvecs arrays.  The input 's' is a
  scaling factor that should be applied to each of these coefficients.
  ----------------------------------------------------------------------------*/

void erkStep_ApplyForcing(ARKodeERKStepMem step_mem, sunrealtype* stage_times,
                          sunrealtype* stage_coefs, int jmax, int* nvec)
{
  sunrealtype tau, taui;
  int j, k;

  /* Shortcuts to step_mem data */
  sunrealtype* vals  = step_mem->cvals;
  N_Vector* vecs     = step_mem->Xvecs;
  sunrealtype tshift = step_mem->tshift;
  sunrealtype tscale = step_mem->tscale;
  int nforcing       = step_mem->nforcing;
  N_Vector* forcing  = step_mem->forcing;

  /* Offset into vals and vecs arrays */
  int offset = *nvec;

  /* Initialize scaling values, set vectors */
  for (k = 0; k < nforcing; k++)
  {
    vals[offset + k] = ZERO;
    vecs[offset + k] = forcing[k];
  }

  for (j = 0; j < jmax; j++)
  {
    tau  = (stage_times[j] - tshift) / tscale;
    taui = ONE;

    for (k = 0; k < nforcing; k++)
    {
      vals[offset + k] += stage_coefs[j] * taui;
      taui *= tau;
    }
  }

  /* Update vector count for linear combination */
  *nvec += nforcing;
}

/*------------------------------------------------------------------------------
  erkStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = f(t,y) + p(t). This function is
  primarily intended for use with multirate integration methods (e.g., MRIStep)
  where ERKStep is used to solve a modified ODE at a fast time scale. The
  polynomial is of the form

  p(t) = sum_{i = 0}^{nvecs - 1} forcing[i] * ((t - tshift) / (tscale))^i

  where tshift and tscale are used to normalize the time t (e.g., with MRIGARK
  methods).
  ----------------------------------------------------------------------------*/

int erkStep_SetInnerForcing(ARKodeMem ark_mem, sunrealtype tshift,
                            sunrealtype tscale, N_Vector* forcing, int nvecs)
{
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (nvecs > 0)
  {
    /* store forcing inputs */
    step_mem->tshift   = tshift;
    step_mem->tscale   = tscale;
    step_mem->forcing  = forcing;
    step_mem->nforcing = nvecs;

    /* If cvals and Xvecs are not allocated then erkStep_Init has not been
       called and the number of stages has not been set yet. These arrays will
       be allocated in erkStep_Init and take into account the value of nforcing.
       On subsequent calls will check if enough space has allocated in case
       nforcing has increased since the original allocation. */
    if (step_mem->cvals != NULL && step_mem->Xvecs != NULL)
    {
      /* check if there are enough reusable arrays for fused operations */
      if ((step_mem->nfusedopvecs - nvecs) < (step_mem->stages + 1))
      {
        /* free current work space */
        if (step_mem->cvals != NULL)
        {
          free(step_mem->cvals);
          ark_mem->lrw -= step_mem->nfusedopvecs;
        }
        if (step_mem->Xvecs != NULL)
        {
          free(step_mem->Xvecs);
          ark_mem->liw -= step_mem->nfusedopvecs;
        }

        /* allocate reusable arrays for fused vector operations */
        step_mem->nfusedopvecs = step_mem->stages + 1 + nvecs;

        step_mem->cvals = NULL;
        step_mem->cvals = (sunrealtype*)calloc(step_mem->nfusedopvecs,
                                               sizeof(sunrealtype));
        if (step_mem->cvals == NULL) { return (ARK_MEM_FAIL); }
        ark_mem->lrw += step_mem->nfusedopvecs;

        step_mem->Xvecs = NULL;
        step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs,
                                            sizeof(N_Vector));
        if (step_mem->Xvecs == NULL) { return (ARK_MEM_FAIL); }
        ark_mem->liw += step_mem->nfusedopvecs;
      }
    }
  }
  else
  {
    /* disable forcing */
    step_mem->tshift   = ZERO;
    step_mem->tscale   = ONE;
    step_mem->forcing  = NULL;
    step_mem->nforcing = 0;
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
