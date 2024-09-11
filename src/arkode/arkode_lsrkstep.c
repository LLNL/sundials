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
 * This is the implementation file for ARKODE's LSRK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_lsrkstep_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

void* LSRKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0,
                     SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  sunbooleantype nvectorOK;
  int retval;

  /* Check that fe is supplied */
  if (fe == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  /* Check that fi is NULL until IMEX module is ready */
  if (fi != NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "\n\nNO IMEX-LSRK support yet, set fi = NULL\n");
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
  nvectorOK = lsrkStep_CheckNVector(y0);
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

  /* Allocate ARKodeLSRKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeLSRKStepMem)malloc(sizeof(struct ARKodeLSRKStepMemRec));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeLSRKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init              = lsrkStep_Init;
  ark_mem->step_fullrhs           = lsrkStep_FullRHS;
  ark_mem->step                   = lsrkStep_TakeStepRKC;
  ark_mem->step_printallstats     = lsrkStep_PrintAllStats;
  ark_mem->step_writeparameters   = lsrkStep_WriteParameters;
  ark_mem->step_resize            = lsrkStep_Resize;
  ark_mem->step_free              = lsrkStep_Free;
  ark_mem->step_printmem          = lsrkStep_PrintMem;
  ark_mem->step_setdefaults       = lsrkStep_SetDefaults;
  ark_mem->step_getestlocalerrors = lsrkStep_GetEstLocalErrors;
  ark_mem->step_mem               = (void*)step_mem;
  ark_mem->step_supports_adaptive = SUNTRUE;

  /* Set default values for optional inputs */
  retval = lsrkStep_SetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Allocate the general LSRK stepper vectors using y0 as a template */

  /* Copy the input parameters into ARKODE state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Update the ARKODE workspace requirements -- UPDATE */
  ark_mem->liw += 0; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
  ark_mem->lrw += 0;

  /* Initialize all the counters */
  step_mem->nfe = 0;
  step_mem->nfi = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Specify preferred interpolation type */
  ark_mem->interp_type = ARK_INTERP_LAGRANGE;

  return ((void*)ark_mem);
}

/*---------------------------------------------------------------
  LSRKStepReInit:

  This routine re-initializes the LSRKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int LSRKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                   N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that fe is supplied */
  if (fe == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (ARK_ILL_INPUT);
  }

  /* Check that fi is NULL until IMEX module is ready */
  if (fi != NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "\n\nNO IMEX-LSRK support yet, set fi = NULL");
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
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return (retval);
  }

  /* Initialize all the counters, flags and stats */
  step_mem->nfe = 0;
  step_mem->nfi = 0;
  step_mem->domeignfe = 0;
  step_mem->ndomeigupdates = 0;
  step_mem->stagemax = 0;
  step_mem->sprmax = 0;
  step_mem->sprmin = 0;  
  step_mem->nstsig = 0;
  step_mem->newdomeig = SUNTRUE;
  step_mem->jacatt    = SUNFALSE;

  return (ARK_SUCCESS);
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_Resize:

  This routine resizes the memory within the LSRKStep module.
  ---------------------------------------------------------------*/
int lsrkStep_Resize(ARKodeMem ark_mem, N_Vector y0,
                    SUNDIALS_MAYBE_UNUSED sunrealtype hscale,
                    SUNDIALS_MAYBE_UNUSED sunrealtype t0, ARKVecResizeFn resize,
                    void* resize_data)
{
  ARKodeLSRKStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Determine change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Resize the internal vector storage */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                    &step_mem->Fe[0]))
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_Free frees all LSRKStep memory.
  ---------------------------------------------------------------*/
void lsrkStep_Free(ARKodeMem ark_mem)
{
  ARKodeLSRKStepMem step_mem;

  /* nothing to do if ark_mem is already NULL */
  if (ark_mem == NULL) { return; }

  /* conditional frees on non-NULL LSRKStep module */
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeLSRKStepMem)ark_mem->step_mem;

    /* free the RHS vectors */
    if (step_mem->Fe != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->Fe[0]);
      free(step_mem->Fe);
      step_mem->Fe = NULL;
      ark_mem->liw -= 1;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL)
    {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= 5;
    }
    if (step_mem->Xvecs != NULL)
    {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= 5;
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }
}

/*---------------------------------------------------------------
  lsrkStep_PrintMem:

  This routine outputs the memory from the LSRKStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
void lsrkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities - required allocations*/

#endif

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output integer quantities */
  fprintf(outfile, "LSRKStep: reqstages      = %i\n", step_mem->reqstages);
  fprintf(outfile, "LSRKStep: nstsig         = %i\n", step_mem->nstsig);
  fprintf(outfile, "LSRKStep: stagemax       = %i\n", step_mem->stagemax);
  fprintf(outfile, "LSRKStep: stagemaxlimit  = %i\n", step_mem->stagemaxlimit);
  fprintf(outfile, "LSRKStep: domeigfreq     = %i\n", step_mem->domeigfreq);

  /* output long integer quantities */
  fprintf(outfile, "LSRKStep: nfe            = %li\n", step_mem->nfe);
  fprintf(outfile, "LSRKStep: domeignfe      = %li\n", step_mem->domeignfe);
  fprintf(outfile, "LSRKStep: ndomeigupdates = %li\n", step_mem->ndomeigupdates);

  /* output sunrealtype quantities */
  fprintf(outfile, "LSRKStep: DomEig        = %f + i%f\n", step_mem->lambdaR,
          step_mem->lambdaI);
  fprintf(outfile, "LSRKStep: sprad         = %f\n", step_mem->sprad);
  fprintf(outfile, "LSRKStep: sprmax        = %f\n", step_mem->sprmax);
  fprintf(outfile, "LSRKStep: sprmin        = %f\n", step_mem->sprmin);
  fprintf(outfile, "LSRKStep: domeigsfty    = %f\n", step_mem->domeigsfty);

  /* output sunbooleantype quantities */
  fprintf(outfile, "LSRKStep: isextDomEig  = %d\n", step_mem->isextDomEig);
  fprintf(outfile, "LSRKStep: newdomeig    = %d\n", step_mem->newdomeig);
  fprintf(outfile, "LSRKStep: jacatt       = %d\n", step_mem->jacatt);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */

#endif
}

/*---------------------------------------------------------------
  lsrkStep_Init:

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
int lsrkStep_Init(ARKodeMem ark_mem, int init_type)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return (ARK_SUCCESS);
  }
  /* enforce use of arkEwtSmallReal if using a fixed step size
     and an internal error weight function */
  if (ark_mem->fixedstep && !ark_mem->user_efun)
  {
    ark_mem->user_efun = SUNFALSE;
    ark_mem->efun      = arkEwtSetSmallReal;
    ark_mem->e_data    = ark_mem;
  }

  /* Store method and embedding orders now that LSRK method choice is finalized */
  step_mem->q = ark_mem->hadapt_mem->q = 2;
  step_mem->p = ark_mem->hadapt_mem->p = 2;

  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate Fe if needed */
  if (step_mem->Fe == NULL)
  {
    step_mem->Fe = (N_Vector*)calloc(1, sizeof(N_Vector));
  }
  if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fe[0])))
  {
    return (ARK_MEM_FAIL);
  }
  ark_mem->liw += 1; /* pointers */

  /* Allocate reusable arrays for fused vector interface */
  if (step_mem->cvals == NULL)
  {
    step_mem->cvals = (sunrealtype*)calloc(5, sizeof(sunrealtype));
    if (step_mem->cvals == NULL) { return (ARK_MEM_FAIL); }
    ark_mem->lrw += 5;
  }
  if (step_mem->Xvecs == NULL)
  {
    step_mem->Xvecs = (N_Vector*)calloc(5, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL) { return (ARK_MEM_FAIL); }
    ark_mem->liw += 5; /* pointers */
  }

  /* Signal to shared arkode module that full RHS evaluations are required */
  ark_mem->call_fullrhs = SUNTRUE;

  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  lsrkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function, f(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If this function is called in ARK_FULLRHS_START or ARK_FULLRHS_END mode and
  evaluating the RHS functions is necessary, we store the vector f(t,y) in Fe[0]
  for reuse in the first stage of the subsequent time step.

  In ARK_FULLRHS_END mode we check if the method is "stiffly accurate" and, if
  appropriate, copy the vector F[stages - 1] to F[0] for reuse in the first
  stage of the subsequent time step.

  ARK_FULLRHS_OTHER mode is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to store the
  intermediate parts so that they do not interfere with the other two modes.
  ----------------------------------------------------------------------------*/
int lsrkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:

    /* compute the RHS */
    if (!(ark_mem->fn_is_current))
    {
      retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* copy RHS vector into output */
    N_VScale(ONE, step_mem->Fe[0], f);

    break;

  case ARK_FULLRHS_END:
    /* No further action is needed since the currently 
    available methods evaluate the RHS at the end of each time step*/

    break;

  case ARK_FULLRHS_OTHER:

    /* call f */
    retval = step_mem->fe(t, y, f, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
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
  lsrkStep_TakeStepRKC:

  This routine serves the primary purpose of the LSRKStepRKC module:
  it performs a single RKC step (with embedding, if possible).

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
int lsrkStep_TakeStepRKC(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  sunrealtype* cvals;
  sunrealtype w0, w1, temp1, temp2, arg, bjm1, bjm2, mus, thjm1, thjm2, zjm1,
    zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2, zj, dzj, d2zj, bj, ajm1, mu, nu, thj;
  sunrealtype onep54 = 1.54, c13 = 13.0, p8 = 0.8, p4 = 0.4;
  N_Vector* Xvecs;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Compute Dominated Eigenvalue and update stats */
  if ((step_mem->newdomeig))
  {
    retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
    if (retval != ARK_SUCCESS) { return (retval); }
    step_mem->ndomeigupdates++;
  }

  /* determine the number of required stages */
  for (int ss = 1; ss < step_mem->stagemaxlimit; ss++)
  {
    if (SUNSQR(ss) >= (onep54 * SUNRabs(ark_mem->h) * step_mem->sprad))
    {
      step_mem->reqstages = SUNMAX(ss, 2);
      break;
    }
  }
  step_mem->stagemax = SUNMAX(step_mem->reqstages, step_mem->stagemax);

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current) && ark_mem->initsetup)
  {
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  /* A tentative solution at t+h is returned in
     y and its slope is evaluated in temp1.  */

  w0 = (ONE + TWO / (c13 * SUNSQR((sunrealtype)(step_mem->reqstages))));

  temp1 = SUNSQR(w0) - ONE;
  temp2 = SUNRsqrt(temp1);
  arg   = step_mem->reqstages * log(w0 + temp2);

  w1 = sinh(arg) * temp1 /
       (cosh(arg) * step_mem->reqstages * temp2 - w0 * sinh(arg));

  bjm1 = ONE / SUNSQR(TWO * w0);
  bjm2 = bjm1;

  /* Evaluate the first stage */
  N_VScale(ONE, ark_mem->yn, ark_mem->tempv1);

  mus = w1 * bjm1;

  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h * mus, ark_mem->fn, ark_mem->tempv2);

  thjm2  = ZERO;
  thjm1  = mus;
  zjm1   = w0;
  zjm2   = ONE;
  dzjm1  = ONE;
  dzjm2  = ZERO;
  d2zjm1 = ZERO;
  d2zjm2 = ZERO;

  /* Evaluate stages j = 2,...,step_mem->reqstages */
  for (int j = 2; j <= step_mem->reqstages; j++)
  {
    zj   = TWO * w0 * zjm1 - zjm2;
    dzj  = TWO * w0 * dzjm1 - dzjm2 + TWO * zjm1;
    d2zj = TWO * w0 * d2zjm1 - d2zjm2 + FOUR * dzjm1;
    bj   = d2zj / SUNSQR(dzj);
    ajm1 = ONE - zjm1 * bjm1;
    mu   = TWO * w0 * bj / bjm1;
    nu   = -bj / bjm2;
    mus  = mu * w1 / w0;

    /* Use the ycur array for temporary storage here */
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h * thjm1, ark_mem->tempv2,
                          ark_mem->ycur, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    cvals[0] = mus * ark_mem->h;
    Xvecs[0] = ark_mem->ycur;
    cvals[1] = nu;
    Xvecs[1] = ark_mem->tempv1;
    cvals[2] = ONE - mu - nu;
    Xvecs[2] = ark_mem->yn;
    cvals[3] = mu;
    Xvecs[3] = ark_mem->tempv2;
    cvals[4] = -mus * ajm1 * ark_mem->h;
    Xvecs[4] = ark_mem->fn;

    retval = N_VLinearCombination(5, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    thj = mu * thjm1 + nu * thjm2 + mus * (ONE - ajm1);

    /* Shift the data for the next stage */
    if (j < step_mem->reqstages)
    {
      /* To avoid two data copies we swap ARKODE's tempv1 and tempv2 pointers*/
      N_Vector* ptrtempv1 = &(ark_mem->tempv1);
      N_Vector* ptrtempv2 = &(ark_mem->tempv2);

      N_Vector temp = *ptrtempv2;
      *ptrtempv2    = *ptrtempv1;
      *ptrtempv1    = temp;

      N_VScale(ONE, ark_mem->ycur, ark_mem->tempv2);

      thjm2  = thjm1;
      thjm1  = thj;
      bjm2   = bjm1;
      bjm1   = bj;
      zjm2   = zjm1;
      zjm1   = zj;
      dzjm2  = dzjm1;
      dzjm1  = dzj;
      d2zjm2 = d2zjm1;
      d2zjm1 = d2zj;
    }
  }
  step_mem->nfe += step_mem->reqstages;

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur,
                          ark_mem->tempv2, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    /* Estimate the local error and compute its weighted RMS norm */
    cvals[0] = p8;
    Xvecs[0] = ark_mem->yn;
    cvals[1] = -p8;
    Xvecs[1] = ark_mem->ycur;
    cvals[2] = p4 * ark_mem->h;
    Xvecs[2] = ark_mem->fn;
    cvals[3] = p4 * ark_mem->h;
    Xvecs[3] = ark_mem->tempv2;

    retval = N_VLinearCombination(4, cvals, Xvecs, ark_mem->tempv1);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
    lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr);
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepRKL:

  This routine serves the primary purpose of the LSRKStepRKL module:
  it performs a single RKL step (with embedding, if possible).

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
int lsrkStep_TakeStepRKL(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  sunrealtype* cvals;
  sunrealtype w1, bjm1, bjm2, mus, bj, ajm1, cjm1, temj, cj, mu, nu;
  sunrealtype p8 = 0.8, p4 = 0.4;
  N_Vector* Xvecs;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Compute Dominated Eigenvalue and update stats */
  if ((step_mem->newdomeig))
  {
    retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
    if (retval != ARK_SUCCESS) { return (retval); }
    step_mem->ndomeigupdates++;
  }

  /* determine the number of required stages */
  for (int ss = 1; ss < step_mem->stagemaxlimit; ss++)
  {
    if ((SUNSQR(ss) + ss - 2) >= 2 * (SUNRabs(ark_mem->h) * step_mem->sprad))
    {
      step_mem->reqstages = SUNMAX(ss, 2);
      break;
    }
  }
  step_mem->stagemax = SUNMAX(step_mem->reqstages, step_mem->stagemax);

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current) && ark_mem->initsetup)
  {
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  /* A tentative solution at t+h is returned in
     y and its slope is evaluated in temp1.  */

  w1 = FOUR / ((step_mem->reqstages + TWO) * (step_mem->reqstages - ONE));

  bjm2 = ONE / THREE;
  bjm1 = bjm2;

  /* Evaluate the first stage */
  N_VScale(ONE, ark_mem->yn, ark_mem->tempv1);

  mus  = w1 * bjm1;
  cjm1 = mus;

  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h * mus, ark_mem->fn, ark_mem->tempv2);

  /* Evaluate stages j = 2,...,step_mem->reqstages */
  for (int j = 2; j <= step_mem->reqstages; j++)
  {
    temj = (j + TWO) * (j - ONE);
    bj   = temj / (TWO * j * (j + ONE));
    ajm1 = ONE - bjm1;
    mu   = (TWO * j - ONE) / j * (bj / bjm1);
    nu   = -(j - ONE) / j * (bj / bjm2);
    mus  = w1 * mu;
    cj   = temj * w1 / FOUR;

    /* Use the ycur array for temporary storage here */
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h * cjm1, ark_mem->tempv2,
                          ark_mem->ycur, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    cvals[0] = mus * ark_mem->h;
    Xvecs[0] = ark_mem->ycur;
    cvals[1] = nu;
    Xvecs[1] = ark_mem->tempv1;
    cvals[2] = ONE - mu - nu;
    Xvecs[2] = ark_mem->yn;
    cvals[3] = mu;
    Xvecs[3] = ark_mem->tempv2;
    cvals[4] = -mus * ajm1 * ark_mem->h;
    Xvecs[4] = ark_mem->fn;

    retval = N_VLinearCombination(5, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    /* Shift the data for the next stage */
    if (j < step_mem->reqstages)
    {
      N_VScale(ONE, ark_mem->tempv2, ark_mem->tempv1);
      N_VScale(ONE, ark_mem->ycur, ark_mem->tempv2);

      cjm1 = cj;
      bjm2 = bjm1;
      bjm1 = bj;
    }
  }
  step_mem->nfe += step_mem->reqstages;

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur,
                          ark_mem->tempv2, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    /* Estimate the local error and compute its weighted RMS norm */
    cvals[0] = p8;
    Xvecs[0] = ark_mem->yn;
    cvals[1] = -p8;
    Xvecs[1] = ark_mem->ycur;
    cvals[2] = p4 * ark_mem->h;
    Xvecs[2] = ark_mem->fn;
    cvals[3] = p4 * ark_mem->h;
    Xvecs[3] = ark_mem->tempv2;

    retval = N_VLinearCombination(4, cvals, Xvecs, ark_mem->tempv1);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
    lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr);
  }
  return (ARK_SUCCESS);
}

// int lsrkStep_TakeStepRKG(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
// {
//   printf("\nRKG is not supported yet! Try RKC or RKL instead.\n");
//   return (ARK_ILL_INPUT);
// }

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs2:

  This routine serves the primary purpose of the LSRKStepSSPs2 module:
  it performs a single SSPs2 step (with embedding).

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
int lsrkStep_TakeStepSSPs2(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;

  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  sunrealtype rs     = (sunrealtype)step_mem->reqstages;
  sunrealtype sm1inv = ONE / (rs - ONE);
  sunrealtype bt1, bt2, bt3;

  if (step_mem->reqstages == 2)
  {
    bt1 = 0.694021459207626;
    bt3 = 1 - 0.694021459207626;
  }
  else
  {
    bt1 = (rs + ONE) / (rs * rs);
    bt2 = ONE / rs;
    bt3 = (rs - ONE) / (rs * rs);
  }

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current) && ark_mem->initsetup)
  {
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  /* A tentative solution at t+h is returned in
     y and its slope is evaluated in temp1.  */

  N_VLinearSum(ONE, ark_mem->yn, sm1inv * ark_mem->h, ark_mem->fn, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->yn, bt1 * ark_mem->h, ark_mem->fn, ark_mem->tempv1);

  /* Evaluate stages j = 2,...,step_mem->reqstages - 1 */
  for (int j = 2; j < step_mem->reqstages; j++)
  {
    retval =
      step_mem->fe(ark_mem->tcur + ((sunrealtype)j - ONE) * sm1inv * ark_mem->h,
                   ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, sm1inv * ark_mem->h, ark_mem->fn,
                 ark_mem->ycur);
    N_VLinearSum(ONE, ark_mem->tempv1, bt2 * ark_mem->h, ark_mem->fn,
                 ark_mem->tempv1);
  }
  /* Evaluate the last stage for j = step_mem->reqstages */
  retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur, ark_mem->fn,
                        ark_mem->user_data);
  if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

  N_VLinearSum(ONE / rs, ark_mem->yn, ONE / (sm1inv * rs), ark_mem->ycur,
               ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur, ark_mem->h / rs, ark_mem->fn, ark_mem->ycur);

  N_VLinearSum(ONE, ark_mem->tempv1, bt3 * ark_mem->h, ark_mem->fn,
               ark_mem->tempv1);

  step_mem->nfe += step_mem->reqstages - 1;

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);

    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }
  if (*dsmPtr <= ONE || ark_mem->fixedstep)
  {
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur,
                          ark_mem->fn, ark_mem->user_data);
    ark_mem->fn_is_current = SUNTRUE;
    step_mem->nfe++;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs3:

  This routine serves the primary purpose of the LSRKStepSSPs3 module:
  it performs a single SSPs3 step (with embedding).

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
int lsrkStep_TakeStepSSPs3(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;

  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  sunrealtype rs  = (sunrealtype)step_mem->reqstages;
  sunrealtype rn  = SUNRsqrt(rs);
  sunrealtype rat = ONE / (rs - rn);
  int in          = (int)rn;

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current) && ark_mem->initsetup)
  {
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  /* A tentative solution at t+h is returned in
     y and its slope is evaluated in temp1.  */

  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h * rat, ark_mem->fn, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h / rs, ark_mem->fn, ark_mem->tempv1);

  /* Evaluate stages j = 2,...,step_mem->reqstages */
  for (int j = 2; j <= (int)((in - 1) * (in - 2) / 2); j++)
  {
    retval =
      step_mem->fe(ark_mem->tcur + ((sunrealtype)j - ONE) * rat * ark_mem->h,
                   ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, ark_mem->h * rat, ark_mem->fn,
                 ark_mem->ycur);
    N_VLinearSum(ONE, ark_mem->tempv1, ark_mem->h / rs, ark_mem->fn,
                 ark_mem->tempv1);
  }
  N_VScale(ONE, ark_mem->ycur, ark_mem->tempv2);
  for (int j = (int)((in - 1) * (in - 2) / 2 + 1);
       j <= (int)(in * (in + 1) / 2 - 1); j++)
  {
    retval =
      step_mem->fe(ark_mem->tcur + ((sunrealtype)j - ONE) * rat * ark_mem->h,
                   ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, ark_mem->h * rat, ark_mem->fn,
                 ark_mem->ycur);
    N_VLinearSum(ONE, ark_mem->tempv1, ark_mem->h / rs, ark_mem->fn,
                 ark_mem->tempv1);
  }
  retval = step_mem->fe(ark_mem->tcur +
                          rat * (rn * (rn + ONE) / TWO - ONE) * ark_mem->h,
                        ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
  if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

  N_VLinearSum(rn / (TWO * rn - ONE), ark_mem->tempv2,
               (rn - ONE) / (TWO * rn - ONE), ark_mem->ycur, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur,
               (rn - ONE) * rat * ark_mem->h / (TWO * rn - ONE), ark_mem->fn,
               ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->tempv1, ark_mem->h / rs, ark_mem->fn,
               ark_mem->tempv1);

  for (int j = (int)(in * (in + 1) / 2 + 1); j <= step_mem->reqstages; j++)
  {
    retval = step_mem->fe(ark_mem->tcur +
                            ((sunrealtype)j - rn - ONE) * rat * ark_mem->h,
                          ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, ark_mem->h * rat, ark_mem->fn,
                 ark_mem->ycur);
    N_VLinearSum(ONE, ark_mem->tempv1, ark_mem->h / rs, ark_mem->fn,
                 ark_mem->tempv1);
  }

  step_mem->nfe += step_mem->reqstages - 1;

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);

    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }
  if (*dsmPtr <= ONE || ark_mem->fixedstep)
  {
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur,
                          ark_mem->fn, ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSP104:

  This routine serves the primary purpose of the LSRKStepSSP104 module:
  it performs a single SSP104 step (with embedding).

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
int lsrkStep_TakeStepSSP104(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;

  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  sunrealtype onesixth = 1.0 / 6.0;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Call the full RHS if needed. If this is the first step then we may need to
     evaluate or copy the RHS values from an  earlier evaluation (e.g., to
     compute h0). For subsequent steps treat this RHS evaluation as an
     evaluation at the end of the just completed step to potentially reuse
     (FSAL methods) RHS evaluations from the end of the last step. */

  if (!(ark_mem->fn_is_current) && ark_mem->initsetup)
  {
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  /* A tentative solution at t+h is returned in
     y and its slope is evaluated in temp1.  */

  N_VScale(ONE, ark_mem->yn, ark_mem->tempv2);

  N_VLinearSum(ONE, ark_mem->yn, onesixth * ark_mem->h, ark_mem->fn,
               ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->yn, 1.0 / 5.0 * ark_mem->h, ark_mem->fn,
               ark_mem->tempv1);

  /* Evaluate stages j = 2,...,step_mem->reqstages */
  for (int j = 2; j <= 5; j++)
  {
    retval = step_mem->fe(ark_mem->tcur +
                            ((sunrealtype)j - ONE) * onesixth * ark_mem->h,
                          ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, onesixth * ark_mem->h, ark_mem->fn,
                 ark_mem->ycur);
    if (j == 4)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, 3.0 / 10.0 * ark_mem->h, ark_mem->fn,
                   ark_mem->tempv1);
    }
  }
  N_VLinearSum(1.0 / 25.0, ark_mem->tempv2, 9.0 / 25.0, ark_mem->ycur,
               ark_mem->tempv2);
  N_VLinearSum(15, ark_mem->tempv2, -5, ark_mem->ycur, ark_mem->ycur);
  for (int j = 6; j <= 9; j++)
  {
    retval = step_mem->fe(ark_mem->tcur +
                            ((sunrealtype)j - 4.0) * onesixth * ark_mem->h,
                          ark_mem->ycur, ark_mem->fn, ark_mem->user_data);
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

    N_VLinearSum(ONE, ark_mem->ycur, onesixth * ark_mem->h, ark_mem->fn,
                 ark_mem->ycur);

    if (j == 7)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, 1.0 / 5.0 * ark_mem->h, ark_mem->fn,
                   ark_mem->tempv1);
    }
    if (j == 9)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, 3.0 / 10.0 * ark_mem->h, ark_mem->fn,
                   ark_mem->tempv1);
    }
  }

  retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur, ark_mem->fn,
                        ark_mem->user_data);
  if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }

  N_VLinearSum(ONE, ark_mem->tempv2, 3.0 / 5.0, ark_mem->ycur, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur, 1.0 / 10.0 * ark_mem->h, ark_mem->fn,
               ark_mem->ycur);

  step_mem->nfe += 9;

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);

    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }
  if (*dsmPtr <= ONE || ark_mem->fixedstep)
  {
    retval = step_mem->fe(ark_mem->tcur + ark_mem->h, ark_mem->ycur,
                          ark_mem->fn, ark_mem->user_data);
    step_mem->nfe++;
    ark_mem->fn_is_current = SUNTRUE;
    if (retval != ARK_SUCCESS) { return (ARK_RHSFUNC_FAIL); }
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_AccessARKODEStepMem:

  Shortcut routine to unpack both ark_mem and step_mem structures
  from void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int lsrkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeLSRKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  /* access ARKodeLSRKStepMem structure */
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_LSRKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeLSRKStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_AccessStepMem:

  Shortcut routine to unpack the step_mem structure from
  ark_mem.  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int lsrkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeLSRKStepMem* step_mem)
{
  /* access ARKodeLSRKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_LSRKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeLSRKStepMem)ark_mem->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype lsrkStep_CheckNVector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvscale == NULL) || (tmpl->ops->nvwrmsnorm == NULL) ||
      (tmpl->ops->nvspace == NULL))
  {
    return (SUNFALSE);
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  lsrkStep_ComputeNewDomEig:

  This routine computes new DomEig and returns SUN_SUCCESS.
  ---------------------------------------------------------------*/

int lsrkStep_ComputeNewDomEig(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem)
{
  int retval = SUN_SUCCESS;

  if ((step_mem->isextDomEig))
  {
    retval = step_mem->extDomEig(ark_mem->tn, ark_mem->ycur, &step_mem->lambdaR,
                                 &step_mem->lambdaI, ark_mem->user_data);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to estimate the dominant eigenvalue");
      return (ARK_DOMEIG_FAIL);
    }

    if (step_mem->lambdaR * ark_mem->h > ZERO)
    {
      printf("\nlambdaR*h must be nonpositive\n");
      return (ARK_ILL_INPUT);
    }
    else if (step_mem->lambdaR == 0 && SUNRabs(step_mem->lambdaI) > 0)
    {
      printf("\nDomEig cannot be purely imaginary\n");
      return (ARK_ILL_INPUT);
    }

    step_mem->lambdaR *= step_mem->domeigsfty;
    step_mem->lambdaI *= step_mem->domeigsfty;
    step_mem->sprad =
      SUNRsqrt(SUNSQR(step_mem->lambdaR) + SUNSQR(step_mem->lambdaI));
  }
  else
  {
    printf("\nInternal DomEig is not supported yet!");
    printf(
      "\nCall LSRKStepSetDomEigFn to provide an external DomEig function\n");

    return (ARK_ILL_INPUT);
  }
  step_mem->jacatt = SUNTRUE;

  step_mem->sprmax = (step_mem->sprad > step_mem->sprmax) ? step_mem->sprad
                                                          : step_mem->sprmax;

  if (step_mem->sprad < step_mem->sprmin || ark_mem->nst == 0)
  {
    step_mem->sprmin = step_mem->sprad;
  }

  step_mem->newdomeig = SUNFALSE;

  return retval;
}

/*---------------------------------------------------------------
  lsrkStep_DomEigUpdateLogic:

  This routine checks if the step is accepted or not and reassigns
  the DomEig update flags accordingly.
  ---------------------------------------------------------------*/

void lsrkStep_DomEigUpdateLogic(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem,
                                sunrealtype dsm)
{
  if (dsm <= ONE || ark_mem->fixedstep)
  {
    N_VScale(ONE, ark_mem->tempv2, ark_mem->fn);
    ark_mem->fn_is_current = SUNTRUE;

    step_mem->jacatt    = (step_mem->constJac == SUNTRUE);
    step_mem->nstsig    = (step_mem->nstsig + 1) % step_mem->domeigfreq;
    step_mem->newdomeig = SUNFALSE;
    if (step_mem->nstsig == 0) { step_mem->newdomeig = !step_mem->jacatt; }
  }
  else { step_mem->newdomeig = !step_mem->jacatt; }
}

/*===============================================================
  EOF
  ===============================================================*/
