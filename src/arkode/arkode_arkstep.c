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
 * This is the implementation file for ARKODE's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "arkode/arkode_butcher.h"
#include "arkode_arkstep_impl.h"
#include "arkode_impl.h"
#include "arkode_interp_impl.h"

#define FIXED_LIN_TOL

/*===============================================================
  ARKStep Exported functions -- Required
  ===============================================================*/

void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0,
                    SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunbooleantype nvectorOK;
  int retval;

  /* Check that at least one of fe, fi is supplied and is to be used */
  if (fe == NULL && fi == NULL)
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
  nvectorOK = arkStep_CheckNVector(y0);
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

  /* Allocate ARKodeARKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeARKStepMem)malloc(sizeof(struct ARKodeARKStepMemRec));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeARKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol   = arkStep_AttachLinsol;
  ark_mem->step_attachmasssol  = arkStep_AttachMasssol;
  ark_mem->step_disablelsetup  = arkStep_DisableLSetup;
  ark_mem->step_disablemsetup  = arkStep_DisableMSetup;
  ark_mem->step_getlinmem      = arkStep_GetLmem;
  ark_mem->step_getmassmem     = arkStep_GetMassMem;
  ark_mem->step_getimplicitrhs = arkStep_GetImplicitRHS;
  ark_mem->step_mmult          = NULL;
  ark_mem->step_getgammas      = arkStep_GetGammas;
  ark_mem->step_init           = arkStep_Init;
  ark_mem->step_fullrhs        = arkStep_FullRHS;
  ark_mem->step                = arkStep_TakeStep_Z;
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

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Allocate the general ARK stepper vectors using y0 as a template */
  /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on
     (based on the number of ARK stages) */

  /* Clone the input vector to create sdata, zpred and zcor */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata)))
  {
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zpred)))
  {
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zcor)))
  {
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Update the ARKODE workspace requirements */
  ark_mem->liw += 41; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
  ark_mem->lrw += 10;

  /* If an implicit component is to be solved, create default Newton NLS object */
  step_mem->ownNLS = SUNFALSE;
  if (step_mem->implicit)
  {
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (NLS == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error creating default Newton solver");
      ARKStepFree((void**)&ark_mem);
      return (NULL);
    }
    retval = ARKStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error attaching default Newton solver");
      ARKStepFree((void**)&ark_mem);
      return (NULL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit       = NULL;
  step_mem->lsetup      = NULL;
  step_mem->lsolve      = NULL;
  step_mem->lfree       = NULL;
  step_mem->lmem        = NULL;
  step_mem->lsolve_type = -1;

  /* Set the mass matrix solver addresses to NULL */
  step_mem->minit       = NULL;
  step_mem->msetup      = NULL;
  step_mem->mmult       = NULL;
  step_mem->msolve      = NULL;
  step_mem->mfree       = NULL;
  step_mem->mass_mem    = NULL;
  step_mem->mass_type   = MASS_IDENTITY;
  step_mem->msolve_type = -1;

  /* Initialize initial error norm  */
  step_mem->eRNrm = ONE;

  /* Initialize all the counters */
  step_mem->nfe       = 0;
  step_mem->nfi       = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;
  step_mem->nls_fails = 0;

  /* Initialize fused op work space */
  step_mem->cvals        = NULL;
  step_mem->Xvecs        = NULL;
  step_mem->nfusedopvecs = 0;

  /* Initialize external polynomial forcing data */
  step_mem->expforcing = SUNFALSE;
  step_mem->impforcing = SUNFALSE;
  step_mem->forcing    = NULL;
  step_mem->nforcing   = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKStepFree((void**)&ark_mem);
    return (NULL);
  }

  return ((void*)ark_mem);
}

/*---------------------------------------------------------------
  ARKStepResize:

  This routine resizes the memory within the ARKStep module.
  It first resizes the main ARKODE infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int ARKStepResize(void* arkode_mem, N_Vector y0, sunrealtype hscale,
                  sunrealtype t0, ARKVecResizeFn resize, void* resize_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKODE infrastructure memory */
  retval = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to resize main ARKODE infrastructure");
    return (retval);
  }

  /* Resize the sdata, zpred and zcor vectors */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                    &step_mem->sdata))
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                    &step_mem->zpred))
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                    &step_mem->zcor))
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  /* Resize the ARKStep vectors */
  /*     Fe */
  if (step_mem->Fe != NULL)
  {
    for (i = 0; i < step_mem->stages; i++)
    {
      if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                        &step_mem->Fe[i]))
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to resize vector");
        return (ARK_MEM_FAIL);
      }
    }
  }
  /*     Fi */
  if (step_mem->Fi != NULL)
  {
    for (i = 0; i < step_mem->stages; i++)
    {
      if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                        &step_mem->Fi[i]))
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to resize vector");
        return (ARK_MEM_FAIL);
      }
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

    /* attach new Newton NLS object to ARKStep */
    retval = ARKStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Error attaching default Newton solver");
      return (ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* reset nonlinear solver counters */
  if (step_mem->NLS != NULL) { step_mem->nsetups = 0; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepReInit:

  This routine re-initializes the ARKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int ARKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                  N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that at least one of fe, fi is supplied and is to be used */
  if (fe == NULL && fi == NULL)
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
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Copy the input parameters into ARKODE state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Initialize initial error norm  */
  step_mem->eRNrm = ONE;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to reinitialize main ARKODE infrastructure");
    return (retval);
  }

  /* Initialize all the counters */
  step_mem->nfe     = 0;
  step_mem->nfi     = 0;
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepReset:

  This routine resets the ARKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int ARKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return (retval);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSStolerances, ARKStepSVtolerances, ARKStepWFtolerances,
  ARKStepResStolerance, ARKStepResVtolerance, ARKStepResFtolerance:

  These routines set integration tolerances (wrappers for general
  ARKODE utility routines)
  ---------------------------------------------------------------*/
int ARKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkSStolerances(ark_mem, reltol, abstol));
}

int ARKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkSVtolerances(ark_mem, reltol, abstol));
}

int ARKStepWFtolerances(void* arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkWFtolerances(ark_mem, efun));
}

int ARKStepResStolerance(void* arkode_mem, sunrealtype rabstol)
{
  /* unpack ark_mem, call arkResStolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkResStolerance(ark_mem, rabstol));
}

int ARKStepResVtolerance(void* arkode_mem, N_Vector rabstol)
{
  /* unpack ark_mem, call arkResVtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkResVtolerance(ark_mem, rabstol));
}

int ARKStepResFtolerance(void* arkode_mem, ARKRwtFn rfun)
{
  /* unpack ark_mem, call arkResFtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkResFtolerance(ark_mem, rfun));
}

/*---------------------------------------------------------------
  ARKStepRootInit:

  Initialize (attach) a rootfinding problem to the stepper
  (wrappers for general ARKODE utility routine)
  ---------------------------------------------------------------*/
int ARKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  return (arkRootInit(ark_mem, nrtfn, g));
}

/*---------------------------------------------------------------
  ARKStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKODE utility routine)
  ---------------------------------------------------------------*/
int ARKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkEvolve(ark_mem, tout, yout, tret, itask);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return (retval);
}

/*---------------------------------------------------------------
  ARKStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKODE utility routine)
  ---------------------------------------------------------------*/
int ARKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkGetDky(ark_mem, t, k, dky);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return (retval);
}

/*---------------------------------------------------------------
  ARKStepComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
int ARKStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepFree frees all ARKStep memory, and then calls an ARKODE
  utility routine to free the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void ARKStepFree(void** arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL) { return; }

  /* conditional frees on non-NULL ARKStep module */
  ark_mem = (ARKodeMem)(*arkode_mem);
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

    /* free the Butcher tables */
    if (step_mem->Be != NULL)
    {
      ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->Be);
      step_mem->Be = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }
    if (step_mem->Bi != NULL)
    {
      ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->Bi);
      step_mem->Bi = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
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

    /* free the mass matrix solver memory */
    if (step_mem->mfree != NULL)
    {
      step_mem->mfree((void*)ark_mem);
      step_mem->mass_mem = NULL;
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
    if (step_mem->Fe != NULL)
    {
      for (j = 0; j < step_mem->stages; j++)
      {
        arkFreeVec(ark_mem, &step_mem->Fe[j]);
      }
      free(step_mem->Fe);
      step_mem->Fe = NULL;
      ark_mem->liw -= step_mem->stages;
    }
    if (step_mem->Fi != NULL)
    {
      for (j = 0; j < step_mem->stages; j++)
      {
        arkFreeVec(ark_mem, &step_mem->Fi[j]);
      }
      free(step_mem->Fi);
      step_mem->Fi = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free stage vectors */
    if (step_mem->z != NULL)
    {
      for (j = 0; j < step_mem->stages; j++)
      {
        arkFreeVec(ark_mem, &step_mem->z[j]);
      }
      free(step_mem->z);
      step_mem->z = NULL;
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

  /* free memory for overall ARKODE infrastructure */
  arkFree(arkode_mem);
}

/*---------------------------------------------------------------
  ARKStepPrintMem:

  This routine outputs the memory from the ARKStep structure and
  the main ARKODE infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void ARKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* if outfile==NULL, set it to stdout */
  if (outfile == NULL) { outfile = stdout; }

  /* output data from main ARKODE infrastructure */
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile, "ARKStep: q = %i\n", step_mem->q);
  fprintf(outfile, "ARKStep: p = %i\n", step_mem->p);
  fprintf(outfile, "ARKStep: istage = %i\n", step_mem->istage);
  fprintf(outfile, "ARKStep: stages = %i\n", step_mem->stages);
  fprintf(outfile, "ARKStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile, "ARKStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile, "ARKStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile, "ARKStep: lsolve_type = %i\n", step_mem->lsolve_type);
  fprintf(outfile, "ARKStep: msolve_type = %i\n", step_mem->msolve_type);
  fprintf(outfile, "ARKStep: convfail = %i\n", step_mem->convfail);

  /* output long integer quantities */
  fprintf(outfile, "ARKStep: nfe = %li\n", step_mem->nfe);
  fprintf(outfile, "ARKStep: nfi = %li\n", step_mem->nfi);
  fprintf(outfile, "ARKStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile, "ARKStep: nstlp = %li\n", step_mem->nstlp);

  /* output boolean quantities */
  fprintf(outfile, "ARKStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile, "ARKStep: user_linear_timedep = %i\n",
          step_mem->linear_timedep);
  fprintf(outfile, "ARKStep: user_explicit = %i\n", step_mem->explicit);
  fprintf(outfile, "ARKStep: user_implicit = %i\n", step_mem->implicit);
  fprintf(outfile, "ARKStep: jcur = %i\n", step_mem->jcur);

  /* output sunrealtype quantities */
  if (step_mem->Be != NULL)
  {
    fprintf(outfile, "ARKStep: explicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Be, outfile);
  }
  if (step_mem->Bi != NULL)
  {
    fprintf(outfile, "ARKStep: implicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Bi, outfile);
  }
  fprintf(outfile, "ARKStep: gamma = %" RSYM "\n", step_mem->gamma);
  fprintf(outfile, "ARKStep: gammap = %" RSYM "\n", step_mem->gammap);
  fprintf(outfile, "ARKStep: gamrat = %" RSYM "\n", step_mem->gamrat);
  fprintf(outfile, "ARKStep: crate = %" RSYM "\n", step_mem->crate);
  fprintf(outfile, "ARKStep: eRNrm = %" RSYM "\n", step_mem->eRNrm);
  fprintf(outfile, "ARKStep: nlscoef = %" RSYM "\n", step_mem->nlscoef);
  fprintf(outfile, "ARKStep: crdown = %" RSYM "\n", step_mem->crdown);
  fprintf(outfile, "ARKStep: rdiv = %" RSYM "\n", step_mem->rdiv);
  fprintf(outfile, "ARKStep: dgmax = %" RSYM "\n", step_mem->dgmax);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "ARKStep: sdata:\n");
  N_VPrintFile(step_mem->sdata, outfile);
  fprintf(outfile, "ARKStep: zpred:\n");
  N_VPrintFile(step_mem->zpred, outfile);
  fprintf(outfile, "ARKStep: zcor:\n");
  N_VPrintFile(step_mem->zcor, outfile);
  if (step_mem->Fe != NULL)
    for (i = 0; i < step_mem->stages; i++)
    {
      fprintf(outfile, "ARKStep: Fe[%i]:\n", i);
      N_VPrintFile(step_mem->Fe[i], outfile);
    }
  if (step_mem->Fi != NULL)
    for (i = 0; i < step_mem->stages; i++)
    {
      fprintf(outfile, "ARKStep: Fi[%i]:\n", i);
      N_VPrintFile(step_mem->Fi[i], outfile);
    }
#endif
}

/*===============================================================
  ARKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKODE
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  ARKStep module.
  ---------------------------------------------------------------*/
int arkStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup, ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type, void* lmem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* free any existing system solver */
  if (step_mem->lfree != NULL) { step_mem->lfree(arkode_mem); }

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit       = linit;
  step_mem->lsetup      = lsetup;
  step_mem->lsolve      = lsolve;
  step_mem->lfree       = lfree;
  step_mem->lmem        = lmem;
  step_mem->lsolve_type = lsolve_type;

  /* Reset all linear solver counters */
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_AttachMasssol:

  This routine attaches the set of mass matrix linear solver
  interface routines, data structure, and solver type to the
  ARKStep module.
  ---------------------------------------------------------------*/
int arkStep_AttachMasssol(void* arkode_mem, ARKMassInitFn minit,
                          ARKMassSetupFn msetup, ARKMassMultFn mmult,
                          ARKMassSolveFn msolve, ARKMassFreeFn mfree,
                          sunbooleantype time_dep,
                          SUNLinearSolver_Type msolve_type, void* mass_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* free any existing mass matrix solver */
  if (step_mem->mfree != NULL) { step_mem->mfree(arkode_mem); }

  /* Attach the provided routines, data structure and solve type */
  step_mem->minit       = minit;
  step_mem->msetup      = msetup;
  step_mem->mmult       = mmult;
  step_mem->msolve      = msolve;
  step_mem->mfree       = mfree;
  step_mem->mass_mem    = mass_mem;
  step_mem->mass_type   = (time_dep) ? MASS_TIMEDEP : MASS_FIXED;
  step_mem->msolve_type = msolve_type;

  /* Attach mmult function pointer to ark_mem as well */
  ark_mem->step_mmult = mmult;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
void arkStep_DisableLSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem == NULL) { return; }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->step_mem == NULL) { return; }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* nullify the lsetup function pointer */
  step_mem->lsetup = NULL;
}

/*---------------------------------------------------------------
  arkStep_DisableMSetup:

  This routine NULLifies the msetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
void arkStep_DisableMSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem == NULL) { return; }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->step_mem == NULL) { return; }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* nullify the msetup function pointer */
  step_mem->msetup = NULL;
}

/*---------------------------------------------------------------
  arkStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.
  ---------------------------------------------------------------*/
void* arkStep_GetLmem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return lmem */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (NULL); }
  return (step_mem->lmem);
}

/*---------------------------------------------------------------
  arkStep_GetMassMem:

  This routine returns the mass matrix solver interface memory
  structure, mass_mem.
  ---------------------------------------------------------------*/
void* arkStep_GetMassMem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return mass_mem */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (NULL); }
  return (step_mem->mass_mem);
}

/*---------------------------------------------------------------
  arkStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
ARKRhsFn arkStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return fi */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (NULL); }
  return (step_mem->fi);
}

/*---------------------------------------------------------------
  arkStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
int arkStep_GetGammas(void* arkode_mem, sunrealtype* gamma, sunrealtype* gamrat,
                      sunbooleantype** jcur, sunbooleantype* dgamma_fail)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set outputs */
  step_mem     = (ARKodeARKStepMem)ark_mem->step_mem;
  *gamma       = step_mem->gamma;
  *gamrat      = step_mem->gamrat;
  *jcur        = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  For all initialization types, this routine sets the relevant
  TakeStep routine based on the current problem configuration.

  With initialization type FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK stages,
    method order, or solver options
  - checks for consistency between the system and mass matrix
    linear solvers (if applicable)
  - initializes and sets up the system and mass matrix linear
    solvers (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)
  - allocates the interpolation data structure (if needed based
    on ARKStep solver options)
  - updates the call_fullrhs flag if necessary

  With initialization type FIRST_INIT or RESIZE_INIT, this routine:
  - sets the relevant TakeStep routine based on the current
    problem configuration
  - checks for consistency between the system and mass matrix
    linear solvers (if applicable)
  - initializes and sets up the system and mass matrix linear
    solvers (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
int arkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int j, retval;
  sunbooleantype reset_efun;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return (ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    /* enforce use of arkEwtSmallReal if using a fixed step size for
       an explicit method, an internal error weight function, and not
       using an iterative mass matrix solver with rwt=ewt */
    reset_efun = SUNTRUE;
    if (step_mem->implicit) { reset_efun = SUNFALSE; }
    if (!ark_mem->fixedstep) { reset_efun = SUNFALSE; }
    if (ark_mem->user_efun) { reset_efun = SUNFALSE; }
    if (ark_mem->rwt_is_ewt &&
        (step_mem->msolve_type == SUNLINEARSOLVER_ITERATIVE))
    {
      reset_efun = SUNFALSE;
    }
    if (ark_mem->rwt_is_ewt &&
        (step_mem->msolve_type == SUNLINEARSOLVER_MATRIX_ITERATIVE))
    {
      reset_efun = SUNFALSE;
    }
    if (reset_efun)
    {
      ark_mem->user_efun = SUNFALSE;
      ark_mem->efun      = arkEwtSetSmallReal;
      ark_mem->e_data    = ark_mem;
    }

    /* Create Butcher tables (if not already set) */
    retval = arkStep_SetButcherTables(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Could not create Butcher table(s)");
      return (ARK_ILL_INPUT);
    }

    /* Check that Butcher tables are OK */
    retval = arkStep_CheckButcherTables(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in Butcher table(s)");
      return (ARK_ILL_INPUT);
    }

    /* Retrieve/store method and embedding orders now that tables are finalized */
    if (step_mem->Bi != NULL)
    {
      step_mem->q = ark_mem->hadapt_mem->q = step_mem->Bi->q;
      step_mem->p = ark_mem->hadapt_mem->p = step_mem->Bi->p;
    }
    else
    {
      step_mem->q = ark_mem->hadapt_mem->q = step_mem->Be->q;
      step_mem->p = ark_mem->hadapt_mem->p = step_mem->Be->p;
    }

    /* Ensure that if adaptivity is enabled, then method includes embedding coefficients */
    if (!ark_mem->fixedstep && (step_mem->p == 0))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Adaptive timestepping cannot be performed without embedding coefficients");
      return (ARK_ILL_INPUT);
    }

    /* Relaxation is incompatible with implicit RHS deduction */
    if (ark_mem->relax_enabled && step_mem->implicit && step_mem->deduce_rhs)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Relaxation cannot be performed when deducing implicit RHS values");
      return ARK_ILL_INPUT;
    }

    /* Allocate ARK RHS vector memory, update storage requirements */
    /*   Allocate Fe[0] ... Fe[stages-1] if needed */
    if (step_mem->explicit)
    {
      if (step_mem->Fe == NULL)
      {
        step_mem->Fe = (N_Vector*)calloc(step_mem->stages, sizeof(N_Vector));
      }
      for (j = 0; j < step_mem->stages; j++)
      {
        if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fe[j])))
        {
          return (ARK_MEM_FAIL);
        }
      }
      ark_mem->liw += step_mem->stages; /* pointers */
    }

    /*   Allocate Fi[0] ... Fi[stages-1] if needed */
    if (step_mem->implicit)
    {
      if (step_mem->Fi == NULL)
      {
        step_mem->Fi = (N_Vector*)calloc(step_mem->stages, sizeof(N_Vector));
      }
      for (j = 0; j < step_mem->stages; j++)
      {
        if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fi[j])))
        {
          return (ARK_MEM_FAIL);
        }
      }
      ark_mem->liw += step_mem->stages; /* pointers */
    }

    /* Allocate stage storage for relaxation with implicit/IMEX methods or if a
       fixed mass matrix is present (since we store f(t,y) not M^{-1} f(t,y)) */
    if (ark_mem->relax_enabled &&
        (step_mem->implicit || step_mem->mass_type == MASS_FIXED))
    {
      if (step_mem->z == NULL)
      {
        step_mem->z = (N_Vector*)calloc(step_mem->stages, sizeof(N_Vector));
      }
      for (j = 0; j < step_mem->stages; j++)
      {
        if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->z[j])))
        {
          return (ARK_MEM_FAIL);
        }
      }
      ark_mem->liw += step_mem->stages; /* pointers */
    }

    /* Allocate reusable arrays for fused vector operations */
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
      step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs,
                                          sizeof(N_Vector));
      if (step_mem->Xvecs == NULL) { return (ARK_MEM_FAIL); }
      ark_mem->liw += step_mem->nfusedopvecs; /* pointers */
    }

    /* Allocate workspace for MRI forcing -- need to allocate here as the
       number of stages may not bet set before this point and we assume
       SetInnerForcing has been called before the first step i.e., methods
       start with a fast integration */
    if (step_mem->expforcing || step_mem->impforcing)
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

  /* set appropriate TakeStep routine based on problem configuration */
  /*    (only one choice for now) */
  ark_mem->step = arkStep_TakeStep_Z;

  /* Check for consistency between mass system and system linear system modules
     (e.g., if lsolve is direct, msolve needs to match) */
  if ((step_mem->mass_type != MASS_IDENTITY) && step_mem->lmem)
  {
    if (step_mem->lsolve_type != step_mem->msolve_type)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Incompatible linear and mass matrix solvers");
      return (ARK_ILL_INPUT);
    }
  }

  /* Perform mass matrix solver initialization and setup (if applicable) */
  if (step_mem->mass_type != MASS_IDENTITY)
  {
    /* Call minit (if it exists) */
    if (step_mem->minit != NULL)
    {
      retval = step_mem->minit((void*)ark_mem);
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_MASSINIT_FAIL, __LINE__, __func__,
                        __FILE__, MSG_ARK_MASSINIT_FAIL);
        return (ARK_MASSINIT_FAIL);
      }
    }

    /* Call msetup (if it exists) */
    if (step_mem->msetup != NULL)
    {
      retval = step_mem->msetup((void*)ark_mem, ark_mem->tcur, ark_mem->tempv1,
                                ark_mem->tempv2, ark_mem->tempv3);
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, __LINE__, __func__,
                        __FILE__, MSG_ARK_MASSSETUP_FAIL);
        return (ARK_MASSSETUP_FAIL);
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
    retval = arkStep_NlsInit(ark_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to initialize SUNNonlinearSolver object");
      return (ARK_NLS_INIT_FAIL);
    }
  }

  /* Signal to shared arkode module that full RHS evaluations are required */
  ark_mem->call_fullrhs = SUNTRUE;

  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  arkStep_FullRHS:

  Rewriting the problem
    My' = fe(t,y) + fi(t,y)
  in the form
    y' = M^{-1}*[ fe(t,y) + fi(t,y) ],
  this routine computes the full right-hand side vector,
    f = M^{-1}*[ fe(t,y) + fi(t,y) ]

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If this function is called in ARK_FULLRHS_START or ARK_FULLRHS_END mode and
  evaluating the RHS functions is necessary, we store the vectors fe(t,y) and
  fi(t,y) in Fe[0] and Fi[0] for possible reuse in the first stage of the
  subsequent time step.

  In ARK_FULLRHS_END mode we check if the method is stiffly accurate and, if
  appropriate, copy the vectors Fe[stages - 1] and Fi[stages - 1] to Fe[0] and
  Fi[0] for possible reuse in the first stage of the subsequent time step.

  ARK_FULLRHS_OTHER mode is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to store the
  intermediate parts so that they do not interfere with the other two modes.
  ----------------------------------------------------------------------------*/
int arkStep_FullRHS(void* arkode_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int nvec, retval;
  sunbooleantype recomputeRHS;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  sunrealtype stage_coefs = ONE;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* setup mass-matrix if required (use output f as a temporary) */
  if ((step_mem->mass_type == MASS_TIMEDEP) && (step_mem->msetup != NULL))
  {
    retval = step_mem->msetup((void*)ark_mem, t, f, ark_mem->tempv2,
                              ark_mem->tempv3);
    if (retval != ARK_SUCCESS) { return (ARK_MASSSETUP_FAIL); }
  }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:

    /* compute the full RHS */
    if (!(ark_mem->fn_is_current))
    {
      /* compute the explicit component */
      if (step_mem->explicit)
      {
        retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
        step_mem->nfe++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
          return (ARK_RHSFUNC_FAIL);
        }

        /* compute and store M(t)^{-1} fe */
        if (step_mem->mass_type == MASS_TIMEDEP)
        {
          retval = step_mem->msolve((void*)ark_mem, step_mem->Fe[0],
                                    step_mem->nlscoef / ark_mem->h);
          if (retval)
          {
            arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                            __FILE__, "Mass matrix solver failure");
            return ARK_MASSSOLVE_FAIL;
          }
        }
      }

      /* compute the implicit component */
      if (step_mem->implicit)
      {
        retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
        step_mem->nfi++;
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
          return (ARK_RHSFUNC_FAIL);
        }

        /* compute and store M(t)^{-1} fi */
        if (step_mem->mass_type == MASS_TIMEDEP)
        {
          retval = step_mem->msolve((void*)ark_mem, step_mem->Fi[0],
                                    step_mem->nlscoef / ark_mem->h);
          if (retval)
          {
            arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                            __FILE__, "Mass matrix solver failure");
            return ARK_MASSSOLVE_FAIL;
          }
        }
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit)
    {
      /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    }
    else if (step_mem->implicit)
    {
      /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    }
    else
    {
      /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    /* compute M^{-1} f for output but do not store */
    if (step_mem->mass_type == MASS_FIXED)
    {
      retval = step_mem->msolve((void*)ark_mem, f,
                                step_mem->nlscoef / ark_mem->h);
      if (retval)
      {
        arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                        __FILE__, "Mass matrix solver failure");
        return ARK_MASSSOLVE_FAIL;
      }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if (step_mem->expforcing || step_mem->impforcing)
    {
      cvals[0] = ONE;
      Xvecs[0] = f;
      nvec     = 1;
      arkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
      N_VLinearCombination(nvec, cvals, Xvecs, f);
    }

    break;

  case ARK_FULLRHS_END:

    /* compute the full RHS */
    if (!(ark_mem->fn_is_current))
    {
      /* determine if RHS functions need to be recomputed */
      recomputeRHS = SUNFALSE;

      if (step_mem->explicit)
      {
        if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Be))
        {
          recomputeRHS = SUNTRUE;
        }
      }

      if (step_mem->implicit)
      {
        if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Bi))
        {
          recomputeRHS = SUNTRUE;
        }
      }

      /* Stiffly Accurate methods are not SA when relaxation is enabled */
      if (ark_mem->relax_enabled) { recomputeRHS = SUNTRUE; }

      /* recompute RHS functions */
      if (recomputeRHS)
      {
        /* compute the explicit component */
        if (step_mem->explicit)
        {
          retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
          step_mem->nfe++;
          if (retval != 0)
          {
            arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                            __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
            return (ARK_RHSFUNC_FAIL);
          }

          /* compute and store M(t)^{-1} fi */
          if (step_mem->mass_type == MASS_TIMEDEP)
          {
            retval = step_mem->msolve((void*)ark_mem, step_mem->Fe[0],
                                      step_mem->nlscoef / ark_mem->h);
            if (retval)
            {
              arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                              __FILE__, "Mass matrix solver failure");
              return ARK_MASSSOLVE_FAIL;
            }
          }
        }

        /* compute the implicit component */
        if (step_mem->implicit)
        {
          retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
          step_mem->nfi++;
          if (retval != 0)
          {
            arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                            __FILE__, MSG_ARK_RHSFUNC_FAILED, t);
            return (ARK_RHSFUNC_FAIL);
          }

          /* compute and store M(t)^{-1} fi */
          if (step_mem->mass_type == MASS_TIMEDEP)
          {
            retval = step_mem->msolve((void*)ark_mem, step_mem->Fi[0],
                                      step_mem->nlscoef / ark_mem->h);
            if (retval)
            {
              arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                              __FILE__, "Mass matrix solver failure");
              return ARK_MASSSOLVE_FAIL;
            }
          }
        }
      }
      else
      {
        if (step_mem->explicit)
        {
          N_VScale(ONE, step_mem->Fe[step_mem->stages - 1], step_mem->Fe[0]);
        }
        if (step_mem->implicit)
        {
          N_VScale(ONE, step_mem->Fi[step_mem->stages - 1], step_mem->Fi[0]);
        }
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit)
    {
      /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    }
    else if (step_mem->implicit)
    {
      /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    }
    else
    {
      /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    /* compute M^{-1} f for output but do not store */
    if (step_mem->mass_type == MASS_FIXED)
    {
      retval = step_mem->msolve((void*)ark_mem, f,
                                step_mem->nlscoef / ark_mem->h);
      if (retval)
      {
        arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                        __FILE__, "Mass matrix solver failure");
        return ARK_MASSSOLVE_FAIL;
      }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if (step_mem->expforcing || step_mem->impforcing)
    {
      cvals[0] = ONE;
      Xvecs[0] = f;
      nvec     = 1;
      arkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
      N_VLinearCombination(nvec, cvals, Xvecs, f);
    }

    break;

  case ARK_FULLRHS_OTHER:

    /* compute the explicit component and store in ark_tempv2 */
    if (step_mem->explicit)
    {
      retval = step_mem->fe(t, y, ark_mem->tempv2, ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* compute the implicit component and store in sdata */
    if (step_mem->implicit)
    {
      retval = step_mem->fi(t, y, step_mem->sdata, ark_mem->user_data);
      step_mem->nfi++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return (ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit)
    { /* ImEx */
      N_VLinearSum(ONE, step_mem->sdata, ONE, ark_mem->tempv2, f);
    }
    else if (step_mem->implicit)
    { /* implicit */
      N_VScale(ONE, step_mem->sdata, f);
    }
    else
    { /* explicit */
      N_VScale(ONE, ark_mem->tempv2, f);
    }

    /* compute M^{-1} f for output but do not store */
    if (step_mem->mass_type != MASS_IDENTITY)
    {
      retval = step_mem->msolve((void*)ark_mem, f,
                                step_mem->nlscoef / ark_mem->h);
      if (retval)
      {
        arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__,
                        __FILE__, "Mass matrix solver failure");
        return ARK_MASSSOLVE_FAIL;
      }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if (step_mem->expforcing || step_mem->impforcing)
    {
      cvals[0] = ONE;
      Xvecs[0] = f;
      nvec     = 1;
      arkStep_ApplyForcing(step_mem, &t, &stage_coefs, 1, &nvec);
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
  arkStep_TakeStep_Z:

  This routine serves the primary purpose of the ARKStep module:
  it performs a single ARK step (with embedding, if possible).
  This version solves for each ARK stage vector, z_i.

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int arkStep_TakeStep_Z(void* arkode_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval, is, is_start, mode;
  sunbooleantype implicit_stage;
  sunbooleantype deduce_stage;
  sunbooleantype save_stages;
  sunbooleantype stiffly_accurate;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  N_Vector zcor0;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if problem will involve no algebraic solvers, initialize nflagPtr to success */
  if ((!step_mem->implicit) && (step_mem->mass_type == MASS_IDENTITY))
  {
    *nflagPtr = ARK_SUCCESS;
  }

  /* call nonlinear solver setup if it exists */
  if (step_mem->NLS)
  {
    if ((step_mem->NLS)->ops->setup)
    {
      zcor0 = ark_mem->tempv3;
      N_VConst(ZERO,
               zcor0); /* set guess to all 0 (since using predictor-corrector form) */
      retval = SUNNonlinSolSetup(step_mem->NLS, zcor0, ark_mem);
      if (retval < 0) { return (ARK_NLS_SETUP_FAIL); }
      if (retval > 0) { return (ARK_NLS_SETUP_RECVR); }
    }
  }

  /* check if we need to store stage values */
  save_stages = SUNFALSE;
  if (ark_mem->relax_enabled &&
      (step_mem->implicit || step_mem->mass_type == MASS_FIXED))
  {
    save_stages = SUNTRUE;
  }

  /* check for implicit method with explicit first stage */
  implicit_stage = SUNFALSE;
  is_start       = 1;
  if (step_mem->implicit)
  {
    if (SUNRabs(step_mem->Bi->A[0][0]) > TINY)
    {
      implicit_stage = SUNTRUE;
      is_start       = 0;
    }
  }

  /* explicit first stage -- store stage if necessary for relaxation */
  if (is_start == 1 && save_stages)
  {
    N_VScale(ONE, ark_mem->yn, step_mem->z[0]);
  }

  /* check if the method is Stiffly Accurate (SA) */
  stiffly_accurate = SUNTRUE;
  if (step_mem->explicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Be))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  if (step_mem->implicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Bi))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  /* Call the full RHS if needed e.g., an explicit first stage. If this is the
     first step then we may need to evaluate or copy the RHS values from an
     earlier evaluation (e.g., to compute h0). For subsequent steps treat this
     RHS evaluation as an evaluation at the end of the just completed step to
     potentially reuse (FSAL methods) or save (stiffly accurate methods with an
     implicit first stage using Hermite interpolation) RHS evaluations from the
     end of the last step. */

  if ((!implicit_stage ||
       (stiffly_accurate && ark_mem->interp_type == ARK_INTERP_HERMITE)) &&
      !(ark_mem->fn_is_current))
  {
    mode   = (ark_mem->initsetup) ? ARK_FULLRHS_START : ARK_FULLRHS_END;
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                   ark_mem->fn, mode);
    if (retval) { return ARK_RHSFUNC_FAIL; }
    ark_mem->fn_is_current = SUNTRUE;
  }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  if (is_start == 1)
  {
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                       "ARKODE::arkStep_TakeStep_Z", "start-stage",
                       "step = %li, stage = %i, implicit = %i, h = %" RSYM
                       ", tcur = %" RSYM,
                       ark_mem->nst, 0, implicit_stage, ark_mem->h,
                       ark_mem->tcur);
#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::arkStep_TakeStep_Z", "explicit stage",
                       "z[%i] =", 0);
    N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
    if (step_mem->implicit)
    {
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "implicit RHS",
                         "Fi[%i] =", 0);
      N_VPrintFile(step_mem->Fi[0], ARK_LOGGER->debug_fp);
    }
    if (step_mem->explicit)
    {
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "explicit RHS",
                         "Fe[%i] =", 0);
      N_VPrintFile(step_mem->Fe[0], ARK_LOGGER->debug_fp);
    }
#endif
  }
#endif

  /* loop over internal stages to the step */
  for (is = is_start; is < step_mem->stages; is++)
  {
    /* store current stage index */
    step_mem->istage = is;

    /* determine whether implicit solve is required */
    implicit_stage = SUNFALSE;
    if (step_mem->implicit)
    {
      if (SUNRabs(step_mem->Bi->A[is][is]) > TINY) { implicit_stage = SUNTRUE; }
    }

    /* determine if the stage RHS will be deduced from the implicit solve */
    deduce_stage = step_mem->deduce_rhs && implicit_stage;

    /* set current stage time(s) */
    if (step_mem->implicit)
    {
      ark_mem->tcur = ark_mem->tn + step_mem->Bi->c[is] * ark_mem->h;
    }
    else { ark_mem->tcur = ark_mem->tn + step_mem->Be->c[is] * ark_mem->h; }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::arkStep_TakeStep_Z", "start-stage",
                       "step = %li, stage = %i, implicit = %i, h = %" RSYM
                       ", tcur = %" RSYM,
                       ark_mem->nst, is, implicit_stage, ark_mem->h,
                       ark_mem->tcur);
#endif

    /* setup time-dependent mass matrix */
    if ((step_mem->mass_type == MASS_TIMEDEP) && (step_mem->msetup != NULL))
    {
      retval = step_mem->msetup((void*)ark_mem, ark_mem->tcur, ark_mem->tempv1,
                                ark_mem->tempv2, ark_mem->tempv3);
      if (retval != ARK_SUCCESS) { return (ARK_MASSSETUP_FAIL); }
    }

    /* if implicit, call built-in and user-supplied predictors
       (results placed in zpred) */
    if (implicit_stage)
    {
      retval = arkStep_Predict(ark_mem, is, step_mem->zpred);
      if (retval != ARK_SUCCESS) { return (retval); }

      /* if a user-supplied predictor routine is provided, call that here.
         Note that arkStep_Predict is *still* called, so this user-supplied
         routine can just 'clean up' the built-in prediction, if desired. */
      if (step_mem->stage_predict)
      {
        retval = step_mem->stage_predict(ark_mem->tcur, step_mem->zpred,
                                         ark_mem->user_data);
        if (retval < 0) { return (ARK_USER_PREDICT_FAIL); }
        if (retval > 0) { return (TRY_AGAIN); }
      }
    }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::arkStep_TakeStep_Z", "predictor", "zpred =", "");
    N_VPrintFile(step_mem->zpred, ARK_LOGGER->debug_fp);
#endif

    /* set up explicit data for evaluation of ARK stage (store in sdata) */
    retval = arkStep_StageSetup(ark_mem, implicit_stage);
    if (retval != ARK_SUCCESS) { return (retval); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
    SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                       "ARKODE::arkStep_TakeStep_Z", "rhs data", "sdata =", "");
    N_VPrintFile(step_mem->sdata, ARK_LOGGER->debug_fp);
#endif

    /* perform implicit solve if required */
    if (implicit_stage)
    {
      /* implicit solve result is stored in ark_mem->ycur;
         return with positive value on anything but success */
      *nflagPtr = arkStep_Nls(ark_mem, *nflagPtr);
      if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "implicit stage",
                         "z[%i] =", is);
      N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

      /* otherwise no implicit solve is needed */
    }
    else
    {
      /* if M is fixed, solve with it to compute update (place back in sdata) */
      if (step_mem->mass_type == MASS_FIXED)
      {
        /* perform solve; return with positive value on anything but success */
        *nflagPtr = step_mem->msolve((void*)ark_mem, step_mem->sdata,
                                     step_mem->nlscoef);
        if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }
      }

      /* set y to be yn + sdata (either computed in arkStep_StageSetup,
         or updated in prev. block) */
      N_VLinearSum(ONE, ark_mem->yn, ONE, step_mem->sdata, ark_mem->ycur);

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "explicit stage",
                         "z[%i] =", is);
      N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    /* NOTE: with internally inconsistent IMEX methods (c_i^E != c_i^I) the value
       of tcur corresponds to the stage time from the implicit table (c_i^I). */
    if (ark_mem->ProcessStage != NULL)
    {
      retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) { return (ARK_POSTPROCESS_STAGE_FAIL); }
    }

    /* successful stage solve */

    /*    store stage (if necessary for relaxation) */
    if (save_stages) { N_VScale(ONE, ark_mem->ycur, step_mem->z[is]); }

    /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
    if (step_mem->implicit)
    {
      if (!deduce_stage)
      {
        retval = step_mem->fi(ark_mem->tcur, ark_mem->ycur, step_mem->Fi[is],
                              ark_mem->user_data);
        step_mem->nfi++;
      }
      else if (step_mem->mass_type == MASS_FIXED)
      {
        retval = step_mem->mmult((void*)ark_mem, step_mem->zcor, ark_mem->tempv1);
        if (retval != ARK_SUCCESS) { return (ARK_MASSMULT_FAIL); }
        N_VLinearSum(ONE / step_mem->gamma, ark_mem->tempv1,
                     -ONE / step_mem->gamma, step_mem->sdata, step_mem->Fi[is]);
      }
      else
      {
        N_VLinearSum(ONE / step_mem->gamma, step_mem->zcor,
                     -ONE / step_mem->gamma, step_mem->sdata, step_mem->Fi[is]);
      }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "implicit RHS",
                         "Fi[%i] =", is);
      N_VPrintFile(step_mem->Fi[is], ARK_LOGGER->debug_fp);
#endif

      if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
      if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }
    }

    /*    store explicit RHS */
    if (step_mem->explicit)
    {
      retval = step_mem->fe(ark_mem->tn + step_mem->Be->c[is] * ark_mem->h,
                            ark_mem->ycur, step_mem->Fe[is], ark_mem->user_data);
      step_mem->nfe++;

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                         "ARKODE::arkStep_TakeStep_Z", "explicit RHS",
                         "Fe[%i] =", is);
      N_VPrintFile(step_mem->Fe[is], ARK_LOGGER->debug_fp);
#endif

      if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
      if (retval > 0) { return (ARK_UNREC_RHSFUNC_ERR); }
    }

    /* if using a time-dependent mass matrix, update Fe[is] and/or Fi[is] with M(t)^{-1} */
    if (step_mem->mass_type == MASS_TIMEDEP)
    {
      /* If the implicit stage was deduced, it already includes M(t)^{-1} */
      if (step_mem->implicit && !deduce_stage)
      {
        *nflagPtr = step_mem->msolve((void*)ark_mem, step_mem->Fi[is],
                                     step_mem->nlscoef);
#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::arkStep_TakeStep_Z", "M^{-1} implicit RHS",
                           "Fi[%i] =", is);
        N_VPrintFile(step_mem->Fi[is], ARK_LOGGER->debug_fp);
#endif
        if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }
      }
      if (step_mem->explicit)
      {
        *nflagPtr = step_mem->msolve((void*)ark_mem, step_mem->Fe[is],
                                     step_mem->nlscoef);
#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
        SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG,
                           "ARKODE::arkStep_TakeStep_Z", "M^{-1} explicit RHS",
                           "Fe[%i] =", is);
        N_VPrintFile(step_mem->Fe[is], ARK_LOGGER->debug_fp);
#endif
        if (*nflagPtr != ARK_SUCCESS) { return (TRY_AGAIN); }
      }
    }

  } /* loop over stages */

  /* compute time-evolved solution (in ark_ycur), error estimate (in dsm).
     This can fail recoverably due to nonconvergence of the mass matrix solve,
     so handle that appropriately. */
  if (step_mem->mass_type == MASS_FIXED)
  {
    *nflagPtr = arkStep_ComputeSolutions_MassFixed(ark_mem, dsmPtr);
  }
  else { *nflagPtr = arkStep_ComputeSolutions(ark_mem, dsmPtr); }
  if (*nflagPtr < 0) { return (*nflagPtr); }
  if (*nflagPtr > 0) { return (TRY_AGAIN); }

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG, "ARKODE::arkStep_TakeStep_Z",
                     "updated solution", "ycur =", "");
  N_VPrintFile(ark_mem->ycur, ARK_LOGGER->debug_fp);
#endif

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO,
                     "ARKODE::arkStep_TakeStep_Z", "end-step",
                     "step = %li, h = %" RSYM ", dsm = %" RSYM ", nflag = %d",
                     ark_mem->nst, ark_mem->h, *dsmPtr, *nflagPtr);
#endif

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int arkStep_AccessStepMem(void* arkode_mem, const char* fname,
                          ARKodeMem* ark_mem, ARKodeARKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeARKStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkStep_CheckNVector(N_Vector tmpl)
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
  arkStep_SetButcherTables

  This routine determines the ERK/DIRK/ARK method to use, based
  on the desired accuracy and information on whether the problem
  is explicit, implicit or imex.
  ---------------------------------------------------------------*/
int arkStep_SetButcherTables(ARKodeMem ark_mem)
{
  int etable, itable;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* if tables have already been specified, just return */
  if ((step_mem->Be != NULL) || (step_mem->Bi != NULL))
  {
    return (ARK_SUCCESS);
  }

  /* initialize table numbers to illegal values */
  etable = itable = -1;

  /**** ImEx methods ****/
  if (step_mem->explicit && step_mem->implicit)
  {
    switch (step_mem->q)
    {
    case (2):
      etable = ARKSTEP_DEFAULT_ARK_ETABLE_2;
      itable = ARKSTEP_DEFAULT_ARK_ITABLE_2;
      break;
    case (3):
      etable = ARKSTEP_DEFAULT_ARK_ETABLE_3;
      itable = ARKSTEP_DEFAULT_ARK_ITABLE_3;
      break;
    case (4):
      etable = ARKSTEP_DEFAULT_ARK_ETABLE_4;
      itable = ARKSTEP_DEFAULT_ARK_ITABLE_4;
      break;
    case (5):
      etable = ARKSTEP_DEFAULT_ARK_ETABLE_5;
      itable = ARKSTEP_DEFAULT_ARK_ITABLE_5;
      break;
    default: /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "No ImEx method at requested order, using q=5.");
      etable = ARKSTEP_DEFAULT_ARK_ETABLE_5;
      itable = ARKSTEP_DEFAULT_ARK_ITABLE_5;
      break;
    }

    /**** implicit methods ****/
  }
  else if (step_mem->implicit)
  {
    switch (step_mem->q)
    {
    case (2): itable = ARKSTEP_DEFAULT_DIRK_2; break;
    case (3): itable = ARKSTEP_DEFAULT_DIRK_3; break;
    case (4): itable = ARKSTEP_DEFAULT_DIRK_4; break;
    case (5): itable = ARKSTEP_DEFAULT_DIRK_5; break;
    default: /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "No implicit method at requested order, using q=5.");
      itable = ARKSTEP_DEFAULT_DIRK_5;
      break;
    }

    /**** explicit methods ****/
  }
  else
  {
    switch (step_mem->q)
    {
    case (2): etable = ARKSTEP_DEFAULT_ERK_2; break;
    case (3): etable = ARKSTEP_DEFAULT_ERK_3; break;
    case (4): etable = ARKSTEP_DEFAULT_ERK_4; break;
    case (5): etable = ARKSTEP_DEFAULT_ERK_5; break;
    case (6): etable = ARKSTEP_DEFAULT_ERK_6; break;
    case (7): etable = ARKSTEP_DEFAULT_ERK_7; break;
    case (8): etable = ARKSTEP_DEFAULT_ERK_8; break;
    case (9): etable = ARKSTEP_DEFAULT_ERK_9; break;
    default: /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "No explicit method at requested order, using q=9.");
      etable = ARKSTEP_DEFAULT_ERK_9;
      break;
    }
  }

  if (etable > -1) { step_mem->Be = ARKodeButcherTable_LoadERK(etable); }
  if (itable > -1) { step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable); }

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* set [redundant] ARK stored values for stage numbers and method orders */
  if (step_mem->Be != NULL)
  {
    step_mem->stages = step_mem->Be->stages;
    step_mem->q      = step_mem->Be->q;
    step_mem->p      = step_mem->Be->p;
  }
  if (step_mem->Bi != NULL)
  {
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q      = step_mem->Bi->q;
    step_mem->p      = step_mem->Bi->p;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_CheckButcherTables

  This routine runs through the explicit and/or implicit Butcher
  tables to ensure that they meet all necessary requirements,
  including:
    strictly lower-triangular (ERK)
    lower-triangular with some nonzeros on diagonal (IRK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int arkStep_CheckButcherTables(ARKodeMem ark_mem)
{
  int i, j;
  sunbooleantype okay;
  ARKodeARKStepMem step_mem;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* check that the expected tables are set */
  if (step_mem->explicit && step_mem->Be == NULL)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "explicit table is NULL!");
    return (ARK_INVALID_TABLE);
  }

  if (step_mem->implicit && step_mem->Bi == NULL)
  {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "implicit table is NULL!");
    return (ARK_INVALID_TABLE);
  }

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
    if (step_mem->implicit)
    {
      if (step_mem->Bi->d == NULL)
      {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__,
                        __FILE__, "no implicit embedding!");
        return (ARK_INVALID_TABLE);
      }
    }
    if (step_mem->explicit)
    {
      if (step_mem->Be->d == NULL)
      {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__,
                        __FILE__, "no explicit embedding!");
        return (ARK_INVALID_TABLE);
      }
    }
  }

  /* check that ERK table is strictly lower triangular */
  if (step_mem->explicit)
  {
    okay = SUNTRUE;
    for (i = 0; i < step_mem->stages; i++)
    {
      for (j = i; j < step_mem->stages; j++)
      {
        if (SUNRabs(step_mem->Be->A[i][j]) > tol) { okay = SUNFALSE; }
      }
    }
    if (!okay)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Ae Butcher table is implicit!");
      return (ARK_INVALID_TABLE);
    }
  }

  /* check that IRK table is implicit and lower triangular */
  if (step_mem->implicit)
  {
    okay = SUNFALSE;
    for (i = 0; i < step_mem->stages; i++)
    {
      if (SUNRabs(step_mem->Bi->A[i][i]) > tol) { okay = SUNTRUE; }
    }
    if (!okay)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Ai Butcher table is explicit!");
      return (ARK_INVALID_TABLE);
    }

    okay = SUNTRUE;
    for (i = 0; i < step_mem->stages; i++)
    {
      for (j = i + 1; j < step_mem->stages; j++)
      {
        if (SUNRabs(step_mem->Bi->A[i][j]) > tol) { okay = SUNFALSE; }
      }
    }
    if (!okay)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "Ai Butcher table has entries above diagonal!");
      return (ARK_INVALID_TABLE);
    }
  }

  /* Check if the method is compatible with relaxation */
  if (ark_mem->relax_enabled)
  {
    if (step_mem->q < 2)
    {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                      "The Butcher table(s) must be at least second order!");
      return ARK_INVALID_TABLE;
    }

    if (step_mem->explicit)
    {
      /* Check if all b values are positive */
      for (i = 0; i < step_mem->stages; i++)
      {
        if (step_mem->Be->b[i] < ZERO)
        {
          arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__,
                          __FILE__,
                          "The explicit Butcher table has a negative b value!");
          return ARK_INVALID_TABLE;
        }
      }
    }

    if (step_mem->implicit)
    {
      /* Check if all b values are positive */
      for (i = 0; i < step_mem->stages; i++)
      {
        if (step_mem->Bi->b[i] < ZERO)
        {
          arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__,
                          __FILE__,
                          "The implicit Butcher table has a negative b value!");
          return ARK_INVALID_TABLE;
        }
      }
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
int arkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess)
{
  int i, retval, jstage, nvec;
  sunrealtype tau;
  sunrealtype h;
  ARKodeARKStepMem step_mem;
  sunrealtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0) &&
      (step_mem->predictor < 4))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Interpolation structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step, use initial condition as guess */
  if (ark_mem->initsetup)
  {
    N_VScale(ONE, ark_mem->yn, yguess);
    return (ARK_SUCCESS);
  }

  /* set evaluation time tau as relative shift from previous successful time */
  tau = step_mem->Bi->c[istage] * ark_mem->h / ark_mem->hold;

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
      jstage = (step_mem->Bi->c[i] != ZERO) ? i : jstage;
    }

    /* if using the trivial predictor, break */
    if (jstage == -1) { break; }

    /* find the "optimal" previous stage to use */
    for (i = 0; i < istage; i++)
    {
      if ((step_mem->Bi->c[i] > step_mem->Bi->c[jstage]) &&
          (step_mem->Bi->c[i] != ZERO))
      {
        jstage = i;
      }
    }

    /* set stage time, stage RHS and interpolation values */
    h    = ark_mem->h * step_mem->Bi->c[jstage];
    tau  = ark_mem->h * step_mem->Bi->c[istage];
    nvec = 0;
    if (step_mem->implicit)
    { /* Implicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fi[jstage];
      nvec += 1;
    }
    if (step_mem->explicit)
    { /* Explicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fe[jstage];
      nvec += 1;
    }

    /* call predictor routine */
    retval = arkPredict_Bootstrap(ark_mem, h, tau, nvec, cvals, Xvecs, yguess);
    if (retval != ARK_ILL_INPUT) { return (retval); }
    break;

  case 5:

    /***** Minimal correction predictor: use all previous stage
           information in this step *****/

    /* set arrays for fused vector operation */
    nvec = 0;
    if (step_mem->explicit)
    { /* Explicit pieces */
      for (jstage = 0; jstage < istage; jstage++)
      {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[istage][jstage];
        Xvecs[nvec] = step_mem->Fe[jstage];
        nvec += 1;
      }
    }
    if (step_mem->implicit)
    { /* Implicit pieces */
      for (jstage = 0; jstage < istage; jstage++)
      {
        cvals[nvec] = ark_mem->h * step_mem->Bi->A[istage][jstage];
        Xvecs[nvec] = step_mem->Fi[jstage];
        nvec += 1;
      }
    }
    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec += 1;

    /* compute predictor */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yguess);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }
    return (ARK_SUCCESS);
    break;
  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->yn, yguess);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_StageSetup

  This routine sets up the stage data for computing the RK
  residual, along with the step- and method-related factors
  gamma, gammap and gamrat.

  The internal behavior of this setup depends on two factors:
  (a) whether the stage is explicit or implicit, and
  (b) the form of mass matrix (I, M, M(t)).

  In each of the following:
  * yn is the previous time step solution,
  * r corresponds to the nonlinear residual,
  * z=zp+zc corresponds to the updated stage solution, where
    zp is the predictor (zp=yn for explicit stages), and zc
    is the corrector,
  * i is the current stage index, and
  * g = h*Ai(i,i) is the implicit coefficient.

  Explicit, I:
      z = yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
     <=>
      zc = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes zc, stored in step_mem->sdata.

  Implicit, I:
      z = yn - h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
             - h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = zc - g*Fi(i) - s
      s = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Explicit, M:
      M*z = M*yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
     <=>
      M*zc = s
      s = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Implicit, M:
      M*z = M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                 + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = M*zc - g*Fi(i) - s
      s = M*(yn - zp) + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Explicit, M(t):
      z = yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j)+Ai(i,j)*Fi(j))
     <=>
      zc = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j)+Ai(i,j)*Fi(j))
    This routine computes zc, stored in step_mem->sdata.

  Implicit, M(t):
      M(t)*z = M(t)*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                       + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = M(t)*(zc - s) - g*Fi(i)
      s = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.


  Thus _internal_ to this routine, we have 3 modes:

  Explicit (any):
    sdata = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
  Implicit, M:
    sdata = M*(yn - zp) + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
  Implicit, I or M(t):
    sdata = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))

  ---------------------------------------------------------------*/
int arkStep_StageSetup(ARKodeMem ark_mem, sunbooleantype implicit)
{
  /* local data */
  ARKodeARKStepMem step_mem;
  int retval, i, j, jmax, nvec;
  sunrealtype* cj;
  sunrealtype** Aij;
  sunrealtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = step_mem->istage;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Update gamma if stage is implicit */
  if (implicit)
  {
    step_mem->gamma = ark_mem->h * step_mem->Bi->A[i][i];
    if (ark_mem->firststage) { step_mem->gammap = step_mem->gamma; }
    step_mem->gamrat = (ark_mem->firststage)
                         ? ONE
                         : step_mem->gamma /
                             step_mem->gammap; /* protect x/x != 1.0 */
  }

  /* If implicit, initialize sdata to yn - zpred (here: zpred = zp), and set
     first entries for eventual N_VLinearCombination call */
  nvec = 0;
  if (implicit)
  {
    N_VLinearSum(ONE, ark_mem->yn, -ONE, step_mem->zpred, step_mem->sdata);
    cvals[0] = ONE;
    Xvecs[0] = step_mem->sdata;
    nvec     = 1;
  }

  /* If implicit with fixed M!=I, update sdata with M*sdata */
  if (implicit && (step_mem->mass_type == MASS_FIXED))
  {
    N_VScale(ONE, step_mem->sdata, ark_mem->tempv1);
    retval = step_mem->mmult((void*)ark_mem, ark_mem->tempv1, step_mem->sdata);
    if (retval != ARK_SUCCESS) { return (ARK_MASSMULT_FAIL); }
  }

  /* Update sdata with prior stage information */
  if (step_mem->explicit)
  { /* Explicit pieces */
    for (j = 0; j < i; j++)
    {
      cvals[nvec] = ark_mem->h * step_mem->Be->A[i][j];
      Xvecs[nvec] = step_mem->Fe[j];
      nvec += 1;
    }
  }
  if (step_mem->implicit)
  { /* Implicit pieces */
    for (j = 0; j < i; j++)
    {
      cvals[nvec] = ark_mem->h * step_mem->Bi->A[i][j];
      Xvecs[nvec] = step_mem->Fi[j];
      nvec += 1;
    }
  }

  /* apply external polynomial (MRI) forcing (M = I required) */
  if (step_mem->expforcing || step_mem->impforcing)
  {
    if (step_mem->expforcing)
    {
      jmax = i;
      Aij  = step_mem->Be->A;
      cj   = step_mem->Be->c;
    }
    else
    {
      jmax = i + 1;
      Aij  = step_mem->Bi->A;
      cj   = step_mem->Bi->c;
    }

    for (j = 0; j < jmax; j++)
    {
      step_mem->stage_times[j] = ark_mem->tn + cj[j] * ark_mem->h;
      step_mem->stage_coefs[j] = ark_mem->h * Aij[i][j];
    }

    arkStep_ApplyForcing(step_mem, step_mem->stage_times, step_mem->stage_coefs,
                         jmax, &nvec);
  }

  /* call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
  if (retval != 0) { return (ARK_VECTOROP_ERR); }

  /* return with success */
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_mem->tempv1, in case
  the calling routine wishes to examine the error locations.

  This version assumes either an identity or time-dependent mass
  matrix (identical steps).
  ---------------------------------------------------------------*/
int arkStep_ComputeSolutions(ARKodeMem ark_mem, sunrealtype* dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  sunrealtype* cj;
  sunrealtype* bj;
  sunrealtype* dj;
  sunbooleantype stiffly_accurate;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;

  /* check if the method is stiffly accurate */
  stiffly_accurate = SUNTRUE;

  if (step_mem->explicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Be))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  if (step_mem->implicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Bi))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  /* If the method is stiffly accurate, ycur is already the new solution */

  if (!stiffly_accurate)
  {
    /* Compute time step solution (if necessary) */
    /*   set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = ark_mem->yn;
    nvec     = 1;
    for (j = 0; j < step_mem->stages; j++)
    {
      if (step_mem->explicit)
      { /* Explicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit)
      { /* Implicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if (step_mem->expforcing || step_mem->impforcing)
    {
      if (step_mem->expforcing)
      {
        cj = step_mem->Be->c;
        bj = step_mem->Be->b;
      }
      else
      {
        cj = step_mem->Bi->c;
        bj = step_mem->Bi->b;
      }

      for (j = 0; j < step_mem->stages; j++)
      {
        step_mem->stage_times[j] = ark_mem->tn + cj[j] * ark_mem->h;
        step_mem->stage_coefs[j] = ark_mem->h * bj[j];
      }

      arkStep_ApplyForcing(step_mem, step_mem->stage_times,
                           step_mem->stage_coefs, step_mem->stages, &nvec);
    }

    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }
  }

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    /* set arrays for fused vector operation */
    nvec = 0;
    for (j = 0; j < step_mem->stages; j++)
    {
      if (step_mem->explicit)
      { /* Explicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit)
      { /* Implicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if (step_mem->expforcing || step_mem->impforcing)
    {
      if (step_mem->expforcing)
      {
        cj = step_mem->Be->c;
        bj = step_mem->Be->b;
        dj = step_mem->Be->d;
      }
      else
      {
        cj = step_mem->Bi->c;
        bj = step_mem->Bi->b;
        dj = step_mem->Bi->d;
      }

      for (j = 0; j < step_mem->stages; j++)
      {
        step_mem->stage_times[j] = ark_mem->tn + cj[j] * ark_mem->h;
        step_mem->stage_coefs[j] = ark_mem->h * (bj[j] - dj[j]);
      }

      arkStep_ApplyForcing(step_mem, step_mem->stage_times,
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

/*---------------------------------------------------------------
  arkStep_ComputeSolutions_MassFixed

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_mem->tempv1, in case
  the calling routine wishes to examine the error locations.

  This version assumes a fixed mass matrix.
  ---------------------------------------------------------------*/
int arkStep_ComputeSolutions_MassFixed(ARKodeMem ark_mem, sunrealtype* dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  sunbooleantype stiffly_accurate;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem)ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;

  /* check if the method is stiffly accurate */
  stiffly_accurate = SUNTRUE;

  if (step_mem->explicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Be))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  if (step_mem->implicit)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(step_mem->Bi))
    {
      stiffly_accurate = SUNFALSE;
    }
  }

  /* If the method is stiffly accurate, ycur is already the new solution */

  if (!stiffly_accurate)
  {
    /* compute y RHS (store in y) */
    /*   set arrays for fused vector operation */
    nvec = 0;
    for (j = 0; j < step_mem->stages; j++)
    {
      if (step_mem->explicit)
      { /* Explicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit)
      { /* Implicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to compute RHS */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    /* solve for y update (stored in y) */
    retval = step_mem->msolve((void*)ark_mem, y, step_mem->nlscoef);
    if (retval < 0)
    {
      *dsmPtr = SUN_RCONST(2.0); /* indicate too much error, step with smaller step */
      N_VScale(ONE, ark_mem->yn, y); /* place old solution into y */
      return (CONV_FAIL);
    }

    /* compute y = yn + update */
    N_VLinearSum(ONE, ark_mem->yn, ONE, y, y);
  }

  /* compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    /* compute yerr RHS vector */
    /*   set arrays for fused vector operation */
    nvec = 0;
    for (j = 0; j < step_mem->stages; j++)
    {
      if (step_mem->explicit)
      { /* Explicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit)
      { /* Implicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to compute yerr RHS */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) { return (ARK_VECTOROP_ERR); }

    /* solve for yerr */
    retval = step_mem->msolve((void*)ark_mem, yerr, step_mem->nlscoef);
    if (retval < 0)
    {
      *dsmPtr = SUN_RCONST(2.0); /* next attempt will reduce step by 'etacf';
                                 insert dsmPtr placeholder here */
      return (CONV_FAIL);
    }
    /* fill error norm */
    *dsmPtr = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  Utility routines for interfacing with MRIStep
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  ARKStepCreateMRIStepInnerStepper

  Wraps an ARKStep memory structure as an MRIStep inner stepper.
  ----------------------------------------------------------------------------*/

int ARKStepCreateMRIStepInnerStepper(void* inner_arkode_mem,
                                     MRIStepInnerStepper* stepper)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  retval = arkStep_AccessStepMem(inner_arkode_mem,
                                 "ARKStepCreateMRIStepInnerStepper", &ark_mem,
                                 &step_mem);
  if (retval)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The ARKStep memory pointer is NULL");
    return ARK_ILL_INPUT;
  }

  retval = MRIStepInnerStepper_Create(ark_mem->sunctx, stepper);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetContent(*stepper, inner_arkode_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetEvolveFn(*stepper, arkStep_MRIStepInnerEvolve);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetFullRhsFn(*stepper,
                                            arkStep_MRIStepInnerFullRhs);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetResetFn(*stepper, arkStep_MRIStepInnerReset);
  if (retval != ARK_SUCCESS) { return (retval); }

  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  arkStep_MRIStepInnerEvolve

  Implementation of MRIStepInnerStepperEvolveFn to advance the inner (fast)
  ODE IVP.
  ----------------------------------------------------------------------------*/

int arkStep_MRIStepInnerEvolve(MRIStepInnerStepper stepper, sunrealtype t0,
                               sunrealtype tout, N_Vector y)
{
  void* arkode_mem;           /* arkode memory             */
  sunrealtype tret;           /* return time               */
  sunrealtype tshift, tscale; /* time normalization values */
  N_Vector* forcing;          /* forcing vectors           */
  int nforcing;               /* number of forcing vectors */
  int retval;                 /* return value              */

  /* extract the ARKODE memory struct */
  retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get the forcing data */
  retval = MRIStepInnerStepper_GetForcingData(stepper, &tshift, &tscale,
                                              &forcing, &nforcing);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set the inner forcing data */
  retval = arkStep_SetInnerForcing(arkode_mem, tshift, tscale, forcing, nforcing);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set the stop time */
  retval = ARKStepSetStopTime(arkode_mem, tout);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* evolve inner ODE */
  retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
  if (retval < 0) { return (retval); }

  /* disable inner forcing */
  retval = arkStep_SetInnerForcing(arkode_mem, ZERO, ONE, NULL, 0);
  if (retval != ARK_SUCCESS) { return (retval); }

  return (ARK_SUCCESS);
}

/*------------------------------------------------------------------------------
  arkStep_MRIStepInnerFullRhs

  Implementation of MRIStepInnerStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/

int arkStep_MRIStepInnerFullRhs(MRIStepInnerStepper stepper, sunrealtype t,
                                N_Vector y, N_Vector f, int mode)
{
  void* arkode_mem;
  int retval;

  /* extract the ARKODE memory struct */
  retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  return (arkStep_FullRHS(arkode_mem, t, y, f, mode));
}

/*------------------------------------------------------------------------------
  arkStep_MRIStepInnerReset

  Implementation of MRIStepInnerStepperResetFn to reset the inner (fast) stepper
  state.
  ----------------------------------------------------------------------------*/

int arkStep_MRIStepInnerReset(MRIStepInnerStepper stepper, sunrealtype tR,
                              N_Vector yR)
{
  void* arkode_mem;
  int retval;

  /* extract the ARKODE memory struct */
  retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  return (ARKStepReset(arkode_mem, tR, yR));
}

/*------------------------------------------------------------------------------
  arkStep_ApplyForcing

  Determines the scaling values and vectors necessary for the MRI polynomial
  forcing terms. This occurs through appending scaling values and N_Vector
  pointers to the underlying cvals and Xvecs arrays in the step_mem structure.

  stage_times -- The times at which to evaluate the forcing.

  stage_coefs -- Scaling factors (method A, b, or b - d coefficients) applied to
  forcing vectors.

  jmax -- the number of values in stage_times and stage_coefs (the stage index
  for explicit methods or the index + 1 for implicit methods).

  nvec -- On input, nvec is the next available entry in the cvals/Xvecs arrays.
  This value is incremented for each value/vector appended to the cvals/Xvecs
  arrays so on return it is the total number of values/vectors in the linear
  combination.
  ----------------------------------------------------------------------------*/

void arkStep_ApplyForcing(ARKodeARKStepMem step_mem, sunrealtype* stage_times,
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
  arkStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = fe(t,y) + fi(t,y) + p(t). This
  function is primarily intended for use with multirate integration methods
  (e.g., MRIStep) where ARKStep is used to solve a modified ODE at a fast time
  scale. The polynomial is of the form

  p(t) = sum_{i = 0}^{nvecs - 1} forcing[i] * ((t - tshift) / (tscale))^i

  where tshift and tscale are used to normalize the time t (e.g., with MRIGARK
  methods).
  ----------------------------------------------------------------------------*/

int arkStep_SetInnerForcing(void* arkode_mem, sunrealtype tshift,
                            sunrealtype tscale, N_Vector* forcing, int nvecs)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (nvecs > 0)
  {
    /* enable forcing */
    if (step_mem->explicit)
    {
      step_mem->expforcing = SUNTRUE;
      step_mem->impforcing = SUNFALSE;
    }
    else
    {
      step_mem->expforcing = SUNFALSE;
      step_mem->impforcing = SUNTRUE;
    }
    step_mem->tshift   = tshift;
    step_mem->tscale   = tscale;
    step_mem->forcing  = forcing;
    step_mem->nforcing = nvecs;

    /* If cvals and Xvecs are not allocated then arkStep_Init has not been
       called and the number of stages has not been set yet. These arrays will
       be allocated in arkStep_Init and take into account the value of nforcing.
       On subsequent calls will check if enough space has allocated in case
       nforcing has increased since the original allocation. */
    if (step_mem->cvals != NULL && step_mem->Xvecs != NULL)
    {
      /* check if there are enough reusable arrays for fused operations */
      if ((step_mem->nfusedopvecs - nvecs) < (2 * step_mem->stages + 2))
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
        step_mem->nfusedopvecs = 2 * step_mem->stages + 2 + nvecs;

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
    step_mem->expforcing = SUNFALSE;
    step_mem->impforcing = SUNFALSE;
    step_mem->tshift     = ZERO;
    step_mem->tscale     = ONE;
    step_mem->forcing    = NULL;
    step_mem->nforcing   = 0;
  }

  return (0);
}

/* -----------------------------------------------------------------------------
 * arkStep_RelaxDeltaE
 *
 * Computes the change in the relaxation functions for use in relaxation methods
 * delta_e = h * sum_i b_i * <relax_jac(z_i), f_i>
 *
 * With implicit and IMEX methods it is necessary to store the method stages
 * (or compute the delta_e estimate along the way) to avoid inconsistencies
 * between z_i, F(z_i), and J_relax(z_i) that arise from reconstructing stages
 * from stored RHS values like with ERK methods. As such the take step function
 * stores the stages along the way but only when there is an implicit RHS. When
 * a fixed mass matrix is present the stages are also stored to avoid additional
 * mass matrix solves in reconstructing the stages for an ERK method.
 *
 * Future: when ARKStep can exploit the structure of low storage methods, it may
 * be necessary to compute the delta_e estimate along the way with explicit
 * methods to avoid storing additional RHS or stage values.
 * ---------------------------------------------------------------------------*/

int arkStep_RelaxDeltaE(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                        long int* num_relax_jac_evals, sunrealtype* delta_e_out)
{
  int i, j, nvec, retval;
  sunrealtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;
  N_Vector z_stage = ark_mem->tempv2;
  N_Vector J_relax = ark_mem->tempv3;
  N_Vector rhs_tmp = NULL;
  sunrealtype bi   = ONE;

  /* Access the stepper memory structure */
  if (!(ark_mem->step_mem))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKSTEP_NO_MEM);
    return ARK_MEM_NULL;
  }
  step_mem = (ARKodeARKStepMem)(ark_mem->step_mem);

  /* Initialize output */
  *delta_e_out = ZERO;

  /* Set arrays for fused vector operation */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  for (i = 0; i < step_mem->stages; i++)
  {
    if (step_mem->implicit || step_mem->mass_type == MASS_FIXED)
    {
      /* Use stored stages */
      z_stage = step_mem->z[i];
    }
    else
    {
      /* Reconstruct explicit stages */
      nvec = 0;

      cvals[nvec] = ONE;
      Xvecs[nvec] = ark_mem->yn;
      nvec++;

      for (j = 0; j < i; j++)
      {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[i][j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec++;
      }

      retval = N_VLinearCombination(nvec, cvals, Xvecs, z_stage);
      if (retval) { return ARK_VECTOROP_ERR; }
    }

    /* Evaluate the Jacobian at z_i */
    retval = relax_jac_fn(z_stage, J_relax, ark_mem->user_data);
    (*num_relax_jac_evals)++;
    if (retval < 0) { return ARK_RELAX_JAC_FAIL; }
    if (retval > 0) { return ARK_RELAX_JAC_RECV; }

    /* Reset temporary RHS alias */
    rhs_tmp = z_stage;

    /* Compute delta_e = h * sum_i b_i * <relax_jac(z_i), f_i> */
    if (step_mem->explicit && step_mem->implicit)
    {
      N_VLinearSum(step_mem->Be->b[i], step_mem->Fe[i], step_mem->Bi->b[i],
                   step_mem->Fi[i], rhs_tmp);
      bi = ONE;
    }
    else if (step_mem->explicit)
    {
      if (step_mem->mass_type == MASS_FIXED)
      {
        N_VScale(ONE, step_mem->Fe[i], rhs_tmp);
      }
      else { rhs_tmp = step_mem->Fe[i]; }
      bi = step_mem->Be->b[i];
    }
    else
    {
      if (step_mem->mass_type == MASS_FIXED)
      {
        N_VScale(ONE, step_mem->Fi[i], rhs_tmp);
      }
      else { rhs_tmp = step_mem->Fi[i]; }
      bi = step_mem->Bi->b[i];
    }

    if (step_mem->mass_type == MASS_FIXED)
    {
      retval = step_mem->msolve((void*)ark_mem, rhs_tmp, step_mem->nlscoef);
      if (retval) { return ARK_MASSSOLVE_FAIL; }
    }

    /* Update estimate of relaxation function change */
    if (J_relax->ops->nvdotprodlocal && J_relax->ops->nvdotprodmultiallreduce)
    {
      *delta_e_out += bi * N_VDotProdLocal(J_relax, rhs_tmp);
    }
    else { *delta_e_out += bi * N_VDotProd(J_relax, rhs_tmp); }
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
 * arkStep_GetOrder
 *
 * Returns the method order
 * ---------------------------------------------------------------------------*/

int arkStep_GetOrder(ARKodeMem ark_mem)
{
  ARKodeARKStepMem step_mem = (ARKodeARKStepMem)(ark_mem->step_mem);
  return step_mem->q;
}

/*===============================================================
  EOF
  ===============================================================*/
