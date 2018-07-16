/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for ARKode's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

#define NO_DEBUG_OUTPUT
/* #define DEBUG_OUTPUT */
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)

#define FIXED_LIN_TOL


/*===============================================================
  ARKStep Exported functions -- Required
  ===============================================================*/

void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  booleantype nvectorOK;
  int iret;

  /* Create ark_mem structure */
  ark_mem = arkCreate();
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Allocate ARKodeARKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeARKStepMem) malloc(sizeof(struct ARKodeARKStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
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
  ark_mem->step                = arkStep_TakeStep;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for ARKStep optional inputs */
  iret = ARKStepSetDefaults((void *)ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKode::ARKStep",
                    "ARKStepCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (!step_mem->implicit && !step_mem->explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = arkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Allocate the general ARK stepper vectors using y0 as a template */
  /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on
     (based on the number of ARK stages) */

  /* Clone the input vector to create sdata, zpred */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata)))
    return(NULL);
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zpred)))
    return(NULL);

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit  = NULL;
  step_mem->lsetup = NULL;
  step_mem->lsolve = NULL;
  step_mem->lfree  = NULL;
  step_mem->lmem   = NULL;
  step_mem->lsolve_type = -1;

  /* Set the mass matrix solver addresses to NULL */
  step_mem->minit    = NULL;
  step_mem->msetup   = NULL;
  step_mem->mmult    = NULL;
  step_mem->msolve   = NULL;
  step_mem->mfree    = NULL;
  step_mem->mass_mem = NULL;
  step_mem->msetuptime = -RCONST(99999999999.0);
  step_mem->msolve_type = -1;

  /* Allocate step adaptivity structure, set default values, note storage */
  step_mem->hadapt_mem = arkAdaptInit();
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepCreate",
                    "Allocation of step adaptivity structure failed");
    return(NULL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;

  /* Initialize initial error norm  */
  step_mem->eRNrm = 1.0;

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->nfi          = 0;
  step_mem->ncfn         = 0;
  step_mem->netf         = 0;
  step_mem->nni          = 0;
  step_mem->nsetups      = 0;
  step_mem->nstlp        = 0;

  /* Initialize main ARKode infrastructure */
  iret = arkInit(ark_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKode::ARKStep", "ARKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
 ARKStepResize:

 This routine resizes the memory within the ARKStep module.
 It first resizes the main ARKode infrastructure memory, and
 then resizes its own data.
---------------------------------------------------------------*/
int ARKStepResize(void *arkode_mem, N_Vector y0, realtype hscale,
                  realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int ier, i, flag;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResize", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResize", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKode infrastructure memory */
  flag = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (flag != ARK_SUCCESS) {
    arkProcessError(ark_mem, flag, "ARKode::ARKStep", "ARKStepResize",
                    "Unable to resize main ARKode infrastructure");
    return(flag);
  }

  /* Resize the sdata and zpred vectors */
  if (step_mem->sdata != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->sdata);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (step_mem->zpred != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->zpred);
    if (ier != ARK_SUCCESS)  return(ier);
  }

  /* Resize the ARKStep vectors */
  /*     Fe */
  if (step_mem->Fe != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, y0, &step_mem->Fe[i]);
      if (ier != ARK_SUCCESS)  return(ier);
    }
  }
  /*     Fi */
  if (step_mem->Fi != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, y0, &step_mem->Fi[i]);
      if (ier != ARK_SUCCESS)  return(ier);
    }
  }

  /* Resize fixed-point solver memory */
  if ( (step_mem->use_fp) && (step_mem->fp_mem == NULL) ) {
    ier = arkResizeFPData(ark_mem, step_mem->fp_mem, resize,
                          resize_data, lrw_diff, liw_diff, y0);
    if (ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepResize", MSG_ARK_MEM_FAIL);
      return(ARK_MEM_FAIL);
    }
  }

  /* Re-initialize the linear solver (assumes that solver
     memory has already been resized appropriately) */
  if (step_mem->implicit) {
    if (step_mem->linit != NULL) {
      ier = step_mem->linit(ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKode::ARKStep",
                        "ARKStepResize", MSG_ARK_LINIT_FAIL);
        return(ARK_LINIT_FAIL);
      }
    }
  }

  /* Re-initialize the mass matrix solver (assumes that solver
     memory has already been resized appropriately) */
  if (step_mem->minit != NULL) {
    ier = step_mem->minit((void *) ark_mem);
    if (ier != 0) {
      arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKode::ARKStep",
                      "ARKStepResize", MSG_ARK_MASSINIT_FAIL);
      return(ARK_MASSINIT_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


int ARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                  ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int iret;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access ARKStep memory */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepReInit",
                    "Missing time step module memory");
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  iret = arkReInit(ark_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKode::ARKStep", "ARKStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(iret);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (!step_mem->implicit && !step_mem->explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Destroy/Reinitialize time step adaptivity structure (if present) */
  if (step_mem->hadapt_mem != NULL) {
    free(step_mem->hadapt_mem);
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepReInit",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }
  /* Initialize initial error norm  */
  step_mem->eRNrm = 1.0;

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->nfi          = 0;
  step_mem->ncfn         = 0;
  step_mem->netf         = 0;
  step_mem->nni          = 0;
  step_mem->nsetups      = 0;
  step_mem->nstlp        = 0;

  return(ARK_SUCCESS);
}


int ARKStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}


int ARKStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}


int ARKStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}


int ARKStepResStolerance(void *arkode_mem, realtype rabstol)
{
  /* unpack ark_mem, call arkResStolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResStolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResStolerance(ark_mem, rabstol));
}


int ARKStepResVtolerance(void *arkode_mem, N_Vector rabstol)
{
  /* unpack ark_mem, call arkResVtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResVtolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResVtolerance(ark_mem, rabstol));
}


int ARKStepResFtolerance(void *arkode_mem, ARKRwtFn rfun)
{
  /* unpack ark_mem, call arkResFtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResFtolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResFtolerance(ark_mem, rfun));
}


int ARKStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


int ARKStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkEvolve(ark_mem, tout, yout, tret, itask));
}


int ARKStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetDky(ark_mem, t, k, dky));
}


/*---------------------------------------------------------------
  ARKStepFree frees all ARKStep memory, and then calls an ARKode
  utility routine to free the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void ARKStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL ARKStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

    /* free the time step adaptivity module */
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
      ark_mem->lrw -= ARK_ADAPT_LRW;
      ark_mem->liw -= ARK_ADAPT_LIW;
    }

    /* free the Butcher tables */
    if (step_mem->Be != NULL) {
      ButcherTableSpace(step_mem->Be, &Bliw, &Blrw);
      FreeButcherTable(step_mem->Be);
      step_mem->Be = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }
    if (step_mem->Bi != NULL) {
      ButcherTableSpace(step_mem->Bi, &Bliw, &Blrw);
      FreeButcherTable(step_mem->Bi);
      step_mem->Bi = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }

    /* free the fixed-point solver memory */
    if (step_mem->fp_mem != NULL) {
      arkFreeFPData(ark_mem, step_mem->fp_mem);
      step_mem->fp_mem = NULL;
    }

    /* free the linear solver memory */
    if (step_mem->lfree != NULL) {
      step_mem->lfree((void *) ark_mem);
      step_mem->lmem = NULL;
    }

    /* free the mass matrix solver memory */
    if (step_mem->mfree != NULL) {
      step_mem->mfree((void *) ark_mem);
      step_mem->mass_mem = NULL;
    }

    /* free the sdata and zpred vectors */
    if (step_mem->sdata != NULL) {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }
    if (step_mem->zpred != NULL) {
      arkFreeVec(ark_mem, &step_mem->zpred);
      step_mem->zpred = NULL;
    }

    /* free the RHS vectors */
    if (step_mem->Fe != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->Fe[j]);
      free(step_mem->Fe);
      step_mem->Fe = NULL;
      ark_mem->liw -= step_mem->stages;
    }
    if (step_mem->Fi != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->Fi[j]);
      free(step_mem->Fi);
      step_mem->Fi = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL) {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= (2*step_mem->stages + 1);
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (2*step_mem->stages + 1);
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;

  }

  /* free memory for overall ARKode infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
 ARKStepPrintMem:

 This routine outputs the memory from the ARKStep structure and
 the main ARKode infrastructure to a specified file pointer
 (useful when debugging).
---------------------------------------------------------------*/
void ARKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int i;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepPrintMem", MSG_ARK_NO_MEM);
    return;
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepPrintMem", MSG_ARKSTEP_NO_MEM);
    return;
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* output data from main ARKode infrastructure */
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"ARKStep: q = %i\n", step_mem->q);
  fprintf(outfile,"ARKStep: p = %i\n", step_mem->p);
  fprintf(outfile,"ARKStep: istage = %i\n", step_mem->istage);
  fprintf(outfile,"ARKStep: stages = %i\n", step_mem->stages);
  fprintf(outfile,"ARKStep: mnewt = %i\n", step_mem->mnewt);
  fprintf(outfile,"ARKStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile,"ARKStep: maxnef = %i\n", step_mem->maxnef);
  fprintf(outfile,"ARKStep: maxncf = %i\n", step_mem->maxncf);
  fprintf(outfile,"ARKStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile,"ARKStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile,"ARKStep: lsolve_type = %i\n", step_mem->lsolve_type);
  fprintf(outfile,"ARKStep: msolve_type = %i\n", step_mem->msolve_type);

  /* output long integer quantities */
  fprintf(outfile,"ARKStep: nst_attempts = %li\n", step_mem->nst_attempts);
  fprintf(outfile,"ARKStep: nfe = %li\n", step_mem->nfe);
  fprintf(outfile,"ARKStep: nfi = %li\n", step_mem->nfi);
  fprintf(outfile,"ARKStep: ncfn = %li\n", step_mem->ncfn);
  fprintf(outfile,"ARKStep: netf = %li\n", step_mem->netf);
  fprintf(outfile,"ARKStep: nni = %li\n", step_mem->nni);
  fprintf(outfile,"ARKStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile,"ARKStep: use_fp = %i\n", step_mem->use_fp);
  if (step_mem->fp_mem != NULL) {
    fprintf(outfile,"ARKStep: fixed-point solver structure:\n");
    fprintf(outfile, "    m = %li\n", step_mem->fp_mem->m);
    if (step_mem->fp_mem->imap != NULL)
      for (i=0; i<step_mem->fp_mem->m; i++)
        fprintf(outfile, "   imap[%i] = %li\n", i, step_mem->fp_mem->imap[i]);
    if (step_mem->fp_mem->R != NULL) {
      fprintf(outfile, "    R =  ");
      for (i=0; i<step_mem->fp_mem->m*step_mem->fp_mem->m; i++)
        fprintf(outfile, "%"RSYM"  ", step_mem->fp_mem->R[i]);
      fprintf(outfile, "\n");
    }
    if (step_mem->fp_mem->gamma != NULL) {
      fprintf(outfile, "    gamma =  ");
      for (i=0; i<step_mem->fp_mem->m; i++)
        fprintf(outfile, "%"RSYM"  ", step_mem->fp_mem->gamma[i]);
      fprintf(outfile, "\n");
    }
  }
  fprintf(outfile,"ARKStep: nstlp = %li\n", step_mem->nstlp);

  /* output boolean quantities */
  fprintf(outfile,"ARKStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile,"ARKStep: user_linear_timedep = %i\n", step_mem->linear_timedep);
  fprintf(outfile,"ARKStep: user_explicit = %i\n", step_mem->explicit);
  fprintf(outfile,"ARKStep: user_implicit = %i\n", step_mem->implicit);
  fprintf(outfile,"ARKStep: hadapt_pq = %i\n", step_mem->hadapt_pq);
  fprintf(outfile,"ARKStep: jcur = %i\n", step_mem->jcur);

  /* output realtype quantities */
  if (step_mem->Be != NULL) {
    fprintf(outfile,"ARKStep: explicit Butcher table:\n");
    WriteButcherTable(step_mem->Be, outfile);
  }
  if (step_mem->Bi != NULL) {
    fprintf(outfile,"ARKStep: implicit Butcher table:\n");
    WriteButcherTable(step_mem->Bi, outfile);
  }
  fprintf(outfile,"ARKStep: gamma = %"RSYM"\n", step_mem->gamma);
  fprintf(outfile,"ARKStep: gammap = %"RSYM"\n", step_mem->gammap);
  fprintf(outfile,"ARKStep: gamrat = %"RSYM"\n", step_mem->gamrat);
  fprintf(outfile,"ARKStep: crate = %"RSYM"\n", step_mem->crate);
  fprintf(outfile,"ARKStep: eRNrm = %"RSYM"\n", step_mem->eRNrm);
  fprintf(outfile,"ARKStep: nlscoef = %"RSYM"\n", step_mem->nlscoef);
  if (step_mem->hadapt_mem != NULL) {
    fprintf(outfile,"ARKStep: timestep adaptivity structure:\n");
    arkPrintAdaptMem(step_mem->hadapt_mem, outfile);
  }
  fprintf(outfile,"ARKStep: crdown = %"RSYM"\n", step_mem->crdown);
  fprintf(outfile,"ARKStep: rdiv = %"RSYM"\n", step_mem->rdiv);
  fprintf(outfile,"ARKStep: dgmax = %"RSYM"\n", step_mem->dgmax);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */
  if (step_mem->sdata != NULL) {
    fprintf(outfile, "ARKStep: sdata:\n");
    N_VPrint_Serial(step_mem->sdata);
  }
  if (step_mem->zpred != NULL) {
    fprintf(outfile, "ARKStep: zpred:\n");
    N_VPrint_Serial(step_mem->zpred);
  }
  if (step_mem->Fe != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fe[%i]:\n", i);
      N_VPrint_Serial(step_mem->Fe[i]);
    }
  if (step_mem->Fi != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fi[%i]:\n", i);
      N_VPrint_Serial(step_mem->Fi[i]);
    }
#endif
}


/*===============================================================
  ARKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
 arkStep_AttachLinsol:

 This routine attaches the various set of system linear solver
 interface routines, data structure, and solver type to the
 ARKStep module.
---------------------------------------------------------------*/
int arkStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup,
                         ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         int lsolve_type, void *lmem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_AttachLinsol", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_AttachLinsol", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* free any existing system solver */
  if (step_mem->lfree != NULL)  step_mem->lfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit       = linit;
  step_mem->lsetup      = lsetup;
  step_mem->lsolve      = lsolve;
  step_mem->lfree       = lfree;
  step_mem->lmem        = lmem;
  step_mem->lsolve_type = lsolve_type;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_AttachMasssol:

 This routine attaches the set of mass matrix linear solver
 interface routines, data structure, and solver type to the
 ARKStep module.
---------------------------------------------------------------*/
int arkStep_AttachMasssol(void* arkode_mem, ARKMassInitFn minit,
                          ARKMassSetupFn msetup,
                          ARKMassMultFn mmult,
                          ARKMassSolveFn msolve,
                          ARKMassFreeFn mfree,
                          int msolve_type, void *mass_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_AttachMasssol", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_AttachMasssol", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* free any existing mass matrix solver */
  if (step_mem->mfree != NULL)  step_mem->mfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->minit       = minit;
  step_mem->msetup      = msetup;
  step_mem->mmult       = mmult;
  step_mem->msolve      = msolve;
  step_mem->mfree       = mfree;
  step_mem->mass_mem    = mass_mem;
  step_mem->msolve_type = msolve_type;

  /* Attach mmult function pointer to ark_mem as well */
  ark_mem->step_mmult = mmult;

  return(ARK_SUCCESS);
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
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

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
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

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

  /* access ARKodeARKStepMem structure, and return lmem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetLmem", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetLmem", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(step_mem->lmem);
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

  /* access ARKodeARKStepMem structure, and return lmem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetMassMem", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetMassMem", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(step_mem->mass_mem);
}


/*---------------------------------------------------------------
 arkStep_GetImplicitRHS:

 This routine returns the implicit RHS function pointer, fi.
---------------------------------------------------------------*/
ARKRhsFn arkStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure, and return fi */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetImplicitRHS", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetImplicitRHS", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(step_mem->fi);
}


/*---------------------------------------------------------------
 arkStep_GetGammas:

 This routine fills the current value of gamma, and states
 whether the gamma ratio fails the dgmax criteria.
---------------------------------------------------------------*/
int arkStep_GetGammas(void* arkode_mem, realtype *gamma,
                      realtype *gamrat, booleantype **jcur,
                      booleantype *dgamma_fail)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure, and return fi */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetGammas", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_GetGammas", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  *gamma  = step_mem->gamma;
  *gamrat = step_mem->gamrat;
  *jcur = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_Init:

 Called from within arkInitialSetup, this routine:
 - sets/checks the ARK Butcher tables to be used
 - allocates any memory that depends on the number of ARK stages,
   method order, or solver options
 - checks for consistency between the system and mass matrix
   linear solvers (if applicable)
 - initializes and sets up the system and mass matrix linear
   solvers (if applicable)
 - allocates the interpolation data structure (if needed based
   on ARKStep solver options)
---------------------------------------------------------------*/
int arkStep_Init(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int ier, j;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Init", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Init", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* destroy adaptivity structure if fixed-stepping is requested */
  if (ark_mem->fixedstep)
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }

  /* Set first step growth factor */
  if (step_mem->hadapt_mem != NULL)
    step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->etamx1;

  /* Create Butcher tables (if not already set) */
  ier = arkStep_SetButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                    "Could not create Butcher table(s)");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher tables are OK */
  ier = arkStep_CheckButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_Init", "Error in Butcher table(s)");
    return(ARK_ILL_INPUT);
  }

  /* note Butcher table space requirements */
  ButcherTableSpace(step_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;
  ButcherTableSpace(step_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate Fe[0] ... Fe[stages-1] if needed */
  if (step_mem->explicit) {
    if (step_mem->Fe == NULL)
      step_mem->Fe = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
    for (j=0; j<step_mem->stages; j++) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fe[j])))
        return(ARK_MEM_FAIL);
    }
    ark_mem->liw += step_mem->stages;  /* pointers */
  }

  /*   Allocate Fi[0] ... Fi[stages-1] if needed */
  if (step_mem->implicit) {
    if (step_mem->Fi == NULL)
      step_mem->Fi = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
    for (j=0; j<step_mem->stages; j++) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fi[j])))
        return(ARK_MEM_FAIL);
    }
    ark_mem->liw += step_mem->stages;  /* pointers */
  }

  /* Allocate reusable arrays for fused vector interface */
  if (step_mem->cvals == NULL) {
    step_mem->cvals = (realtype *) calloc(2*step_mem->stages+1, sizeof(realtype));
    if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
    ark_mem->lrw += (2*step_mem->stages + 1);
  }
  if (step_mem->Xvecs == NULL) {
    step_mem->Xvecs = (N_Vector *) calloc(2*step_mem->stages+1, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
    ark_mem->liw += (2*step_mem->stages + 1);   /* pointers */
  }

  /* Check for consistency between linear system modules
     (if lsolve is direct, msolve needs to match) */
  if (step_mem->mass_mem != NULL) {  /* M != I */
    if (step_mem->lsolve_type != step_mem->msolve_type) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                      "Incompatible linear and mass matrix solvers");
      return(ARK_ILL_INPUT);
    }
  }

  /* Perform mass matrix solver initialization and setup (if applicable) */
  if (step_mem->mass_mem != NULL) {

    /* Call minit (if it exists) */
    if (step_mem->minit != NULL) {
      ier = step_mem->minit((void *) ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_MASSINIT_FAIL);
        return(ARK_MASSINIT_FAIL);
      }
    }

    /* Call msetup (if it exists) -- use ewt, ark_tempv2 and sdata as temp vectors */
    if (step_mem->msetup != NULL) {
      ier = step_mem->msetup((void *) ark_mem, ark_mem->ewt,
                                ark_mem->tempv2, step_mem->sdata);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_MASSSETUP_FAIL);
        return(ARK_MASSSETUP_FAIL);
        step_mem->msetuptime = ark_mem->tcur;
      }
    }
  }

  /* Check if lsolve function exists and call linit (if it exists) */
  if (step_mem->implicit && !step_mem->use_fp) {
    if (step_mem->lsolve == NULL) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_Init", MSG_ARK_LSOLVE_NULL);
      return(ARK_ILL_INPUT);
    }
    if (step_mem->linit) {
      ier = step_mem->linit(ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_LINIT_FAIL);
        return(ARK_LINIT_FAIL);
      }
    }
  }

  /* Allocate fixed-point solver memory */
  if ( (step_mem->use_fp) && (step_mem->fp_mem == NULL) ){
    step_mem->fp_mem = arkAllocFPData(ark_mem, step_mem->fp_m);
    if (step_mem->fp_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "arkStep_Init",
                      "Unable to allocate fixed-point solver memory");
      return(ARK_MEM_FAIL);
    }
  }

  /* Allocate interpolation memory (if unallocated, and needed) */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0)) {
    ark_mem->interp = arkInterpCreate(ark_mem);
    if (ark_mem->interp == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "arkStep_Init",
                      "Unable to allocate interpolation structure");
      return(ARK_MEM_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_FullRHS:

 Rewriting the problem
    My' = fe(t,y) + fi(t,y)
 in the form
    y' = M^{-1}*[ fe(t,y) + fi(t,y) ],
 this routine computes the full right-hand side vector,
    f = M^{-1}*[ fe(t,y) + fi(t,y) ]

 This will be called in one of three 'modes':
     0 -> called at the beginning of a simulation
     1 -> called at the end of a successful step
     2 -> called elsewhere (e.g. for dense output)

 If it is called in mode 0, we store the vectors fe(t,y) and
 fi(t,y) in Fe[0] and Fi[0] for possible reuse in the first
 stage of the subsequent time step.

 If it is called in mode 1 and the ARK method coefficients
 support it, we may just copy vectors Fe[stages] and Fi[stages]
 to fill f instead of calling fe() and fi().

 Mode 2 is only called for dense output in-between steps, so we
 strive to store the intermediate parts so that they do not
 interfere with the other two modes.
---------------------------------------------------------------*/
int arkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;
  booleantype recomputeRHS;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_FullRHS", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_FullRHS", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if the problem involves a non-identity mass matrix and setup is
     required, do so here (use f, ark_tempv2 and sdata as a temporaries) */
  if ( (step_mem->mass_mem != NULL) && (step_mem->msetup != NULL) )
    if (SUNRabs(step_mem->msetuptime - t) > FUZZ_FACTOR*ark_mem->uround) {
      retval = step_mem->msetup((void *) ark_mem, f, ark_mem->tempv2,
                                   step_mem->sdata);
      if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
      step_mem->msetuptime = t;
    }

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vectors fe(t,y) and fi(t,y) in Fe[0] and Fi[0] for
     possible reuse in the first stage of the subsequent time step */
  case 0:

    /* call fe if the problem has an explicit component */
    if (step_mem->explicit) {
      retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fi if the problem has an implicit component */
    if (step_mem->implicit) {
      retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
      step_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    } else {                                              /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    break;


  /* Mode 1: called at the end of a successful step
     If the ARK method coefficients support it, we just copy the last stage RHS vectors
     to fill f instead of calling fe() and fi().
     Copy the results to Fe[0] and Fi[0] if the ARK coefficients support it. */
  case 1:

    /* determine if explicit/implicit RHS functions need to be recomputed */
    recomputeRHS = SUNFALSE;
    if ( step_mem->explicit && (SUNRabs(step_mem->Be->c[step_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;
    if ( step_mem->implicit && (SUNRabs(step_mem->Bi->c[step_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;

    /* base RHS calls on recomputeRHS argument */
    if (recomputeRHS) {

      /* call fe if the problem has an explicit component */
      if (step_mem->explicit) {
        retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
        step_mem->nfe++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
          return(ARK_RHSFUNC_FAIL);
        }
      }

      /* call fi if the problem has an implicit component */
      if (step_mem->implicit) {
        retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
        step_mem->nfi++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
          return(ARK_RHSFUNC_FAIL);
        }
      }
    } else {
      if (step_mem->explicit)
        N_VScale(ONE, step_mem->Fe[step_mem->stages-1], step_mem->Fe[0]);
      if (step_mem->implicit)
        N_VScale(ONE, step_mem->Fi[step_mem->stages-1], step_mem->Fi[0]);
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    } else {                                              /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    break;

  /*  Mode 2: called for dense output in-between steps
      store the intermediate calculations in such a way as to not
      interfere with the other two modes */
  default:

    /* call fe if the problem has an explicit component (store in ark_tempv2) */
    if (step_mem->explicit) {
      retval = step_mem->fe(t, y, ark_mem->tempv2, ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fi if the problem has an implicit component (store in sdata) */
    if (step_mem->implicit) {
      retval = step_mem->fi(t, y, step_mem->sdata, ark_mem->user_data);
      step_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->sdata, ONE, ark_mem->tempv2, f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->sdata, f);
    } else {                                              /* explicit */
      N_VScale(ONE, ark_mem->tempv2, f);
    }

    break;
  }


  /* if M != I, then update f = M^{-1}*f
     MAY WANT TO MOVE THIS TO BE CONDITIONALLY CALLED INSIDE THE PRECEDING CASE BLOCKS */
  if (step_mem->mass_mem != NULL) {
    /* scale RHS */
    N_VScale(ark_mem->h, f, f);
    /* update fnew = M^{-1}*fnew */
    retval = step_mem->msolve((void *) ark_mem, f, step_mem->nlscoef);
    /* scale result */
    N_VScale(ONE/ark_mem->h, f, f);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKode::ARKStep",
                      "arkStep_FullRHS", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_TakeStep:

 This routine serves the primary purpose of the ARKStep module:
 it performs a single successful ARK step (with embedding, if
 possible).  Multiple attempts may be taken in this process --
 once a step completes with successful (non)linear solves at
 each stage and passes the error estimate, the routine returns
 successfully.  If it cannot do so, it returns with an
 appropriate error flag.
---------------------------------------------------------------*/
int arkStep_TakeStep(void* arkode_mem)
{
  realtype dsm;
  int retval, ncf, nef, is, nflag, kflag, eflag;
  booleantype implicit_stage;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_TakeStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_TakeStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  ncf = nef = 0;
  nflag = FIRST_CALL;
  eflag = ARK_SUCCESS;
  kflag = SOLVE_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {

    /* increment attempt counter */
    step_mem->nst_attempts++;

    /* Loop over internal stages to the step */
    for (is=0; is<step_mem->stages; is++) {

      /* store current stage index */
      step_mem->istage = is;

      /* Set current stage time(s) */
      if (step_mem->implicit)
        ark_mem->tcur = ark_mem->tn + step_mem->Bi->c[is]*ark_mem->h;
      else
        ark_mem->tcur = ark_mem->tn + step_mem->Be->c[is]*ark_mem->h;

#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
         ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

      /* determine whether implicit solve is required */
      implicit_stage = SUNFALSE;
      if (step_mem->implicit)
        if (SUNRabs(step_mem->Bi->A[is][is]) > TINY)
          implicit_stage = SUNTRUE;

      /* Call predictor for current stage solution (result placed in zpred) */
      if (implicit_stage) {
        eflag = arkStep_Predict(ark_mem, is, step_mem->zpred);
        if (eflag != ARK_SUCCESS)  return (eflag);
      } else {
        N_VScale(ONE, ark_mem->yn, step_mem->zpred);
      }

#ifdef DEBUG_OUTPUT
 printf("predictor:\n");
 N_VPrint_Serial(step_mem->zpred);
#endif

      /* Set up data for evaluation of ARK stage residual (data stored in sdata) */
      eflag = arkStep_StageSetup(ark_mem);
      if (eflag != ARK_SUCCESS)  return (eflag);

#ifdef DEBUG_OUTPUT
 printf("rhs data:\n");
 N_VPrint_Serial(step_mem->sdata);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report)
        fprintf(ark_mem->diagfp, "step  %li  %"RSYM"  %i  %"RSYM"\n",
                ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

      /* perform implicit solve if required */
      if (implicit_stage) {

        /* perform implicit solve (result is stored in ark_mem->ycur) */
        nflag = arkStep_Nls(ark_mem, nflag);

#ifdef DEBUG_OUTPUT
 printf("nonlinear solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

        /* check for convergence (on failure, h will have been modified) */
        kflag = arkStep_HandleNFlag(ark_mem, &nflag, &ncf);

        /* If fixed time-stepping is used, then anything other than a
           successful solve must result in an error */
        if (ark_mem->fixedstep && (kflag != SOLVE_SUCCESS))
          return(kflag);

        /* If h reduced and step needs to be retried, break loop */
        if (kflag == PREDICT_AGAIN) break;

        /* Return if nonlinear solve failed and recovery not possible. */
        if (kflag != SOLVE_SUCCESS) return(kflag);

      /* otherwise no implicit solve is needed */
      } else {

        /* if M!=I, solve with M to compute update (place back in sdata) */
        if (step_mem->mass_mem != NULL) {

          /* /\* call msetup if needed -- NEED TEMPVEC1 and TEMPVEC2 *\/ */
          /* if (step_mem->msetup != NULL) { */
          /*   if (SUNRabs(step_mem->msetuptime - ark_mem->tcur) > */
          /*       FUZZ_FACTOR*ark_mem->uround) { */
          /*     eflag = step_mem->msetup((void *) ark_mem, ark_mem->tcur, */
          /*                                 ark_mem->ycur, TEMPVEC1, TEMPVEC2); */
          /*     if (eflag != ARK_SUCCESS) { */
          /*       arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKode::ARKStep", */
          /*                       "arkStep_TakeStep", MSG_ARK_MASSSSETUP_FAIL); */
          /*       return(ARK_MASSSETUP_FAIL); */
          /*     } */
          /*     step_mem->msetuptime = tend; */
          /*   } */
          /* } */

          /* perform mass matrix solve */
          nflag = step_mem->msolve((void *) ark_mem, step_mem->sdata,
                                      step_mem->nlscoef);

          /* check for convergence (on failure, h will have been modified) */
          kflag = arkStep_HandleNFlag(ark_mem, &nflag, &ncf);

          /* If fixed time-stepping is used, then anything other than a
             successful solve must result in an error */
          if (ark_mem->fixedstep && (kflag != SOLVE_SUCCESS))
            return(kflag);

          /* If h reduced and step needs to be retried, break loop */
          if (kflag == PREDICT_AGAIN) break;

          /* Return if solve failed and recovery not possible. */
          if (kflag != SOLVE_SUCCESS) return(kflag);

        /* if M==I, set y to be zpred + RHS data computed in arkStep_StageSetup */
        } else {
          N_VLinearSum(ONE, step_mem->sdata, ONE,
                       step_mem->zpred, ark_mem->ycur);
        }

#ifdef DEBUG_OUTPUT
 printf("explicit solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

      }

      /* successful stage solve */
      /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
      if (step_mem->implicit) {
        retval = step_mem->fi(ark_mem->tcur, ark_mem->ycur,
                                 step_mem->Fi[is], ark_mem->user_data);
        step_mem->nfi++;
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

      /*    store explicit RHS if necessary
            (already computed at first stage of purely explicit runs) */
      if (step_mem->explicit) {
          retval = step_mem->fe(ark_mem->tn + step_mem->Be->c[is]*ark_mem->h,
                                   ark_mem->ycur, step_mem->Fe[is], ark_mem->user_data);
          step_mem->nfe++;
          if (retval < 0)  return(ARK_RHSFUNC_FAIL);
          if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

    } /* loop over stages */

    /* if h has changed due to convergence failure and a new
       prediction is needed, continue to next attempt at step
       (cannot occur if fixed time stepping is enabled) */
    if (kflag == PREDICT_AGAIN)  continue;

    /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
    retval = arkStep_ComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);

#ifdef DEBUG_OUTPUT
 printf("error estimate = %"RSYM"\n", dsm);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "  etest  %li  %"RSYM"  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, dsm);

    /* Perform time accuracy error test (if failure, updates h for next try) */
    if (!ark_mem->fixedstep)
      eflag = arkStep_DoErrorTest(ark_mem, &nflag, &nef, dsm);

#ifdef DEBUG_OUTPUT
 printf("error test flag = %i\n", eflag);
#endif

    /* Restart step attempt (recompute all stages) if error test fails recoverably */
    if (eflag == TRY_AGAIN)  continue;

    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS)  return(eflag);

    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  } /* loop over step attempts */


  /* The step has completed successfully, clean up and
     consider change of step size */
  retval = arkStep_PrepareNextStep(ark_mem, dsm);
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
 arkStep_CheckNVector:

 This routine checks if all required vector operations are
 present.  If any of them is missing it returns SUNFALSE.
---------------------------------------------------------------*/
booleantype arkStep_CheckNVector(N_Vector tmpl)
{
  if ( (tmpl->ops->nvclone     == NULL) ||
       (tmpl->ops->nvdestroy   == NULL) ||
       (tmpl->ops->nvlinearsum == NULL) ||
       (tmpl->ops->nvconst     == NULL) ||
       (tmpl->ops->nvscale     == NULL) ||
       (tmpl->ops->nvwrmsnorm  == NULL) )
    return(SUNFALSE);
  return(SUNTRUE);
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

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_SetButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if tables have already been specified, just return */
  if ( (step_mem->Be != NULL) || (step_mem->Bi != NULL) )
    return(ARK_SUCCESS);

  /* initialize table numbers to illegal values */
  etable = itable = -1;

  /**** ImEx methods ****/
  if (step_mem->explicit && step_mem->implicit) {

    switch (step_mem->q) {

    case(2):
    case(3):
      etable = DEFAULT_ARK_ETABLE_3;
      itable = DEFAULT_ARK_ITABLE_3;
      break;
    case(4):
      etable = DEFAULT_ARK_ETABLE_4;
      itable = DEFAULT_ARK_ITABLE_4;
      break;
    case(5):
      etable = DEFAULT_ARK_ETABLE_5;
      itable = DEFAULT_ARK_ITABLE_5;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No ImEx method at requested order, using q=5.");
      etable = DEFAULT_ARK_ETABLE_5;
      itable = DEFAULT_ARK_ITABLE_5;
      break;
    }

  /**** implicit methods ****/
  } else if (step_mem->implicit) {

    switch (step_mem->q) {
    case(2):
      itable = DEFAULT_DIRK_2;
      break;
    case(3):
      itable = DEFAULT_DIRK_3;
      break;
    case(4):
      itable = DEFAULT_DIRK_4;
      break;
    case(5):
      itable = DEFAULT_DIRK_5;;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No implicit method at requested order, using q=5.");
      itable = DEFAULT_DIRK_5;
      break;
    }

  /**** explicit methods ****/
  } else {

    switch (step_mem->q) {
    case(2):
      etable = DEFAULT_ERK_2;
      break;
    case(3):
      etable = DEFAULT_ERK_3;
      break;
    case(4):
      etable = DEFAULT_ERK_4;
      break;
    case(5):
      etable = DEFAULT_ERK_5;
      break;
    case(6):
      etable = DEFAULT_ERK_6;
      break;
    case(7):
    case(8):
      etable = DEFAULT_ERK_8;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No explicit method at requested order, using q=6.");
      etable = DEFAULT_ERK_6;
      break;
    }

  }

  if (etable > -1)
    step_mem->Be = ARKodeLoadButcherTable_ERK(etable);
  if (itable > -1)
    step_mem->Bi = ARKodeLoadButcherTable_DIRK(itable);

  /* set [redundant] ARK stored values for stage numbers and method orders */
  if (step_mem->Be != NULL) {
    step_mem->stages = step_mem->Be->stages;
    step_mem->q = step_mem->Be->q;
    step_mem->p = step_mem->Be->p;
  }
  if (step_mem->Bi != NULL) {
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q = step_mem->Bi->q;
    step_mem->p = step_mem->Bi->p;
  }

  return(ARK_SUCCESS);
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

Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
---------------------------------------------------------------*/
int arkStep_CheckButcherTables(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeARKStepMem step_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check that ERK table is strictly lower triangular */
  if (step_mem->explicit) {
    okay = SUNTRUE;
    for (i=0; i<step_mem->stages; i++)
      for (j=i; j<step_mem->stages; j++)
        if (SUNRabs(step_mem->Be->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ae Butcher table is implicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is implicit and lower triangular */
  if (step_mem->implicit) {
    okay = SUNFALSE;
    for (i=0; i<step_mem->stages; i++)
      if (SUNRabs(step_mem->Bi->A[i][i]) > tol)
        okay = SUNTRUE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table is explicit!");
      return(ARK_ILL_INPUT);
    }

    okay = SUNTRUE;
    for (i=0; i<step_mem->stages; i++)
      for (j=i+1; j<step_mem->stages; j++)
        if (SUNRabs(step_mem->Bi->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table has entries above diagonal!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_Predict

 This routine computes the prediction for a specific internal
 stage solution, storing the result in ypred.  The
 prediction is done using the interpolation structure in
 extrapolation mode, hence stages "far" from the previous time
 interval are predicted using lower order polynomials than the
 "nearby" stages.
---------------------------------------------------------------*/
int arkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess)
{
  int i, retval, ord, jstage, nvec;
  realtype tau;
  realtype h, a0, a1, a2;
  ARKodeARKStepMem step_mem;
  realtype tau_tol = 0.5;
  realtype tau_tol2 = 0.75;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Predict", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Predict",
                    "Interpolation structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* local shortcuts to fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->nst == 0 || ark_mem->resized) {
    N_VScale(ONE, ark_mem->yn, yguess);
    return(ARK_SUCCESS);
  }

  /* set evaluation time tau relative shift from previous successful time */
  tau = step_mem->Bi->c[istage]*ark_mem->h/ark_mem->hold;

  /* use requested predictor formula */
  switch (step_mem->predictor) {

  case 1:

    /***** Interpolatory Predictor 1 -- all to max order *****/
    retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 2:

    /***** Interpolatory Predictor 2 -- decrease order w/ increasing level of extrapolation *****/
    retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 3:

    /***** Cutoff predictor: max order interpolatory output for stages "close"
           to previous step, first-order predictor for subsequent stages *****/
    retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i,
           construct a quadratic Hermite interpolant for prediction; otherwise
           use the trivial predictor.  The actual calculations are performed in
           arkPredict_Bootstrap, but here we need to determine the appropriate
           stage, c_j, to use. *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (step_mem->mass_mem != NULL)  {
      N_VScale(ONE, ark_mem->yn, yguess);
      break;
    }

    /* determine if any previous stages in step meet criteria */
    jstage = -1;
    for (i=0; i<istage; i++)
      jstage = (step_mem->Bi->c[i] != ZERO) ? i : jstage;

    /* if using the trivial predictor, break */
    if (jstage == -1)  break;

    /* find the "optimal" previous stage to use */
    for (i=0; i<istage; i++)
      if ( (step_mem->Bi->c[i] > step_mem->Bi->c[jstage]) &&
           (step_mem->Bi->c[i] != ZERO) )
        jstage = i;

    /* set stage time, stage RHS and interpolation values */
    h = ark_mem->h * step_mem->Bi->c[jstage];
    tau = ark_mem->h * step_mem->Bi->c[istage];
    nvec = 0;
    if (step_mem->implicit) {    /* Implicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fi[jstage];
      nvec += 1;
    }
    if (step_mem->explicit) {    /* Explicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fe[jstage];
      nvec += 1;
    }
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yguess);
    if (retval != ARK_SUCCESS)  break;

    /* call predictor routine */
    retval = arkPredict_Bootstrap(ark_mem, h, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 5:

    /***** Minimal correction predictor: use all previous stage
           information in this step *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (step_mem->mass_mem != NULL)  {
      N_VScale(ONE, ark_mem->yn, yguess);
      break;
    }

    /* set arrays for fused vector operation */
    nvec = 0;
    if (step_mem->explicit) {       /* Explicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[istage][jstage];
        Xvecs[nvec] = step_mem->Fe[jstage];
        nvec += 1;
      }
    }
    if (step_mem->implicit) {      /* Implicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
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
    if (retval != 0) return(ARK_VECTOROP_ERR);

    return(ARK_SUCCESS);
    break;

  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->yn, yguess);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_StageSetup

 This routine sets up the stage data for computing the RK
 residual, along with the step- and method-related factors
 gamma, gammap and gamrat.

 At the ith stage, we compute the residual vector:
    r = -M*zi + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                     + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = -M*zp - M*zc + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                            + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = (-M*zc + gamma*Fi(zi)) + (M*yn - M*zp + data)
 where zi = zp + zc.  In the above form of the residual,
 the first group corresponds to the current solution
 correction, and the second group corresponds to existing data.
 This routine computes this existing data, (M*yn - M*zp + data)
 and stores in step_mem->sdata.
---------------------------------------------------------------*/
int arkStep_StageSetup(ARKodeMem ark_mem)
{
  /* local data */
  ARKodeARKStepMem step_mem;
  int retval, i, j, nvec;
  N_Vector tmp = ark_mem->tempv1;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_StageSetup", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = step_mem->istage;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* If predictor==5, then sdata=0, otherwise set sdata appropriately */
  if ( (step_mem->predictor == 5) && (step_mem->mass_mem == NULL) ) {

    N_VConst(ZERO, step_mem->sdata);

  } else {

    /* Initialize sdata to ycur - zpred (here: ycur = yn and zpred = zp) */
    N_VLinearSum(ONE, ark_mem->yn, -ONE, step_mem->zpred,
                 step_mem->sdata);

    /* If M!=I, replace sdata with M*sdata, so that sdata = M*(yn-zpred) */
    if (step_mem->mass_mem != NULL) {
      N_VScale(ONE, step_mem->sdata, tmp);
      retval = step_mem->mmult((void *) ark_mem, tmp, step_mem->sdata);
      if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    }

    /* Update rhs with prior stage information */
    /*   set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = step_mem->sdata;
    nvec = 1;
    if (step_mem->explicit)    /* Explicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[i][j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
    if (step_mem->implicit)    /* Implicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * step_mem->Bi->A[i][j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }

    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
    if (retval != 0) return(ARK_VECTOROP_ERR);

  }

  /* Update gamma (if the method contains an implicit component) */
  if (step_mem->implicit) {
    step_mem->gamma = ark_mem->h * step_mem->Bi->A[i][i];
    if (ark_mem->firststage)
      step_mem->gammap = step_mem->gamma;
    step_mem->gamrat = (ark_mem->firststage) ?
      ONE : step_mem->gamma / step_mem->gammap;  /* protect x/x != 1.0 */
  }

  /* return with success */
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_Nls

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 It calls one of arkStep_NlsAccelFP, arkStep_Ls or
 arkStep_NlsNewton to do the work.

 Upon entry, the predicted solution is held in step_mem->zpred;
 this array is never changed throughout this routine.  If an
 initial attempt at solving the nonlinear system fails (e.g. due
 to a stale Jacobian), this allows for new attempts at the
 solution.

 Upon a successful solve, the solution is held in ark_mem->ycur.
---------------------------------------------------------------*/
int arkStep_Nls(ARKodeMem ark_mem, int nflag)
{
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Nls", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* call the appropriate solver */

  /*   fixed point */
  if (step_mem->use_fp)
    return(arkStep_NlsAccelFP(ark_mem, nflag));

  /*   linearly implicit (one Newton iteration) */
  if (step_mem->linear)
    return(arkStep_Ls(ark_mem, nflag));

  /*   Newton */
  return(arkStep_NlsNewton(ark_mem, nflag));

}


/*---------------------------------------------------------------
 arkStep_NlsResid

 This routine evaluates the negative nonlinear residual for the
 additive Runge-Kutta method.  It assumes that any data from
 previous time steps/stages is contained in step_mem, and
 merely combines this old data with the current implicit ODE RHS
 vector, fy = fi(t,y), to compute the nonlinear residual r.

 At the ith stage, we compute the residual vector:
    r = -M*zi + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                     + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = -M*zp - M*zc + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                            + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = (-M*zc + gamma*Fi(zi)) + (M*yn - M*zp + data)
 where zi = zp + zc, and where
    zc is stored in the input z
    Fi(zi) is stored in fz
    (M*yn-M*zp+data) is stored in step_mem->sdata,
 so we really just compute
    r = -M*z + gamma*fz + step_mem->sdata

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
int arkStep_NlsResid(ARKodeMem ark_mem, N_Vector z, N_Vector fz,
                     N_Vector r)
{

  /* temporary variables */
  int retval;
  realtype c[3];
  N_Vector X[3];
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_NlsResid", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* put M*z in r */
  if (step_mem->mass_mem != NULL) {
    retval = step_mem->mmult((void *) ark_mem, z, r);
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    X[0] = r;
  } else {
    X[0] = z;
  }

  /* update with My, sdata and gamma*fy */
  c[0] = -ONE;
  c[1] = ONE;
  X[1] = step_mem->sdata;
  c[2] = step_mem->gamma;
  X[2] = fz;
  retval = N_VLinearCombination(3, c, X, r);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_NlsNewton

 This routine handles the Newton iteration for implicit portions
 of the ARK and DIRK methods. It calls lsetup if indicated,
 performs a modified Newton iteration (using a fixed
 Jacobian / precond), and retries a failed attempt at Newton
 iteration if that is indicated.

 Upon entry, the predicted solution is held in step_mem->zpred;
 this array is never changed throughout this routine.  If an
 initial attempt at solving the nonlinear system fails (e.g. due to
 a stale Jacobian), this allows for new attempts that first revert
 ark_mem->ycur back to this initial guess before trying again.

 Upon a successful solve, the solution is held in ark_mem->ycur.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+
   ARK_LSETUP_FAIL    |-> halt the integration
   ARK_LSOLVE_FAIL   -+

   CONV_FAIL        -+
   RHSFUNC_RECVR    -+-> predict again or stop if too many
---------------------------------------------------------------*/
int arkStep_NlsNewton(ARKodeMem ark_mem, int nflag)
{
  N_Vector b, z, zcor;
  int convfail, retval, ier, m, is;
  booleantype callSetup;
  realtype del, delp, dcon;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_NlsNewton", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  b    = ark_mem->tempv1;  /* rename tempv1 as b for readability */
  z    = ark_mem->ycur;    /* rename ycur as z for readability */
  zcor = ark_mem->tempv2;  /* rename tempv2 as zcor for readability */
  is   = step_mem->istage;

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (step_mem->lsetup) {
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (ark_mem->firststage) || (step_mem->msbp < 0) ||
      (ark_mem->nst >= step_mem->nstlp + abs(step_mem->msbp)) ||
      (SUNRabs(step_mem->gamrat-ONE) > step_mem->dgmax);
  } else {
    step_mem->crate = ONE;
    callSetup = SUNFALSE;
  }

  /* Looping point for attempts at solution of the nonlinear system:
     Evaluate f at predicted z, store result in Fi[is].
     Call lsetup if indicated, setting statistics and gamma factors.
     Zero out the correction array (zcor).
     Copy the prediction (zpred) into the output (z).
     Performs the modified Newton iteration using the existing lsetup.
     Repeat process if a recoverable failure occurred (convergence
     failure with stale Jacobian). */
  for(;;) {

    retval = step_mem->fi(ark_mem->tcur, step_mem->zpred,
                             step_mem->Fi[is], ark_mem->user_data);
    step_mem->nfi++;
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

    /* update system matrix/factorization if necessary */
    if (callSetup) {

      /* Solver diagnostics reporting */
      if (ark_mem->report)  fprintf(ark_mem->diagfp, "  lsetup\n");

      /* call lsetup, using zcor, z and b as temporary vectors */
      ier = step_mem->lsetup(ark_mem, convfail, ark_mem->tcur,
                                step_mem->zpred, step_mem->Fi[is],
                                &step_mem->jcur, zcor, z, b);
      step_mem->nsetups++;
      callSetup = SUNFALSE;
      ark_mem->firststage = SUNFALSE;
      step_mem->gamrat = step_mem->crate = ONE;
      step_mem->gammap = step_mem->gamma;
      step_mem->nstlp  = ark_mem->nst;

      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set zcor to zero and load prediction into z vector */
    N_VConst(ZERO, zcor);
    N_VScale(ONE, step_mem->zpred, z);


    /***********************************
     Do the modified Newton iteration:
     - Reuses a single call to ark_mem->lsetup (called by the enclosing loop).
     - Upon entry, the predicted solution is held in step_mem->zpred.
     - Here, zcor accumulates the difference between the predictor
       and the final stage solution.
     - Upon a successful solve, the solution is held in z (ark_mem->ycur).
    */

    /* Initialize temporary variables for use in iteration */
    step_mem->mnewt = m = 0;
    del = delp = ZERO;

    /* Reset the stored residual norm (for iterative linear solvers) */
    step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;

    /* Looping point for Newton iteration */
    for(;;) {

      /* Evaluate the nonlinear system residual, put result into b */
      retval = arkStep_NlsResid(ark_mem, zcor, step_mem->Fi[is], b);
      if (retval != ARK_SUCCESS) {
        ier = ARK_RHSFUNC_FAIL;
        break;
      }

#ifdef DEBUG_OUTPUT
 printf("residual\n");
 N_VPrint_Serial(b);
#endif

      /* Call the lsolve function, overwriting b with solution */
      retval = step_mem->lsolve(ark_mem, b, ark_mem->tcur,
                                   z, step_mem->Fi[is],
                                   step_mem->eRNrm, step_mem->mnewt);
      step_mem->nni++;
      if (retval < 0) {
        ier = ARK_LSOLVE_FAIL;
        break;
      }

#ifdef DEBUG_OUTPUT
 printf("linear system solution\n");
 N_VPrint_Serial(b);
#endif

      /* If lsolve had a recoverable failure and Jacobian data is
         not current, signal to try the solution again */
      if (retval > 0) {
        if ((!step_mem->jcur) && (step_mem->lsetup))
          ier = TRY_AGAIN;
        else
          ier = CONV_FAIL;
        break;
      }

      /* Get WRMS norm of correction; add to zcor and y */
      del = N_VWrmsNorm(b, ark_mem->ewt);
      N_VLinearSum(ONE, zcor, ONE, b, zcor);
      N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

#ifdef DEBUG_OUTPUT
 printf("corrected solution\n");
 N_VPrint_Serial(z);
#endif

      /* Compute the nonlinear error estimate.  If m > 0, an estimate of the convergence
         rate constant is stored in crate, and used in the subsequent estimates */
      if (m > 0)
        step_mem->crate = SUNMAX(step_mem->crdown*step_mem->crate, del/delp);
      dcon = SUNMIN(step_mem->crate, ONE) * del / step_mem->nlscoef;

      /* compute the forcing term for linear solver tolerance */
      step_mem->eRNrm = SUNMIN(step_mem->crate, ONE) * del
                         * RCONST(0.1) * step_mem->nlscoef;
#ifdef FIXED_LIN_TOL
      /* reset if a fixed linear solver tolerance is desired */
      step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;
#endif

#ifdef DEBUG_OUTPUT
 printf("Newton iter %i,  del = %"RSYM",  crate = %"RSYM"\n", m, del, step_mem->crate);
 printf("   dcon = %"RSYM"\n", dcon);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report)
        fprintf(ark_mem->diagfp, "    newt  %i  %"RSYM"  %"RSYM"\n", m, del, dcon);

      if (dcon <= ONE) {
        step_mem->jcur = SUNFALSE;
        ier = ARK_SUCCESS;
        break;
      }

      /* update Newton iteration counter */
      step_mem->mnewt = ++m;

      /* Stop at maxcor iterations or if iteration seems to be diverging.
         If still not converged and Jacobian data is not current, signal
         to try the solution again */
      if ( (m == step_mem->maxcor) ||
           ((m >= 2) && (del > step_mem->rdiv*delp)) ) {
        if ((!step_mem->jcur) && (step_mem->lsetup))
          ier = TRY_AGAIN;
        else
          ier = CONV_FAIL;
        break;
      }

      /* Save norm of correction, evaluate fi, and loop again */
      delp = del;
      retval = step_mem->fi(ark_mem->tcur, z,
                               step_mem->Fi[step_mem->istage],
                               ark_mem->user_data);
      step_mem->nfi++;
      if (retval < 0) {
        ier = ARK_RHSFUNC_FAIL;
        break;
      }
      if (retval > 0) {
        if ((!step_mem->jcur) && (step_mem->lsetup))
          ier = TRY_AGAIN;
        else
          ier = RHSFUNC_RECVR;
        break;
      }

    }
    /* end modified Newton iteration */
    /*********************************/

    /* If there is a convergence failure and the Jacobian-related
       data appears not to be current, loop again with a call to lsetup
       in which convfail=ARK_FAIL_BAD_J.  Otherwise return. */
    if (ier != TRY_AGAIN) return(ier);

    callSetup = SUNTRUE;
    convfail = ARK_FAIL_BAD_J;
  }
}


/*---------------------------------------------------------------
 arkStep_NlsAccelFP

 This routine handles the Anderson-accelerated fixed-point
 iteration for implicit portions of the ARK and DIRK methods. It
 performs an accelerated fixed-point iteration (without need for
 a Jacobian or preconditioner.

 Upon entry, the predicted solution is held in step_mem->zpred.

 Upon a successful solve, the solution is held in ark_mem->ycur.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  ---> halt the integration

   CONV_FAIL         -+
   RHSFUNC_RECVR     -+-> predict again or stop if too many
---------------------------------------------------------------*/
int arkStep_NlsAccelFP(ARKodeMem ark_mem, int nflag)
{
  /* local variables */
  int retval;
  realtype del, delp, dcon, tcur, *R, *gamma;
  long int maa;
  void *udata;
  N_Vector zold, z, zpred, ftemp, fval, tempv;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_NlsAccelFP", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* local shortcut variables */
  tcur = ark_mem->tcur;
  R  = step_mem->fp_mem->R;
  gamma = step_mem->fp_mem->gamma;
  maa = step_mem->fp_mem->m;
  udata = ark_mem->user_data;
  zold  = ark_mem->tempv2;
  z     = ark_mem->ycur;
  zpred = step_mem->zpred;
  ftemp = step_mem->Fi[step_mem->istage];
  fval  = step_mem->fp_mem->fval;
  tempv = ark_mem->tempv1;

  /* Initialize temporary variables for use in iteration */
  del = delp = ZERO;
  step_mem->crate = ONE;

  /* Initialize iteration counter */
  step_mem->mnewt = 0;

  /* Load prediction into z vector, zold */
  N_VScale(ONE, zpred, z);
  N_VScale(ONE, zpred, zold);

  /* Looping point for attempts at solution of the nonlinear system:
       Evaluate f at z, store result in ftemp.
       Performs the accelerated fixed-point iteration.
       Repeat process if a recoverable failure occurred (convergence
          failure with stale Jacobian). */
  for(;;) {

    /* evaluate implicit ODE RHS function, store in ftemp */
    retval = step_mem->fi(tcur, z, ftemp, udata);
    step_mem->nfi++;
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

    /* store difference between current guess and prediction in tempv */
    N_VLinearSum(ONE, z, -ONE, zpred, tempv);

    /* evaluate nonlinear residual, store in fval */
    retval = arkStep_NlsResid(ark_mem, tempv, ftemp, fval);
    if (retval != ARK_SUCCESS) return(ARK_RHSFUNC_FAIL);

    /* convert nonlinear residual result to a fixed-point function result */
    /* NOTE: AS IMPLEMENTED, DOES NOT WORK WITH NON-IDENTITY MASS MATRIX */
    N_VLinearSum(ONE, z, ONE, fval, fval);

    /* perform fixed point update */
    if (maa == 0) {
      /* plain fixed-point solver, copy residual into z */
      N_VScale(ONE, fval, z);
    } else {
      /* Anderson-accelerated solver */
      N_VScale(ONE, zold, z);   /* update guess */
      retval = arkStep_AndersonAcc(ark_mem, fval, tempv, z,  /* accelerate guess */
                                   zold, (int)(step_mem->mnewt), R, gamma);
    }
    step_mem->nni++;

    /* compute the nonlinear error estimate.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the subsequent estimates */
    N_VLinearSum(ONE, z, -ONE, zold, tempv);
    del = N_VWrmsNorm(tempv, ark_mem->ewt);
    if (step_mem->mnewt > 0)
      step_mem->crate = SUNMAX(step_mem->crdown*step_mem->crate, del/delp);
    dcon = SUNMIN(step_mem->crate, ONE) * del / step_mem->nlscoef;

#ifdef DEBUG_OUTPUT
 printf("FP iter %i,  del = %"RSYM",  crate = %"RSYM"\n", step_mem->mnewt, del, step_mem->crate);
 printf("   dcon = %"RSYM"\n", dcon);
 printf("Fixed-point correction:\n");
 N_VPrint_Serial(tempv);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "    fp  %i  %"RSYM"  %"RSYM"\n", step_mem->mnewt, del, dcon);

    /* update iteration counter */
    step_mem->mnewt++;

    if (dcon <= ONE)  return(ARK_SUCCESS);

    /* Stop at maxcor iterations or if iteration seems to be diverging */
    if ((step_mem->mnewt == step_mem->maxcor) ||
        ((step_mem->mnewt >= 2) && (del > step_mem->rdiv*delp)))
      return(CONV_FAIL);

    /* Update current solution guess */
    N_VScale(ONE, z, zold);

    /* Save norm of correction and loop again */
    delp = del;

  }
}


/*---------------------------------------------------------------
 arkStep_AndersonAcc

 This routine computes the Anderson-accelerated fixed point
 iterate itself, as used by the nonlinear solver
 arkStep_NlsAccelFP().

 Upon entry, the predicted solution is held in xold;
 this array is never changed throughout this routine.

 The result of the routine is held in x.

 Possible return values:
   ARK_SUCCESS   ---> successful completion

---------------------------------------------------------------*/
int arkStep_AndersonAcc(ARKodeMem ark_mem, N_Vector gval,
                        N_Vector fv, N_Vector x, N_Vector xold,
                        int iter, realtype *R, realtype *gamma)
{
  /* local variables */
  int      nvec, retval;
  long int i_pt, i, j, lAA, *ipt_map, maa;
  realtype alfa, a, b, temp, c, s;
  N_Vector vtemp2, gold, fold, *df, *dg, *Q;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_AndersonAcc", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* local shortcut variables */
  vtemp2 = ark_mem->ycur;     /* rename ycur as vtemp2 for readability */
  ipt_map = step_mem->fp_mem->imap;
  maa = step_mem->fp_mem->m;
  gold = step_mem->fp_mem->gold;
  fold = step_mem->fp_mem->fold;
  df = step_mem->fp_mem->df;
  dg = step_mem->fp_mem->dg;
  Q = step_mem->fp_mem->q;
  cvals = step_mem->fp_mem->cvals;
  Xvecs = step_mem->fp_mem->Xvecs;

  for (i=0; i<maa; i++)  ipt_map[i]=0;
  i_pt = iter-1 - ((iter-1)/maa)*maa;
  N_VLinearSum(ONE, gval, -ONE, xold, fv);
  if (iter > 0) {
    /* compute dg_new = gval - gval_old*/
    N_VLinearSum(ONE, gval, -ONE, gold, dg[i_pt]);
    /* compute df_new = fval - fval_old */
    N_VLinearSum(ONE, fv, -ONE, fold, df[i_pt]);
  }

  N_VScale(ONE, gval, gold);
  N_VScale(ONE, fv, fold);

  if (iter == 0) {
    N_VScale(ONE, gval, x);
  } else {
    if (iter == 1) {
      R[0] = SUNRsqrt(N_VDotProd(df[i_pt], df[i_pt]));
      alfa = ONE/R[0];
      N_VScale(alfa, df[i_pt], Q[i_pt]);
      ipt_map[0] = 0;
    } else if (iter <= maa) {
      N_VScale(ONE, df[i_pt], vtemp2);
      for (j=0; j<(iter-1); j++) {
        ipt_map[j] = j;
        R[(iter-1)*maa+j] = N_VDotProd(Q[j], vtemp2);
        N_VLinearSum(ONE, vtemp2, -R[(iter-1)*maa+j], Q[j], vtemp2);
      }
      R[(iter-1)*maa+iter-1] = SUNRsqrt(N_VDotProd(vtemp2, vtemp2));
      if (R[(iter-1)*maa+iter-1] == ZERO) {
        N_VScale(ZERO, vtemp2, Q[i_pt]);
      } else {
        N_VScale((ONE/R[(iter-1)*maa+iter-1]), vtemp2, Q[i_pt]);
      }
      ipt_map[iter-1] = iter-1;
    } else {
      /* Delete left-most column vector from QR factorization */
      for (i=0; i < maa-1; i++) {
        a = R[(i+1)*maa + i];
        b = R[(i+1)*maa + i+1];
        temp = SUNRsqrt(a*a + b*b);
        c = a / temp;
        s = b / temp;
        R[(i+1)*maa + i] = temp;
        R[(i+1)*maa + i+1] = 0.0;
        /* OK to re-use temp */
        if (i < maa-1) {
          for (j = i+2; j < maa; j++) {
            a = R[j*maa + i];
            b = R[j*maa + i+1];
            temp = c * a + s * b;
            R[j*maa + i+1] = -s*a + c*b;
            R[j*maa + i] = temp;
          }
        }
        N_VLinearSum(c, Q[i], s, Q[i+1], vtemp2);
        N_VLinearSum(-s, Q[i], c, Q[i+1], Q[i+1]);
        N_VScale(ONE, vtemp2, Q[i]);
      }

      /* Shift R to the left by one. */
      for (i = 1; i < maa; i++) {
        for (j = 0; j < maa-1; j++) {
          R[(i-1)*maa + j] = R[i*maa + j];
        }
      }

      /* Add the new df vector */
      N_VScale(ONE,df[i_pt],vtemp2);
      for (j=0; j < (maa-1); j++) {
        R[(maa-1)*maa+j] = N_VDotProd(Q[j],vtemp2);
        N_VLinearSum(ONE,vtemp2,-R[(maa-1)*maa+j],Q[j],vtemp2);
      }
      R[(maa-1)*maa+maa-1] = SUNRsqrt(N_VDotProd(vtemp2,vtemp2));
      N_VScale((1/R[(maa-1)*maa+maa-1]),vtemp2,Q[maa-1]);

      /* Update the iteration map */
      j = 0;
      for (i=i_pt+1; i<maa; i++)
        ipt_map[j++] = i;
      for (i=0; i<(i_pt+1); i++)
        ipt_map[j++] = i;
    }

    /* Solve least squares problem and update solution */
    lAA = iter;
    if (maa < iter)  lAA = maa;
    retval = N_VDotProdMulti(lAA, fv, Q, gamma);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = gval;
    nvec = 1;
    for (i=lAA-1; i>-1; i--) {
      for (j=i+1; j<lAA; j++)
        gamma[i] = gamma[i] - R[j*maa+i]*gamma[j];
      if (gamma[i] == ZERO) {
        gamma[i] = ZERO;
      } else {
        gamma[i] = gamma[i]/R[i*maa+i];
      }
      cvals[nvec] = -gamma[i];
      Xvecs[nvec] = dg[ipt_map[i]];
      nvec += 1;
    }
    /* update solution */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, x);
    if (retval != 0)  return(ARK_VECTOROP_ERR);
  }

  return 0;
}


/*---------------------------------------------------------------
 arkStep_Ls

 This routine attempts to solve the linear system associated
 with a single implicit step of the ARK method.  This should
 only be called if the user has specified that the implicit
 problem is linear.  In this routine, we assume that the problem
 depends linearly on the solution.  Additionally, if the Jacobian
 is not time dependent we only call lsetup on changes to gamma;
 otherwise we call lsetup at every call.  In all cases, we then
 call the user-specified linear solver (with a tight tolerance)
 to compute the time-evolved solution.

 Upon entry, the predicted solution is held in step_mem->zpred.

 Upon a successful solve, the solution is held in ark_mem->ycur.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+
   ARK_LSETUP_FAIL    |-> halt the integration
   ARK_LSOLVE_FAIL   -+

   RHSFUNC_RECVR     --> predict again or stop if too many
---------------------------------------------------------------*/
int arkStep_Ls(ARKodeMem ark_mem, int nflag)
{
  N_Vector b, z, zcor;
  int convfail, retval, ier, is;
  booleantype callSetup;
  realtype del;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Ls", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  b    = ark_mem->tempv1;   /* rename tempv1 as b for readability */
  z    = ark_mem->ycur;     /* rename ycur as z for readability */
  zcor = ark_mem->tempv2;   /* rename tempv2 as zcor for readability */
  is   = step_mem->istage;

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = (nflag == FIRST_CALL) ? ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (step_mem->lsetup) {
    callSetup = (ark_mem->firststage) ||
      (step_mem->linear_timedep) || (step_mem->msbp < 0) ||
      (SUNRabs(step_mem->gamrat-ONE) > step_mem->dgmax);
  } else {
    callSetup = SUNFALSE;
  }

  /* update implicit RHS, store in ark_mem->Fi[is] */
  retval = step_mem->fi(ark_mem->tcur, step_mem->zpred,
                           step_mem->Fi[is], ark_mem->user_data);
  step_mem->nfi++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* update system matrix/factorization if necessary */
  if (callSetup) {

    /* Solver diagnostics reporting */
    if (ark_mem->report)  fprintf(ark_mem->diagfp, "  lsetup\n");

    /* call lsetup, using zcor, z and b as temporary vectors */
    ier = step_mem->lsetup(ark_mem, convfail, ark_mem->tcur,
                              step_mem->zpred, step_mem->Fi[is],
                              &step_mem->jcur, zcor, z, b);
    step_mem->nsetups++;
    callSetup = SUNFALSE;
    ark_mem->firststage = SUNFALSE;
    step_mem->gamrat = step_mem->crate = ONE;
    step_mem->gammap = step_mem->gamma;
    step_mem->nstlp  = ark_mem->nst;

    /* Return if lsetup failed */
    if (ier < 0) return(ARK_LSETUP_FAIL);
    if (ier > 0) return(CONV_FAIL);
  }

  /* Set zcor to zero and load prediction into y vector */
  N_VConst(ZERO, zcor);
  N_VScale(ONE, step_mem->zpred, z);


  /* Do a single Newton iteration */

  /*   Initialize temporary variables for use in iteration */
  step_mem->mnewt = 0;
  del = ZERO;

  /*   Set the stored residual norm to force an "accurate" initial linear solve */
  step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;

  /*   Evaluate the nonlinear system residual, put result into b */
  retval = arkStep_NlsResid(ark_mem, zcor, step_mem->Fi[is], b);
  if (retval != ARK_SUCCESS)  return (ARK_RHSFUNC_FAIL);

  /*   Call the lsolve function */
  retval = step_mem->lsolve(ark_mem, b, ark_mem->tcur, z, step_mem->Fi[is],
                               step_mem->eRNrm, step_mem->mnewt);
  step_mem->nni++;
  if (retval != 0)  return (ARK_LSOLVE_FAIL);

  /*   Get WRMS norm of correction; add correction to zcor and z */
  del = N_VWrmsNorm(b, ark_mem->ewt);
  N_VLinearSum(ONE, zcor, ONE, b, zcor);
  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  /*   Solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "    newt  %i  %"RSYM"  %g\n", 0, del, 0.0);

  /* clean up and return */
  step_mem->jcur = SUNFALSE;
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_HandleNFlag

 This routine takes action on the return value nflag = *nflagPtr
 returned by arkNls, as follows:

 If arkStep_Nls succeeded in solving the nonlinear system, then
 arkHandleNFlag returns the constant SOLVE_SUCCESS, which tells
 arkStep it is safe to continue with other stage solves, or to
 perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value ARK_LSETUP_FAIL.

 If it failed due to an unrecoverable failure in solve, then we return
 the value ARK_LSOLVE_FAIL.

 If it failed due to an unrecoverable failure in rhs, then we return
 the value ARK_RHSFUNC_FAIL.

 Otherwise, a recoverable failure occurred when solving the
 nonlinear system (arkNls returned nflag == CONV_FAIL or RHSFUNC_RECVR).
 In this case, if using fixed time step sizes, or if ncf is now equal
 to maxncf, or if |h| = hmin, then we return the value ARK_CONV_FAILURE
 (if nflag=CONV_FAIL) or ARK_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
 If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
 PREDICT_AGAIN, telling arkStep to reattempt the step.
---------------------------------------------------------------*/
int arkStep_HandleNFlag(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr)
{
  int nflag;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_HandleNFlag", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  nflag = *nflagPtr;

  if (nflag == ARK_SUCCESS) return(SOLVE_SUCCESS);

  /* The nonlinear soln. failed; increment ncfn */
  step_mem->ncfn++;

  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->fixedstep)    return(ARK_CONV_FAILURE);

  /* Otherwise, access adaptivity structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep", "arkStep_HandleNFlag",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);

  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  hadapt_mem->etamax = ONE;

  /* If we had maxncf failures, or if |h| = hmin,
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((*ncfPtr == step_mem->maxncf) ||
      (SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM)) {
    if (nflag == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);
  }

  /* Reduce step size; return to reattempt the step */
  ark_mem->eta = SUNMAX(hadapt_mem->etacf,
                        ark_mem->hmin / SUNRabs(ark_mem->h));
  ark_mem->h *= ark_mem->eta;
  ark_mem->next_h = ark_mem->h;
  *nflagPtr = PREV_CONV_FAIL;

  return(PREDICT_AGAIN);
}


/*---------------------------------------------------------------
 arkStep_ComputeSolutions

 This routine calculates the final RK solution using the existing
 data.  This solution is placed directly in ark_ycur.  This routine
 also computes the error estimate ||y-ytilde||_WRMS, where ytilde
 is the embedded solution, and the norm weights come from
 ark_ewt.  This norm value is returned.  The vector form of this
 estimated error (y-ytilde) is stored in ark_tempv1, in case the
 calling routine wishes to examine the error locations.

 Note: at this point in the step, the vectors ark_tempv1,
 sdata and zpred may all be used as temporary vectors.
---------------------------------------------------------------*/
int arkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm)
{
  /* local data */
  realtype tend;
  int ier, j, nvec;
  N_Vector y, yerr;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_ComputeSolutions", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;
  tend = ark_mem->tn + ark_mem->h;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsm = ZERO;

  /* Compute updated solution and error estimate based on whether
     a non-identity mass matrix is present */
  if (step_mem->mass_mem != NULL) {   /* M != I */

    /* setup mass matrix, using y, sdata, zpred as temporaries */
    if (step_mem->msetup != NULL)
      if (SUNRabs(step_mem->msetuptime - tend) > FUZZ_FACTOR*ark_mem->uround) {
        ier = step_mem->msetup((void *) ark_mem, y,
                                  step_mem->sdata, step_mem->zpred);
        if (ier != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
        step_mem->msetuptime = tend;
      }

    /* compute y RHS (store in y) */
    /*   set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<step_mem->stages; j++) {
      if (step_mem->explicit) {      /* Explicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit) {      /* Implicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to compute RHS */
    ier = N_VLinearCombination(nvec, cvals, Xvecs, y);
    if (ier != 0) return(ARK_VECTOROP_ERR);

    /* solve for y update (stored in y) */
    ier = step_mem->msolve((void *) ark_mem, y, step_mem->nlscoef);
    if (ier < 0) {
      *dsm = 2.0;         /* indicate too much error, step with smaller step */
      N_VScale(ONE, ark_mem->yn, y);      /* place old solution into y */
      return(CONV_FAIL);
    }

    /* compute y = yn + update */
    N_VLinearSum(ONE, ark_mem->yn, ONE, y, y);


    /* compute yerr (if step adaptivity enabled) */
    if (!ark_mem->fixedstep) {

      /* compute yerr RHS vector (store in yerr) */
      /*   set arrays for fused vector operation */
      nvec = 0;
      for (j=0; j<step_mem->stages; j++) {
        if (step_mem->explicit) {        /* Explicit pieces */
          cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
          Xvecs[nvec] = step_mem->Fe[j];
          nvec += 1;
        }
        if (step_mem->implicit) {        /* Implicit pieces */
          cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
          Xvecs[nvec] = step_mem->Fi[j];
          nvec += 1;
        }
      }

      /*   call fused vector operation to compute yerr RHS */
      ier = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
      if (ier != 0) return(ARK_VECTOROP_ERR);

      /* solve for yerr */
      ier = step_mem->msolve((void *) ark_mem, yerr, step_mem->nlscoef);
      if (ier < 0) {
        *dsm = 2.0;         /* indicate too much error, step with smaller step */
        return(CONV_FAIL);
      }
      /* fill error norm */
      *dsm = N_VWrmsNorm(yerr, ark_mem->ewt);
    }

  } else {                          /* M == I */

    /* Compute time step solution */
    /*   set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = ark_mem->yn;
    nvec = 1;
    for (j=0; j<step_mem->stages; j++) {
      if (step_mem->explicit) {      /* Explicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit) {      /* Implicit pieces */
        cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to do the work */
    ier = N_VLinearCombination(nvec, cvals, Xvecs, y);
    if (ier != 0) return(ARK_VECTOROP_ERR);

    /* Compute yerr (if step adaptivity enabled) */
    if (!ark_mem->fixedstep) {

      /* set arrays for fused vector operation */
      nvec = 0;
      for (j=0; j<step_mem->stages; j++) {
        if (step_mem->explicit) {        /* Explicit pieces */
          cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
          Xvecs[nvec] = step_mem->Fe[j];
          nvec += 1;
        }
        if (step_mem->implicit) {        /* Implicit pieces */
          cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
          Xvecs[nvec] = step_mem->Fi[j];
          nvec += 1;
        }
      }

      /* call fused vector operation to do the work */
      ier = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
      if (ier != 0) return(ARK_VECTOROP_ERR);

      /* fill error norm */
      *dsm = N_VWrmsNorm(yerr, ark_mem->ewt);
    }

  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_DoErrorTest

 This routine performs the local error test for the ARK method.
 The weighted local error norm dsm is passed in, and
 the test dsm ?<= 1 is made.

 If the test passes, arkDoErrorTest returns ARK_SUCCESS.

 If the test fails, we revert to the last successful solution
 time, and:
   - if maxnef error test failures have occurred or if
     SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
   - otherwise: update time step factor eta based on local error
     estimate and reduce h.  Then set *nflagPtr to PREV_ERR_FAIL,
     and return TRY_AGAIN.
---------------------------------------------------------------*/
int arkStep_DoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                        int *nefPtr, realtype dsm)
{
  realtype ehist2, hhist2;
  int retval, k;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_DoErrorTest", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep", "arkDoErrorTest",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */
  if (dsm <= ONE) return(ARK_SUCCESS);

  /* Test failed; increment counters, set nflag */
  (*nefPtr)++;
  step_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  /* At |h| = hmin or maxnef failures, return ARK_ERR_FAILURE */
  if ((SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM) ||
      (*nefPtr == step_mem->maxnef))
    return(ARK_ERR_FAILURE);

  /* Set etamax=1 to prevent step size increase at end of this step */
  hadapt_mem->etamax = ONE;

  /* Temporarily update error history array for recomputation of h */
  ehist2 = hadapt_mem->ehist[2];
  hadapt_mem->ehist[2] = hadapt_mem->ehist[1];
  hadapt_mem->ehist[1] = hadapt_mem->ehist[0];
  hadapt_mem->ehist[0] = dsm*hadapt_mem->bias;

  /* Temporarily update step history array for recomputation of h */
  hhist2 = hadapt_mem->hhist[2];
  hadapt_mem->hhist[2] = hadapt_mem->hhist[1];
  hadapt_mem->hhist[1] = hadapt_mem->hhist[0];
  hadapt_mem->hhist[0] = ark_mem->h;

  /* Compute accuracy-based time step estimate (updated ark_eta) */
  k = (step_mem->hadapt_pq) ? step_mem->q : step_mem->p;
  retval = arkAdapt((void*) ark_mem, step_mem->hadapt_mem, ark_mem->ycur,
                    ark_mem->tcur, ark_mem->h, k, ark_mem->nst);
  if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);

  /* Revert error history array */
  hadapt_mem->ehist[0] = hadapt_mem->ehist[1];
  hadapt_mem->ehist[1] = hadapt_mem->ehist[2];
  hadapt_mem->ehist[2] = ehist2;

  /* Revert step history array */
  hadapt_mem->hhist[0] = hadapt_mem->hhist[1];
  hadapt_mem->hhist[1] = hadapt_mem->hhist[2];
  hadapt_mem->hhist[2] = hhist2;

  /* Enforce failure bounds on eta, update h, and return for retry of step */
  if (*nefPtr >= hadapt_mem->small_nef)
    ark_mem->eta = SUNMIN(ark_mem->eta, hadapt_mem->etamxf);
  ark_mem->h *= ark_mem->eta;
  ark_mem->next_h = ark_mem->h;
  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
 arkStep_PrepareNextStep

 This routine handles ARK-specific updates following a successful
 step: copying the ARK result to the current solution vector,
 updating the error/step history arrays, and setting the
 prospective step size, hprime, for the next step.  Along with
 hprime, it sets the ratio eta=hprime/h.  It also updates other
 state variables related to a change of step size.
---------------------------------------------------------------*/
int arkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval, k;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_PrepareNextStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Update step size and error history arrays */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->ehist[2] = step_mem->hadapt_mem->ehist[1];
    step_mem->hadapt_mem->ehist[1] = step_mem->hadapt_mem->ehist[0];
    step_mem->hadapt_mem->ehist[0] = dsm*step_mem->hadapt_mem->bias;
    step_mem->hadapt_mem->hhist[2] = step_mem->hadapt_mem->hhist[1];
    step_mem->hadapt_mem->hhist[1] = step_mem->hadapt_mem->hhist[0];
    step_mem->hadapt_mem->hhist[0] = ark_mem->h;
  }

  /* If fixed time-stepping requested, defer
     step size changes until next step */
  if (ark_mem->fixedstep){
    ark_mem->hprime = ark_mem->h;
    ark_mem->eta = ONE;
    return(ARK_SUCCESS);
  }

  /* If etamax = 1, defer step size changes until next step,
     and reset etamax */
  if (step_mem->hadapt_mem != NULL)
    if (step_mem->hadapt_mem->etamax == ONE) {
      ark_mem->hprime = ark_mem->h;
      ark_mem->eta = ONE;
      step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->growth;
      return(ARK_SUCCESS);
    }

  /* Adjust ark_eta in arkAdapt */
  if (step_mem->hadapt_mem != NULL) {
    k = (step_mem->hadapt_pq) ? step_mem->q : step_mem->p;
    retval = arkAdapt((void*) ark_mem, step_mem->hadapt_mem,
                      ark_mem->ycur, ark_mem->tn + ark_mem->h,
                      ark_mem->h, k, ark_mem->nst+1);
    if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);
  }

  /* Set hprime value for next step size */
  ark_mem->hprime = ark_mem->h * ark_mem->eta;

  /* Reset growth factor for subsequent time step */
  if (step_mem->hadapt_mem != NULL)
    step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->growth;

  return(ARK_SUCCESS);
}


/*===============================================================
   EOF
===============================================================*/
