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

int ARKStepCreate(void* arkode_mem, ARKRhsFn fe,
                  ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  booleantype nvectorOK;
  int iret;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "ARKStepCreate", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKStepCreate", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Allocate ARKodeARKStepMem structure, and initialize to zero */
  arkstep_mem = NULL;
  arkstep_mem = (ARKodeARKStepMem) malloc(sizeof(struct ARKodeARKStepMemRec));
  if (arkstep_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep",
                    "ARKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  memset(arkstep_mem, 0, sizeof(struct ARKodeARKStepMemRec));

  /* free any existing time step module attached to ARKode */
  if (ark_mem->step_free != NULL)
    ark_mem->step_free(ark_mem);
  
  /* Attach arkstep_mem structure and function pointers to ark_mem */
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
  ark_mem->step_resize         = arkStep_Resize;
  ark_mem->step_fullrhs        = arkStep_FullRHS;
  if (ark_mem->fixedstep) {
    ark_mem->step              = arkStep_FixedStep;
  } else {
    ark_mem->step              = arkStep_AdaptiveStep;
  }
  ark_mem->step_print          = arkStep_PrintMem;
  ark_mem->step_free           = arkStep_Free;
  ark_mem->step_mem            = (void*) arkstep_mem;

  /* Set default values for ARKStep optional inputs */
  iret = ARKStepSetDefaults((void *)ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ARKStep",
                    "ARKStepCreate", 
                    "Error setting default solver options");
    return(iret);
  }
  
  /* Set implicit/explicit problem based on function pointers */
  arkstep_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  arkstep_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (!arkstep_mem->implicit && !arkstep_mem->explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = arkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Allocate the general ARK stepper vectors using y0 as a template */
  /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on 
     (based on the number of ARK stages) */

  /* Clone the input vector to create sdata, zpred */
  if (!arkAllocVec(ark_mem, y0, &(arkstep_mem->sdata)))
    return(ARK_MEM_FAIL);
  if (!arkAllocVec(ark_mem, y0, &(arkstep_mem->zpred)))
    return(ARK_MEM_FAIL);
  
  /* Copy the input parameters into ARKODE state */
  arkstep_mem->fe = fe;
  arkstep_mem->fi = fi;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;
  
  /* Set the linear solver addresses to NULL (we check != NULL later) */
  arkstep_mem->linit  = NULL;
  arkstep_mem->lsetup = NULL;
  arkstep_mem->lsolve = NULL;
  arkstep_mem->lfree  = NULL;
  arkstep_mem->lmem   = NULL;
  arkstep_mem->lsolve_type = -1;

  /* Set the mass matrix solver addresses to NULL */
  arkstep_mem->minit    = NULL;
  arkstep_mem->msetup   = NULL;
  arkstep_mem->mmult    = NULL;
  arkstep_mem->msolve   = NULL;
  arkstep_mem->mfree    = NULL;
  arkstep_mem->mass_mem = NULL;
  arkstep_mem->msetuptime = -RCONST(99999999999.0);
  arkstep_mem->msolve_type = -1;

  /* Allocate step adaptivity structure, set default values, note storage */
  arkstep_mem->hadapt_mem = arkAdaptInit();
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep", "ARKStepCreate", 
                    "Allocation of step adaptivity structure failed");
    return(ARK_MEM_FAIL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;
  
  /* Initialize initial error norm  */
  arkstep_mem->eRNrm = 1.0;
  
  /* Initialize all the counters */
  arkstep_mem->nst_attempts = 0;
  arkstep_mem->nfe          = 0;
  arkstep_mem->nfi          = 0;
  arkstep_mem->ncfn         = 0;
  arkstep_mem->netf         = 0;
  arkstep_mem->nni          = 0;
  arkstep_mem->nsetups      = 0;
  arkstep_mem->nstlp        = 0;

  /* Initialize main ARKode infrastructure */
  iret = arkodeInit(ark_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ARKStep", "ARKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(iret);
  }
  
  return(ARK_SUCCESS);
}



int ARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                  ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  int iret;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "ARKStepReInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access ARKStep memory */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "ARKStepReInit",
                    "Missing time step module memory");
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  iret = arkodeReInit(ark_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ARKStep", "ARKStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(iret);
  }
  
  /* Set implicit/explicit problem based on function pointers */
  arkstep_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  arkstep_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (!arkstep_mem->implicit && !arkstep_mem->explicit) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  arkstep_mem->fe = fe;
  arkstep_mem->fi = fi;

  /* Destroy/Reinitialize time step adaptivity structure (if present) */
  if (arkstep_mem->hadapt_mem != NULL) {
    free(arkstep_mem->hadapt_mem);
    arkstep_mem->hadapt_mem = arkAdaptInit();
    if (arkstep_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep", "ARKStepReInit", 
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }
  /* Initialize initial error norm  */
  arkstep_mem->eRNrm = 1.0;
  
  /* Initialize all the counters */
  arkstep_mem->nst_attempts = 0;
  arkstep_mem->nfe          = 0;
  arkstep_mem->nfi          = 0;
  arkstep_mem->ncfn         = 0;
  arkstep_mem->netf         = 0;
  arkstep_mem->nni          = 0;
  arkstep_mem->nsetups      = 0;
  arkstep_mem->nstlp        = 0;

  return(ARK_SUCCESS);
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AttachLinsol", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AttachLinsol", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* free any existing system solver */
  if (arkstep_mem->lfree != NULL)  arkstep_mem->lfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  arkstep_mem->linit       = linit;
  arkstep_mem->lsetup      = lsetup;
  arkstep_mem->lsolve      = lsolve;
  arkstep_mem->lfree       = lfree;
  arkstep_mem->lmem        = lmem;
  arkstep_mem->lsolve_type = lsolve_type;
  
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AttachMasssol", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AttachMasssol", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* free any existing mass matrix solver */
  if (arkstep_mem->mfree != NULL)  arkstep_mem->mfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  arkstep_mem->minit       = minit;
  arkstep_mem->msetup      = msetup;
  arkstep_mem->mmult       = mmult;
  arkstep_mem->msolve      = msolve;
  arkstep_mem->mfree       = mfree;
  arkstep_mem->mass_mem    = mass_mem;
  arkstep_mem->msolve_type = msolve_type;

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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* nullify the lsetup function pointer */
  arkstep_mem->lsetup = NULL;
}


/*---------------------------------------------------------------
 arkStep_DisableMSetup:

 This routine NULLifies the msetup function pointer in the
 ARKStep module.
---------------------------------------------------------------*/
void arkStep_DisableMSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* nullify the msetup function pointer */
  arkstep_mem->msetup = NULL;
}


/*---------------------------------------------------------------
 arkStep_GetLmem:

 This routine returns the system linear solver interface memory 
 structure, lmem.
---------------------------------------------------------------*/
void* arkStep_GetLmem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure, and return lmem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetLmem", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetLmem", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(arkstep_mem->lmem);
}


/*---------------------------------------------------------------
 arkStep_GetMassMem:

 This routine returns the mass matrix solver interface memory 
 structure, mass_mem.
---------------------------------------------------------------*/
void* arkStep_GetMassMem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure, and return lmem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetMassMem", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetMassMem", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(arkstep_mem->mass_mem);
}


/*---------------------------------------------------------------
 arkStep_GetImplicitRHS:

 This routine returns the implicit RHS function pointer, fi.
---------------------------------------------------------------*/
ARKRhsFn arkStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure, and return fi */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetImplicitRHS", MSG_ARK_NO_MEM);
    return(NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetImplicitRHS", MSG_ARKSTEP_NO_MEM);
    return(NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  return(arkstep_mem->fi);
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure, and return fi */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetGammas", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_GetGammas", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  *gamma  = arkstep_mem->gamma;
  *gamrat = arkstep_mem->gamrat;
  *jcur = &arkstep_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= arkstep_mem->dgmax);

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
  ARKodeARKStepMem arkstep_mem;
  sunindextype Blrw, Bliw;
  int ier, j;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Init", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Init", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* destroy adaptivity structure if fixed-stepping is requested */
  if (ark_mem->fixedstep) 
    if (arkstep_mem->hadapt_mem != NULL) {
      free(arkstep_mem->hadapt_mem);
      arkstep_mem->hadapt_mem = NULL;
    }
  
  /* Set first step growth factor */
  if (arkstep_mem->hadapt_mem != NULL)
    arkstep_mem->hadapt_mem->etamax = arkstep_mem->hadapt_mem->etamx1;

  /* Create Butcher tables (if not already set) */
  ier = arkStep_SetButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", "arkStep_Init", 
                    "Could not create Butcher table(s)");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher tables are OK */
  ier = arkStep_CheckButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "arkStep_Init", "Error in Butcher table(s)");
    return(ARK_ILL_INPUT);
  }

  /* note Butcher table space requirements */
  ButcherTableSpace(arkstep_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;
  ButcherTableSpace(arkstep_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;
  
  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate Fe[0] ... Fe[stages-1] if needed */
  if (arkstep_mem->explicit) {
    if (arkstep_mem->Fe == NULL) 
      arkstep_mem->Fe = (N_Vector *) calloc(arkstep_mem->stages, sizeof(N_Vector));
    for (j=0; j<arkstep_mem->stages; j++) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(arkstep_mem->Fe[j])))
        return(ARK_MEM_FAIL);
    }
    ark_mem->liw += arkstep_mem->stages;  /* pointers */
  }

  /*   Allocate Fi[0] ... Fi[stages-1] if needed */
  if (arkstep_mem->implicit) {
    if (arkstep_mem->Fi == NULL) 
      arkstep_mem->Fi = (N_Vector *) calloc(arkstep_mem->stages, sizeof(N_Vector));
    for (j=0; j<arkstep_mem->stages; j++) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(arkstep_mem->Fi[j])))
        return(ARK_MEM_FAIL);
    }
    ark_mem->liw += arkstep_mem->stages;  /* pointers */
  }

  /* Allocate reusable arrays for fused vector interface */
  if (arkstep_mem->cvals == NULL) {
    arkstep_mem->cvals = (realtype *) calloc(2*arkstep_mem->stages+1, sizeof(realtype));
    if (arkstep_mem->cvals == NULL)  return(ARK_MEM_FAIL);
    ark_mem->lrw += (2*arkstep_mem->stages + 1);
  }
  if (arkstep_mem->Xvecs == NULL) {
    arkstep_mem->Xvecs = (N_Vector *) calloc(2*arkstep_mem->stages+1, sizeof(N_Vector));
    if (arkstep_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
    ark_mem->liw += (2*arkstep_mem->stages + 1);   /* pointers */
  }
  
  /* Check for consistency between linear system modules 
     (if lsolve is direct, msolve needs to match) */
  if (arkstep_mem->mass_mem != NULL) {  /* M != I */
    if (arkstep_mem->lsolve_type != arkstep_mem->msolve_type) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", "arkStep_Init", 
                      "Incompatible linear and mass matrix solvers");
      return(ARK_ILL_INPUT);
    }
  }

  /* Perform mass matrix solver initialization and setup (if applicable) */
  if (arkstep_mem->mass_mem != NULL) {

    /* Call minit (if it exists) */
    if (arkstep_mem->minit != NULL) {
      ier = arkstep_mem->minit((void *) ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKODE::ARKStep", 
                        "arkStep_Init", MSG_ARK_MASSINIT_FAIL);
        return(ARK_MASSINIT_FAIL);
      }
    }

    /* Call msetup (if it exists) -- use ewt, ark_tempv2 and sdata as temp vectors */
    if (arkstep_mem->msetup != NULL) {
      ier = arkstep_mem->msetup((void *) ark_mem, ark_mem->ewt,
                                ark_mem->tempv2, arkstep_mem->sdata);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKODE::ARKStep", 
                        "arkStep_Init", MSG_ARK_MASSSETUP_FAIL);
        return(ARK_MASSSETUP_FAIL);
        arkstep_mem->msetuptime = ark_mem->tcur;
      }
    }
  }

  /* Check if lsolve function exists and call linit (if it exists -- 
     LIKELY NEED TO MODIFY ARK_MEM INPUT) */
  if (arkstep_mem->implicit && !arkstep_mem->use_fp) {
    if (arkstep_mem->lsolve == NULL) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_Init", MSG_ARK_LSOLVE_NULL);
      return(ARK_ILL_INPUT);
    }
    if (arkstep_mem->linit) {
      ier = arkstep_mem->linit(ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE::ARKStep", 
                        "arkStep_Init", MSG_ARK_LINIT_FAIL);
        return(ARK_LINIT_FAIL);
      }
    }
  }

  /* Allocate fixed-point solver memory (LIKELY NEED TO MODIFY ARK_MEM INPUT) */
  if ( (arkstep_mem->use_fp) && (arkstep_mem->fp_mem == NULL) ){
    arkstep_mem->fp_mem = arkAllocFPData(ark_mem, arkstep_mem->fp_m);
    if (arkstep_mem->fp_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep", "arkStep_Init",
                      "Unable to allocate fixed-point solver memory");
      return(ARK_MEM_FAIL);
    }
  }

  /* Allocate interpolation memory (if unallocated, and needed) */
  if ((ark_mem->interp == NULL) && (arkstep_mem->predictor > 0)) {
    ark_mem->interp = arkInterpCreate(ark_mem);
    if (ark_mem->interp == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep", "arkStep_Init",
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
  ARKodeARKStepMem arkstep_mem;
  int retval;
  booleantype recomputeRHS;
  
  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_FullRHS", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_FullRHS", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if the problem involves a non-identity mass matrix and setup is
     required, do so here (use f, ark_tempv2 and sdata as a temporaries) */
  if ( (arkstep_mem->mass_mem != NULL) && (arkstep_mem->msetup != NULL) )
    if (SUNRabs(arkstep_mem->msetuptime - t) > FUZZ_FACTOR*ark_mem->uround) {
      retval = arkstep_mem->msetup((void *) ark_mem, f, ark_mem->tempv2,
                                   arkstep_mem->sdata);
      if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
      arkstep_mem->msetuptime = t;
    }

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vectors fe(t,y) and fi(t,y) in Fe[0] and Fi[0] for 
     possible reuse in the first stage of the subsequent time step */
  case 0:

    /* call fe if the problem has an explicit component */
    if (arkstep_mem->explicit) {
      retval = arkstep_mem->fe(t, y, arkstep_mem->Fe[0], ark_mem->user_data);
      arkstep_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fi if the problem has an implicit component */
    if (arkstep_mem->implicit) {
      retval = arkstep_mem->fi(t, y, arkstep_mem->Fi[0], ark_mem->user_data);
      arkstep_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep", 
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vector(s) into output */
    if (arkstep_mem->explicit && arkstep_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, arkstep_mem->Fi[0], ONE, arkstep_mem->Fe[0], f);
    } else if (arkstep_mem->implicit) {                   /* implicit */
      N_VScale(ONE, arkstep_mem->Fi[0], f);
    } else {                                              /* explicit */
      N_VScale(ONE, arkstep_mem->Fe[0], f);
    }

    break;


  /* Mode 1: called at the end of a successful step
     If the ARK method coefficients support it, we just copy the last stage RHS vectors 
     to fill f instead of calling fe() and fi().  
     Copy the results to Fe[0] and Fi[0] if the ARK coefficients support it. */
  case 1:

    /* determine if explicit/implicit RHS functions need to be recomputed */
    recomputeRHS = SUNFALSE;
    if ( arkstep_mem->explicit && (SUNRabs(arkstep_mem->Be->c[arkstep_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;
    if ( arkstep_mem->implicit && (SUNRabs(arkstep_mem->Bi->c[arkstep_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;

    /* base RHS calls on recomputeRHS argument */
    if (recomputeRHS) {
      
      /* call fe if the problem has an explicit component */
      if (arkstep_mem->explicit) {
        retval = arkstep_mem->fe(t, y, arkstep_mem->Fe[0], ark_mem->user_data);
        arkstep_mem->nfe++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep",
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
          return(ARK_RHSFUNC_FAIL);
        }
      }
      
      /* call fi if the problem has an implicit component */
      if (arkstep_mem->implicit) {
        retval = arkstep_mem->fi(t, y, arkstep_mem->Fi[0], ark_mem->user_data);
        arkstep_mem->nfi++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep", 
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
          return(ARK_RHSFUNC_FAIL);
        }
      }
    } else {
      if (arkstep_mem->explicit) 
        N_VScale(ONE, arkstep_mem->Fe[arkstep_mem->stages-1], arkstep_mem->Fe[0]);
      if (arkstep_mem->implicit) 
        N_VScale(ONE, arkstep_mem->Fi[arkstep_mem->stages-1], arkstep_mem->Fi[0]);
    }

    /* combine RHS vector(s) into output */
    if (arkstep_mem->explicit && arkstep_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, arkstep_mem->Fi[0], ONE, arkstep_mem->Fe[0], f);
    } else if (arkstep_mem->implicit) {                   /* implicit */
      N_VScale(ONE, arkstep_mem->Fi[0], f);
    } else {                                              /* explicit */
      N_VScale(ONE, arkstep_mem->Fe[0], f);
    }
    
    break;

  /*  Mode 2: called for dense output in-between steps
      store the intermediate calculations in such a way as to not 
      interfere with the other two modes */
  default:

    /* call fe if the problem has an explicit component (store in ark_tempv2) */
    if (arkstep_mem->explicit) {
      retval = arkstep_mem->fe(t, y, ark_mem->tempv2, ark_mem->user_data);
      arkstep_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fi if the problem has an implicit component (store in sdata) */
    if (arkstep_mem->implicit) {
      retval = arkstep_mem->fi(t, y, arkstep_mem->sdata, ark_mem->user_data);
      arkstep_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ARKStep", 
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, y);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* combine RHS vector(s) into output */
    if (arkstep_mem->explicit && arkstep_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, arkstep_mem->sdata, ONE, ark_mem->tempv2, f);
    } else if (arkstep_mem->implicit) {                   /* implicit */
      N_VScale(ONE, arkstep_mem->sdata, f);
    } else {                                              /* explicit */
      N_VScale(ONE, ark_mem->tempv2, f);
    }
    
    break;
  }

  
  /* if M != I, then update f = M^{-1}*f 
     MAY WANT TO MOVE THIS TO BE CONDITIONALLY CALLED INSIDE THE PRECEDING CASE BLOCKS */
  if (arkstep_mem->mass_mem != NULL) {
    /* scale RHS */
    N_VScale(ark_mem->h, f, f);
    /* update fnew = M^{-1}*fnew */
    retval = arkstep_mem->msolve((void *) ark_mem, f, arkstep_mem->nlscoef); 
    /* scale result */
    N_VScale(ONE/ark_mem->h, f, f);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKODE::ARKStep", 
                      "arkStep_FullRHS", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_AdaptiveStep:

 This routine serves the primary purpose of the ARKStep module: 
 it performs a single successful embedded ARK step (if possible).  
 Multiple attempts may be taken in this process -- once a step 
 completes with successful (non)linear solves at each stage and 
 passes the error estimate, the routine returns successfully.  
 If it cannot do so, it returns with an appropriate error flag.
---------------------------------------------------------------*/
int arkStep_AdaptiveStep(void* arkode_mem)
{
  realtype dsm;
  int retval, ncf, nef, is, nflag, kflag, eflag;
  booleantype implicit_stage;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AdaptiveStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AdaptiveStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  ncf = nef = 0;
  nflag = FIRST_CALL;
  eflag = ARK_SUCCESS;
  kflag = SOLVE_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {  

    /* increment attempt counter */
    arkstep_mem->nst_attempts++;

    /* Loop over internal stages to the step */
    for (is=0; is<arkstep_mem->stages; is++) {

      /* store current stage index */
      arkstep_mem->istage = is;

      /* Set current stage time(s) */
      if (arkstep_mem->implicit)
        ark_mem->tcur = ark_mem->tn + arkstep_mem->Bi->c[is]*ark_mem->h;
      else
        ark_mem->tcur = ark_mem->tn + arkstep_mem->Be->c[is]*ark_mem->h;
        
#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n", 
         ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif
      
      /* determine whether implicit solve is required */
      implicit_stage = SUNFALSE;
      if (arkstep_mem->implicit) 
        if (SUNRabs(arkstep_mem->Bi->A[is][is]) > TINY)
          implicit_stage = SUNTRUE;

      /* Call predictor for current stage solution (result placed in zpred) */
      if (implicit_stage) {
        eflag = arkStep_Predict(ark_mem, is, arkstep_mem->zpred);
        if (eflag != ARK_SUCCESS)  return (eflag);
      } else {
        N_VScale(ONE, ark_mem->yn, arkstep_mem->zpred);
      }

#ifdef DEBUG_OUTPUT
 printf("predictor:\n");
 N_VPrint_Serial(arkstep_mem->zpred);
#endif
      
      /* Set up data for evaluation of ARK stage residual (data stored in sdata) */
      eflag = arkStep_StageSetup(ark_mem);
      if (eflag != ARK_SUCCESS)  return (eflag);

#ifdef DEBUG_OUTPUT
 printf("rhs data:\n");
 N_VPrint_Serial(arkstep_mem->sdata);
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
        if (arkstep_mem->mass_mem != NULL) {

          /* /\* call msetup if needed -- NEED TEMPVEC1 and TEMPVEC2 *\/ */
          /* if (arkstep_mem->msetup != NULL) { */
          /*   if (SUNRabs(arkstep_mem->msetuptime - ark_mem->tcur) > */
          /*       FUZZ_FACTOR*ark_mem->uround) { */
          /*     eflag = arkstep_mem->msetup((void *) ark_mem, ark_mem->tcur, */
          /*                                 ark_mem->ycur, TEMPVEC1, TEMPVEC2); */
          /*     if (eflag != ARK_SUCCESS) { */
          /*       arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKODE::ARKStep", */
          /*                       "arkStep_AdaptiveStep", MSG_ARK_MASSSSETUP_FAIL); */
          /*       return(ARK_MASSSETUP_FAIL); */
          /*     } */
          /*     arkstep_mem->msetuptime = tend; */
          /*   } */
          /* } */

          /* perform mass matrix solve */
          nflag = arkstep_mem->msolve((void *) ark_mem, arkstep_mem->sdata,
                                      arkstep_mem->nlscoef); 
          
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
          N_VLinearSum(ONE, arkstep_mem->sdata, ONE,
                       arkstep_mem->zpred, ark_mem->ycur);
        }

#ifdef DEBUG_OUTPUT
 printf("explicit solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

      }

      /* successful stage solve */
      /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
      if (arkstep_mem->implicit) {
        retval = arkstep_mem->fi(ark_mem->tcur, ark_mem->ycur,
                                 arkstep_mem->Fi[is], ark_mem->user_data);
        arkstep_mem->nfi++; 
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

      /*    store explicit RHS if necessary
            (already computed at first stage of purely explicit runs) */
      if (arkstep_mem->explicit) {
          retval = arkstep_mem->fe(ark_mem->tn + arkstep_mem->Be->c[is]*ark_mem->h,
                                   ark_mem->ycur, arkstep_mem->Fe[is], ark_mem->user_data);
          arkstep_mem->nfe++;
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
 arkStep_FixedStep:

 This routine serves the primary purpose of the ARKStep module: 
 it performs a single non-embedded ARK step (if possible).  
 Multiple attempts may be taken in this process -- once a step 
 completes with successful (non)linear solves at each stage, 
 the routine returns successfully.  If it cannot do so, it 
 returns with an appropriate error flag.
---------------------------------------------------------------*/
int arkStep_FixedStep(void* arkode_mem)
{
  realtype dsm;
  int retval, ncf, nef, is, nflag, kflag;
  booleantype implicit_stage;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_FixedStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_FixedStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  ncf = nef = 0;
  nflag = FIRST_CALL;
  kflag = SOLVE_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {  

    /* increment attempt counter */
    arkstep_mem->nst_attempts++;

    /* Loop over internal stages to the step */
    for (is=0; is<arkstep_mem->stages; is++) {

      /* store current stage index */
      arkstep_mem->istage = is;

      /* Set current stage time(s) */
      if (arkstep_mem->implicit)
        ark_mem->tcur = ark_mem->tn + arkstep_mem->Bi->c[is]*ark_mem->h;
      else
        ark_mem->tcur = ark_mem->tn + arkstep_mem->Be->c[is]*ark_mem->h;
        
#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n", 
         ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif
      
      /* determine whether implicit solve is required */
      implicit_stage = SUNFALSE;
      if (arkstep_mem->implicit) 
        if (SUNRabs(arkstep_mem->Bi->A[is][is]) > TINY)
          implicit_stage = SUNTRUE;

      /* Call predictor for current stage solution (result placed in zpred) */
      if (implicit_stage) {
        if (arkStep_Predict(ark_mem, is, arkstep_mem->zpred) != ARK_SUCCESS)
          return (ARK_MASSSOLVE_FAIL);
      } else {
        N_VScale(ONE, ark_mem->yn, arkstep_mem->zpred);
      }

#ifdef DEBUG_OUTPUT
 printf("predictor:\n");
 N_VPrint_Serial(arkstep_mem->zpred);
#endif
      
      /* Set up data for evaluation of ARK stage residual (data stored in sdata) */
      if (arkStep_StageSetup(ark_mem) != ARK_SUCCESS)
        return (ARK_MASSMULT_FAIL);

#ifdef DEBUG_OUTPUT
 printf("rhs data:\n");
 N_VPrint_Serial(arkstep_mem->sdata);
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
        if (kflag != SOLVE_SUCCESS) return(kflag);

        /* If h reduced and step needs to be retried, break loop */
        if (kflag == PREDICT_AGAIN) break;

        /* Return if nonlinear solve failed and recovery not possible. */
        if (kflag != SOLVE_SUCCESS) return(kflag);

      /* otherwise no implicit solve is needed */
      } else {

        /* if M!=I, solve with M to compute update (place back in sdata) */
        if (arkstep_mem->mass_mem != NULL) {

          /* /\* call msetup if needed -- NEED TEMPVEC1 and TEMPVEC2 *\/ */
          /* if (arkstep_mem->msetup != NULL) { */
          /*   if (SUNRabs(arkstep_mem->msetuptime - ark_mem->tcur) > */
          /*       FUZZ_FACTOR*ark_mem->uround) { */
          /*     kflag = arkstep_mem->msetup((void *) ark_mem, ark_mem->tcur, */
          /*                                 ark_mem->ycur, TEMPVEC1, TEMPVEC2); */
          /*     if (kflag != ARK_SUCCESS) { */
          /*       arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKODE::ARKStep", */
          /*                       "arkStep_FixedStep", MSG_ARK_MASSSSETUP_FAIL); */
          /*       return(ARK_MASSSETUP_FAIL); */
          /*     } */
          /*     arkstep_mem->msetuptime = tend; */
          /*   } */
          /* } */

          /* perform mass matrix solve */
          nflag = arkstep_mem->msolve((void *) ark_mem, arkstep_mem->sdata,
                                      arkstep_mem->nlscoef); 
          
          /* check for convergence (on failure, h will have been modified) */
          kflag = arkStep_HandleNFlag(ark_mem, &nflag, &ncf);

          /* If fixed time-stepping is used, then anything other than a 
             successful solve must result in an error */
          if (kflag != SOLVE_SUCCESS) return(kflag);

          /* If h reduced and step needs to be retried, break loop */
          if (kflag == PREDICT_AGAIN) break;

          /* Return if solve failed and recovery not possible. */
          if (kflag != SOLVE_SUCCESS) return(kflag);

        /* if M==I, set y to be zpred + RHS data computed in arkStep_StageSetup */
        } else {
          N_VLinearSum(ONE, arkstep_mem->sdata, ONE,
                       arkstep_mem->zpred, ark_mem->ycur);
        }

#ifdef DEBUG_OUTPUT
 printf("explicit solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

      }

      /* successful stage solve */
      /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
      if (arkstep_mem->implicit) {
        retval = arkstep_mem->fi(ark_mem->tcur, ark_mem->ycur,
                                 arkstep_mem->Fi[is], ark_mem->user_data);
        arkstep_mem->nfi++; 
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

      /*    store explicit RHS if necessary
            (already computed at first stage of purely explicit runs) */
      if (arkstep_mem->explicit) {
          retval = arkstep_mem->fe(ark_mem->tn + arkstep_mem->Be->c[is]*ark_mem->h,
                                   ark_mem->ycur, arkstep_mem->Fe[is], ark_mem->user_data);
          arkstep_mem->nfe++;
          if (retval < 0)  return(ARK_RHSFUNC_FAIL);
          if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

    } /* loop over stages */

    /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
    retval = arkStep_ComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);

    break;

  } /* loop over step attempts */


  /* The step has completed successfully, clean up and 
     consider change of step size */
  retval = arkStep_PrepareNextStep(ark_mem, dsm); 
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_Resize:

 This routine resizes the memory within the ARKStep module.
---------------------------------------------------------------*/
int arkStep_Resize(void* arkode_mem, ARKVecResizeFn resize,
                   void* resize_data, sunindextype lrw_diff,
                   sunindextype liw_diff, N_Vector tmpl)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  int ier, i;
  
  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Resize", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Resize", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Resize the sdata and zpred vectors */
  if (arkstep_mem->sdata != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, tmpl, &arkstep_mem->sdata);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  if (arkstep_mem->zpred != NULL) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, tmpl, &arkstep_mem->zpred);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  
  /* Resize the ARKStep vectors */
  /*     Fe */
  if (arkstep_mem->Fe != NULL) {
    for (i=0; i<arkstep_mem->stages; i++) {
      ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, tmpl, &arkstep_mem->Fe[i]);
      if (ier != ARK_SUCCESS)  return(ier);
    }
  }
  /*     Fi */
  if (arkstep_mem->Fi != NULL) {
    for (i=0; i<arkstep_mem->stages; i++) {
      ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, tmpl, &arkstep_mem->Fi[i]);
      if (ier != ARK_SUCCESS)  return(ier);
    }
  }

  /* Resize fixed-point solver memory */
  if ( (arkstep_mem->use_fp) && (arkstep_mem->fp_mem == NULL) ) {
    ier = arkResizeFPData(ark_mem, arkstep_mem->fp_mem, resize, 
                          resize_data, lrw_diff, liw_diff, tmpl);
    if (ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ARKStep",
                      "arkStep_Resize", MSG_ARK_MEM_FAIL);
      return(ARK_MEM_FAIL);
    }
  }

  /* Re-initialize the linear solver (assumes that solver 
     memory has already been resized appropriately) */
  if (arkstep_mem->implicit) {
    if (arkstep_mem->linit != NULL) {
      ier = arkstep_mem->linit(ark_mem);
      if (ier != 0) {
        arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE::ARKStep", 
                        "arkStep_Resize", MSG_ARK_LINIT_FAIL);
        return(ARK_LINIT_FAIL);
      }
    }
  }

  /* Re-initialize the mass matrix solver (assumes that solver 
     memory has already been resized appropriately) */
  if (arkstep_mem->minit != NULL) {
    ier = arkstep_mem->minit((void *) ark_mem);
    if (ier != 0) {
      arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKODE::ARKStep", 
                      "arkStep_Resize", MSG_ARK_MASSINIT_FAIL);
      return(ARK_MASSINIT_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 arkStep_PrintMem:

 This routine outputs the ARKStep structure to a specified file 
 pointer (useful when debugging).
---------------------------------------------------------------*/
void arkStep_PrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  int i;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_PrintMem", MSG_ARK_NO_MEM);
    return;
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_PrintMem", MSG_ARKSTEP_NO_MEM);
    return;
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  /* output integer quantities */
  fprintf(outfile,"ARKStep: q = %i\n", arkstep_mem->q);
  fprintf(outfile,"ARKStep: p = %i\n", arkstep_mem->p);
  fprintf(outfile,"ARKStep: istage = %i\n", arkstep_mem->istage);
  fprintf(outfile,"ARKStep: stages = %i\n", arkstep_mem->stages);
  fprintf(outfile,"ARKStep: mnewt = %i\n", arkstep_mem->mnewt);
  fprintf(outfile,"ARKStep: maxcor = %i\n", arkstep_mem->maxcor);
  fprintf(outfile,"ARKStep: maxnef = %i\n", arkstep_mem->maxnef);
  fprintf(outfile,"ARKStep: maxncf = %i\n", arkstep_mem->maxncf);
  fprintf(outfile,"ARKStep: msbp = %i\n", arkstep_mem->msbp);
  fprintf(outfile,"ARKStep: predictor = %i\n", arkstep_mem->predictor);
  fprintf(outfile,"ARKStep: lsolve_type = %i\n", arkstep_mem->lsolve_type);
  fprintf(outfile,"ARKStep: msolve_type = %i\n", arkstep_mem->msolve_type);

  /* output long integer quantities */
  fprintf(outfile,"ARKStep: nst_attempts = %li\n", arkstep_mem->nst_attempts);
  fprintf(outfile,"ARKStep: nfe = %li\n", arkstep_mem->nfe);
  fprintf(outfile,"ARKStep: nfi = %li\n", arkstep_mem->nfi);
  fprintf(outfile,"ARKStep: ncfn = %li\n", arkstep_mem->ncfn);
  fprintf(outfile,"ARKStep: netf = %li\n", arkstep_mem->netf);
  fprintf(outfile,"ARKStep: nni = %li\n", arkstep_mem->nni);
  fprintf(outfile,"ARKStep: nsetups = %li\n", arkstep_mem->nsetups);
  fprintf(outfile,"ARKStep: use_fp = %i\n", arkstep_mem->use_fp);
  if (arkstep_mem->fp_mem != NULL) {
    fprintf(outfile,"ARKStep: fixed-point solver structure:\n");
    fprintf(outfile, "    m = %li\n", arkstep_mem->fp_mem->m);
    if (arkstep_mem->fp_mem->imap != NULL)
      for (i=0; i<arkstep_mem->fp_mem->m; i++)
        fprintf(outfile, "   imap[%i] = %li\n", i, arkstep_mem->fp_mem->imap[i]);
    if (arkstep_mem->fp_mem->R != NULL) {
      fprintf(outfile, "    R =  ");
      for (i=0; i<arkstep_mem->fp_mem->m*arkstep_mem->fp_mem->m; i++)
        fprintf(outfile, "%"RSYM"  ", arkstep_mem->fp_mem->R[i]);
      fprintf(outfile, "\n");
    }
    if (arkstep_mem->fp_mem->gamma != NULL) {
      fprintf(outfile, "    gamma =  ");
      for (i=0; i<arkstep_mem->fp_mem->m; i++)
        fprintf(outfile, "%"RSYM"  ", arkstep_mem->fp_mem->gamma[i]);
      fprintf(outfile, "\n");
    }
  }
  fprintf(outfile,"ARKStep: nstlp = %li\n", arkstep_mem->nstlp);

  /* output boolean quantities */
  fprintf(outfile,"ARKStep: user_linear = %i\n", arkstep_mem->linear);
  fprintf(outfile,"ARKStep: user_linear_timedep = %i\n", arkstep_mem->linear_timedep);
  fprintf(outfile,"ARKStep: user_explicit = %i\n", arkstep_mem->explicit);
  fprintf(outfile,"ARKStep: user_implicit = %i\n", arkstep_mem->implicit);
  fprintf(outfile,"ARKStep: hadapt_pq = %i\n", arkstep_mem->hadapt_pq);
  fprintf(outfile,"ARKStep: jcur = %i\n", arkstep_mem->jcur);

  /* output realtype quantities */
  if (arkstep_mem->Be != NULL) {
    fprintf(outfile,"ARKStep: explicit Butcher table:\n");
    WriteButcherTable(arkstep_mem->Be, outfile);
  }
  if (arkstep_mem->Bi != NULL) {
    fprintf(outfile,"ARKStep: implicit Butcher table:\n");
    WriteButcherTable(arkstep_mem->Bi, outfile);
  }
  fprintf(outfile,"ARKStep: gamma = %"RSYM"\n", arkstep_mem->gamma);
  fprintf(outfile,"ARKStep: gammap = %"RSYM"\n", arkstep_mem->gammap);
  fprintf(outfile,"ARKStep: gamrat = %"RSYM"\n", arkstep_mem->gamrat);
  fprintf(outfile,"ARKStep: crate = %"RSYM"\n", arkstep_mem->crate);
  fprintf(outfile,"ARKStep: eRNrm = %"RSYM"\n", arkstep_mem->eRNrm);
  fprintf(outfile,"ARKStep: nlscoef = %"RSYM"\n", arkstep_mem->nlscoef);
  if (arkstep_mem->hadapt_mem != NULL) {
    fprintf(outfile,"ARKStep: timestep adaptivity structure:\n");
    arkPrintAdaptMem(arkstep_mem->hadapt_mem, outfile);
  }
  fprintf(outfile,"ARKStep: crdown = %"RSYM"\n", arkstep_mem->crdown);
  fprintf(outfile,"ARKStep: rdiv = %"RSYM"\n", arkstep_mem->rdiv);
  fprintf(outfile,"ARKStep: dgmax = %"RSYM"\n", arkstep_mem->dgmax);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */  
  if (arkstep_mem->sdata != NULL) {
    fprintf(outfile, "ARKStep: sdata:\n");
    N_VPrint_Serial(arkstep_mem->sdata);
  }
  if (arkstep_mem->zpred != NULL) {
    fprintf(outfile, "ARKStep: zpred:\n");
    N_VPrint_Serial(arkstep_mem->zpred);
  }
  if (arkstep_mem->Fe != NULL) 
    for (i=0; i<arkstep_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fe[%i]:\n", i);
      N_VPrint_Serial(arkstep_mem->Fe[i]);
    }
  if (arkstep_mem->Fi != NULL) 
    for (i=0; i<arkstep_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fi[%i]:\n", i);
      N_VPrint_Serial(arkstep_mem->Fi[i]);
    }
#endif
}


/*---------------------------------------------------------------
 arkStep_Free:

 This routine frees all memory internal to the ARKStep module, 
 as well as the module itself.
---------------------------------------------------------------*/
int arkStep_Free(void* arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* trivially return if the ARKStep module is not allocated */
  if (arkode_mem == NULL)  return(ARK_SUCCESS);
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return(ARK_SUCCESS);
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* free the time step adaptivity module */
  if (arkstep_mem->hadapt_mem != NULL) {
    free(arkstep_mem->hadapt_mem);
    arkstep_mem->hadapt_mem = NULL;
    ark_mem->lrw -= ARK_ADAPT_LRW;
    ark_mem->liw -= ARK_ADAPT_LIW;
  }

  /* free the Butcher tables */
  if (arkstep_mem->Be != NULL) {
    ButcherTableSpace(arkstep_mem->Be, &Bliw, &Blrw);
    FreeButcherTable(arkstep_mem->Be);
    arkstep_mem->Be = NULL;
    ark_mem->liw -= Bliw;
    ark_mem->lrw -= Blrw;
  }
  if (arkstep_mem->Bi != NULL) {
    ButcherTableSpace(arkstep_mem->Bi, &Bliw, &Blrw);
    FreeButcherTable(arkstep_mem->Bi);
    arkstep_mem->Bi = NULL;
    ark_mem->liw -= Bliw;
    ark_mem->lrw -= Blrw;
  }

  /* free the fixed-point solver memory */
  if (arkstep_mem->fp_mem != NULL) {
    arkFreeFPData(ark_mem, arkstep_mem->fp_mem);
    arkstep_mem->fp_mem = NULL;
  }

  /* free the linear solver memory */
  if (arkstep_mem->lfree != NULL) {
    arkstep_mem->lfree((void *) ark_mem);
    arkstep_mem->lmem = NULL;
  }

  /* free the mass matrix solver memory */
  if (arkstep_mem->mfree != NULL) {
    arkstep_mem->mfree((void *) ark_mem);
    arkstep_mem->mass_mem = NULL;
  }

  /* free the sdata and zpred vectors */
  if (arkstep_mem->sdata != NULL) {
    arkFreeVec(ark_mem, &arkstep_mem->sdata);
    arkstep_mem->sdata = NULL;
  }
  if (arkstep_mem->zpred != NULL) {
    arkFreeVec(ark_mem, &arkstep_mem->zpred);
    arkstep_mem->zpred = NULL;
  }
  
  /* free the RHS vectors */
  if (arkstep_mem->Fe != NULL) {
    for(j=0; j<arkstep_mem->stages; j++) 
      arkFreeVec(ark_mem, &arkstep_mem->Fe[j]);
    free(arkstep_mem->Fe);
    arkstep_mem->Fe = NULL;
    ark_mem->liw -= arkstep_mem->stages;
  }
  if (arkstep_mem->Fi != NULL) {
    for(j=0; j<arkstep_mem->stages; j++) 
      arkFreeVec(ark_mem, &arkstep_mem->Fi[j]);
    free(arkstep_mem->Fi);
    arkstep_mem->Fi = NULL;
    ark_mem->liw -= arkstep_mem->stages;
  }

  /* free the reusable arrays for fused vector interface */
  if (arkstep_mem->cvals != NULL) {
    free(arkstep_mem->cvals);
    arkstep_mem->cvals = NULL;
    ark_mem->lrw -= (2*arkstep_mem->stages + 1);
  }
  if (arkstep_mem->Xvecs != NULL) {
    free(arkstep_mem->Xvecs);
    arkstep_mem->Xvecs = NULL;
    ark_mem->liw -= (2*arkstep_mem->stages + 1);
  }

  /* free the time stepper module itself */
  free(ark_mem->step_mem);
  ark_mem->step_mem = NULL;
  
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_SetButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if tables have already been specified, just return */
  if ( (arkstep_mem->Be != NULL) || (arkstep_mem->Bi != NULL) )
    return(ARK_SUCCESS);

  /* initialize table numbers to illegal values */
  etable = itable = -1;
  
  /**** ImEx methods ****/
  if (arkstep_mem->explicit && arkstep_mem->implicit) {
    
    switch (arkstep_mem->q) {
      
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
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_SetButcherTables", 
                      "No ImEx method at requested order, using q=5.");
      etable = DEFAULT_ARK_ETABLE_5;
      itable = DEFAULT_ARK_ITABLE_5;
      break;
    }
    
  /**** implicit methods ****/
  } else if (arkstep_mem->implicit) {

    switch (arkstep_mem->q) {
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
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_SetButcherTables", 
                      "No implicit method at requested order, using q=5.");
      itable = DEFAULT_DIRK_5;
      break;
    }

  /**** explicit methods ****/
  } else {

    switch (arkstep_mem->q) {
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
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_SetButcherTables", 
                      "No explicit method at requested order, using q=6.");
      etable = DEFAULT_ERK_6;
      break;
    }

  }

  if (etable > -1)
    arkstep_mem->Be = ARKodeLoadButcherTable(etable);
  if (itable > -1)
    arkstep_mem->Bi = ARKodeLoadButcherTable(itable);

  /* set [redundant] ARK stored values for stage numbers and method orders */
  if (arkstep_mem->Be != NULL) {
    arkstep_mem->stages = arkstep_mem->Be->stages;
    arkstep_mem->q = arkstep_mem->Be->q;
    arkstep_mem->p = arkstep_mem->Be->p;
  }
  if (arkstep_mem->Bi != NULL) {
    arkstep_mem->stages = arkstep_mem->Bi->stages;
    arkstep_mem->q = arkstep_mem->Bi->q;
    arkstep_mem->p = arkstep_mem->Bi->p;
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
  ARKodeARKStepMem arkstep_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_CheckButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check that ERK table is strictly lower triangular */
  if (arkstep_mem->explicit) {
    okay = SUNTRUE;
    for (i=0; i<arkstep_mem->stages; i++)
      for (j=i; j<arkstep_mem->stages; j++)
        if (SUNRabs(arkstep_mem->Be->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_CheckButcherTables",
                      "Ae Butcher table is implicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is implicit and lower triangular */
  if (arkstep_mem->implicit) {
    okay = SUNFALSE;
    for (i=0; i<arkstep_mem->stages; i++)
      if (SUNRabs(arkstep_mem->Bi->A[i][i]) > tol)
        okay = SUNTRUE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table is explicit!");
      return(ARK_ILL_INPUT);
    }

    okay = SUNTRUE;
    for (i=0; i<arkstep_mem->stages; i++)
      for (j=i+1; j<arkstep_mem->stages; j++)
        if (SUNRabs(arkstep_mem->Bi->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table has entries above diagonal!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that method order q > 0 */
  if (arkstep_mem->q < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "arkStep_CheckButcherTables",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((arkstep_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "arkStep_CheckButcherTables",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that stages > 0 */
  if (arkstep_mem->stages < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
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
  ARKodeARKStepMem arkstep_mem;
  realtype tau_tol = 0.5;
  realtype tau_tol2 = 0.75;
  realtype* cvals;
  N_Vector* Xvecs;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (arkstep_mem->predictor > 0)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Predict", 
                    "Interpolation structure is NULL");
    return(ARK_MEM_NULL);
  }
  
  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Predict", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* local shortcuts to fused vector operations */
  cvals = arkstep_mem->cvals;
  Xvecs = arkstep_mem->Xvecs;
  
  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->nst == 0 || ark_mem->resized) {
    N_VScale(ONE, ark_mem->yn, yguess);
    return(ARK_SUCCESS);
  }

  /* set evaluation time tau relative shift from previous successful time */
  tau = arkstep_mem->Bi->c[istage]*ark_mem->h/ark_mem->hold;

  /* use requested predictor formula */
  switch (arkstep_mem->predictor) {

  case 1:

    /***** Interpolatory Predictor 1 -- all to max order *****/
    retval = arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                               0, ark_mem->dense_q, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 2:

    /***** Interpolatory Predictor 2 -- decrease order w/ increasing stage *****/
    if (tau <= tau_tol) {
      ord = 3;
    } else if (tau <= tau_tol2) {
      ord = 2;
    } else {
      ord = 1;
    }
    retval = arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                               0, ord, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 3:

    /***** Cutoff predictor: max order interpolatory output for stages "close" 
           to previous step, first-order predictor for subsequent stages *****/
    if (tau <= tau_tol) {
      retval = arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                                 0, ark_mem->dense_q, yguess);
    } else {
      retval = arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                                 0, 1, yguess);
    }
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i, 
           construct a quadratic Hermite interpolant for prediction; otherwise 
           use the trivial predictor. *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (arkstep_mem->mass_mem != NULL)  {
      N_VScale(ONE, ark_mem->yn, yguess);
      break;
    }

    /* determine if any previous stages in step meet criteria */
    jstage = -1;
    for (i=0; i<istage; i++)
      jstage = (arkstep_mem->Bi->c[i] != ZERO) ? i : jstage;

    /* if using the trivial predictor, break */
    if (jstage == -1)  break;

    /* find the "optimal" previous stage to use */
    for (i=0; i<istage; i++) 
      if ( (arkstep_mem->Bi->c[i] > arkstep_mem->Bi->c[jstage]) &&
           (arkstep_mem->Bi->c[i] != ZERO) )
        jstage = i;

    /* evaluate the quadratic Hermite interpolant for the prediction */
    h = ark_mem->h * arkstep_mem->Bi->c[jstage];
    tau = ark_mem->h * arkstep_mem->Bi->c[istage];
    a0 = ONE;
    a2 = tau*tau/TWO/h;
    a1 = tau - a2;

    /* set arrays for fused vector operation */
    cvals[0] = a0;
    Xvecs[0] = ark_mem->yn;
    cvals[1] = a1;
    Xvecs[1] = ark_mem->interp->fnew;
    nvec = 2;
    if (arkstep_mem->implicit) {    /* Implicit piece */
      cvals[nvec] = a2;
      Xvecs[nvec] = arkstep_mem->Fi[jstage];
      nvec += 1;
    }
    if (arkstep_mem->explicit) {    /* Explicit piece */
      cvals[nvec] = a2;
      Xvecs[nvec] = arkstep_mem->Fe[jstage];
      nvec += 1;
    }

    /* compute predictor */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yguess);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    
    return(ARK_SUCCESS);
    break;

  case 5:

    /***** Minimal correction predictor: use all previous stage 
           information in this step *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (arkstep_mem->mass_mem != NULL)  {
      N_VScale(ONE, ark_mem->yn, yguess);
      break;
    }

    /* set arrays for fused vector operation */
    nvec = 0;
    if (arkstep_mem->explicit) {       /* Explicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
        cvals[nvec] = ark_mem->h * arkstep_mem->Be->A[istage][jstage];
        Xvecs[nvec] = arkstep_mem->Fe[jstage];
        nvec += 1;
      }
    }
    if (arkstep_mem->implicit) {      /* Implicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
        cvals[nvec] = ark_mem->h * arkstep_mem->Bi->A[istage][jstage];
        Xvecs[nvec] = arkstep_mem->Fi[jstage];
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
 and stores in arkstep_mem->sdata.
---------------------------------------------------------------*/
int arkStep_StageSetup(ARKodeMem ark_mem)
{
  /* local data */
  ARKodeARKStepMem arkstep_mem;
  int retval, i, j, nvec;
  N_Vector tmp = ark_mem->tempv1;
  realtype* cvals;
  N_Vector* Xvecs;
  
  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_StageSetup", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = arkstep_mem->istage;

  /* local shortcuts for fused vector operations */
  cvals = arkstep_mem->cvals;
  Xvecs = arkstep_mem->Xvecs;
  
  /* If predictor==5, then sdata=0, otherwise set sdata appropriately */
  if ( (arkstep_mem->predictor == 5) && (arkstep_mem->mass_mem == NULL) ) {

    N_VConst(ZERO, arkstep_mem->sdata);

  } else {

    /* Initialize sdata to ycur - zpred (here: ycur = yn and zpred = zp) */
    N_VLinearSum(ONE, ark_mem->yn, -ONE, arkstep_mem->zpred,
                 arkstep_mem->sdata);

    /* If M!=I, replace sdata with M*sdata, so that sdata = M*(yn-zpred) */
    if (arkstep_mem->mass_mem != NULL) {
      N_VScale(ONE, arkstep_mem->sdata, tmp);
      retval = arkstep_mem->mmult((void *) ark_mem, tmp, arkstep_mem->sdata);
      if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    }

    /* Update rhs with prior stage information */
    /*   set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = arkstep_mem->sdata;
    nvec = 1;
    if (arkstep_mem->explicit)    /* Explicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * arkstep_mem->Be->A[i][j];
        Xvecs[nvec] = arkstep_mem->Fe[j];
        nvec += 1;
      }
    if (arkstep_mem->implicit)    /* Implicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * arkstep_mem->Bi->A[i][j];
        Xvecs[nvec] = arkstep_mem->Fi[j];
        nvec += 1;
      }
 
    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, arkstep_mem->sdata);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    
  }

  /* Update gamma (if the method contains an implicit component) */
  if (arkstep_mem->implicit) {
    arkstep_mem->gamma = ark_mem->h * arkstep_mem->Bi->A[i][i];
    if (ark_mem->firststage)  
      arkstep_mem->gammap = arkstep_mem->gamma;
    arkstep_mem->gamrat = (ark_mem->firststage) ? 
      ONE : arkstep_mem->gamma / arkstep_mem->gammap;  /* protect x/x != 1.0 */
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

 Upon entry, the predicted solution is held in arkstep_mem->zpred; 
 this array is never changed throughout this routine.  If an 
 initial attempt at solving the nonlinear system fails (e.g. due 
 to a stale Jacobian), this allows for new attempts at the 
 solution.

 Upon a successful solve, the solution is held in ark_mem->ycur.
---------------------------------------------------------------*/
int arkStep_Nls(ARKodeMem ark_mem, int nflag)
{
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Nls", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* call the appropriate solver */

  /*   fixed point */
  if (arkstep_mem->use_fp)
    return(arkStep_NlsAccelFP(ark_mem, nflag));

  /*   linearly implicit (one Newton iteration) */
  if (arkstep_mem->linear)
    return(arkStep_Ls(ark_mem, nflag));
  
  /*   Newton */
  return(arkStep_NlsNewton(ark_mem, nflag));

}


/*---------------------------------------------------------------
 arkStep_NlsResid

 This routine evaluates the negative nonlinear residual for the 
 additive Runge-Kutta method.  It assumes that any data from 
 previous time steps/stages is contained in arkstep_mem, and 
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
    (M*yn-M*zp+data) is stored in arkstep_mem->sdata,
 so we really just compute
    r = -M*z + gamma*fz + arkstep_mem->sdata

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
int arkStep_NlsResid(ARKodeMem ark_mem, N_Vector z, N_Vector fz,
                     N_Vector r)
{

  /* temporary variables */
  int retval;
  realtype c[3];
  N_Vector X[3];
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_NlsResid", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* put M*z in r */
  if (arkstep_mem->mass_mem != NULL) {
    retval = arkstep_mem->mmult((void *) ark_mem, z, r);
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    X[0] = r;
  } else {
    X[0] = z;
  }

  /* update with My, sdata and gamma*fy */
  c[0] = -ONE;
  c[1] = ONE;
  X[1] = arkstep_mem->sdata;
  c[2] = arkstep_mem->gamma;
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

 Upon entry, the predicted solution is held in arkstep_mem->zpred; 
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_NlsNewton", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  b    = ark_mem->tempv1;  /* rename tempv1 as b for readability */
  z    = ark_mem->ycur;    /* rename ycur as z for readability */
  zcor = ark_mem->tempv2;  /* rename tempv2 as zcor for readability */
  is   = arkstep_mem->istage;

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (arkstep_mem->lsetup) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (ark_mem->firststage) || (arkstep_mem->msbp < 0) ||
      (ark_mem->nst >= arkstep_mem->nstlp + abs(arkstep_mem->msbp)) || 
      (SUNRabs(arkstep_mem->gamrat-ONE) > arkstep_mem->dgmax);
  } else {  
    arkstep_mem->crate = ONE;
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

    retval = arkstep_mem->fi(ark_mem->tcur, arkstep_mem->zpred, 
                             arkstep_mem->Fi[is], ark_mem->user_data);
    arkstep_mem->nfi++; 
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);
    
    /* update system matrix/factorization if necessary */
    if (callSetup) {

      /* Solver diagnostics reporting */
      if (ark_mem->report)  fprintf(ark_mem->diagfp, "  lsetup\n");

      /* call lsetup, using zcor, z and b as temporary vectors */
      ier = arkstep_mem->lsetup(ark_mem, convfail, ark_mem->tcur,
                                arkstep_mem->zpred, arkstep_mem->Fi[is],
                                &arkstep_mem->jcur, zcor, z, b);
      arkstep_mem->nsetups++;
      callSetup = SUNFALSE;
      ark_mem->firststage = SUNFALSE;
      arkstep_mem->gamrat = arkstep_mem->crate = ONE; 
      arkstep_mem->gammap = arkstep_mem->gamma;
      arkstep_mem->nstlp  = ark_mem->nst;

      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set zcor to zero and load prediction into z vector */
    N_VConst(ZERO, zcor);
    N_VScale(ONE, arkstep_mem->zpred, z);


    /***********************************
     Do the modified Newton iteration:
     - Reuses a single call to ark_mem->lsetup (called by the enclosing loop). 
     - Upon entry, the predicted solution is held in arkstep_mem->zpred.
     - Here, zcor accumulates the difference between the predictor 
       and the final stage solution.
     - Upon a successful solve, the solution is held in z (ark_mem->ycur). 
    */

    /* Initialize temporary variables for use in iteration */
    arkstep_mem->mnewt = m = 0;
    del = delp = ZERO;

    /* Reset the stored residual norm (for iterative linear solvers) */
    arkstep_mem->eRNrm = RCONST(0.1) * arkstep_mem->nlscoef;

    /* Looping point for Newton iteration */
    for(;;) {

      /* Evaluate the nonlinear system residual, put result into b */
      retval = arkStep_NlsResid(ark_mem, zcor, arkstep_mem->Fi[is], b);
      if (retval != ARK_SUCCESS) {
        ier = ARK_RHSFUNC_FAIL;
        break;
      }

#ifdef DEBUG_OUTPUT
 printf("residual\n");
 N_VPrint_Serial(b);
#endif

      /* Call the lsolve function, overwriting b with solution */
      retval = arkstep_mem->lsolve(ark_mem, b, ark_mem->tcur,
                                   z, arkstep_mem->Fi[is],
                                   arkstep_mem->eRNrm, arkstep_mem->mnewt); 
      arkstep_mem->nni++;
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
        if ((!arkstep_mem->jcur) && (arkstep_mem->lsetup)) 
          ier = TRY_AGAIN;
        else 
          ier = CONV_FAIL;
        break;
      }

      /* Get WRMS norm of correction; add to zcor and y */
      del = N_VWrmsNorm(b, ark_mem->ewt);
      N_VLinearSum(ONE, zcor, ONE, b, zcor);
      N_VLinearSum(ONE, arkstep_mem->zpred, ONE, zcor, z);

#ifdef DEBUG_OUTPUT
 printf("corrected solution\n");
 N_VPrint_Serial(z);
#endif
      
      /* Compute the nonlinear error estimate.  If m > 0, an estimate of the convergence
         rate constant is stored in crate, and used in the subsequent estimates */
      if (m > 0) 
        arkstep_mem->crate = SUNMAX(arkstep_mem->crdown*arkstep_mem->crate, del/delp);
      dcon = SUNMIN(arkstep_mem->crate, ONE) * del / arkstep_mem->nlscoef;

      /* compute the forcing term for linear solver tolerance */
      arkstep_mem->eRNrm = SUNMIN(arkstep_mem->crate, ONE) * del
                         * RCONST(0.1) * arkstep_mem->nlscoef;
#ifdef FIXED_LIN_TOL
      /* reset if a fixed linear solver tolerance is desired */
      arkstep_mem->eRNrm = RCONST(0.1) * arkstep_mem->nlscoef;
#endif

#ifdef DEBUG_OUTPUT
 printf("Newton iter %i,  del = %"RSYM",  crate = %"RSYM"\n", m, del, arkstep_mem->crate);
 printf("   dcon = %"RSYM"\n", dcon);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report) 
        fprintf(ark_mem->diagfp, "    newt  %i  %"RSYM"  %"RSYM"\n", m, del, dcon);
    
      if (dcon <= ONE) {
        arkstep_mem->jcur = SUNFALSE;
        ier = ARK_SUCCESS;
        break;
      }

      /* update Newton iteration counter */
      arkstep_mem->mnewt = ++m;
    
      /* Stop at maxcor iterations or if iteration seems to be diverging.
         If still not converged and Jacobian data is not current, signal 
         to try the solution again */
      if ( (m == arkstep_mem->maxcor) || 
           ((m >= 2) && (del > arkstep_mem->rdiv*delp)) ) {
        if ((!arkstep_mem->jcur) && (arkstep_mem->lsetup)) 
          ier = TRY_AGAIN;
        else
          ier = CONV_FAIL;
        break;
      }
    
      /* Save norm of correction, evaluate fi, and loop again */
      delp = del;
      retval = arkstep_mem->fi(ark_mem->tcur, z, 
                               arkstep_mem->Fi[arkstep_mem->istage], 
                               ark_mem->user_data);
      arkstep_mem->nfi++;
      if (retval < 0) {
        ier = ARK_RHSFUNC_FAIL;
        break;
      }
      if (retval > 0) {
        if ((!arkstep_mem->jcur) && (arkstep_mem->lsetup)) 
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

 Upon entry, the predicted solution is held in arkstep_mem->zpred.

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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_NlsAccelFP", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  /* local shortcut variables */
  tcur = ark_mem->tcur;
  R  = arkstep_mem->fp_mem->R;
  gamma = arkstep_mem->fp_mem->gamma;
  maa = arkstep_mem->fp_mem->m;
  udata = ark_mem->user_data;
  zold  = ark_mem->tempv2;
  z     = ark_mem->ycur;
  zpred = arkstep_mem->zpred;
  ftemp = arkstep_mem->Fi[arkstep_mem->istage];
  fval  = arkstep_mem->fp_mem->fval;
  tempv = ark_mem->tempv1;

  /* Initialize temporary variables for use in iteration */
  del = delp = ZERO;
  arkstep_mem->crate = ONE;

  /* Initialize iteration counter */
  arkstep_mem->mnewt = 0;

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
    retval = arkstep_mem->fi(tcur, z, ftemp, udata);
    arkstep_mem->nfi++;
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
                                   zold, (int)(arkstep_mem->mnewt), R, gamma);
    }
    arkstep_mem->nni++;

    /* compute the nonlinear error estimate.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the subsequent estimates */
    N_VLinearSum(ONE, z, -ONE, zold, tempv);
    del = N_VWrmsNorm(tempv, ark_mem->ewt);
    if (arkstep_mem->mnewt > 0)
      arkstep_mem->crate = SUNMAX(arkstep_mem->crdown*arkstep_mem->crate, del/delp);
    dcon = SUNMIN(arkstep_mem->crate, ONE) * del / arkstep_mem->nlscoef;

#ifdef DEBUG_OUTPUT
 printf("FP iter %i,  del = %"RSYM",  crate = %"RSYM"\n", arkstep_mem->mnewt, del, arkstep_mem->crate);
 printf("   dcon = %"RSYM"\n", dcon);
 printf("Fixed-point correction:\n");
 N_VPrint_Serial(tempv);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "    fp  %i  %"RSYM"  %"RSYM"\n", arkstep_mem->mnewt, del, dcon);

    /* update iteration counter */
    arkstep_mem->mnewt++;

    if (dcon <= ONE)  return(ARK_SUCCESS);

    /* Stop at maxcor iterations or if iteration seems to be diverging */
    if ((arkstep_mem->mnewt == arkstep_mem->maxcor) ||
        ((arkstep_mem->mnewt >= 2) && (del > arkstep_mem->rdiv*delp))) 
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_AndersonAcc", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* local shortcut variables */
  vtemp2 = ark_mem->ycur;     /* rename ycur as vtemp2 for readability */
  ipt_map = arkstep_mem->fp_mem->imap;
  maa = arkstep_mem->fp_mem->m;
  gold = arkstep_mem->fp_mem->gold;
  fold = arkstep_mem->fp_mem->fold;
  df = arkstep_mem->fp_mem->df;
  dg = arkstep_mem->fp_mem->dg;
  Q = arkstep_mem->fp_mem->q;
  cvals = arkstep_mem->fp_mem->cvals;
  Xvecs = arkstep_mem->fp_mem->Xvecs;
  
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

 Upon entry, the predicted solution is held in arkstep_mem->zpred.

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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_Ls", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  b    = ark_mem->tempv1;   /* rename tempv1 as b for readability */
  z    = ark_mem->ycur;     /* rename ycur as z for readability */
  zcor = ark_mem->tempv2;   /* rename tempv2 as zcor for readability */
  is   = arkstep_mem->istage;

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = (nflag == FIRST_CALL) ? ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (arkstep_mem->lsetup) {    
    callSetup = (ark_mem->firststage) || 
      (arkstep_mem->linear_timedep) || (arkstep_mem->msbp < 0) ||
      (SUNRabs(arkstep_mem->gamrat-ONE) > arkstep_mem->dgmax);
  } else {  
    callSetup = SUNFALSE;
  }
  
  /* update implicit RHS, store in ark_mem->Fi[is] */
  retval = arkstep_mem->fi(ark_mem->tcur, arkstep_mem->zpred, 
                           arkstep_mem->Fi[is], ark_mem->user_data);
  arkstep_mem->nfi++; 
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);
  
  /* update system matrix/factorization if necessary */
  if (callSetup) {

    /* Solver diagnostics reporting */
    if (ark_mem->report)  fprintf(ark_mem->diagfp, "  lsetup\n");

    /* call lsetup, using zcor, z and b as temporary vectors */
    ier = arkstep_mem->lsetup(ark_mem, convfail, ark_mem->tcur,
                              arkstep_mem->zpred, arkstep_mem->Fi[is],
                              &arkstep_mem->jcur, zcor, z, b);
    arkstep_mem->nsetups++;
    callSetup = SUNFALSE;
    ark_mem->firststage = SUNFALSE;
    arkstep_mem->gamrat = arkstep_mem->crate = ONE; 
    arkstep_mem->gammap = arkstep_mem->gamma;
    arkstep_mem->nstlp  = ark_mem->nst;

    /* Return if lsetup failed */
    if (ier < 0) return(ARK_LSETUP_FAIL);
    if (ier > 0) return(CONV_FAIL);
  }

  /* Set zcor to zero and load prediction into y vector */
  N_VConst(ZERO, zcor);
  N_VScale(ONE, arkstep_mem->zpred, z);
  

  /* Do a single Newton iteration */

  /*   Initialize temporary variables for use in iteration */
  arkstep_mem->mnewt = 0;
  del = ZERO;

  /*   Set the stored residual norm to force an "accurate" initial linear solve */
  arkstep_mem->eRNrm = RCONST(0.1) * arkstep_mem->nlscoef;

  /*   Evaluate the nonlinear system residual, put result into b */
  retval = arkStep_NlsResid(ark_mem, zcor, arkstep_mem->Fi[is], b);
  if (retval != ARK_SUCCESS)  return (ARK_RHSFUNC_FAIL);

  /*   Call the lsolve function */
  retval = arkstep_mem->lsolve(ark_mem, b, ark_mem->tcur, z, arkstep_mem->Fi[is], 
                               arkstep_mem->eRNrm, arkstep_mem->mnewt); 
  arkstep_mem->nni++;
  if (retval != 0)  return (ARK_LSOLVE_FAIL);
    
  /*   Get WRMS norm of correction; add correction to zcor and z */
  del = N_VWrmsNorm(b, ark_mem->ewt);
  N_VLinearSum(ONE, zcor, ONE, b, zcor);
  N_VLinearSum(ONE, arkstep_mem->zpred, ONE, zcor, z);

  /*   Solver diagnostics reporting */
  if (ark_mem->report) 
    fprintf(ark_mem->diagfp, "    newt  %i  %"RSYM"  %g\n", 0, del, 0.0);

  /* clean up and return */ 
  arkstep_mem->jcur = SUNFALSE;
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_HandleNFlag", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  nflag = *nflagPtr;
  
  if (nflag == ARK_SUCCESS) return(SOLVE_SUCCESS);

  /* The nonlinear soln. failed; increment ncfn */
  arkstep_mem->ncfn++;
  
  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->fixedstep)    return(ARK_CONV_FAILURE);

  /* Otherwise, access adaptivity structure */
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", "arkStep_HandleNFlag",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);
  
  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  hadapt_mem->etamax = ONE;

  /* If we had maxncf failures, or if |h| = hmin,
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((*ncfPtr == arkstep_mem->maxncf) ||
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_ComputeSolutions", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;
  tend = ark_mem->tn + ark_mem->h;

  /* local shortcuts for fused vector operations */
  cvals = arkstep_mem->cvals;
  Xvecs = arkstep_mem->Xvecs;
  
  /* initialize output */
  *dsm = ZERO;

  /* Compute updated solution and error estimate based on whether
     a non-identity mass matrix is present */
  if (arkstep_mem->mass_mem != NULL) {   /* M != I */

    /* setup mass matrix, using y, sdata, zpred as temporaries */
    if (arkstep_mem->msetup != NULL)
      if (SUNRabs(arkstep_mem->msetuptime - tend) > FUZZ_FACTOR*ark_mem->uround) {
        ier = arkstep_mem->msetup((void *) ark_mem, y,
                                  arkstep_mem->sdata, arkstep_mem->zpred);
        if (ier != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
        arkstep_mem->msetuptime = tend;
      }

    /* compute y RHS (store in y) */
    /*   set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<arkstep_mem->stages; j++) {
      if (arkstep_mem->explicit) {      /* Explicit pieces */
        cvals[nvec] = ark_mem->h * arkstep_mem->Be->b[j];
        Xvecs[nvec] = arkstep_mem->Fe[j];
        nvec += 1;
      }
      if (arkstep_mem->implicit) {      /* Implicit pieces */
        cvals[nvec] = ark_mem->h * arkstep_mem->Bi->b[j];
        Xvecs[nvec] = arkstep_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to compute RHS */
    ier = N_VLinearCombination(nvec, cvals, Xvecs, y);
    if (ier != 0) return(ARK_VECTOROP_ERR);
    
    /* solve for y update (stored in y) */
    ier = arkstep_mem->msolve((void *) ark_mem, y, arkstep_mem->nlscoef); 
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
      for (j=0; j<arkstep_mem->stages; j++) {
        if (arkstep_mem->explicit) {        /* Explicit pieces */
          cvals[nvec] = ark_mem->h * (arkstep_mem->Be->b[j] - arkstep_mem->Be->d[j]);
          Xvecs[nvec] = arkstep_mem->Fe[j];
          nvec += 1;
        }
        if (arkstep_mem->implicit) {        /* Implicit pieces */
          cvals[nvec] = ark_mem->h * (arkstep_mem->Bi->b[j] - arkstep_mem->Bi->d[j]);
          Xvecs[nvec] = arkstep_mem->Fi[j];
          nvec += 1;
        }
      }
 
      /*   call fused vector operation to compute yerr RHS */
      ier = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
      if (ier != 0) return(ARK_VECTOROP_ERR);
      
      /* solve for yerr */
      ier = arkstep_mem->msolve((void *) ark_mem, yerr, arkstep_mem->nlscoef); 
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
    for (j=0; j<arkstep_mem->stages; j++) {
      if (arkstep_mem->explicit) {      /* Explicit pieces */
        cvals[nvec] = ark_mem->h * arkstep_mem->Be->b[j];
        Xvecs[nvec] = arkstep_mem->Fe[j];
        nvec += 1;
      }
      if (arkstep_mem->implicit) {      /* Implicit pieces */
        cvals[nvec] = ark_mem->h * arkstep_mem->Bi->b[j];
        Xvecs[nvec] = arkstep_mem->Fi[j];
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
      for (j=0; j<arkstep_mem->stages; j++) {
        if (arkstep_mem->explicit) {        /* Explicit pieces */
          cvals[nvec] = ark_mem->h * (arkstep_mem->Be->b[j] - arkstep_mem->Be->d[j]);
          Xvecs[nvec] = arkstep_mem->Fe[j];
          nvec += 1;
        }
        if (arkstep_mem->implicit) {        /* Implicit pieces */
          cvals[nvec] = ark_mem->h * (arkstep_mem->Bi->b[j] - arkstep_mem->Bi->d[j]);
          Xvecs[nvec] = arkstep_mem->Fi[j];
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_DoErrorTest", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", "arkDoErrorTest",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters, set nflag */
  (*nefPtr)++;
  arkstep_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  /* At |h| = hmin or maxnef failures, return ARK_ERR_FAILURE */
  if ((SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM) ||
      (*nefPtr == arkstep_mem->maxnef))
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
  k = (arkstep_mem->hadapt_pq) ? arkstep_mem->q : arkstep_mem->p;
  retval = arkAdapt((void*) ark_mem, arkstep_mem->hadapt_mem, ark_mem->ycur,
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
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep",
                    "arkStep_PrepareNextStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Update step size and error history arrays */
  if (arkstep_mem->hadapt_mem != NULL) {
    arkstep_mem->hadapt_mem->ehist[2] = arkstep_mem->hadapt_mem->ehist[1];
    arkstep_mem->hadapt_mem->ehist[1] = arkstep_mem->hadapt_mem->ehist[0];
    arkstep_mem->hadapt_mem->ehist[0] = dsm*arkstep_mem->hadapt_mem->bias;
    arkstep_mem->hadapt_mem->hhist[2] = arkstep_mem->hadapt_mem->hhist[1];
    arkstep_mem->hadapt_mem->hhist[1] = arkstep_mem->hadapt_mem->hhist[0];
    arkstep_mem->hadapt_mem->hhist[0] = ark_mem->h;
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
  if (arkstep_mem->hadapt_mem != NULL)
    if (arkstep_mem->hadapt_mem->etamax == ONE) {
      ark_mem->hprime = ark_mem->h;
      ark_mem->eta = ONE;
      arkstep_mem->hadapt_mem->etamax = arkstep_mem->hadapt_mem->growth;
      return(ARK_SUCCESS);
    }

  /* Adjust ark_eta in arkAdapt */
  if (arkstep_mem->hadapt_mem != NULL) {
    k = (arkstep_mem->hadapt_pq) ? arkstep_mem->q : arkstep_mem->p;
    retval = arkAdapt((void*) ark_mem, arkstep_mem->hadapt_mem,
                      ark_mem->ycur, ark_mem->tn + ark_mem->h,
                      ark_mem->h, k, ark_mem->nst+1);
    if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);
  }
  
  /* Set hprime value for next step size */
  ark_mem->hprime = ark_mem->h * ark_mem->eta;
  
  /* Reset growth factor for subsequent time step */
  if (arkstep_mem->hadapt_mem != NULL)
    arkstep_mem->hadapt_mem->etamax = arkstep_mem->hadapt_mem->growth;

  return(ARK_SUCCESS);
}


/*===============================================================
   EOF
===============================================================*/
