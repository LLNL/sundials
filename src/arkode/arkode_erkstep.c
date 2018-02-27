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
#include "arkode_erkstep_impl.h"
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



/*===============================================================
  ERKStep Exported functions -- Required
  ===============================================================*/

int ERKStepCreate(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  booleantype nvectorOK;
  int iret;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepCreate", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "ERKStepCreate", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Allocate ARKodeERKStepMem structure, and initialize to zero */
  erkstep_mem = NULL;
  erkstep_mem = (ARKodeERKStepMem) malloc(sizeof(struct ARKodeERKStepMemRec));
  if (erkstep_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep",
                    "ERKStepCreate", MSGARK_ARKMEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  memset(erkstep_mem, 0, sizeof(struct ARKodeERKStepMemRec));

  /* free any existing time step module attached to ARKode */
  if (ark_mem->ark_step_free != NULL)
    ark_mem->ark_step_free(ark_mem);
  
  /* Attach erkstep_mem structure and function pointers to ark_mem */
  ark_mem->ark_step_init    = erkStep_Init;
  ark_mem->ark_step_resize  = erkStep_Resize;
  ark_mem->ark_step_fullrhs = erkStep_FullRHS;
  if (ark_mem->ark_fixedstep) {
    ark_mem->ark_step       = erkStep_FixedStep;
  } else {
    ark_mem->ark_step       = erkStep_AdaptiveStep;
  }
  ark_mem->ark_step_print   = erkStep_PrintMem;
  ark_mem->ark_step_free    = erkStep_Free;
  ark_mem->ark_step_mem     = (void*) erkstep_mem;

  /* Set default values for ERKStep optional inputs */
  iret = ERKStepSetDefaults((void *) ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ERKStep",
                    "ERKStepCreate", 
                    "Error setting default solver options");
    return(iret);
  }
  
  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepCreate", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = erkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "ERKStepCreate", MSGARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Allocate the general ERK stepper vectors using y0 as a template */
  /* NOTE: F will be allocated later on.  Once we determine which 
     N_Vectors are no longer needed in the main ARKode memory structure, 
     they should be moved here */ 

  /* Copy the input parameters into ARKODE state */
  erkstep_mem->f = f;

  /* Update the ARKode workspace requirements -- UPDATE */
  ark_mem->ark_liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->ark_lrw += 10;
  
  /* Allocate step adaptivity structure, set default values, note storage */
  erkstep_mem->hadapt_mem = arkAdaptInit();
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep", "ERKStepCreate", 
                    "Allocation of step adaptivity structure failed");
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_lrw += ARK_ADAPT_LRW;
  ark_mem->ark_liw += ARK_ADAPT_LIW;
  
  /* Initialize all the counters */
  erkstep_mem->nst_attempts = 0;
  erkstep_mem->nfe          = 0;
  erkstep_mem->netf         = 0;

  /* Initialize main ARKode infrastructure */
  iret = arkodeInit(arkode_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ERKStep", "ERKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(iret);
  }
  
  return(ARK_SUCCESS);
}



int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  int iret;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepReInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access ERKStep memory */
  if (ark_mem->ark_step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "ERKStepReInit",
                    "Missing time step module memory");
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;
  
  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "ERKStepReInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  iret = arkodeReInit(arkode_mem, t0, y0);
  if (iret != ARK_SUCCESS) {
    arkProcessError(ark_mem, iret, "ARKODE::ERKStep", "ERKStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(iret);
  }
  
  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep",
                    "ERKStepReInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  erkstep_mem->f = f;

  /* Destroy/Reinitialize time step adaptivity structure (if present) */
  if (erkstep_mem->hadapt_mem != NULL) {
    free(erkstep_mem->hadapt_mem);
    erkstep_mem->hadapt_mem = arkAdaptInit();
    if (erkstep_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE::ERKStep", "ERKStepReInit", 
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }
  
  /* Initialize all the counters */
  erkstep_mem->nst_attempts = 0;
  erkstep_mem->nfe          = 0;
  erkstep_mem->netf         = 0;

  return(ARK_SUCCESS);
}


/*===============================================================
  ERKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
 erkStep_Init:

 Called from within arkInitialSetup, this routine:
 - sets/checks the ARK Butcher tables to be used
 - allocates any memory that depends on the number of ARK stages,
   method order, or solver options
 - checks for consistency between the system and mass matrix 
   linear solvers (if applicable)
 - initializes and sets up the system and mass matrix linear 
   solvers (if applicable)
---------------------------------------------------------------*/
int erkStep_Init(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  sunindextype Blrw, Bliw;
  int ier, j;

  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_Init", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_Init", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* destroy adaptivity structure if fixed-stepping is requested */
  if (ark_mem->ark_fixedstep) 
    if (erkstep_mem->hadapt_mem != NULL) {
      free(erkstep_mem->hadapt_mem);
      erkstep_mem->hadapt_mem = NULL;
    }
  
  /* Set first step growth factor */
  if (erkstep_mem->hadapt_mem != NULL)
    erkstep_mem->hadapt_mem->etamax = erkstep_mem->hadapt_mem->etamx1;

  /* Create Butcher table (if not already set) */
  ier = erkStep_SetButcherTable(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", "erkStep_Init", 
                    "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher table are OK */
  ier = erkStep_CheckButcherTable(ark_mem);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_Init", "Error in Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* note Butcher table space requirements */
  ButcherTableSpace(erkstep_mem->B, &Bliw, &Blrw);
  ark_mem->ark_liw += Bliw;
  ark_mem->ark_lrw += Blrw;
  
  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate F[0] ... F[stages-1] if needed */
  if (erkstep_mem->F == NULL) 
    erkstep_mem->F = (N_Vector *) calloc(erkstep_mem->stages, sizeof(N_Vector));
  for (j=0; j<erkstep_mem->stages; j++) {
    if (!arkAllocVec(ark_mem, ark_mem->ark_ewt, &(erkstep_mem->F[j])))
      return(ARK_MEM_FAIL);
  }
  ark_mem->ark_liw += erkstep_mem->stages;  /* pointers */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_FullRHS:

 This is just a wrapper to call the user-supplied RHS function, 
 f(t,y).
 
 This will be called in one of three 'modes':
     0 -> called at the beginning of a simulation
     1 -> called at the end of a successful step
     2 -> called elsewhere (e.g. for dense output)

 If it is called in mode 0, we store the vectors f(t,y) in F[0] 
 for possible reuse in the first stage of the subsequent time step.

 If it is called in mode 1 and the method coefficients 
 support it, we may just copy vectors F[stages] to fill f instead 
 of calling f().

 Mode 2 is only called for dense output in-between steps, so we 
 strive to store the intermediate parts so that they do not 
 interfere with the other two modes.
---------------------------------------------------------------*/
int erkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  int i, s, retval;
  booleantype recomputeRHS;
  
  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_FullRHS", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_FullRHS", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vectors f(t,y) in F[0] for possible reuse 
     in the first stage of the subsequent time step */
  case 0:

    /* call f */
    retval = erkstep_mem->f(t, y, erkstep_mem->F[0], ark_mem->ark_user_data);
    erkstep_mem->nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                      "erkStep_FullRHS", MSGARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* copy RHS vector into output */
    N_VScale(ONE, erkstep_mem->F[0], f);

    break;


  /* Mode 1: called at the end of a successful step
     If the method coefficients support it, we just copy the last stage RHS vectors 
     to fill f instead of calling f(t,y).  
     Copy the results to F[0] if the coefficients support it. */
  case 1:

    /* determine if explicit/implicit RHS functions need to be recomputed */
    recomputeRHS = SUNFALSE;
    s = erkstep_mem->B->stages;
    for (i=0; i<s; i++) 
      if (SUNRabs(erkstep_mem->B->b[i] - erkstep_mem->B->A[s-1][i])>TINY)
        recomputeRHS = SUNTRUE;

    /* base RHS calls on recomputeRHS argument */
    if (recomputeRHS) {
      
      /* call f */
      retval = erkstep_mem->f(t, y, erkstep_mem->F[0], ark_mem->ark_user_data);
      erkstep_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                        "erkStep_FullRHS", MSGARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }

    } else {
      N_VScale(ONE, erkstep_mem->F[erkstep_mem->stages-1], erkstep_mem->F[0]);
    }

    /* copy RHS vector into output */
    N_VScale(ONE, erkstep_mem->F[0], f);
    
    break;

  /*  Mode 2: called for dense output in-between steps
      store the intermediate calculations in such a way as to not 
      interfere with the other two modes */
  default:

    /* call f */
    retval = erkstep_mem->f(t, y, ark_mem->ark_tempv2, ark_mem->ark_user_data);
    erkstep_mem->nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE::ERKStep",
                      "erkStep_FullRHS", MSGARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* copy RHS vector into output */
    N_VScale(ONE, ark_mem->ark_tempv2, f);
    
    break;
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_AdaptiveStep:

 This routine serves the primary purpose of the ERKStep module: 
 it performs a single successful embedded ERK step (if possible).  
 Multiple attempts may be taken in this process -- once a step 
 completes with successful (non)linear solves at each stage and 
 passes the error estimate, the routine returns successfully.  
 If it cannot do so, it returns with an appropriate error flag.
---------------------------------------------------------------*/
int erkStep_AdaptiveStep(void* arkode_mem)
{
  realtype dsm, hA, tstage;
  int retval, nef, is, eflag, js;
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_AdaptiveStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_AdaptiveStep", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;
  
  nef = 0;
  eflag = ARK_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {  

    /* increment attempt counter */
    erkstep_mem->nst_attempts++;

#ifdef DEBUG_OUTPUT
 printf("stage 0 RHS:\n");
 N_VPrint_Serial(erkstep_mem->F[0]);
#endif
    
    /* Loop over internal stages to the step; since the method is explicit
       the first stage RHS is just the full RHS from the start of the step */
    for (is=1; is<erkstep_mem->stages; is++) {

      /* Set current stage time(s) */
      tstage = ark_mem->ark_tnew + erkstep_mem->B->c[is]*ark_mem->ark_h;
        
#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n", 
         ark_mem->ark_nst, is, ark_mem->ark_h, tstage);
#endif
      
      /* Initialize y to ycur */
      N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);

      /* Iterate over each prior stage updating rhs */
      for (js=0; js<is; js++) {
        hA = ark_mem->ark_h * erkstep_mem->B->A[is][js];
        N_VLinearSum(hA, erkstep_mem->F[js], ONE, ark_mem->ark_y, ark_mem->ark_y);
      }

      /* compute updated RHS */
      retval = erkstep_mem->f(tstage, ark_mem->ark_y, erkstep_mem->F[is],
                              ark_mem->ark_user_data);
      erkstep_mem->nfe++;
      if (retval < 0)  return(ARK_RHSFUNC_FAIL);
      if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef DEBUG_OUTPUT
 printf("RHS:\n");
 N_VPrint_Serial(erkstep_mem->F[is]);
#endif

    } /* loop over stages */

    /* compute time-evolved solution (in ark_y), error estimate (in dsm) */
    retval = erkStep_ComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);    /* msetup failure */

#ifdef DEBUG_OUTPUT
 printf("error estimate = %"RSYM"\n", dsm);
 printf("updated solution:\n");
 N_VPrint_Serial(ark_mem->ark_y);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->ark_report) 
      fprintf(ark_mem->ark_diagfp, "  etest  %li  %"RSYM"  %"RSYM"\n", 
              ark_mem->ark_nst, ark_mem->ark_h, dsm);

    /* Perform time accuracy error test (if failure, updates h for next try) */
    if (!ark_mem->ark_fixedstep) 
      eflag = erkStep_DoErrorTest(ark_mem, &nef, dsm);
        
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
  retval = erkStep_PrepareNextStep(ark_mem, dsm); 
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_FixedStep:

 This routine serves the primary purpose of the ERKStep module: 
 it performs a single non-embedded ERK step.
---------------------------------------------------------------*/
int erkStep_FixedStep(void* arkode_mem)
{
  realtype dsm, hA, tstage;
  int retval, nef, is, eflag, js;
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_AdaptiveStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_AdaptiveStep", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;
  
  nef = 0;
  eflag = ARK_SUCCESS;

  /* increment attempt counter */
  erkstep_mem->nst_attempts++;

  /* Loop over internal stages to the step; since the method is explicit
     the first stage RHS is just the full RHS from the start of the step */
  for (is=1; is<erkstep_mem->stages; is++) {

    /* Set current stage time(s) */
    tstage = ark_mem->ark_tnew + erkstep_mem->B->c[is]*ark_mem->ark_h;
        
#ifdef DEBUG_OUTPUT
    printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n", 
           ark_mem->ark_nst, is, ark_mem->ark_h, tstage);
#endif
      
    /* Initialize y to ycur */
    N_VScale(ONE, ark_mem->ark_ycur, ark_mem->ark_y);

    /* Iterate over each prior stage updating rhs */
    for (js=0; js<is; js++) {
      hA = ark_mem->ark_h * erkstep_mem->B->A[is][js];
      N_VLinearSum(hA, erkstep_mem->F[js], ONE, ark_mem->ark_y, ark_mem->ark_y);
    }

    /* compute updated RHS */
    retval = erkstep_mem->f(tstage, ark_mem->ark_y, erkstep_mem->F[is],
                            ark_mem->ark_user_data);
    erkstep_mem->nfe++;
    if (retval < 0)  return(ARK_RHSFUNC_FAIL);
    if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
    
  } /* loop over stages */

  /* compute time-evolved solution (in ark_y), error estimate (in dsm) */
  retval = erkStep_ComputeSolutions(ark_mem, &dsm);
  if (retval < 0)  return(retval);

  /* clean up prepare for next step */
  retval = erkStep_PrepareNextStep(ark_mem, dsm); 
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_Resize:

 This routine resizes the memory within the ERKStep module.
---------------------------------------------------------------*/
int erkStep_Resize(void* arkode_mem, ARKVecResizeFn resize,
                   void* resize_data, sunindextype lrw_diff,
                   sunindextype liw_diff, N_Vector tmpl)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  int ier, i;
  
  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_Resize", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_Resize", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;
  
  /* Resize the RHS vectors */
  for (i=0; i<erkstep_mem->stages; i++) {
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, tmpl, &erkstep_mem->F[i]);
    if (ier != ARK_SUCCESS)  return(ier);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_PrintMem:

 This routine outputs the ERKStep structure to a specified file 
 pointer (useful when debugging).
---------------------------------------------------------------*/
void erkStep_PrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  int i;

  /* access ARKodeERKStepMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_PrintMem", MSGARK_NO_MEM);
    return;
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_PrintMem", MSGERKSTEP_NO_MEM);
    return;
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;
  
  /* output integer quantities */
  fprintf(outfile,"ERKStep: q = %i\n", erkstep_mem->q);
  fprintf(outfile,"ERKStep: p = %i\n", erkstep_mem->p);
  fprintf(outfile,"ERKStep: stages = %i\n", erkstep_mem->stages);
  fprintf(outfile,"ERKStep: maxnef = %i\n", erkstep_mem->maxnef);

  /* output long integer quantities */
  fprintf(outfile,"ERKStep: nst_attempts = %li\n", erkstep_mem->nst_attempts);
  fprintf(outfile,"ERKStep: nfe = %li\n", erkstep_mem->nfe);
  fprintf(outfile,"ERKStep: netf = %li\n", erkstep_mem->netf);

  /* output boolean quantities */
  fprintf(outfile,"ERKStep: hadapt_pq = %i\n", erkstep_mem->hadapt_pq);

  /* output realtype quantities */
  fprintf(outfile,"ERKStep: Butcher table:\n");
  WriteButcherTable(erkstep_mem->B, outfile);
  if (erkstep_mem->hadapt_mem != NULL) {
    fprintf(outfile,"ERKStep: timestep adaptivity structure:\n");
    arkPrintAdaptMem(erkstep_mem->hadapt_mem, outfile);
  }

#ifdef DEBUG_OUTPUT
  /* output vector quantities */  
  for (i=0; i<erkstep_mem->stages; i++) {
    fprintf(outfile,"ERKStep: F[%i]:\n", i);
    N_VPrint_Serial(erkstep_mem->F[i]);
  }
#endif
}


/*---------------------------------------------------------------
 erkStep_Free:

 This routine frees all memory internal to the ERKStep module, 
 as well as the module itself.
---------------------------------------------------------------*/
int erkStep_Free(void* arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* trivially return if the ERKStep module is not allocated */
  if (arkode_mem == NULL)  return(ARK_SUCCESS);
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_step_mem==NULL)  return(ARK_SUCCESS);
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* free the time step adaptivity module */
  if (erkstep_mem->hadapt_mem != NULL) {
    free(erkstep_mem->hadapt_mem);
    erkstep_mem->hadapt_mem = NULL;
    ark_mem->ark_lrw -= ARK_ADAPT_LRW;
    ark_mem->ark_liw -= ARK_ADAPT_LIW;
  }

  /* free the Butcher table */
  if (erkstep_mem->B != NULL) {
    ButcherTableSpace(erkstep_mem->B, &Bliw, &Blrw);
    FreeButcherTable(erkstep_mem->B);
    erkstep_mem->B = NULL;
    ark_mem->ark_liw -= Bliw;
    ark_mem->ark_lrw -= Blrw;
  }

  /* free the RHS vectors */
  if (erkstep_mem->F != NULL) {
    for(j=0; j<erkstep_mem->stages; j++) 
      arkFreeVec(ark_mem, &erkstep_mem->F[j]);
    free(erkstep_mem->F);
    erkstep_mem->F = NULL;
  }

  /* free the time stepper module itself */
  free(ark_mem->ark_step_mem);
  ark_mem->ark_step_mem = NULL;
  
  return(ARK_SUCCESS);
}



/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
 erkStep_CheckNVector:

 This routine checks if all required vector operations are 
 present.  If any of them is missing it returns SUNFALSE.
---------------------------------------------------------------*/
booleantype erkStep_CheckNVector(N_Vector tmpl)
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
 erkStep_SetButcherTable

 This routine determines the ERK method to use, based on the 
 desired accuracy.
---------------------------------------------------------------*/
int erkStep_SetButcherTable(ARKodeMem ark_mem)
{
  int etable;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_SetButcherTable", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* if table has already been specified, just return */
  if (erkstep_mem->B != NULL)
    return(ARK_SUCCESS);

  /* initialize table number to illegal values */
  etable = -1;

  /* select method based on order */
  switch (erkstep_mem->q) {
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
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_SetButcherTable", 
                    "No explicit method at requested order, using q=6.");
    etable = DEFAULT_ERK_6;
    break;
  }

  if (etable > -1)
    erkstep_mem->B = ARKodeLoadButcherTable(etable);

  /* set [redundant] stored values for stage numbers and method orders */
  if (erkstep_mem->B != NULL) {
    erkstep_mem->stages = erkstep_mem->B->stages;
    erkstep_mem->q = erkstep_mem->B->q;
    erkstep_mem->p = erkstep_mem->B->p;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_CheckButcherTable

 This routine runs through the explicit Butcher table to ensure 
 that it meets all necessary requirements, including:
     strictly lower-triangular (ERK)
     method order q > 0 (all)
     embedding order q > 0 (all -- if adaptive time-stepping enabled)
     stages > 0 (all)

Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
---------------------------------------------------------------*/
int erkStep_CheckButcherTable(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeERKStepMem erkstep_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeERKStepMem structure */
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_CheckButcherTable", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* check that ERK table is strictly lower triangular */
  okay = SUNTRUE;
  for (i=0; i<erkstep_mem->stages; i++)
    for (j=i; j<erkstep_mem->stages; j++)
      if (SUNRabs(erkstep_mem->B->A[i][j]) > tol)
        okay = SUNFALSE;
  if (!okay) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_CheckButcherTable",
                    "Ae Butcher table is implicit!");
    return(ARK_ILL_INPUT);
  }

  /* check that method order q > 0 */
  if (erkstep_mem->q < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_CheckButcherTable",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((erkstep_mem->p < 1) && (!ark_mem->ark_fixedstep)) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_CheckButcherTable",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that stages > 0 */
  if (erkstep_mem->stages < 1) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ERKStep", 
                    "erkStep_CheckButcherTable",
                    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_ComputeSolutions

 This routine calculates the final RK solution using the existing 
 data.  This solution is placed directly in ark_y.  This routine 
 also computes the error estimate ||y-ytilde||_WRMS, where ytilde 
 is the embedded solution, and the norm weights come from 
 ark_ewt.  This norm value is returned.  The vector form of this 
 estimated error (y-ytilde) is stored in ark_tempv1, in case the 
 calling routine wishes to examine the error locations.

 Note: at this point in the step, the vector ark_tempv1 may be 
 used as a temporary vector.
---------------------------------------------------------------*/
int erkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm)
{
  /* local data */
  realtype hb;
  int ier, j;
  N_Vector y, yerr;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_ComputeSolutions", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* set N_Vector shortcuts */
  y    = ark_mem->ark_y;
  yerr = ark_mem->ark_tempv1;

  /* initialize output */
  *dsm = ZERO;

  /* Compute updated solution and error estimate  */
  /*    Initialize solution to yn, error estimate to zero */
  N_VScale(ONE, ark_mem->ark_ycur, y);
  N_VConst(ZERO, yerr);

  /* Compute time step solution */
  for (j=0; j<erkstep_mem->stages; j++) {
    hb = ark_mem->ark_h * erkstep_mem->B->b[j];
    N_VLinearSum(hb, erkstep_mem->F[j], ONE, y, y);
  }

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->ark_fixedstep) {
    for (j=0; j<erkstep_mem->stages; j++) {
      hb = ark_mem->ark_h * (erkstep_mem->B->b[j] - erkstep_mem->B->d[j]);
      N_VLinearSum(hb, erkstep_mem->F[j], ONE, yerr, yerr);
    }

    /* fill error norm */
    *dsm = N_VWrmsNorm(yerr, ark_mem->ark_ewt);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 erkStep_DoErrorTest

 This routine performs the local error test for the ARK method. 
 The weighted local error norm dsm is passed in, and 
 the test dsm ?<= 1 is made.

 If the test passes, arkDoErrorTest returns ARK_SUCCESS. 

 If the test fails, we revert to the last successful solution 
 time, and:
   - if maxnef error test failures have occurred or if 
     SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
   - otherwise: update time step factor eta based on local error 
     estimate and reduce h.  
---------------------------------------------------------------*/
int erkStep_DoErrorTest(ARKodeMem ark_mem, int *nefPtr, realtype dsm)
{
  realtype ehist2, hhist2;
  int retval, k;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_DoErrorTest", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep", "arkDoErrorTest",
                    MSGARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters */
  (*nefPtr)++;
  erkstep_mem->netf++;

  /* At |h| = hmin or maxnef failures, return ARK_ERR_FAILURE */
  if ((SUNRabs(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) ||
      (*nefPtr == erkstep_mem->maxnef))
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
  hadapt_mem->hhist[0] = ark_mem->ark_h;

  /* Compute accuracy-based time step estimate (updated ark_eta) */
  k = (erkstep_mem->hadapt_pq) ? erkstep_mem->q : erkstep_mem->p;
  retval = arkAdapt((void*) ark_mem, erkstep_mem->hadapt_mem, ark_mem->ark_y,
                    ark_mem->ark_tn, ark_mem->ark_h, k, ark_mem->ark_nst);
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
    ark_mem->ark_eta = SUNMIN(ark_mem->ark_eta, hadapt_mem->etamxf);
  ark_mem->ark_h *= ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
 erkStep_PrepareNextStep

 This routine handles ARK-specific updates following a successful 
 step: copying the ARK result to the current solution vector, 
 updating the error/step history arrays, and setting the 
 prospective step size, hprime, for the next step.  Along with 
 hprime, it sets the ratio eta=hprime/h.  It also updates other 
 state variables related to a change of step size.
---------------------------------------------------------------*/
int erkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval, k;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->ark_step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ERKStep",
                    "erkStep_PrepareNextStep", MSGERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->ark_step_mem;

  /* Update step size and error history arrays */
  if (erkstep_mem->hadapt_mem != NULL) {
    erkstep_mem->hadapt_mem->ehist[2] = erkstep_mem->hadapt_mem->ehist[1];
    erkstep_mem->hadapt_mem->ehist[1] = erkstep_mem->hadapt_mem->ehist[0];
    erkstep_mem->hadapt_mem->ehist[0] = dsm*erkstep_mem->hadapt_mem->bias;
    erkstep_mem->hadapt_mem->hhist[2] = erkstep_mem->hadapt_mem->hhist[1];
    erkstep_mem->hadapt_mem->hhist[1] = erkstep_mem->hadapt_mem->hhist[0];
    erkstep_mem->hadapt_mem->hhist[0] = ark_mem->ark_h;
  }
  
  /* If fixed time-stepping requested, defer 
     step size changes until next step */
  if (ark_mem->ark_fixedstep){
    ark_mem->ark_hprime = ark_mem->ark_h;
    ark_mem->ark_eta = ONE;
    return(ARK_SUCCESS);
  }

  /* If etamax = 1, defer step size changes until next step, 
     and reset etamax */
  if (erkstep_mem->hadapt_mem != NULL)
    if (erkstep_mem->hadapt_mem->etamax == ONE) {
      ark_mem->ark_hprime = ark_mem->ark_h;
      ark_mem->ark_eta = ONE;
      erkstep_mem->hadapt_mem->etamax = erkstep_mem->hadapt_mem->growth;
      return(ARK_SUCCESS);
    }

  /* Adjust ark_eta in arkAdapt */
  if (erkstep_mem->hadapt_mem != NULL) {
    k = (erkstep_mem->hadapt_pq) ? erkstep_mem->q : erkstep_mem->p;
    retval = arkAdapt((void*) ark_mem, erkstep_mem->hadapt_mem,
                      ark_mem->ark_y, ark_mem->ark_tn + ark_mem->ark_h,
                      ark_mem->ark_h, k, ark_mem->ark_nst+1);
    if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);
  }
  
  /* Set hprime value for next step size */
  ark_mem->ark_hprime = ark_mem->ark_h * ark_mem->ark_eta;
  
  /* Reset growth factor for subsequent time step */
  if (erkstep_mem->hadapt_mem != NULL)
    erkstep_mem->hadapt_mem->etamax = erkstep_mem->hadapt_mem->growth;

  return(ARK_SUCCESS);
}


/*===============================================================
   EOF
===============================================================*/
