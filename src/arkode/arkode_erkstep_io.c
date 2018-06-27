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
 * This is the implementation file for the optional input and 
 * output functions for the ARKODE ERKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here, 
 * with slightly different names.  The code transition will be 
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_erkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  ERKStep Optional input functions
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetDefaults:

  Resets all ERKStep optional inputs to their default values.  
  Does not change problem-defining function pointers or 
  user_data pointer.  Also leaves alone any data 
  structures/options related to the ARKode infrastructure itself
  (e.g. root-finding).
  ---------------------------------------------------------------*/
int ERKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetDefaults", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  
  /* Set default values for integrator optional inputs */
  erkstep_mem->q         = Q_DEFAULT;      /* method order */
  erkstep_mem->p         = 0;              /* embedding order */
  erkstep_mem->hadapt_pq = SUNFALSE;       /* use embedding order */
  if (erkstep_mem->hadapt_mem != NULL) {
    erkstep_mem->hadapt_mem->etamx1      = ETAMX1;       /* max change on first step */
    erkstep_mem->hadapt_mem->etamxf      = RCONST(0.3);  /* max change on error-failed step */
    erkstep_mem->hadapt_mem->small_nef   = SMALL_NEF ;   /* num error fails before ETAMXF enforced */
    erkstep_mem->hadapt_mem->etacf       = ETACF;        /* max change on convergence failure */
    erkstep_mem->hadapt_mem->HAdapt      = NULL;         /* step adaptivity fn */
    erkstep_mem->hadapt_mem->HAdapt_data = NULL;         /* step adaptivity data */
    erkstep_mem->hadapt_mem->imethod     = 1;            /* PI controller */
    erkstep_mem->hadapt_mem->cfl         = CFLFAC;       /* explicit stability factor */
    erkstep_mem->hadapt_mem->safety      = RCONST(0.99); /* step adaptivity safety factor  */
    erkstep_mem->hadapt_mem->bias        = RCONST(1.2);  /* step adaptivity error bias */
    erkstep_mem->hadapt_mem->growth      = RCONST(25.0); /* step adaptivity growth factor */
    erkstep_mem->hadapt_mem->lbound      = HFIXED_LB;    /* step adaptivity no-change lower bound */
    erkstep_mem->hadapt_mem->ubound      = HFIXED_UB;    /* step adaptivity no-change upper bound */
    erkstep_mem->hadapt_mem->k1          = RCONST(0.8);  /* step adaptivity parameter */
    erkstep_mem->hadapt_mem->k2          = RCONST(0.31); /* step adaptivity parameter */
    erkstep_mem->hadapt_mem->k3          = AD0_K3;       /* step adaptivity parameter */
  }
  erkstep_mem->maxnef = MAXNEF;         /* max error test fails */
  erkstep_mem->stages = 0;              /* no stages */
  erkstep_mem->B      = NULL;           /* no Butcher table */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along 
  with ERKStepSetERKTable or ERKStepSetERKTableNum.  This 
  routine is used to specify a desired method order using 
  default Butcher tables, whereas any user-supplied table will 
  have their own order associated with them.
  ---------------------------------------------------------------*/
int ERKStepSetOrder(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetOrder",  MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    erkstep_mem->q = Q_DEFAULT;
  } else {
    erkstep_mem->q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  erkstep_mem->stages = 0;
  erkstep_mem->p = 0;
  FreeButcherTable(erkstep_mem->B);  erkstep_mem->B = NULL;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetERKTable:

  Specifies to use a customized Butcher table for the explicit 
  portion of the system.
 
  If d==NULL, then the method is automatically flagged as a 
  fixed-step method; a user MUST also call either 
  ARKodeSetFixedStep or ARKodeSetInitStep to set the desired time 
  step size.
  ---------------------------------------------------------------*/
int ERKStepSetERKTable(void *arkode_mem, int s, int q, int p,
                       realtype *c, realtype *A, realtype *b, 
                       realtype *d)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetERKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetERKTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* check for legal inputs */
  if ((c == NULL) || (A == NULL) || (b == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetERKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  erkstep_mem->stages = 0;
  erkstep_mem->q = 0;
  erkstep_mem->p = 0;
  FreeButcherTable(erkstep_mem->B);  erkstep_mem->B = NULL;

  /* set the relevant parameters */
  erkstep_mem->stages = s;
  erkstep_mem->q = q;
  erkstep_mem->p = p;
  erkstep_mem->B = AllocButcherTable(s, (d != NULL));
  erkstep_mem->B->q = q;
  erkstep_mem->B->p = p;
  for (i=0; i<s; i++) {
    erkstep_mem->B->c[i] = c[i];
    erkstep_mem->B->b[i] = b[i];
    for (j=0; j<s; j++) {
      erkstep_mem->B->A[i][j] = A[i*s + j];
    }
  }

  /* set embedding (if applicable), otherwise set as fixed-step method */
  if (d == NULL) {
    ark_mem->fixedstep = SUNTRUE;
    if (erkstep_mem->hadapt_mem != NULL) {
      free(erkstep_mem->hadapt_mem);
      erkstep_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) 
      erkstep_mem->B->d[i] = d[i];
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetERKTableNum:

  Specifies to use a pre-existing Butcher table for the problem, 
  based on the integer flag held in ARKodeLoadButcherTable() 
  within the file arkode_butcher.c.
  ---------------------------------------------------------------*/
int ERKStepSetERKTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetERKTableNum", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetERKTableNum", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* check that argument specifies an explicit table */
  if (itable<MIN_ERK_NUM || itable>MAX_ERK_NUM) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetERKTableNum", 
                    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  erkstep_mem->stages = 0;
  erkstep_mem->q = 0;
  erkstep_mem->p = 0;
  FreeButcherTable(erkstep_mem->B);  erkstep_mem->B = NULL;

  /* fill in table based on argument */
  erkstep_mem->B = ARKodeLoadButcherTable(itable);
  if (erkstep_mem->B == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetERKTableNum", 
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  erkstep_mem->stages = erkstep_mem->B->stages;
  erkstep_mem->q = erkstep_mem->B->q;
  erkstep_mem->p = erkstep_mem->B->p;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ERKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetCFLFraction", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetCFLFraction", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetCFLFraction",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ERKStepSetCFLFraction", "Illegal CFL fraction");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (cfl_frac <= ZERO) {
    hadapt_mem->cfl = CFLFAC;
  } else {
    hadapt_mem->cfl = cfl_frac;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted 
  time step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to the 
  default value.
  ---------------------------------------------------------------*/
int ERKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetSafetyFactor", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetSafetyFactor", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetSafetyFactoy",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ERKStepSetSafetyFactor", "Illegal safety factor");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (safety <= ZERO) {
    hadapt_mem->safety = SAFETY;
  } else {
    hadapt_mem->safety = safety;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int ERKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetErrorBias", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetErrorBias", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetErrorBias", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (bias < 1.0) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses 
  a separate maximum growth factor.  Allowable values must be 
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int ERKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetMaxGrowth", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxGrowth", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (mx_growth == ZERO) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step 
  size will remain unchanged.  Allowable values must enclose the 
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int ERKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetFixedStepBounds", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetFixedStepBounds", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetFixedStepBounds", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* set allowable interval, otherwise set defaults */
  if ((lb <= 1.0) && (ub >= 1.0)) {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  } else {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and 
  optionally, its associated parameters) to use.  All parameters 
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityMethod(void *arkode_mem, int imethod, 
                               int idefault, int pq, 
                               realtype *adapt_params)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetAdaptivityMethod", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetAdaptivityMethod", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetAdaptivityMethod", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ERKStepSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p or q */
  erkstep_mem->hadapt_pq = (pq != 0);

  /* set method parameters */
  if (idefault == 1) {
    switch (hadapt_mem->imethod) {
    case (0):
      hadapt_mem->k1 = AD0_K1; 
      hadapt_mem->k2 = AD0_K2;
      hadapt_mem->k3 = AD0_K3; break;
    case (1):
      hadapt_mem->k1 = AD1_K1;
      hadapt_mem->k2 = AD1_K2; break;
    case (2):
      hadapt_mem->k1 = AD2_K1; break;
    case (3):
      hadapt_mem->k1 = AD3_K1;
      hadapt_mem->k2 = AD3_K2; break;
    case (4):
      hadapt_mem->k1 = AD4_K1;
      hadapt_mem->k2 = AD4_K2; break;
    case (5):
      hadapt_mem->k1 = AD5_K1;
      hadapt_mem->k2 = AD5_K2; 
      hadapt_mem->k3 = AD5_K3; break;
    }
  } else {
    hadapt_mem->k1 = adapt_params[0];
    hadapt_mem->k2 = adapt_params[1];
    hadapt_mem->k3 = adapt_params[2];
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
                           void *h_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetAdaptivityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetAdaptivityFn", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetAdaptivityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;

  /* NULL hfun sets default, otherwise set inputs */
  if (hfun == NULL) {
    hadapt_mem->HAdapt      = NULL;
    hadapt_mem->HAdapt_data = NULL;
    hadapt_mem->imethod     = 0;
  } else {
    hadapt_mem->HAdapt      = hfun;
    hadapt_mem->HAdapt_data = h_data;
    hadapt_mem->imethod     = -1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamx1.  Legal values are greater than 1.0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ERKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetMaxFirstGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxFirstGrowth", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxFirstGrowth",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamxf. Legal values are in the interval (0,1].  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ERKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetMaxEFailGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxEFailGrowth", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxEFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ERKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetSmallNumEFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetSmallNumEFails", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetSmallNumEFails", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability 
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int ERKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
                          void *estab_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeERKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetStabilityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetStabilityFn", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;
  if (erkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetStabilityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = erkstep_mem->hadapt_mem;
  
  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    hadapt_mem->expstab    = arkExpStab;
    hadapt_mem->estab_data = ark_mem;
  } else {
    hadapt_mem->expstab    = EStab;
    hadapt_mem->estab_data = estab_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ERKStepSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepSetMaxErrTestFails", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepSetMaxErrTestFails", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    erkstep_mem->maxnef = MAXNEF;
  } else {
    erkstep_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE optional output functions
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int ERKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetNumExpSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetNumExpSteps", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (erkstep_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = erkstep_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int ERKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetNumAccSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetNumAccSteps", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (erkstep_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = erkstep_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int ERKStepGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetNumStepAttempts", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetNumStepAttempts", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* get value from erkstep_mem */
  *nsteps = erkstep_mem->nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int ERKStepGetNumRhsEvals(void *arkode_mem, long int *fevals)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetNumRhsEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetNumRhsEvals", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* get values from erkstep_mem */
  *fevals = erkstep_mem->nfe;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int ERKStepGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetNumErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetNumErrTestFails", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* get value from erkstep_mem */
  *netfails = erkstep_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetCurrentButcherTable:

  Sets pointers to the Butcher table currently in use.
  ---------------------------------------------------------------*/
int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                  ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetCurrentButcherTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetCurrentButcherTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* get tables from erkstep_mem */
  *B = erkstep_mem->B;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetEstLocalErrors: (updated to the correct vector, but 
  need to verify that it is unchanged between filling the 
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ERKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetEstLocalErrors", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetEstLocalErrors", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ERKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                               long int *accsteps, long int *attempts,
                               long int *fevals, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepGetTimestepperStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepGetTimestepperStats", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, 
     just set expsteps and accsteps to 0 */
  if (erkstep_mem->hadapt_mem == NULL) {
    *expsteps = 0;
    *accsteps = 0;
  } else {
    *expsteps = erkstep_mem->hadapt_mem->nst_exp;
    *accsteps = erkstep_mem->hadapt_mem->nst_acc;
  }

  /* set remaining outputs from erkstep_mem */
  *attempts = erkstep_mem->nst_attempts;
  *fevals   = erkstep_mem->nfe;
  *netfails = erkstep_mem->netf;

  return(ARK_SUCCESS);
}


/*===============================================================
  ERKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ERKStepWriteParameters", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ERKStepWriteParameters", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ERKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",erkstep_mem->q);
  if (erkstep_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            erkstep_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            erkstep_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            erkstep_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            erkstep_mem->hadapt_mem->etacf); 
    fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
            erkstep_mem->hadapt_mem->cfl);
    if (erkstep_mem->hadapt_mem->HAdapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", erkstep_mem->hadapt_mem->imethod);
      fprintf(fp, "     Safety factor = %"RSYM"\n", erkstep_mem->hadapt_mem->safety);
      fprintf(fp, "     Bias factor = %"RSYM"\n", erkstep_mem->hadapt_mem->bias);
      fprintf(fp, "     Growth factor = %"RSYM"\n", erkstep_mem->hadapt_mem->growth);
      fprintf(fp, "     Step growth lower bound = %"RSYM"\n", erkstep_mem->hadapt_mem->lbound);
      fprintf(fp, "     Step growth upper bound = %"RSYM"\n", erkstep_mem->hadapt_mem->ubound);
      fprintf(fp, "     k1 = %"RSYM"\n", erkstep_mem->hadapt_mem->k1);
      fprintf(fp, "     k2 = %"RSYM"\n", erkstep_mem->hadapt_mem->k2);
      fprintf(fp, "     k3 = %"RSYM"\n", erkstep_mem->hadapt_mem->k3);
      if (erkstep_mem->hadapt_mem->expstab == arkExpStab) {
        fprintf(fp, "  Default explicit stability function\n");
      } else {
        fprintf(fp, "  User provided explicit stability function\n");
      }
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",erkstep_mem->maxnef);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeERKStepMem erkstep_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKodeWriteButcher", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKodeWriteButcher", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  erkstep_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (erkstep_mem->B == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeWriteButcher", 
                    "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }
  
  /* print Butcher table to file */
  fprintf(fp, "\nERKStep Butcher table (stages = %i):\n", erkstep_mem->stages);
  for (i=0; i<erkstep_mem->stages; i++) {
    fprintf(fp, "     %"RSYM"",erkstep_mem->B->c[i]);
    for (j=0; j<erkstep_mem->stages; j++) 
      fprintf(fp, " %"RSYM"",erkstep_mem->B->A[i][j]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "            ");
  for (j=0; j<erkstep_mem->stages; j++) 
    fprintf(fp, " %"RSYM"",erkstep_mem->B->b[j]);
  fprintf(fp,"\n");
  fprintf(fp, "            ");
  if (erkstep_mem->B->d != NULL) {
    for (j=0; j<erkstep_mem->stages; j++) 
      fprintf(fp, " %"RSYM"",erkstep_mem->B->d[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
