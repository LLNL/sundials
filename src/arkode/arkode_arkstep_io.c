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
 * output functions for the ARKODE ARKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here, 
 * with slightly different names.  The code transition will be 
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_arkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  ARKStep Optional input functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepSetDefaults:

  Resets all ARKStep optional inputs to their default values.  
  Does not change problem-defining function pointers or 
  user_data pointer.  Also leaves alone any data 
  structures/options related to the ARKode infrastructure itself
  (e.g. root-finding).
  ---------------------------------------------------------------*/
int ARKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetDefaults", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  
  /* Set default values for integrator optional inputs */
  arkstep_mem->q                = Q_DEFAULT;      /* method order */
  arkstep_mem->p                = 0;              /* embedding order */
  arkstep_mem->hadapt_pq        = SUNFALSE;       /* use embedding order */
  arkstep_mem->predictor        = 0;              /* trivial predictor */
  arkstep_mem->linear           = SUNFALSE;       /* nonlinear problem */
  arkstep_mem->linear_timedep   = SUNTRUE;        /* dfi/dy depends on t */
  arkstep_mem->explicit         = SUNTRUE;        /* fe(t,y) will be used */
  arkstep_mem->implicit         = SUNTRUE;        /* fi(t,y) will be used */
  if (arkstep_mem->hadapt_mem != NULL) {
    arkstep_mem->hadapt_mem->etamx1      = ETAMX1;     /* max change on first step */
    arkstep_mem->hadapt_mem->etamxf      = ETAMXF;     /* max change on error-failed step */
    arkstep_mem->hadapt_mem->small_nef   = SMALL_NEF;  /* num error fails before ETAMXF enforced */
    arkstep_mem->hadapt_mem->etacf       = ETACF;      /* max change on convergence failure */
    arkstep_mem->hadapt_mem->HAdapt      = NULL;       /* step adaptivity fn */
    arkstep_mem->hadapt_mem->HAdapt_data = NULL;       /* step adaptivity data */
    arkstep_mem->hadapt_mem->imethod     = 0;          /* PID controller */
    arkstep_mem->hadapt_mem->cfl         = CFLFAC;     /* explicit stability factor */
    arkstep_mem->hadapt_mem->safety      = SAFETY;     /* step adaptivity safety factor  */
    arkstep_mem->hadapt_mem->bias        = BIAS;       /* step adaptivity error bias */
    arkstep_mem->hadapt_mem->growth      = GROWTH;     /* step adaptivity growth factor */
    arkstep_mem->hadapt_mem->lbound      = HFIXED_LB;  /* step adaptivity no-change lower bound */
    arkstep_mem->hadapt_mem->ubound      = HFIXED_UB;  /* step adaptivity no-change upper bound */
    arkstep_mem->hadapt_mem->k1          = AD0_K1;     /* step adaptivity parameter */
    arkstep_mem->hadapt_mem->k2          = AD0_K2;     /* step adaptivity parameter */
    arkstep_mem->hadapt_mem->k3          = AD0_K3;     /* step adaptivity parameter */
  }
  arkstep_mem->maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  arkstep_mem->maxnef           = MAXNEF;         /* max error test fails */
  arkstep_mem->maxncf           = MAXNCF;         /* max convergence fails */
  arkstep_mem->nlscoef          = NLSCOEF;        /* nonlinear tolerance coefficient */
  arkstep_mem->crdown           = CRDOWN;         /* nonlinear convergence estimate coeff. */
  arkstep_mem->rdiv             = RDIV;           /* nonlinear divergence tolerance */
  arkstep_mem->dgmax            = DGMAX;          /* max step change before recomputing J or P */
  arkstep_mem->msbp             = MSBP;           /* max steps between updates to J or P */
  arkstep_mem->use_fp           = SUNFALSE;       /* use Newton solver */
  arkstep_mem->fp_m             = FP_ACCEL_M;     /* num Anderson acceleration vectors */
  arkstep_mem->fp_mem           = NULL;
  arkstep_mem->stages           = 0;              /* no stages */
  arkstep_mem->istage           = 0;              /* current stage */
  arkstep_mem->Be               = NULL;           /* no Butcher tables */
  arkstep_mem->Bi               = NULL;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetOptimalParams:

  Sets all adaptivity and solver parameters to our 'best guess' 
  values, for a given ARKStep integration method (ERK, DIRK, ARK), 
  a given method order, and a given nonlinear solver type.  Should
  only be called after the method order, solver, and integration
  method have been set, and only if time step adaptivity is 
  enabled.
  ---------------------------------------------------------------*/
int ARKStepSetOptimalParams(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetOptimalParams", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetOptimalParams",
                    "Missing ARKStep time stepper module");
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetOptimalParams",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* Choose values based on method, order */

  /*    explicit */
  if (arkstep_mem->explicit && !arkstep_mem->implicit) {
    hadapt_mem->imethod = 1;
    hadapt_mem->safety  = RCONST(0.99);
    hadapt_mem->bias    = RCONST(1.2);
    hadapt_mem->growth  = RCONST(25.0);
    hadapt_mem->k1      = RCONST(0.8);
    hadapt_mem->k2      = RCONST(0.31);
    hadapt_mem->etamxf  = RCONST(0.3);

    /*    implicit */
  } else if (arkstep_mem->implicit && !arkstep_mem->explicit) {
    switch (arkstep_mem->q) {
    case 2:   /* just use standard defaults since better ones unknown */
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = SAFETY;
      hadapt_mem->bias      = BIAS;
      hadapt_mem->growth    = GROWTH;
      hadapt_mem->etamxf    = ETAMXF;
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.001);
      arkstep_mem->maxcor   = 5;
      arkstep_mem->crdown   = CRDOWN;
      arkstep_mem->rdiv     = RDIV;
      arkstep_mem->dgmax    = DGMAX;
      arkstep_mem->msbp     = MSBP;
      break;
    case 3:
      hadapt_mem->imethod   = 2;
      hadapt_mem->safety    = RCONST(0.957);
      hadapt_mem->bias      = RCONST(1.9);
      hadapt_mem->growth    = RCONST(17.6);
      hadapt_mem->etamxf    = RCONST(0.45);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.22);
      arkstep_mem->crdown   = RCONST(0.17);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.19);
      arkstep_mem->msbp     = 60;
      break;
    case 4:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.988);
      hadapt_mem->bias      = RCONST(1.2);
      hadapt_mem->growth    = RCONST(31.5);
      hadapt_mem->k1        = RCONST(0.535);
      hadapt_mem->k2        = RCONST(0.209);
      hadapt_mem->k3        = RCONST(0.148);
      hadapt_mem->etamxf    = RCONST(0.33);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.24);
      arkstep_mem->crdown   = RCONST(0.26);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.16);
      arkstep_mem->msbp     = 31;
      break;
    case 5:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.937);
      hadapt_mem->bias      = RCONST(3.3);
      hadapt_mem->growth    = RCONST(22.0);
      hadapt_mem->k1        = RCONST(0.56);
      hadapt_mem->k2        = RCONST(0.338);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.44);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.25);
      arkstep_mem->crdown   = RCONST(0.4);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.32);
      arkstep_mem->msbp     = 31;
      break;
    }

    /* newton vs fixed-point */
    if (arkstep_mem->use_fp)  arkstep_mem->maxcor = 10;

    /*    imex */
  } else {
    switch (arkstep_mem->q) {
    case 3:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.965);
      hadapt_mem->bias      = RCONST(1.42);
      hadapt_mem->growth    = RCONST(28.7);
      hadapt_mem->k1        = RCONST(0.54);
      hadapt_mem->k2        = RCONST(0.36);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.46);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.22);
      arkstep_mem->crdown   = RCONST(0.17);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.19);
      arkstep_mem->msbp     = 60;
      break;
    case 4:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.97);
      hadapt_mem->bias      = RCONST(1.35);
      hadapt_mem->growth    = RCONST(25.0);
      hadapt_mem->k1        = RCONST(0.543);
      hadapt_mem->k2        = RCONST(0.297);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.47);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.24);
      arkstep_mem->crdown   = RCONST(0.26);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.16);
      arkstep_mem->msbp     = 31;
      break;
    case 5:
      hadapt_mem->imethod   = 1;
      hadapt_mem->safety    = RCONST(0.993);
      hadapt_mem->bias      = RCONST(1.15);
      hadapt_mem->growth    = RCONST(28.5);
      hadapt_mem->k1        = RCONST(0.8);
      hadapt_mem->k2        = RCONST(0.35);
      hadapt_mem->etamxf    = RCONST(0.3);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      arkstep_mem->nlscoef  = RCONST(0.25);
      arkstep_mem->crdown   = RCONST(0.4);
      arkstep_mem->rdiv     = RCONST(2.3);
      arkstep_mem->dgmax    = RCONST(0.32);
      arkstep_mem->msbp     = 31;
      break;
    }

    /* newton vs fixed-point */
    if (arkstep_mem->use_fp)  arkstep_mem->maxcor = 10;

  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along 
  with ARKStepSetERKTable, ARKStepSetIRKTable, ARKStepSetARKTable, 
  ARKStepSetERKTableNum, ARKStepSetIRKTableNum or 
  ARKStepSetARKTableNum.  This routine is used to specify a 
  desired method order using default Butcher tables, whereas 
  any user-supplied table will have their own order associated 
  with them.
  ---------------------------------------------------------------*/
int ARKStepSetOrder(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetOrder",  MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    arkstep_mem->q = Q_DEFAULT;
  } else {
    arkstep_mem->q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  arkstep_mem->stages = 0;
  arkstep_mem->istage = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetLinear:

  Specifies that the implicit portion of the problem is linear, 
  and to tighten the linear solver tolerances while taking only 
  one Newton iteration.  Not useful when used in combination with
  the fixed-point solver.  Automatically tightens DeltaGammaMax 
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the 
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when 
  using an iterative linear solver this flag denotes time 
  dependence of the preconditioner. 
  ---------------------------------------------------------------*/
int ARKStepSetLinear(void *arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetLinear", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetLinear", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set parameters */
  arkstep_mem->linear = SUNTRUE;
  arkstep_mem->linear_timedep = (timedepend == 1);
  arkstep_mem->dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to ARKStepSetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetNonlinear", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetNonlinear", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set parameters */
  arkstep_mem->linear = SUNFALSE;
  arkstep_mem->linear_timedep = SUNTRUE;
  arkstep_mem->dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetFixedPoint:

  Specifies that the implicit portion of the problem should be 
  solved using the accelerated fixed-point solver instead of the 
  Newton iteration.  Allowed values for the dimension of the 
  acceleration space, fp_m, must be non-negative.  Illegal 
  values imply to use the default.
  ---------------------------------------------------------------*/
int ARKStepSetFixedPoint(void *arkode_mem, long int fp_m)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetFixedPoint", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetFixedPoint", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set parameters */
  arkstep_mem->use_fp = SUNTRUE;
  if (fp_m < 0) {
    arkstep_mem->fp_m = FP_ACCEL_M;
  } else {
    arkstep_mem->fp_m = fp_m;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNewton:

  Specifies that the implicit portion of the problem should be 
  solved using the modified Newton solver.  Used to undo a 
  previous call to ARKStepSetFixedPoint.
  ---------------------------------------------------------------*/
int ARKStepSetNewton(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetNewton", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetNewton", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set parameters */
  arkstep_mem->use_fp = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetExplicit:

  Specifies that the implicit portion of the problem is disabled, 
  and to use an explicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetExplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetExplicit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetExplicit", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* ensure that fe is defined */
  if (arkstep_mem->fe == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetExplicit", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  arkstep_mem->explicit = SUNTRUE;
  arkstep_mem->implicit = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetImplicit:

  Specifies that the explicit portion of the problem is disabled, 
  and to use an implicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetImplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetImplicit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetImplicit", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* ensure that fi is defined */
  if (arkstep_mem->fi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetImplicit", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  arkstep_mem->implicit = SUNTRUE;
  arkstep_mem->explicit = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetImEx:

  Specifies that the specifies that problem has both implicit and
  explicit parts, and to use an ARK method (this is the default).
  ---------------------------------------------------------------*/
int ARKStepSetImEx(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetImEx", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetImEx", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* ensure that fe and fi are defined */
  if (arkstep_mem->fe == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetImEx", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }
  if (arkstep_mem->fi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetImEx", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  arkstep_mem->explicit = SUNTRUE;
  arkstep_mem->implicit = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetERKTable:

  Specifies to use a customized Butcher table for the explicit 
  portion of the system (automatically calls ARKStepSetExplicit).
 
  If d==NULL, then the method is automatically flagged as a 
  fixed-step method; a user MUST also call either 
  ARKodeSetFixedStep or ARKodeSetInitStep to set the desired time 
  step size.
  ---------------------------------------------------------------*/
int ARKStepSetERKTable(void *arkode_mem, int s, int q, int p,
                       realtype *c, realtype *A, realtype *b, 
                       realtype *d)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetERKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetERKTable", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check for legal inputs */
  if ((c == NULL) || (A == NULL) || (b == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetERKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* set the relevant parameters */
  arkstep_mem->stages = s;
  arkstep_mem->q = q;
  arkstep_mem->p = p;
  arkstep_mem->Be = AllocButcherTable(s, (d != NULL));
  arkstep_mem->Be->q = q;
  arkstep_mem->Be->p = p;
  for (i=0; i<s; i++) {
    arkstep_mem->Be->c[i] = c[i];
    arkstep_mem->Be->b[i] = b[i];
    for (j=0; j<s; j++) {
      arkstep_mem->Be->A[i][j] = A[i*s + j];
    }
  }

  /* set embedding (if applicable), otherwise set as fixed-step method */
  if (d == NULL) {
    ark_mem->fixedstep = SUNTRUE;
    if (arkstep_mem->hadapt_mem != NULL) {
      free(arkstep_mem->hadapt_mem);
      arkstep_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) 
      arkstep_mem->Be->d[i] = d[i];
  }
  
  /* set method as purely explicit */
  if (ARKStepSetExplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetERKTable", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetIRKTable:

  Specifies to use a customized Butcher table for the implicit 
  portion of the system (automatically calls ARKStepSetImplicit).

  If d==NULL, then the method is automatically flagged as a 
  fixed-step method; a user MUST also call either 
  ARKodeSetFixedStep or ARKodeSetInitStep to set the desired time 
  step size.
  ---------------------------------------------------------------*/
int ARKStepSetIRKTable(void *arkode_mem, int s, int q, int p,
                       realtype *c, realtype *A, realtype *b, 
                       realtype *d)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetIRKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetIRKTable", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check for legal inputs */
  if ((c == NULL) || (A == NULL) || (b == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetIRKTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* set the relevant parameters */
  arkstep_mem->stages = s;
  arkstep_mem->q = q;
  arkstep_mem->p = p;
  arkstep_mem->Bi = AllocButcherTable(s, (d != NULL));
  arkstep_mem->Bi->q = q;
  arkstep_mem->Bi->p = p;
  for (i=0; i<s; i++) {
    arkstep_mem->Bi->c[i] = c[i];
    arkstep_mem->Bi->b[i] = b[i];
    for (j=0; j<s; j++) {
      arkstep_mem->Bi->A[i][j] = A[i*s + j];
    }
  }

  /* set embedding (if applicable), otherwise set as fixed-step method */
  if (d == NULL) {
    ark_mem->fixedstep = SUNTRUE;
    if (arkstep_mem->hadapt_mem != NULL) {
      free(arkstep_mem->hadapt_mem);
      arkstep_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) 
      arkstep_mem->Bi->d[i] = d[i];
  }
  
  /* set method as purely implicit */
  if (ARKStepSetImplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetIRKTable", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetARKTables:

  Specifies to use customized Butcher tables for the ImEx system
  (automatically calls ARKStepSetImEx).

  If either de==NULL or di==NULL, then the method is 
  automatically flagged as a fixed-step method; a user MUST also 
  call either ARKodeSetFixedStep or ARKodeSetInitStep to set the 
  desired time step size.
  ---------------------------------------------------------------*/
int ARKStepSetARKTables(void *arkode_mem, int s, int q, int p,
                        realtype *ci, realtype *ce, 
                        realtype *Ai, realtype *Ae, 
                        realtype *bi, realtype *be, 
                        realtype *di, realtype *de)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetARKTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetARKTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check for legal inputs */
  if ((ci == NULL) || (ce == NULL) || 
      (Ai == NULL) || (Ae == NULL) || 
      (bi == NULL) || (be == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetARKTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* set the relevant parameters */
  arkstep_mem->stages = s;
  arkstep_mem->q = q;
  arkstep_mem->p = p;
  arkstep_mem->Be = AllocButcherTable(s, (de != NULL));
  arkstep_mem->Bi = AllocButcherTable(s, (di != NULL));
  arkstep_mem->Be->p = p;
  arkstep_mem->Bi->q = q;
  arkstep_mem->Be->q = q;
  arkstep_mem->Bi->p = p;
  for (i=0; i<s; i++) {
    arkstep_mem->Be->c[i] = ce[i];
    arkstep_mem->Bi->c[i] = ci[i];
    arkstep_mem->Be->b[i] = be[i];
    arkstep_mem->Bi->b[i] = bi[i];
    for (j=0; j<s; j++) {
      arkstep_mem->Be->A[i][j] = Ae[i*s + j];
      arkstep_mem->Bi->A[i][j] = Ai[i*s + j];
    }
  }

  /* set embeddings (if applicable), otherwise set as fixed-step method */
  if ((de == NULL) || (di == NULL)) {
    ark_mem->fixedstep = SUNTRUE;
    if (arkstep_mem->hadapt_mem != NULL) {
      free(arkstep_mem->hadapt_mem);
      arkstep_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) {
      arkstep_mem->Be->d[i] = de[i];
      arkstep_mem->Bi->d[i] = di[i];
    }
  }
  
  /* set method as ImEx */
  if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetARKTables", MSG_ARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetERKTableNum:

  Specifies to use a pre-existing Butcher table for the explicit 
  portion of the problem, based on the integer flag held in 
  ARKodeLoadButcherTable() within the file arkode_butcher.c 
  (automatically calls ARKStepSetExplicit).
  ---------------------------------------------------------------*/
int ARKStepSetERKTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetERKTableNum", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetERKTableNum", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check that argument specifies an explicit table */
  if (itable<MIN_ERK_NUM || itable>MAX_ERK_NUM) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetERKTableNum", 
                    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* fill in table based on argument */
  arkstep_mem->Be = ARKodeLoadButcherTable(itable);
  if (arkstep_mem->Be == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetERKTableNum", 
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  arkstep_mem->stages = arkstep_mem->Be->stages;
  arkstep_mem->q = arkstep_mem->Be->q;
  arkstep_mem->p = arkstep_mem->Be->p;

  /* set method as purely explicit */
  if (ARKStepSetExplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetERKTableNum", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetIRKTableNum:

  Specifies to use a pre-existing Butcher table for the implicit 
  portion of the problem, based on the integer flag held in 
  ARKodeLoadButcherTable() within the file arkode_butcher.c
  (automatically calls ARKStepSetImplicit).
  ---------------------------------------------------------------*/
int ARKStepSetIRKTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetIRKTableNum", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetIRKTableNum", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check that argument specifies an implicit table */
  if (itable<MIN_DIRK_NUM || itable>MAX_DIRK_NUM) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetIRKTableNum", 
                    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* fill in table based on argument */
  arkstep_mem->Bi = ARKodeLoadButcherTable(itable);
  if (arkstep_mem->Bi == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetIRKTableNum", 
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  arkstep_mem->stages = arkstep_mem->Bi->stages;
  arkstep_mem->q = arkstep_mem->Bi->q;
  arkstep_mem->p = arkstep_mem->Bi->p;

  /* set method as purely implicit */
  if (ARKStepSetImplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetIRKTableNum", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetARKTableNum:

  Specifies to use pre-existing Butcher tables for the ImEx system,
  based on the integer flags held in ARKodeLoadButcherTable() 
  within the file arkode_butcher.c (automatically calls ARKStepSetImEx).
  ---------------------------------------------------------------*/
int ARKStepSetARKTableNum(void *arkode_mem, int itable, int etable)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetARKTableNum", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetARKTableNum", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* ensure that tables match */
  if ( !((etable == ARK324L2SA_ERK_4_2_3) && (itable == ARK324L2SA_DIRK_4_2_3)) &&
       !((etable == ARK436L2SA_ERK_6_3_4) && (itable == ARK436L2SA_DIRK_6_3_4)) &&
       !((etable == ARK548L2SA_ERK_8_4_5) && (itable == ARK548L2SA_DIRK_8_4_5)) ) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetARKTableNum", 
                    "Incompatible Butcher tables for ARK method");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  arkstep_mem->stages = 0;
  arkstep_mem->q = 0;
  arkstep_mem->p = 0;
  FreeButcherTable(arkstep_mem->Be);  arkstep_mem->Be = NULL;
  FreeButcherTable(arkstep_mem->Bi);  arkstep_mem->Bi = NULL;

  /* fill in tables based on arguments */
  arkstep_mem->Bi = ARKodeLoadButcherTable(itable);
  arkstep_mem->Be = ARKodeLoadButcherTable(etable);
  if (arkstep_mem->Bi == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetARKTableNum", 
                    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }
  if (arkstep_mem->Be == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetARKTableNum", 
                    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }
  arkstep_mem->stages = arkstep_mem->Bi->stages;
  arkstep_mem->q = arkstep_mem->Bi->q;
  arkstep_mem->p = arkstep_mem->Bi->p;

  /* set method as ImEx */
  if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetARKTableNum", MSG_ARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetCFLFraction", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetCFLFraction", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetCFLFraction",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetCFLFraction", "Illegal CFL fraction");
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
  ARKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted 
  time step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to the 
  default value.
  ---------------------------------------------------------------*/
int ARKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetSafetyFactor", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetSafetyFactor", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetSafetyFactoy",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetSafetyFactor", "Illegal safety factor");
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
  ARKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetErrorBias", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetErrorBias", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetErrorBias", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (bias < 1.0) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses 
  a separate maximum growth factor.  Allowable values must be 
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxGrowth", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxGrowth", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (mx_growth == ZERO) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step 
  size will remain unchanged.  Allowable values must enclose the 
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetFixedStepBounds", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetFixedStepBounds", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetFixedStepBounds", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

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
  ARKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and 
  optionally, its associated parameters) to use.  All parameters 
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int ARKStepSetAdaptivityMethod(void *arkode_mem, int imethod, 
                               int idefault, int pq, 
                               realtype *adapt_params)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetAdaptivityMethod", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetAdaptivityMethod", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetAdaptivityMethod", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
                    "ARKStepSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p or q */
  arkstep_mem->hadapt_pq = (pq != 0);

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
  ARKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int ARKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
                           void *h_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetAdaptivityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetAdaptivityFn", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetAdaptivityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;

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
  ARKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamx1.  Legal values are greater than 1.0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ARKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxFirstGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxFirstGrowth", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxFirstGrowth",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamxf. Legal values are in the interval (0,1].  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ARKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxEFailGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxEFailGrowth", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxEFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ARKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetSmallNumEFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetSmallNumEFails", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetSmallNumEFails", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int ARKStepSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxCFailGrowth", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxCFailGrowth", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxCFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values 
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinCRDown(void *arkode_mem, realtype crdown)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetNonlinCRDown", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetNonlinCRDown", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) {
    arkstep_mem->crdown = CRDOWN;
  } else {
    arkstep_mem->crdown = crdown;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values 
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinRDiv(void *arkode_mem, realtype rdiv)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetNonlinRDiv", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetNonlinRDiv", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) {
    arkstep_mem->rdiv = RDIV;
  } else {
    arkstep_mem->rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply 
  a reset to the default. 
  ---------------------------------------------------------------*/
int ARKStepSetDeltaGammaMax(void *arkode_mem, realtype dgmax)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetDeltaGammaMax", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetDeltaGammaMax", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) {
    arkstep_mem->dgmax = DGMAX;
  } else {
    arkstep_mem->dgmax = dgmax;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxStepsBetweenLSet:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the number of time steps to wait 
  before calling lsetup; negative values imply recomputation of 
  lsetup at each Newton iteration; a zero value implies a reset 
  to the default. 
  ---------------------------------------------------------------*/
int ARKStepSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxStepsBetweenLSet", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxStepsBetweenLSet", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    arkstep_mem->msbp = MSBP;
  } else {
    arkstep_mem->msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.  
  Non-default choices are {1,2,3,4}, all others will use default 
  (trivial) predictor.
  ---------------------------------------------------------------*/
int ARKStepSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetPredictorMethod", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetPredictorMethod", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set parameters */
  arkstep_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability 
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int ARKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
                          void *estab_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeARKStepMem and ARKodeHAdaptMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetStabilityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetStabilityFn", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  if (arkstep_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetStabilityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = arkstep_mem->hadapt_mem;
  
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
  ARKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxErrTestFails", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxErrTestFails", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    arkstep_mem->maxnef = MAXNEF;
  } else {
    arkstep_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures 
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxConvFails", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxConvFails", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    arkstep_mem->maxncf = MAXNCF;
  } else {
    arkstep_mem->maxncf = maxncf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the 
  default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetMaxNonlinIters", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetMaxNonlinIters", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    arkstep_mem->maxcor = MAXCOR;
  } else {
    arkstep_mem->maxcor = maxcor;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepSetNonlinConvCoef", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepSetNonlinConvCoef", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    arkstep_mem->nlscoef = NLSCOEF;
  } else {
    arkstep_mem->nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE optional output functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int ARKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumExpSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumExpSteps", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (arkstep_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = arkstep_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int ARKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumAccSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumAccSteps", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (arkstep_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = arkstep_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int ARKStepGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumStepAttempts", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumStepAttempts", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* get value from arkstep_mem */
  *nsteps = arkstep_mem->nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int ARKStepGetNumRhsEvals(void *arkode_mem, long int *fe_evals,
                          long int *fi_evals)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumRhsEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumRhsEvals", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* get values from arkstep_mem */
  *fe_evals = arkstep_mem->nfe;
  *fi_evals = arkstep_mem->nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int ARKStepGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumLinSolvSetups", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumLinSolvSetups", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* get value from arkstep_mem */
  *nlinsetups = arkstep_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int ARKStepGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumErrTestFails", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* get value from arkstep_mem */
  *netfails = arkstep_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetCurrentButcherTables:

  Sets pointers to the explicit and implicit Butcher tables 
  currently in use.
  ---------------------------------------------------------------*/
int ARKStepGetCurrentButcherTables(void *arkode_mem,
                                   ARKodeButcherTable *Bi,
                                   ARKodeButcherTable *Be)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetCurrentButcherTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetCurrentButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* get tables from arkstep_mem */
  *Bi = arkstep_mem->Bi;
  *Be = arkstep_mem->Be;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetEstLocalErrors: (updated to the correct vector, but 
  need to verify that it is unchanged between filling the 
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ARKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetEstLocalErrors", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetEstLocalErrors", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ARKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                               long int *accsteps, long int *step_attempts,
                               long int *fe_evals, long int *fi_evals,
                               long int *nlinsetups, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetTimestepperStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetTimestepperStats", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if step adaptivity structure not allocated, 
     just set expsteps and accsteps to 0 */
  if (arkstep_mem->hadapt_mem == NULL) {
    *expsteps = 0;
    *accsteps = 0;
  } else {
    *expsteps = arkstep_mem->hadapt_mem->nst_exp;
    *accsteps = arkstep_mem->hadapt_mem->nst_acc;
  }

  /* set remaining outputs from arkstep_mem */
  *step_attempts = arkstep_mem->nst_attempts;
  *fe_evals      = arkstep_mem->nfe;
  *fi_evals      = arkstep_mem->nfi;
  *nlinsetups    = arkstep_mem->nsetups;
  *netfails      = arkstep_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations 
  ---------------------------------------------------------------*/
int ARKStepGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumNonlinSolvIters", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumNonlinSolvIters", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set output from arkstep_mem */
  *nniters = arkstep_mem->nni;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int ARKStepGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNumNonlinSolvConvFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNumNonlinSolvConvFails", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set output from arkstep_mem */
  *nncfails = arkstep_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int ARKStepGetNonlinSolvStats(void *arkode_mem, long int *nniters, 
                              long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepGetNonlinSolvStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepGetNonlinSolvStats", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set outputs from arkstep_mem */
  *nniters = arkstep_mem->nni;
  *nncfails = arkstep_mem->ncfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ARKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKStepWriteParameters", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKStepWriteParameters", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ARKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",arkstep_mem->q);
  if (arkstep_mem->linear) {
    fprintf(fp, "  Linear implicit problem");
    if (arkstep_mem->linear_timedep) {
      fprintf(fp, " (time-dependent Jacobian)\n");
    } else {
      fprintf(fp, " (time-independent Jacobian)\n");
    }
  }
  if (arkstep_mem->explicit && arkstep_mem->implicit) {
    fprintf(fp, "  ImEx integrator\n");
  } else if (arkstep_mem->implicit) {
    fprintf(fp, "  Implicit integrator\n");
  } else {
    fprintf(fp, "  Explicit integrator\n");
  }
  if (arkstep_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            arkstep_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            arkstep_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            arkstep_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            arkstep_mem->hadapt_mem->etacf); 
    if (arkstep_mem->explicit) 
      fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
              arkstep_mem->hadapt_mem->cfl);
    if (arkstep_mem->hadapt_mem->HAdapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", arkstep_mem->hadapt_mem->imethod);
      fprintf(fp, "     Safety factor = %"RSYM"\n", arkstep_mem->hadapt_mem->safety);
      fprintf(fp, "     Bias factor = %"RSYM"\n", arkstep_mem->hadapt_mem->bias);
      fprintf(fp, "     Growth factor = %"RSYM"\n", arkstep_mem->hadapt_mem->growth);
      fprintf(fp, "     Step growth lower bound = %"RSYM"\n", arkstep_mem->hadapt_mem->lbound);
      fprintf(fp, "     Step growth upper bound = %"RSYM"\n", arkstep_mem->hadapt_mem->ubound);
      fprintf(fp, "     k1 = %"RSYM"\n", arkstep_mem->hadapt_mem->k1);
      fprintf(fp, "     k2 = %"RSYM"\n", arkstep_mem->hadapt_mem->k2);
      fprintf(fp, "     k3 = %"RSYM"\n", arkstep_mem->hadapt_mem->k3);
      if (arkstep_mem->hadapt_mem->expstab == arkExpStab) {
        fprintf(fp, "  Default explicit stability function\n");
      } else {
        fprintf(fp, "  User provided explicit stability function\n");
      }
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",arkstep_mem->maxnef);

  if (arkstep_mem->implicit) {
    fprintf(fp, "  Maximum number of convergence test failures = %i\n",arkstep_mem->maxncf);
    fprintf(fp, "  Implicit predictor method = %i\n",arkstep_mem->predictor);
    fprintf(fp, "  Implicit solver tolerance coefficient = %"RSYM"\n",arkstep_mem->nlscoef);
    fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",arkstep_mem->maxcor);
    fprintf(fp, "  Nonlinear convergence rate constant = %"RSYM"\n",arkstep_mem->crdown);
    fprintf(fp, "  Nonlinear divergence tolerance = %"RSYM"\n",arkstep_mem->rdiv);
    fprintf(fp, "  Gamma factor LSetup tolerance = %"RSYM"\n",arkstep_mem->dgmax);
    fprintf(fp, "  Number of steps between LSetup calls = %i\n",arkstep_mem->msbp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int ARKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int i, j;
  ARKodeMem ark_mem;
  ARKodeARKStepMem arkstep_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKodeWriteButcher", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKodeWriteButcher", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkstep_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* print Butcher tables to file */
  fprintf(fp, "\nARKStep Butcher tables (stages = %i):\n", arkstep_mem->stages);
  if (arkstep_mem->explicit) {
    fprintf(fp, "  Explicit Butcher table:\n");
    for (i=0; i<arkstep_mem->stages; i++) {
      fprintf(fp, "     %"RSYM"",arkstep_mem->Be->c[i]);
      for (j=0; j<arkstep_mem->stages; j++) 
        fprintf(fp, " %"RSYM"",arkstep_mem->Be->A[i][j]);
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<arkstep_mem->stages; j++) 
      fprintf(fp, " %"RSYM"",arkstep_mem->Be->b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    if (arkstep_mem->Be->d != NULL) {
      for (j=0; j<arkstep_mem->stages; j++) 
        fprintf(fp, " %"RSYM"",arkstep_mem->Be->d[j]);
      fprintf(fp,"\n");
    }
  }
  if (arkstep_mem->implicit) {
    fprintf(fp, "  Implicit Butcher table:\n");
    for (i=0; i<arkstep_mem->stages; i++) {
      fprintf(fp, "     %"RSYM"",arkstep_mem->Bi->c[i]);
      for (j=0; j<arkstep_mem->stages; j++) 
        fprintf(fp, " %"RSYM"",arkstep_mem->Bi->A[i][j]);
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<arkstep_mem->stages; j++) 
      fprintf(fp, " %"RSYM"",arkstep_mem->Bi->b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    if (arkstep_mem->Bi->d != NULL) {
      for (j=0; j<arkstep_mem->stages; j++) 
        fprintf(fp, " %"RSYM"",arkstep_mem->Bi->d[j]);
      fprintf(fp,"\n");
    }
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
