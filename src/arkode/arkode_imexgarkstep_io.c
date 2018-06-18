/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on arkode_arkstep_io.c written by Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the implementation file for the optional input and 
 * output functions for the ARKODE ARKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here, 
 * with slightly different names.  The code transition will be 
 * minimal, but the documentation changes will be significant.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_imexgarkstep_impl.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_types.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#define RSYME "11.4e"
#endif


/* =============================================================================
 * IMEXGARKStep Optional input functions
 * ===========================================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepSetDefaults:

  Resets all IMEXGARKStep optional inputs to their default values.  
  Does not change problem-defining function pointers or 
  user_data pointer.  Also leaves alone any data 
  structures/options related to the ARKode infrastructure itself
  (e.g. root-finding).
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDefaults(void* mem)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetDefaults", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  
  /* Set default values for integrator optional inputs */
  step_mem->q                = Q_DEFAULT;      /* method order */
  step_mem->p                = 0;              /* embedding order */
  step_mem->hadapt_pq        = SUNFALSE;       /* use embedding order */
  step_mem->predictor        = 0;              /* trivial predictor */
  step_mem->linear           = SUNFALSE;       /* nonlinear problem */
  step_mem->linear_timedep   = SUNTRUE;        /* dfi/dy depends on t */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->etamx1      = ETAMX1;     /* max change on first step */
    step_mem->hadapt_mem->etamxf      = ETAMXF;     /* max change on error-failed step */
    step_mem->hadapt_mem->small_nef   = SMALL_NEF;  /* num error fails before ETAMXF enforced */
    step_mem->hadapt_mem->etacf       = ETACF;      /* max change on convergence failure */
    step_mem->hadapt_mem->HAdapt      = NULL;       /* step adaptivity fn */
    step_mem->hadapt_mem->HAdapt_data = NULL;       /* step adaptivity data */
    step_mem->hadapt_mem->imethod     = 0;          /* PID controller */
    step_mem->hadapt_mem->cfl         = CFLFAC;     /* explicit stability factor */
    step_mem->hadapt_mem->safety      = SAFETY;     /* step adaptivity safety factor  */
    step_mem->hadapt_mem->bias        = BIAS;       /* step adaptivity error bias */
    step_mem->hadapt_mem->growth      = GROWTH;     /* step adaptivity growth factor */
    step_mem->hadapt_mem->lbound      = HFIXED_LB;  /* step adaptivity no-change lower bound */
    step_mem->hadapt_mem->ubound      = HFIXED_UB;  /* step adaptivity no-change upper bound */
    step_mem->hadapt_mem->k1          = AD0_K1;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k2          = AD0_K2;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k3          = AD0_K3;     /* step adaptivity parameter */
  }
  step_mem->maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  step_mem->maxnef           = MAXNEF;         /* max error test fails */
  step_mem->maxncf           = MAXNCF;         /* max convergence fails */
  step_mem->nlscoef          = NLSCOEF;        /* nonlinear tolerance coefficient */
  step_mem->crdown           = CRDOWN;         /* nonlinear convergence estimate coeff. */
  step_mem->rdiv             = RDIV;           /* nonlinear divergence tolerance */
  step_mem->dgmax            = DGMAX;          /* max step change before recomputing J or P */
  step_mem->msbp             = MSBP;           /* max steps between updates to J or P */
  step_mem->use_fp           = SUNFALSE;       /* use Newton solver */
  step_mem->fp_m             = FP_ACCEL_M;     /* num Anderson acceleration vectors */
  step_mem->fp_mem           = NULL;
  step_mem->stages           = 0;              /* no stages */
  step_mem->istage           = 0;              /* current stage */
  step_mem->Bee              = NULL;           /* no Butcher tables */
  step_mem->Bii              = NULL;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetLinear:

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
int IMEXGARKStepSetLinear(void *mem, int timedepend)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetLinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetLinear", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set parameters */
  step_mem->linear = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to IMEXGARKStepSetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinear(void *mem)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetNonlinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetNonlinear", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set parameters */
  step_mem->linear = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetFixedPoint:

  Specifies that the implicit portion of the problem should be 
  solved using the accelerated fixed-point solver instead of the 
  Newton iteration.  Allowed values for the dimension of the 
  acceleration space, fp_m, must be non-negative.  Illegal 
  values imply to use the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetFixedPoint(void *mem, long int fp_m)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetFixedPoint", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetFixedPoint", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set parameters */
  step_mem->use_fp = SUNTRUE;
  if (fp_m < 0) {
    step_mem->fp_m = FP_ACCEL_M;
  } else {
    step_mem->fp_m = fp_m;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNewton:

  Specifies that the implicit portion of the problem should be 
  solved using the modified Newton solver.  Used to undo a 
  previous call to IMEXGARKStepSetFixedPoint.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNewton(void *mem)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetNewton", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetNewton", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set parameters */
  step_mem->use_fp = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetARKTables:

  Specifies to use customized Butcher tables for the ImEx system
  (automatically calls IMEXGARKStepSetImEx).

  If either b2e == NULL or b2i == NULL, then the method is 
  automatically flagged as a fixed-step method; a user MUST also 
  call either ARKodeSetFixedStep or ARKodeSetInitStep to set the 
  desired time step size.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IMEXGARKStepSetButcherTables(void *mem, int s, 
                                                 int q, int p, 
                                                 realtype *ce, realtype *ci,
                                                 realtype *Aee, realtype *Aei,
                                                 realtype *Aie, realtype *Aii,
                                                 realtype *be, realtype *bi,
                                                 realtype *de, realtype *di)
{
  int i, j;
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetButcherTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetButcherTables", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* check for legal inputs */
  if ((ci  == NULL) || (ce  == NULL) ||
      (Aee == NULL) || (Aei == NULL) ||
      (Aie == NULL) || (Aii == NULL) ||
      (bi  == NULL) || (be  == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetButcherTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  FreeButcherTable(step_mem->Bee); step_mem->Bee = NULL;
  FreeButcherTable(step_mem->Bei); step_mem->Bei = NULL; /* don't need a full Butcher table (only A) */
  FreeButcherTable(step_mem->Bie); step_mem->Bie = NULL; /* don't need a full Butcher table (only A) */
  FreeButcherTable(step_mem->Bii); step_mem->Bii = NULL;

  /* set the relevant parameters */
  step_mem->stages = s;
  step_mem->q = q;
  step_mem->p = p;
  step_mem->Bee = AllocButcherTable(s, (de != NULL));
  step_mem->Bei = AllocButcherTable(s, (de != NULL)); /* don't need a full Butcher table (only A) */
  step_mem->Bie = AllocButcherTable(s, (de != NULL)); /* don't need a full Butcher table (only A) */
  step_mem->Bii = AllocButcherTable(s, (di != NULL));
  step_mem->Bee->p = p;
  step_mem->Bii->q = q;
  step_mem->Bee->q = q;
  step_mem->Bii->p = p;
  for (i=0; i<s; i++) {
    step_mem->Bee->c[i] = ce[i];
    step_mem->Bii->c[i] = ci[i];
    step_mem->Bee->b[i] = be[i];
    step_mem->Bii->b[i] = bi[i];
    for (j=0; j<s; j++) {
      step_mem->Bee->A[i][j] = Aee[i*s + j];
      step_mem->Bei->A[i][j] = Aei[i*s + j]; /* not necessarily square */
      step_mem->Bie->A[i][j] = Aie[i*s + j]; /* not necessarily square */
      step_mem->Bii->A[i][j] = Aii[i*s + j];
    }
  }

  /* set embeddings (if applicable), otherwise set as fixed-step method */
  if ((de == NULL) || (di == NULL)) {
    arkode_mem->fixedstep = SUNTRUE;
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }
  } else {
    for (i=0; i<s; i++) {
      step_mem->Bee->d[i] = de[i];
      step_mem->Bii->d[i] = di[i];
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetCFLFraction(void *mem, realtype cfl_frac)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetCFLFraction", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetCFLFraction", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetCFLFraction",
                    MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(arkode_mem, ARK_ILL_INPUT, "ARKODE", 
                    "IMEXGARKStepSetCFLFraction", "Illegal CFL fraction");
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
  IMEXGARKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted 
  time step size.  Allowable values must be within the open 
  interval (0,1).  A non-positive input implies a reset to the 
  default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetSafetyFactor(void *mem, realtype safety)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetSafetyFactor", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetSafetyFactor", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetSafetyFactoy",MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(arkode_mem, ARK_ILL_INPUT, "ARKODE", 
                    "IMEXGARKStepSetSafetyFactor", "Illegal safety factor");
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
  IMEXGARKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetErrorBias(void *mem, realtype bias)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetErrorBias", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetErrorBias", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetErrorBias", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (bias < 1.0) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses 
  a separate maximum growth factor.  Allowable values must be 
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxGrowth(void *mem, realtype mx_growth)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxGrowth", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxGrowth", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxGrowth", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (mx_growth == ZERO) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step 
  size will remain unchanged.  Allowable values must enclose the 
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetFixedStepBounds(void *mem, realtype lb, realtype ub)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetFixedStepBounds", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetFixedStepBounds", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetFixedStepBounds", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

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
  IMEXGARKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and 
  optionally, its associated parameters) to use.  All parameters 
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetAdaptivityMethod(void *mem, int imethod, 
                               int idefault, int pq, 
                               realtype *adapt_params)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetAdaptivityMethod", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetAdaptivityMethod", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetAdaptivityMethod", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(arkode_mem, ARK_ILL_INPUT, "ARKODE", 
                    "IMEXGARKStepSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p or q */
  step_mem->hadapt_pq = (pq != 0);

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
  IMEXGARKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetAdaptivityFn(void *mem, ARKAdaptFn hfun,
                           void *h_data)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetAdaptivityFn", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetAdaptivityFn", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

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
  IMEXGARKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamx1.  Legal values are greater than 1.0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxFirstGrowth(void *mem, realtype etamx1)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxFirstGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxFirstGrowth", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxFirstGrowth",MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant 
  etamxf. Legal values are in the interval (0,1].  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxEFailGrowth(void *mem, realtype etamxf)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxEFailGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxEFailGrowth", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxEFailGrowth", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetSmallNumEFails(void *mem, int small_nef)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetSmallNumEFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetSmallNumEFails", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetSmallNumEFails", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values 
  imply a reset to the default value. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxCFailGrowth(void *mem, realtype etacf)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxCFailGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxCFailGrowth", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxCFailGrowth", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;
  
  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values 
  imply a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinCRDown(void *mem, realtype crdown)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetNonlinCRDown", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetNonlinCRDown", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) {
    step_mem->crdown = CRDOWN;
  } else {
    step_mem->crdown = crdown;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values 
  imply a reset to the default.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinRDiv(void *mem, realtype rdiv)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetNonlinRDiv", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetNonlinRDiv", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) {
    step_mem->rdiv = RDIV;
  } else {
    step_mem->rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply 
  a reset to the default. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetDeltaGammaMax(void *mem, realtype dgmax)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetDeltaGammaMax", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetDeltaGammaMax", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) {
    step_mem->dgmax = DGMAX;
  } else {
    step_mem->dgmax = dgmax;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxStepsBetweenLSet:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the number of time steps to wait 
  before calling lsetup; negative values imply recomputation of 
  lsetup at each Newton iteration; a zero value implies a reset 
  to the default. 
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxStepsBetweenLSet(void *mem, int msbp)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxStepsBetweenLSet", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxStepsBetweenLSet", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    step_mem->msbp = MSBP;
  } else {
    step_mem->msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.  
  Non-default choices are {1,2,3,4}, all others will use default 
  (trivial) predictor.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetPredictorMethod(void *mem, int pred_method)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetPredictorMethod", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetPredictorMethod", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set parameters */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability 
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int IMEXGARKStepSetStabilityFn(void *mem, ARKExpStabFn EStab,
                          void *estab_data)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;

  /* access ARKodeMem, ARKodeIMEXGARKStepMem and ARKodeHAdaptMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetStabilityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetStabilityFn", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetStabilityFn", MSG_IMEXGARK_NO_ADAPT_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;
  
  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    hadapt_mem->expstab    = arkExpStab;
    hadapt_mem->estab_data = arkode_mem;
  } else {
    hadapt_mem->expstab    = EStab;
    hadapt_mem->estab_data = estab_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxErrTestFails(void *mem, int maxnef)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxErrTestFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxErrTestFails", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    step_mem->maxnef = MAXNEF;
  } else {
    step_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures 
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxConvFails(void *mem, int maxncf)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxConvFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxConvFails", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    step_mem->maxncf = MAXNCF;
  } else {
    step_mem->maxncf = maxncf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the 
  default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetMaxNonlinIters(void *mem, int maxcor)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetMaxNonlinIters", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetMaxNonlinIters", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    step_mem->maxcor = MAXCOR;
  } else {
    step_mem->maxcor = maxcor;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int IMEXGARKStepSetNonlinConvCoef(void *mem, realtype nlscoef)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepSetNonlinConvCoef", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepSetNonlinConvCoef", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    step_mem->nlscoef = NLSCOEF;
  } else {
    step_mem->nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE optional output functions
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumExpSteps(void *mem, long int *nsteps)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumExpSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumExpSteps", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumAccSteps(void *mem, long int *nsteps)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumAccSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumAccSteps", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumStepAttempts(void *mem, long int *nsteps)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumStepAttempts", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumStepAttempts", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* get value from step_mem */
  *nsteps = step_mem->nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumRhsEvals(void *mem, long int *fe_evals,
                          long int *fi_evals)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumRhsEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumRhsEvals", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* get values from step_mem */
  *fe_evals = step_mem->nfe;
  *fi_evals = step_mem->nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumLinSolvSetups(void *mem, long int *nlinsetups)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumLinSolvSetups", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumLinSolvSetups", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumErrTestFails(void *mem, long int *netfails)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumErrTestFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumErrTestFails", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* get value from step_mem */
  *netfails = step_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetCurrentButcherTables:

  Sets pointers to the explicit and implicit Butcher tables 
  currently in use.
  ---------------------------------------------------------------*/
int IMEXGARKStepGetCurrentButcherTables(void *mem,
                                        ARKodeButcherTable *Bee,
                                        ARKodeButcherTable *Bei,
                                        ARKodeButcherTable *Bie,
                                        ARKodeButcherTable *Bii)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetCurrentButcherTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetCurrentButcherTables", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* get tables from step_mem */
  *Bii = step_mem->Bii;
  *Bee = step_mem->Bee;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetEstLocalErrors: (updated to the correct vector, but 
  need to verify that it is unchanged between filling the 
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int IMEXGARKStepGetEstLocalErrors(void *mem, N_Vector ele)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetEstLocalErrors", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetEstLocalErrors", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* copy vector to output */
  N_VScale(ONE, arkode_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int IMEXGARKStepGetTimestepperStats(void *mem, long int *expsteps,
                               long int *accsteps, long int *step_attempts,
                               long int *fe_evals, long int *fi_evals,
                               long int *nlinsetups, long int *netfails)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetTimestepperStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetTimestepperStats", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* if step adaptivity structure not allocated, 
     just set expsteps and accsteps to 0 */
  if (step_mem->hadapt_mem == NULL) {
    *expsteps = 0;
    *accsteps = 0;
  } else {
    *expsteps = step_mem->hadapt_mem->nst_exp;
    *accsteps = step_mem->hadapt_mem->nst_acc;
  }

  /* set remaining outputs from step_mem */
  *step_attempts = step_mem->nst_attempts;
  *fe_evals      = step_mem->nfe;
  *fi_evals      = step_mem->nfi;
  *nlinsetups    = step_mem->nsetups;
  *netfails      = step_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations 
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumNonlinSolvIters(void *mem, long int *nniters)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumNonlinSolvIters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumNonlinSolvIters", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set output from step_mem */
  *nniters = step_mem->nni;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNumNonlinSolvConvFails(void *mem, long int *nncfails)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNumNonlinSolvConvFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNumNonlinSolvConvFails", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set output from step_mem */
  *nncfails = step_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int IMEXGARKStepGetNonlinSolvStats(void *mem, long int *nniters, 
                              long int *nncfails)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepGetNonlinSolvStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepGetNonlinSolvStats", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* set outputs from step_mem */
  *nniters = step_mem->nni;
  *nncfails = step_mem->ncfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  IMEXGARKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  IMEXGARKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int IMEXGARKStepWriteParameters(void *mem, FILE *fp)
{
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "IMEXGARKStepWriteParameters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "IMEXGARKStepWriteParameters", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* print integrator parameters to file */
  fprintf(fp, "IMEXGARKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",step_mem->q);
  if (step_mem->linear) {
    fprintf(fp, "  Linear implicit problem");
    if (step_mem->linear_timedep) {
      fprintf(fp, " (time-dependent Jacobian)\n");
    } else {
      fprintf(fp, " (time-independent Jacobian)\n");
    }
  }
  fprintf(fp, "  ImEx integrator\n");
  if (step_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            step_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            step_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            step_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            step_mem->hadapt_mem->etacf); 
    fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
            step_mem->hadapt_mem->cfl);
    if (step_mem->hadapt_mem->HAdapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", step_mem->hadapt_mem->imethod);
      fprintf(fp, "     Safety factor = %"RSYM"\n", step_mem->hadapt_mem->safety);
      fprintf(fp, "     Bias factor = %"RSYM"\n", step_mem->hadapt_mem->bias);
      fprintf(fp, "     Growth factor = %"RSYM"\n", step_mem->hadapt_mem->growth);
      fprintf(fp, "     Step growth lower bound = %"RSYM"\n", step_mem->hadapt_mem->lbound);
      fprintf(fp, "     Step growth upper bound = %"RSYM"\n", step_mem->hadapt_mem->ubound);
      fprintf(fp, "     k1 = %"RSYM"\n", step_mem->hadapt_mem->k1);
      fprintf(fp, "     k2 = %"RSYM"\n", step_mem->hadapt_mem->k2);
      fprintf(fp, "     k3 = %"RSYM"\n", step_mem->hadapt_mem->k3);
      if (step_mem->hadapt_mem->expstab == arkExpStab) {
        fprintf(fp, "  Default explicit stability function\n");
      } else {
        fprintf(fp, "  User provided explicit stability function\n");
      }
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",step_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",step_mem->maxncf);
  fprintf(fp, "  Implicit predictor method = %i\n",step_mem->predictor);
  fprintf(fp, "  Implicit solver tolerance coefficient = %"RSYM"\n",step_mem->nlscoef);
  fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",step_mem->maxcor);
  fprintf(fp, "  Nonlinear convergence rate constant = %"RSYM"\n",step_mem->crdown);
  fprintf(fp, "  Nonlinear divergence tolerance = %"RSYM"\n",step_mem->rdiv);
  fprintf(fp, "  Gamma factor LSetup tolerance = %"RSYM"\n",step_mem->dgmax);
  fprintf(fp, "  Number of steps between LSetup calls = %i\n",step_mem->msbp);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  IMEXGARKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int IMEXGARKStepWriteButcher(void *mem, FILE *fp)
{
  int i, j;
  ARKodeMem arkode_mem;
  ARKodeIMEXGARKStepMem step_mem;

  /* access ARKodeMem and ARKodeIMEXGARKStepMem structures */
  if (mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
                    "ARKodeWriteButcher", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  arkode_mem = (ARKodeMem) mem;
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "ARKodeWriteButcher", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  if ( (step_mem->Bee == NULL) || (step_mem->Bii == NULL) ) {
    arkProcessError(arkode_mem, ARK_ILL_INPUT, "ARKODE",
                    "IMEXGARKStepWriteButcher",
                    "Butcher tables must both be non-NULL");
    return(ARK_ILL_INPUT);
  }

  /* print Butcher tables to file */
  fprintf(fp, "\nIMEXGARKStep Butcher table (stages = %i):\n\n", step_mem->stages);

  for (i=0; i<step_mem->stages; i++) {
    fprintf(fp, "     %"RSYME"",step_mem->Bee->c[i]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bee->A[i][j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bei->A[i][j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  fprintf(fp, "      ");
  for (j=0; j < (step_mem->stages*2)+1; j++) 
    fprintf(fp, "-----------");
  fprintf(fp,"\n");

  for (i=0; i<step_mem->stages; i++) {
    fprintf(fp, "     %"RSYME"",step_mem->Bii->c[i]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bie->A[i][j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bii->A[i][j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  fprintf(fp, "      ");
  for (j=0; j < (step_mem->stages*2)+1; j++) 
    fprintf(fp, "-----------");
  fprintf(fp,"\n");

  fprintf(fp, "                ");
  fprintf(fp, " | ");
  for (j=0; j<step_mem->stages; j++) 
    fprintf(fp, " %"RSYME"",step_mem->Bee->b[j]);
  fprintf(fp, " | ");
  for (j=0; j<step_mem->stages; j++) 
    fprintf(fp, " %"RSYME"",step_mem->Bii->b[j]);
  fprintf(fp, " | ");
  fprintf(fp,"\n");

  fprintf(fp, "                ");
  fprintf(fp, " | ");
  if (step_mem->Bee->d != NULL) {
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bee->d[j]);
    fprintf(fp, " | ");
    for (j=0; j<step_mem->stages; j++) 
      fprintf(fp, " %"RSYME"",step_mem->Bii->d[j]);
    fprintf(fp, " | ");
    fprintf(fp,"\n");
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
