/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKODE infrastructure; these routines
 * should not be called directly by the user; instead they are
 * provided as utility routines for ARKODE time-step modules
 * to use.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_user_controller.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>


/*===============================================================
  ARKODE optional input utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkSetDefaults:

  Resets all optional inputs to ARKODE default values.  Does not
  change problem-defining function pointers fe and fi or
  user_data pointer.  Also leaves alone any data
  structures/options related to root-finding (those can be reset
  using ARKodeRootInit) or post-processing a step (ProcessStep).
  ---------------------------------------------------------------*/
int arkSetDefaults(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->use_compensated_sums    = SUNFALSE;
  ark_mem->fixedstep               = SUNFALSE;       /* default to use adaptive steps */
  ark_mem->reltol                  = RCONST(1.e-4);  /* relative tolerance */
  ark_mem->itol                    = ARK_SS;         /* scalar-scalar solution tolerances */
  ark_mem->ritol                   = ARK_SS;         /* scalar-scalar residual tolerances */
  ark_mem->Sabstol                 = RCONST(1.e-9);  /* solution absolute tolerance */
  ark_mem->atolmin0                = SUNFALSE;       /* min(abstol) > 0 */
  ark_mem->SRabstol                = RCONST(1.e-9);  /* residual absolute tolerance */
  ark_mem->Ratolmin0               = SUNFALSE;       /* min(Rabstol) > 0 */
  ark_mem->user_efun               = SUNFALSE;       /* no user-supplied ewt function */
  ark_mem->efun                    = arkEwtSetSS;    /* built-in scalar-scalar ewt function */
  ark_mem->e_data                  = ark_mem;        /* ewt function data */
  ark_mem->user_rfun               = SUNFALSE;       /* no user-supplied rwt function */
  ark_mem->rfun                    = arkRwtSet;      /* built-in rwt function */
  ark_mem->r_data                  = ark_mem;        /* rwt function data */
  ark_mem->ehfun                   = arkErrHandler;  /* default error handler fn */
  ark_mem->eh_data                 = ark_mem;        /* error handler data */
  ark_mem->errfp                   = stderr;         /* output stream for errors */
#if SUNDIALS_LOGGING_LEVEL > 0
  ark_mem->errfp                   = (ARK_LOGGER->error_fp) ? ARK_LOGGER->error_fp : stderr;
#endif
  ark_mem->mxstep                  = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->mxhnil                  = MXHNIL;         /* max warns of t+h==t */
  ark_mem->maxnef                  = MAXNEF;         /* max error test fails */
  ark_mem->maxncf                  = MAXNCF;         /* max convergence fails */
  ark_mem->maxconstrfails          = MAXCONSTRFAILS; /* max number of constraint fails */
  ark_mem->hin                     = ZERO;           /* determine initial step on-the-fly */
  ark_mem->hmin                    = ZERO;           /* no minimum step size */
  ark_mem->hmax_inv                = ZERO;           /* no maximum step size */
  ark_mem->tstopset                = SUNFALSE;       /* no stop time set */
  ark_mem->tstopinterp             = SUNFALSE;       /* copy at stop time */
  ark_mem->tstop                   = ZERO;           /* no fixed stop time */
  ark_mem->hadapt_mem->etamx1      = ETAMX1;         /* max change on first step */
  ark_mem->hadapt_mem->etamxf      = ETAMXF;         /* max change on error-failed step */
  ark_mem->hadapt_mem->etamin      = ETAMIN;         /* min bound on time step reduction */
  ark_mem->hadapt_mem->small_nef   = SMALL_NEF;      /* num error fails before ETAMXF enforced */
  ark_mem->hadapt_mem->etacf       = ETACF;          /* max change on convergence failure */
  ark_mem->hadapt_mem->cfl         = CFLFAC;         /* explicit stability factor */
  ark_mem->hadapt_mem->safety      = SAFETY;         /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->growth      = GROWTH;         /* step adaptivity growth factor */
  ark_mem->hadapt_mem->lbound      = HFIXED_LB;      /* step adaptivity no-change lower bound */
  ark_mem->hadapt_mem->ubound      = HFIXED_UB;      /* step adaptivity no-change upper bound */
  ark_mem->hadapt_mem->expstab     = arkExpStab;     /* internal explicit stability fn */
  ark_mem->hadapt_mem->estab_data  = NULL;           /* no explicit stability fn data */
  ark_mem->hadapt_mem->pq          = PQ;             /* embedding order */
  ark_mem->hadapt_mem->p           = 0;              /* no default embedding order */
  ark_mem->hadapt_mem->q           = 0;              /* no default method order */
  ark_mem->hadapt_mem->adjust      = ADJUST;         /* controller order adjustment */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInterpolantType:

  Specifies use of the Lagrange or Hermite interpolation modules.
    itype == ARK_INTERP_HERMITE specifies the Hermite (nonstiff)
      interpolation module.
    itype == ARK_INTERP_LAGRANGE specifies the Lagrange (stiff)
      interpolation module.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input.
     ARK_MEM_FAIL if the interpolation module cannot be allocated.
     ARK_ILL_INPUT if the itype argument is not recognized.
  ---------------------------------------------------------------*/
int arkSetInterpolantType(void *arkode_mem, int itype)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetInterpolantType", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal itype input */
  if ((itype != ARK_INTERP_HERMITE) && (itype != ARK_INTERP_LAGRANGE)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetInterpolantType",
                    "Illegal interpolation type input.");
    return(ARK_ILL_INPUT);
  }

  /* do not change type once the module has been initialized */
  if (ark_mem->initialized) {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKODE",
                    "arkSetInterpolantType",
                    "Type cannot be specified after module initialization.");
    return(ARK_ILL_INPUT);
  }

  /* delete any existing interpolation module */
  if (ark_mem->interp != NULL) {
    arkInterpFree(ark_mem, ark_mem->interp);
    ark_mem->interp = NULL;
  }

  /* create requested interpolation module, initially specifying
     the maximum possible interpolant degree. */
  if (itype == ARK_INTERP_HERMITE)
  {
    ark_mem->interp = arkInterpCreate_Hermite(arkode_mem, ARK_INTERP_MAX_DEGREE);
    ark_mem->interp_type = ARK_INTERP_HERMITE;
  }
  else if (itype == ARK_INTERP_LAGRANGE)
  {
    ark_mem->interp = arkInterpCreate_Lagrange(arkode_mem, ARK_INTERP_MAX_DEGREE);
    ark_mem->interp_type = ARK_INTERP_LAGRANGE;
  }
  else
  {
    ark_mem->interp = NULL;
    ark_mem->interp_type = -1;
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetInterpolantType",
                    "Unable to allocate interpolation structure");
    return(ARK_MEM_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInterpolantDegree:

  Specifies the polynomial degree for the dense output
  interpolation module.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input or nonexistent
       interpolation module.
     ARK_INTERP_FAIL if the interpolation module is already
       initialized.
     ARK_ILL_INPUT if the degree is illegal.
  ---------------------------------------------------------------*/
int arkSetInterpolantDegree(void *arkode_mem, int degree)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetInterpolantDegree", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE",
                    "arkSetInterpolantDegree",
                    "Interpolation module is not yet allocated");
    return(ARK_MEM_NULL);
  }

  /* do not change degree once the module has been initialized */
  if (ark_mem->initialized) {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKODE",
                    "arkSetInterpolantType",
                    "Degree cannot be specified after module initialization.");
    return(ARK_ILL_INPUT);
  }

  /* pass 'degree' to interpolation module, returning its value */
  return(arkInterpSetDegree(ark_mem, ark_mem->interp, degree));
}


/*---------------------------------------------------------------
  arkSetErrHandlerFn:

  Specifies the error handler function
  ---------------------------------------------------------------*/
int arkSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                       void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided values, or defaults, depending on argument */
  if (ehfun == NULL) {
    ark_mem->ehfun   = arkErrHandler;
    ark_mem->eh_data = ark_mem;
  } else {
    ark_mem->ehfun   = ehfun;
    ark_mem->eh_data = eh_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetErrFile:

  Specifies the FILE pointer for output (NULL means no messages)
  ---------------------------------------------------------------*/
int arkSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->errfp = errfp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetUserData:

  Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int arkSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->user_data = user_data;

  /* Set data for efun */
  if (ark_mem->user_efun)
    ark_mem->e_data = user_data;

  /* Set data for rfun */
  if (ark_mem->user_rfun)
    ark_mem->r_data = user_data;

  /* Set data for root finding */
  if (ark_mem->root_mem != NULL)
    ark_mem->root_mem->root_data = user_data;

  /* Set data for post-processing a step */
  if (ark_mem->ProcessStep != NULL)
    ark_mem->ps_data = user_data;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkSetAdaptController:

  Specifies a non-default SUNAdaptController time step controller
  object. If a NULL-valued SUNAdaptController is input, the
  default will be re-enabled.
  ---------------------------------------------------------------*/
int arkSetAdaptController(void *arkode_mem, SUNAdaptController C)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetAdaptController", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  if (ark_mem->hadapt_mem->owncontroller) {
    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUNADAPTCONTROLLER_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptController",
                      "SUNAdaptController_Destroy failure");
      return(ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* On NULL-valued input, create default SUNAdaptController object */
  if (C == NULL) {
    C = SUNAdaptController_PID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptController",
                      "SUNAdaptControllerPID allocation failure");
      return(ARK_MEM_FAIL);
    }
    ark_mem->hadapt_mem->owncontroller = SUNTRUE;
  } else {
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
  }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller = C;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkSetMaxNumSteps:

  Specifies the maximum number of integration steps
  ---------------------------------------------------------------*/
int arkSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0)
    ark_mem->mxstep = MXSTEP_DEFAULT;
  else
    ark_mem->mxstep = mxsteps;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxHnilWarns:

  Specifies the maximum number of warnings for small h
  ---------------------------------------------------------------*/
int arkSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxhnil=0 sets the default, otherwise use input. */
  if (mxhnil == 0) {
    ark_mem->mxhnil = 10;
  } else {
    ark_mem->mxhnil = mxhnil;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInitStep:

  Specifies the initial step size
  ---------------------------------------------------------------*/
int arkSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;
  int retval;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) {
    ark_mem->hin = ZERO;
  } else {
    ark_mem->hin = hin;
  }

  /* Clear previous initial step */
  ark_mem->h0u = ZERO;

  /* Reset error controller (e.g., error and step size history) */
  retval = SUNAdaptController_Reset(ark_mem->hadapt_mem->hcontroller);
  if (retval != SUNADAPTCONTROLLER_SUCCESS) { return(ARK_CONTROLLER_ERR); }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMinStep:

  Specifies the minimum step size
  ---------------------------------------------------------------*/
int arkSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmin = 0 */
  if (hmin <= ZERO) {
    ark_mem->hmin = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * ark_mem->hmax_inv > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetMinStep", MSG_ARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmin = hmin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxStep:

  Specifies the maximum step size
  ---------------------------------------------------------------*/
int arkSetMaxStep(void *arkode_mem, realtype hmax)
{
  realtype hmax_inv;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= ZERO) {
    ark_mem->hmax_inv = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmax and hmin are agreeable */
  hmax_inv = ONE/hmax;
  if (hmax_inv * ark_mem->hmin > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetMaxStep", MSG_ARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmax_inv = hmax_inv;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetStopTime:

  Specifies the time beyond which the integration is not to proceed.
  ---------------------------------------------------------------*/
int arkSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetStopTime", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If ARKODE was called at least once, test if tstop is legal
     (i.e. if it was not already passed).
     If arkSetStopTime is called before the first call to ARKODE,
     tstop will be checked in ARKODE. */
  if (ark_mem->nst > 0) {
    if ( (tstop - ark_mem->tcur) * ark_mem->h < ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                      "arkSetStopTime", MSG_ARK_BAD_TSTOP,
                      tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  ark_mem->tstop    = tstop;
  ark_mem->tstopset = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInterpolateStopTime:

  Specifies to use interpolation to fill the solution output at
  the stop time (instead of a copy).
  ---------------------------------------------------------------*/
int arkSetInterpolateStopTime(void *arkode_mem, booleantype interp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetInterpolateStopTime", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->tstopinterp = interp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkClearStopTime:

  Disable the stop time.
  ---------------------------------------------------------------*/
int arkClearStopTime(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkClearStopTime", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->tstopset = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetFixedStep:

  Specifies to use a fixed time step size instead of performing
  any form of temporal adaptivity.  ARKODE will use this step size
  for all steps (unless tstop is set, in which case it may need to
  modify that last step approaching tstop.  If any solver failure
  occurs in the timestepping module, ARKODE will typically
  immediately return with an error message indicating that the
  selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int arkSetFixedStep(void *arkode_mem, realtype hfixed)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetFixedStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* re-attach internal error weight functions if necessary */
  if ((hfixed == ZERO) && (!ark_mem->user_efun)) {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
      retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    else
      retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    if (retval != ARK_SUCCESS) return(retval);
  }

  /* set ark_mem "fixedstep" entry */
  if (hfixed != ZERO) {
    ark_mem->fixedstep = SUNTRUE;
    ark_mem->hin = hfixed;
  } else {
    ark_mem->fixedstep = SUNFALSE;
  }

  /* Notify ARKODE to use hfixed as the initial step size, and return */
  retval = arkSetInitStep(arkode_mem, hfixed);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetRootDirection:

  Specifies the direction of zero-crossings to be monitored.
  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int arkSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  int i;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE",
                    "arkSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  if (ark_root_mem->nrtfn == 0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetRootDirection", MSG_ARK_NO_ROOT);
    return(ARK_ILL_INPUT);
  }
  for(i=0; i<ark_root_mem->nrtfn; i++)
    ark_root_mem->rootdir[i] = rootdir[i];
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetNoInactiveRootWarn:

  Disables issuing a warning if some root function appears
  to be identically zero at the beginning of the integration
  ---------------------------------------------------------------*/
int arkSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE",
                    "arkSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  ark_root_mem->mxgnull = 0;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetPostprocessStepFn:

  Specifies a user-provided step postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables step
  postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int arkSetPostprocessStepFn(void *arkode_mem,
                            ARKPostProcessFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStep = ProcessStep;
  ark_mem->ps_data     = ark_mem->user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetPostprocessStageFn:

  Specifies a user-provided stage postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables
  stage postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.

  While it is possible to perform postprocessing when
  ARKStepSetDeduceImplicitRhs is enabled, this can cause implicit
  RHS evaluations to be inconsistent with the postprocessed stage
  values.  It is strongly recommended to disable
  ARKStepSetDeduceImplicitRhs in order to guarantee
  postprocessing constraints are enforced.
  ---------------------------------------------------------------*/
int arkSetPostprocessStageFn(void *arkode_mem,
                             ARKPostProcessFn ProcessStage)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetPostprocessStageFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStage = ProcessStage;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetConstraints:

  Activates or Deactivates inequality constraint checking.
  ---------------------------------------------------------------*/
int arkSetConstraints(void *arkode_mem, N_Vector constraints)
{
  realtype temptest;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetConstraints", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If there are no constraints, destroy data structures */
  if (constraints == NULL) {
    arkFreeVec(ark_mem, &ark_mem->constraints);
    ark_mem->constraintsSet = SUNFALSE;
    return(ARK_SUCCESS);
  }

  /* Test if required vector ops. are defined */
  if (constraints->ops->nvdiv         == NULL ||
      constraints->ops->nvmaxnorm     == NULL ||
      constraints->ops->nvcompare     == NULL ||
      constraints->ops->nvconstrmask  == NULL ||
      constraints->ops->nvminquotient == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetConstraints", MSG_ARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Check the constraints vector */
  temptest = N_VMaxNorm(constraints);
  if ((temptest > RCONST(2.5)) || (temptest < HALF)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetConstraints", MSG_ARK_BAD_CONSTR);
    return(ARK_ILL_INPUT);
  }

  /* Allocate the internal constrains vector (if necessary) */
  if (!arkAllocVec(ark_mem, constraints, &ark_mem->constraints))
    return(ARK_MEM_FAIL);

  /* Load the constraints vector */
  N_VScale(ONE, constraints, ark_mem->constraints);
  ark_mem->constraintsSet = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxNumConstrFails:

  Set max number of allowed constraint failures in a step before
  returning an error
  ---------------------------------------------------------------*/
int arkSetMaxNumConstrFails(void *arkode_mem, int maxfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxNumConstrFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing maxfails = 0 sets the default, otherwise set to input */
  if (maxfails <= 0)
    ark_mem->maxconstrfails = MAXCONSTRFAILS;
  else
    ark_mem->maxconstrfails = maxfails;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetAdaptivityMethod:  ***DEPRECATED***

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.

  Users should transition to constructing non-default SUNAdaptController
  objects directly, and providing those directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
int arkSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault,
                           int pq, realtype adapt_params[3])
{
  int retval;
  long int lenrw, leniw;
  realtype k1, k2, k3;
  ARKodeMem ark_mem;
  SUNAdaptController C;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetController", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for illegal inputs */
  if ((idefault != 1) && (adapt_params == NULL)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetAdaptivityMethod",
                    "NULL-valued adapt_params provided");
    return(ARK_ILL_INPUT);
  }

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  if (ark_mem->hadapt_mem->owncontroller) {
    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUNADAPTCONTROLLER_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_Destroy failure");
      return(ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* set adaptivity parameters from inputs */
  k1 = k2 = k3 = ZERO;
  if (idefault != 1) {
    k1 = adapt_params[0];
    k2 = adapt_params[1];
    k3 = adapt_params[2];
  }
  ark_mem->hadapt_mem->pq = pq;

  /* Create new SUNAdaptController object based on "imethod" input, optionally setting
     the specified controller parameters */
  C = NULL;
  switch (imethod) {
  case (ARK_ADAPT_PID):
    C = SUNAdaptController_PID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_PID allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_PID(C, k1, -k2, k3);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_PID failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_PI):
    C = SUNAdaptController_PI(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_PI allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_PI(C, k1, -k2);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_PI failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_I):
    C = SUNAdaptController_I(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_I allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_I(C, k1);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_I failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_EXP_GUS):
    C = SUNAdaptController_ExpGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_ExpGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_ExpGus(C, k1, k2);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_ExpGus failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_IMP_GUS):
    C = SUNAdaptController_ImpGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_ImpGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_ImpGus(C, k1, k2);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_ImpGus failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_IMEX_GUS):
    C = SUNAdaptController_ImExGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNAdaptController_ImExGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    if (idefault != 1) {
      retval = SUNAdaptController_SetParams_ImExGus(C, k1, k2, k3, k3);
      if (retval != SUNADAPTCONTROLLER_SUCCESS) {
        (void) SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                        "arkSetAdaptivityMethod", "SUNAdaptController_SetParams_ImExGus failure");
        return(ARK_CONTROLLER_ERR);
      }
    }
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller = C;
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetAdaptivityFn:  ***DEPRECATED***

  Specifies the user-provided time step adaptivity function to use.
  If 'hfun' is NULL-valued, then the default PID controller will
  be used instead.

  Users should transition to constructing a custom SUNAdaptController
  object, and providing this directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
int arkSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  SUNAdaptController C;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetAdaptivityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  if (ark_mem->hadapt_mem->owncontroller) {
    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUNADAPTCONTROLLER_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityFn",
                      "SUNAdaptController_Destroy failure");
      return(ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* Create new SUNAdaptController object depending on NULL-ity of 'hfun' */
  C = NULL;
  if (hfun == NULL) {
    C = SUNAdaptController_PID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityFn",
                      "SUNAdaptController_PID allocation failure");
      return(ARK_MEM_FAIL);
    }
  } else {
    C = ARKUserControl(ark_mem->sunctx, arkode_mem, hfun, h_data);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityFn",
                      "ARKUserControl allocation failure");
      return(ARK_MEM_FAIL);
    }
  }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUNADAPTCONTROLLER_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller = C;
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetCFLFraction",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for allowable parameters */
  if (cfl_frac >= ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetCFLFraction", "Illegal CFL fraction");
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
  arkSetAdaptivityAdjustment:

  Adjusts the method order supplied to the temporal adaptivity
  controller.  For example, if the user expects order reduction
  due to problem stiffness, they may request that the controller
  assume a reduced order of accuracy for the method by specifying
  a value adjust < 0.
  ---------------------------------------------------------------*/
int arkSetAdaptivityAdjustment(void *arkode_mem, int adjust)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetAdaptivityAdjustment",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* store requested adjustment */
  hadapt_mem->adjust = adjust;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int arkSetSafetyFactor(void *arkode_mem, realtype safety)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetSafetyFactor",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for allowable parameters */
  if (safety >= ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetSafetyFactor", "Illegal safety factor");
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
  arkSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetErrorBias(void *arkode_mem, realtype bias)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetErrorBias",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowed value, otherwise set default */
  if (bias < ONE) {
    retval = SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller, -1.0);
  } else {
    retval = SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller, bias);
  }
  if (retval != SUNADAPTCONTROLLER_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetErrorBias", "SUNAdaptController_SetErrorBias failure");
      return(ARK_CONTROLLER_ERR);
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int arkSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowed value, otherwise set default */
  if (mx_growth <= ONE) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMinReduction:

  Specifies the minimum possible step size reduction factor to be
  allowed between successive integration steps. Allowable values
  must be > 0.0 and < 1.0. Any illegal value implies a reset to
  the default.
  ---------------------------------------------------------------*/
int arkSetMinReduction(void *arkode_mem, realtype eta_min)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMinReduction",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set allowed value, otherwise set default */
  if (eta_min >= ONE || eta_min <= ZERO) {
    hadapt_mem->etamin = ETAMIN;
  } else {
    hadapt_mem->etamin = eta_min;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int arkSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetFixedStepBounds",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowable interval, otherwise set defaults */
  if ((lb <= ONE) && (ub >= ONE)) {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  } else {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxFirstGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxEFailGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetSmallNumEFails",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxCFailGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int arkSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetStabilityFn",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

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
  arkSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    ark_mem->maxnef = MAXNEF;
  } else {
    ark_mem->maxnef = maxnef;
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetMaxConvFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    ark_mem->maxncf = MAXNCF;
  } else {
    ark_mem->maxncf = maxncf;
  }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkSetUseCompensatedSums:

  Specifies that ARKODE should use compensated (Kahan) summation
  where relevant to mitigate roundoff error.
  ---------------------------------------------------------------*/
int arkSetUseCompensatedSums(void *arkode_mem, sunbooleantype onoff)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetUseCompensatedSums", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (onoff) {
    ark_mem->use_compensated_sums = SUNTRUE;
  } else {
    ark_mem->use_compensated_sums = SUNFALSE;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE optional output utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkGetNumStepAttempts:

   Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int arkGetNumStepAttempts(void *arkode_mem, long int *nstep_attempts)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumStepAttempts", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nstep_attempts = ark_mem->nst_attempts;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumSteps:

  Returns the current number of integration steps
  ---------------------------------------------------------------*/
int arkGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->nst;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetActualInitStep:

  Returns the step size used on the first step
  ---------------------------------------------------------------*/
int arkGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hinused = ark_mem->h0u;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetLastStep:

  Returns the step size used on the last successful step
  ---------------------------------------------------------------*/
int arkGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hlast = ark_mem->hold;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentStep:

  Returns the step size to be attempted on the next step
  ---------------------------------------------------------------*/
int arkGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hcur = ark_mem->next_h;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentState:

  Returns the current solution (before or after as step) or
  stage value (during step solve).
  ---------------------------------------------------------------*/
int arkGetCurrentState(void *arkode_mem, N_Vector *state)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetCurrentState", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *state = ark_mem->ycur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentTime:

  Returns the current value of the independent variable
  ---------------------------------------------------------------*/
int arkGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tcur = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetTolScaleFactor:

  Returns a suggested factor for scaling tolerances
  ---------------------------------------------------------------*/
int arkGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tolsfact = ark_mem->tolsf;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetErrWeights:

  This routine returns the current error weight vector.
  ---------------------------------------------------------------*/
int arkGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ewt, eweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetResWeights:

  This routine returns the current residual weight vector.
  ---------------------------------------------------------------*/
int arkGetResWeights(void *arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetResWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->rwt, rweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetWorkSpace:

  Returns integrator work space requirements
  ---------------------------------------------------------------*/
int arkGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *leniw = ark_mem->liw;
  *lenrw = ark_mem->lrw;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumGEvals:

  Returns the current number of calls to g (for rootfinding)
  ---------------------------------------------------------------*/
int arkGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  *ngevals = ark_root_mem->nge;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetRootInfo:

  Returns pointer to array rootsfound showing roots found
  ---------------------------------------------------------------*/
int arkGetRootInfo(void *arkode_mem, int *rootsfound)
{
  int i;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKODE",
                    "arkGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  for (i=0; i<ark_root_mem->nrtfn; i++)
    rootsfound[i] = ark_root_mem->iroots[i];
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetStepStats:

  Returns step statistics
  ---------------------------------------------------------------*/
int arkGetStepStats(void *arkode_mem, long int *nsteps,
                    realtype *hinused, realtype *hlast,
                    realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps  = ark_mem->nst;
  *hinused = ark_mem->h0u;
  *hlast   = ark_mem->hold;
  *hcur    = ark_mem->next_h;
  *tcur    = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumConstrFails:

  Returns the current number of constraint fails
  ---------------------------------------------------------------*/
int arkGetNumConstrFails(void *arkode_mem, long int *nconstrfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumConstrFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nconstrfails = ark_mem->nconstrfails;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int arkGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumExpSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_exp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int arkGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumAccSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_acc;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int arkGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *netfails = ark_mem->netf;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumStepSolveFails:

  Returns the current number of failed steps due to an algebraic
  solver convergence failure.
  ---------------------------------------------------------------*/
int arkGetNumStepSolveFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetNumStepSolveFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nncfails = ark_mem->ncfn;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetUserData:

  Returns the user data pointer
  ---------------------------------------------------------------*/
int arkGetUserData(void *arkode_mem, void** user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *user_data = ark_mem->user_data;

  return(ARK_SUCCESS);
}


/*-----------------------------------------------------------------
  arkPrintAllStats

  Prints the current value of all statistics
  ---------------------------------------------------------------*/

int arkPrintAllStats(void *arkode_mem, FILE *outfile, SUNOutputFormat fmt)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkPrintAllStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  switch(fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Current time                 = %"RSYM"\n", ark_mem->tcur);
    fprintf(outfile, "Steps                        = %ld\n", ark_mem->nst);
    fprintf(outfile, "Step attempts                = %ld\n", ark_mem->nst_attempts);
    fprintf(outfile, "Stability limited steps      = %ld\n", ark_mem->hadapt_mem->nst_exp);
    fprintf(outfile, "Accuracy limited steps       = %ld\n", ark_mem->hadapt_mem->nst_acc);
    fprintf(outfile, "Error test fails             = %ld\n", ark_mem->netf);
    fprintf(outfile, "NLS step fails               = %ld\n", ark_mem->ncfn);
    fprintf(outfile, "Inequality constraint fails  = %ld\n", ark_mem->nconstrfails);
    fprintf(outfile, "Initial step size            = %"RSYM"\n", ark_mem->h0u);
    fprintf(outfile, "Last step size               = %"RSYM"\n", ark_mem->hold);
    fprintf(outfile, "Current step size            = %"RSYM"\n", ark_mem->next_h);
    if (ark_mem->root_mem)
    {
      ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
      fprintf(outfile, "Root fn evals                = %ld\n", ark_root_mem->nge);
    }
    break;
  case SUN_OUTPUTFORMAT_CSV:
    fprintf(outfile, "Time,%"RSYM, ark_mem->tcur);
    fprintf(outfile, ",Steps,%ld", ark_mem->nst);
    fprintf(outfile, ",Step attempts,%ld", ark_mem->nst_attempts);
    fprintf(outfile, ",Stability limited steps,%ld", ark_mem->hadapt_mem->nst_exp);
    fprintf(outfile, ",Accuracy limited steps,%ld", ark_mem->hadapt_mem->nst_acc);
    fprintf(outfile, ",Error test fails,%ld", ark_mem->netf);
    fprintf(outfile, ",NLS step fails,%ld", ark_mem->ncfn);
    fprintf(outfile, ",Inequality constraint fails,%ld", ark_mem->nconstrfails);
    fprintf(outfile, ",Initial step size,%"RSYM, ark_mem->h0u);
    fprintf(outfile, ",Last step size,%"RSYM, ark_mem->hold);
    fprintf(outfile, ",Current step size,%"RSYM, ark_mem->next_h);
    if (ark_mem->root_mem)
    {
      ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
      fprintf(outfile, ",Roof fn evals,%ld", ark_root_mem->nge);
    }
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkPrintAllStats",
                    "Invalid formatting option.");
    return(ARK_ILL_INPUT);
  }

  /* Print relaxation stats */
  if (ark_mem->relax_enabled)
  {
    retval = arkRelaxPrintAllStats(arkode_mem, outfile, fmt);
    if (retval != ARK_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*-----------------------------------------------------------------*/

char *arkGetReturnFlagName(long int flag)
{
  char *name;
  name = (char *)malloc(27*sizeof(char));

  switch(flag) {
  case ARK_SUCCESS:
    sprintf(name,"ARK_SUCCESS");
    break;
  case ARK_TSTOP_RETURN:
    sprintf(name,"ARK_TSTOP_RETURN");
    break;
  case ARK_ROOT_RETURN:
    sprintf(name,"ARK_ROOT_RETURN");
    break;
  case ARK_WARNING:
    sprintf(name,"ARK_WARNING");
    break;
  case ARK_TOO_MUCH_WORK:
    sprintf(name,"ARK_TOO_MUCH_WORK");
    break;
  case ARK_TOO_MUCH_ACC:
    sprintf(name,"ARK_TOO_MUCH_ACC");
    break;
  case ARK_ERR_FAILURE:
    sprintf(name,"ARK_ERR_FAILURE");
    break;
  case ARK_CONV_FAILURE:
    sprintf(name,"ARK_CONV_FAILURE");
    break;
  case ARK_LINIT_FAIL:
    sprintf(name,"ARK_LINIT_FAIL");
    break;
  case ARK_LSETUP_FAIL:
    sprintf(name,"ARK_LSETUP_FAIL");
    break;
  case ARK_LSOLVE_FAIL:
    sprintf(name,"ARK_LSOLVE_FAIL");
    break;
  case ARK_RHSFUNC_FAIL:
    sprintf(name,"ARK_RHSFUNC_FAIL");
    break;
  case ARK_FIRST_RHSFUNC_ERR:
    sprintf(name,"ARK_FIRST_RHSFUNC_ERR");
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    sprintf(name,"ARK_REPTD_RHSFUNC_ERR");
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    sprintf(name,"ARK_UNREC_RHSFUNC_ERR");
    break;
  case ARK_RTFUNC_FAIL:
    sprintf(name,"ARK_RTFUNC_FAIL");
    break;
  case ARK_LFREE_FAIL:
    sprintf(name,"ARK_LFREE_FAIL");
    break;
  case ARK_MASSINIT_FAIL:
    sprintf(name,"ARK_MASSINIT_FAIL");
    break;
  case ARK_MASSSETUP_FAIL:
    sprintf(name,"ARK_MASSSETUP_FAIL");
    break;
  case ARK_MASSSOLVE_FAIL:
    sprintf(name,"ARK_MASSSOLVE_FAIL");
    break;
  case ARK_MASSFREE_FAIL:
    sprintf(name,"ARK_MASSFREE_FAIL");
    break;
  case ARK_MASSMULT_FAIL:
    sprintf(name,"ARK_MASSMULT_FAIL");
    break;
  case ARK_CONSTR_FAIL:
    sprintf(name,"ARK_CONSTR_FAIL");
    break;
  case ARK_MEM_FAIL:
    sprintf(name,"ARK_MEM_FAIL");
    break;
  case ARK_MEM_NULL:
    sprintf(name,"ARK_MEM_NULL");
    break;
  case ARK_ILL_INPUT:
    sprintf(name,"ARK_ILL_INPUT");
    break;
  case ARK_NO_MALLOC:
    sprintf(name,"ARK_NO_MALLOC");
    break;
  case ARK_BAD_K:
    sprintf(name,"ARK_BAD_K");
    break;
  case ARK_BAD_T:
    sprintf(name,"ARK_BAD_T");
    break;
  case ARK_BAD_DKY:
    sprintf(name,"ARK_BAD_DKY");
    break;
  case ARK_TOO_CLOSE:
    sprintf(name,"ARK_TOO_CLOSE");
    break;
  case ARK_VECTOROP_ERR:
    sprintf(name,"ARK_VECTOROP_ERR");
    break;
  case ARK_NLS_INIT_FAIL:
    sprintf(name,"ARK_NLS_INIT_FAIL");
    break;
  case ARK_NLS_SETUP_FAIL:
    sprintf(name,"ARK_NLS_SETUP_FAIL");
    break;
  case ARK_NLS_SETUP_RECVR:
    sprintf(name,"ARK_NLS_SETUP_RECVR");
    break;
  case ARK_NLS_OP_ERR:
    sprintf(name,"ARK_NLS_OP_ERR");
    break;
  case ARK_INNERSTEP_ATTACH_ERR:
    sprintf(name,"ARK_INNERSTEP_ATTACH_ERR");
    break;
  case ARK_INNERSTEP_FAIL:
    sprintf(name,"ARK_INNERSTEP_FAIL");
    break;
  case ARK_OUTERTOINNER_FAIL:
    sprintf(name,"ARK_OUTERTOINNER_FAIL");
    break;
  case ARK_INNERTOOUTER_FAIL:
    sprintf(name,"ARK_INNERTOOUTER_FAIL");
    break;
  case ARK_POSTPROCESS_STEP_FAIL:
    sprintf(name,"ARK_POSTPROCESS_STEP_FAIL");
    break;
  case ARK_POSTPROCESS_STAGE_FAIL:
    sprintf(name,"ARK_POSTPROCESS_STAGE_FAIL");
    break;
  case ARK_USER_PREDICT_FAIL:
    sprintf(name,"ARK_USER_PREDICT_FAIL");
    break;
  case ARK_INTERP_FAIL:
    sprintf(name,"ARK_INTERP_FAIL");
    break;
  case ARK_INVALID_TABLE:
    sprintf(name,"ARK_INVALID_TABLE");
    break;
  case ARK_CONTEXT_ERR:
    sprintf(name,"ARK_CONTEXT_ERR");
    break;
  case ARK_CONTROLLER_ERR:
    sprintf(name,"ARK_CONTROLLER_ERR");
    break;
  case ARK_UNRECOGNIZED_ERROR:
    sprintf(name,"ARK_UNRECOGNIZED_ERROR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}



/*===============================================================
  ARKODE parameter output utility routine
  ===============================================================*/

/*---------------------------------------------------------------
  arkodeWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int arkWriteParameters(ARKodeMem ark_mem, FILE *fp)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkWriteParameters", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* print integrator parameters to file */
  fprintf(fp, "ARKODE solver parameters:\n");
  if (ark_mem->hmin != ZERO)
    fprintf(fp, "  Minimum step size = %" RSYM"\n",ark_mem->hmin);
  if (ark_mem->hmax_inv != ZERO)
    fprintf(fp, "  Maximum step size = %" RSYM"\n",ONE/ark_mem->hmax_inv);
  if (ark_mem->fixedstep)
    fprintf(fp, "  Fixed time-stepping enabled\n");
  if (ark_mem->itol == ARK_WF) {
    fprintf(fp, "  User provided error weight function\n");
  } else {
    fprintf(fp, "  Solver relative tolerance = %" RSYM"\n", ark_mem->reltol);
    if (ark_mem->itol == ARK_SS) {
      fprintf(fp, "  Solver absolute tolerance = %" RSYM"\n", ark_mem->Sabstol);
    } else {
      fprintf(fp, "  Vector-valued solver absolute tolerance\n");
    }
  }
  if (!ark_mem->rwt_is_ewt) {
    if (ark_mem->ritol == ARK_WF) {
      fprintf(fp, "  User provided residual weight function\n");
    } else {
      if (ark_mem->ritol == ARK_SS) {
        fprintf(fp, "  Absolute residual tolerance = %" RSYM"\n", ark_mem->SRabstol);
      } else {
        fprintf(fp, "  Vector-valued residual absolute tolerance\n");
      }
    }
  }
  if (ark_mem->hin != ZERO)
    fprintf(fp, "  Initial step size = %" RSYM"\n",ark_mem->hin);
  fprintf(fp, "\n");
  fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
          ark_mem->hadapt_mem->etamx1);
  fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
          ark_mem->hadapt_mem->etamxf);
  fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
          ark_mem->hadapt_mem->small_nef);
  fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
          ark_mem->hadapt_mem->etacf);
  fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
          ark_mem->hadapt_mem->cfl);
  fprintf(fp, "  Safety factor = %"RSYM"\n", ark_mem->hadapt_mem->safety);
  fprintf(fp, "  Growth factor = %"RSYM"\n", ark_mem->hadapt_mem->growth);
  fprintf(fp, "  Step growth lower bound = %"RSYM"\n", ark_mem->hadapt_mem->lbound);
  fprintf(fp, "  Step growth upper bound = %"RSYM"\n", ark_mem->hadapt_mem->ubound);
  if (ark_mem->hadapt_mem->expstab == arkExpStab) {
    fprintf(fp, "  Default explicit stability function\n");
   } else {
    fprintf(fp, "  User provided explicit stability function\n");
  }
  (void) SUNAdaptController_Write(ark_mem->hadapt_mem->hcontroller, fp);

  fprintf(fp, "  Maximum number of error test failures = %i\n",ark_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",ark_mem->maxncf);

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE + XBraid interface utility functions
  ===============================================================*/


/*---------------------------------------------------------------
  arkSetForcePass:

  Ignore the value of kflag after the temporal error test and
  force the step to pass.
  ---------------------------------------------------------------*/
int arkSetForcePass(void *arkode_mem, booleantype force_pass)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetForcePass", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->force_pass = force_pass;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetLastKFlag:

  The last kflag value retured by the temporal error test.
  ---------------------------------------------------------------*/
int arkGetLastKFlag(void *arkode_mem, int *last_kflag)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkGetLastKFlag", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *last_kflag = ark_mem->last_kflag;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
