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
#include <suncontrol/suncontrol_pid.h>
#include <suncontrol/suncontrol_pi.h>
#include <suncontrol/suncontrol_i.h>
#include <suncontrol/suncontrol_expgus.h>
#include <suncontrol/suncontrol_impgus.h>
#include <suncontrol/suncontrol_imexgus.h>
#include <sunheuristics/sunheuristics_default.h>


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
  int retval;
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
  ark_mem->tstopset                = SUNFALSE;       /* no stop time set */
  ark_mem->tstopinterp             = SUNFALSE;       /* copy at stop time */
  ark_mem->tstop                   = ZERO;           /* no fixed stop time */
  ark_mem->diagfp                  = NULL;           /* no solver diagnostics file */
  ark_mem->report                  = SUNFALSE;       /* don't report solver diagnostics */

  /* Set default values for controller object */
  retval = SUNControlSetDefaults(ark_mem->hcontroller);
  if (retval != SUNCONTROL_SUCCESS) { return(ARK_CONTROLLER_ERR); }

  /* Set default values for heuristics object */
  retval = SUNHeuristicsSetDefaults(ark_mem->hconstraints);
  if (retval != SUNHEURISTICS_SUCCESS) { return(ARK_HEURISTICS_ERR); }

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
  if (itype == ARK_INTERP_HERMITE) {
    ark_mem->interp = arkInterpCreate_Hermite(arkode_mem, ARK_INTERP_MAX_DEGREE);
  } else if (itype == ARK_INTERP_LAGRANGE) {
    ark_mem->interp = arkInterpCreate_Lagrange(arkode_mem, ARK_INTERP_MAX_DEGREE);
  } else {
    ark_mem->interp = NULL;
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
  arkSetDiagnostics:

  Specifies to enable solver diagnostics, and specifies the FILE
  pointer for output (diagfp==NULL disables output)
  ---------------------------------------------------------------*/
int arkSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->diagfp = diagfp;
  if (diagfp != NULL) {
    ark_mem->report = SUNTRUE;
  } else {
    ark_mem->report = SUNFALSE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetController:

  Specifies a non-default SUNControl time step controller object.
  If a NULL-valued SUNControl is input, the default will be
  re-enabled.
  ---------------------------------------------------------------*/
int arkSetController(void *arkode_mem, SUNControl C)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetController", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNControl object */
  retval = SUNControlSpace(ark_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNControlDestroy(ark_mem->hcontroller);
  ark_mem->hcontroller = NULL;

  /* On NULL-valued input, create default SUNControl object */
  if (C == NULL) {
    C = SUNControlPID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetController",
                      "SUNControlPID allocation failure");
      return(ARK_MEM_FAIL);
    }
  }

  /* Attach new SUNControl object */
  retval = SUNControlSpace(C, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hcontroller = C;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetHeuristics:

  Specifies a non-default SUNHeuristics time step constraints
  object. If a NULL-valued SUNHeuristics is input, the default
  will be re-enabled.
  ---------------------------------------------------------------*/
int arkSetHeuristics(void *arkode_mem, SUNHeuristics H)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetHeuristics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNHeuristics object */
  retval = SUNHeuristicsSpace(ark_mem->hconstraints, &lenrw, &leniw);
  if (retval == SUNHEURISTICS_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNHeuristicsDestroy(ark_mem->hconstraints);
  ark_mem->hconstraints = NULL;

  /* On NULL-valued input, create default SUNHeuristics object */
  if (H == NULL) {
    H = SUNHeuristicsDefault(ark_mem->sunctx);
    if (H == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetHeuristics",
                      "SUNHeuristicsDefault allocation failure");
      return(ARK_MEM_FAIL);
    }
  }

  /* Attach new SUNHeuristics object */
  retval = SUNHeuristicsSpace(H, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hconstraints = H;

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
  retval = SUNControlReset(ark_mem->hcontroller);
  if (retval != SUNCONTROL_SUCCESS) { return(ARK_CONTROLLER_ERR); }

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
  any form of temporal adaptivity.  Any nonzero argument will
  result in the use of that fixed step size; an argument of 0
  will re-enable temporal adaptivity.  ARKODE will use this step
  size for all steps (unless tstop is set, in which case it may
  need to modify that last step approaching tstop.  If any solver
  failure occurs in the timestepping module, ARKODE will
  typically immediately return with an error message indicating
  that the selected step size cannot be used.
  ---------------------------------------------------------------*/
int arkSetFixedStep(void *arkode_mem, realtype hfixed)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetFixedStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNControl object */
  retval = SUNControlSpace(ark_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNControlDestroy(ark_mem->hcontroller);
  ark_mem->hcontroller = NULL;

  /* If re-enabling time adaptivity, create default PID controller
     and attach object to ARKODE */
  SUNControl C = NULL;
  if (hfixed == 0)
  {
    C = SUNControlPID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetFixedStep",
                      "SUNControlPID allocation failure");
      return(ARK_MEM_FAIL);
    }

    retval = SUNControlSpace(C, &lenrw, &leniw);
    if (retval == SUNCONTROL_SUCCESS) {
      ark_mem->liw += leniw;
      ark_mem->lrw += lenrw;
    }
    ark_mem->hcontroller = C;
  }

  /* re-attach internal error weight functions if necessary */
  if ((hfixed == ZERO) && (!ark_mem->user_efun)) {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
      retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    else
      retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    if (retval != ARK_SUCCESS) return(retval);
  }

  /* set ark_mem "fixedstep" entry, and set
     stepsize bounds into heuristics module */
  if (hfixed != ZERO) {
    ark_mem->fixedstep = SUNTRUE;
    retval = SUNHeuristicsSetMaxStep(ark_mem->hconstraints, hfixed);
    if (retval != ARK_SUCCESS) return(ARK_HEURISTICS_ERR);
    retval = SUNHeuristicsSetMinStep(ark_mem->hconstraints, hfixed);
    if (retval != ARK_SUCCESS) return(ARK_HEURISTICS_ERR);
  } else {
    ark_mem->fixedstep = SUNFALSE;
  }

  /* Notify ARKODE to use hfixed as the initial step size, and return */
  retval = arkSetInitStep(arkode_mem, hfixed);
  return(retval);
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
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep",
                    "ARKStepSetConstraints", MSG_ARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Check the constraints vector */
  temptest = N_VMaxNorm(constraints);
  if ((temptest > RCONST(2.5)) || (temptest < HALF)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep",
                    "ARKStepSetConstraints", MSG_ARK_BAD_CONSTR);
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

  Users should transition to constructing non-default SUNControl
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
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetController", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNControl object */
  retval = SUNControlSpace(ark_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNControlDestroy(ark_mem->hcontroller);
  ark_mem->hcontroller = NULL;

  /* set adaptivity parameters from inputs or signal use of defaults */
  if (idefault == 1) {
    k1 = -RCONST(1.0);
    k2 = -RCONST(1.0);
    k3 = -RCONST(1.0);
  } else {
    k1 = adapt_params[0];
    k2 = adapt_params[1];
    k3 = adapt_params[2];
  }

  /* Create new SUNControl object based on "imethod" input, optionally setting
     the specified controller parameters */
  SUNControl C = NULL;
  switch (imethod) {
  case (ARK_ADAPT_PID):
    C = SUNControlPID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlPID allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlPID_SetParams(C, pq, k1, k2, k3);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlPID_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  case (ARK_ADAPT_PI):
    C = SUNControlPI(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlPI allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlPI_SetParams(C, pq, k1, k2);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlPI_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  case (ARK_ADAPT_I):
    C = SUNControlI(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlI allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlI_SetParams(C, pq, k1);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlI_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  case (ARK_ADAPT_EXP_GUS):
    C = SUNControlExpGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlExpGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlExpGus_SetParams(C, pq, k1, k2);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlExpGus_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  case (ARK_ADAPT_IMP_GUS):
    C = SUNControlImpGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlImpGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlImpGus_SetParams(C, pq, k1, k2);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlImpGus_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  case (ARK_ADAPT_IMEX_GUS):
    C = SUNControlImExGus(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityMethod",
                      "SUNControlImExGus allocation failure");
      return(ARK_MEM_FAIL);
    }
    retval = SUNControlImExGus_SetParams(C, pq, k1, k2, k3, k3);
    if (retval != SUNCONTROL_SUCCESS) {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, "ARKODE",
                      "arkSetAdaptivityMethod", "SUNControlImExGus_SetParams failure");
      return(ARK_CONTROLLER_ERR);
    }
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
                    "arkSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* Attach new SUNControl object */
  retval = SUNControlSpace(C, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hcontroller = C;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetAdaptivityFn:  ***DEPRECATED***

  Specifies the user-provided time step adaptivity function to use.
  If 'hfun' is NULL-valued, then the default PID controller will
  be used instead.

  Users should transition to constructing a custom SUNControl
  object, and providing this directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
int arkSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetAdaptivityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Remove current SUNControl object */
  retval = SUNControlSpace(ark_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  SUNControlDestroy(ark_mem->hcontroller);
  ark_mem->hcontroller = NULL;

  /* Create new SUNControl object depending on NULL-ity of 'hfun' */
  SUNControl C = NULL;
  if (hfun == NULL) {
    C = SUNControlPID(ark_mem->sunctx);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetAdaptivityFn",
                      "SUNControlPID allocation failure");
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

  /* Attach new SUNControl object */
  retval = SUNControlSpace(C, &lenrw, &leniw);
  if (retval == SUNCONTROL_SUCCESS) {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hcontroller = C;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetStabilityFn:  ***DEPRECATED***

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int arkSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetStabilityFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default. */
  if (EStab == NULL) {
    retval = SUNHeuristicsSetExpStabFn(ark_mem->hconstraints, NULL, NULL);
    if (retval != SUNHEURISTICS_SUCCESS) {
      arkProcessError(ark_mem, ARK_HEURISTICS_ERR, "ARKODE",
                      "arkSetStabilityFn", "SUNHeuristicsSetExpStabFn failure");
      return(ARK_HEURISTICS_ERR);
    }
  }

  /* Otherwise create wrapper structure */
  void* wrapper_data = ARKUserStability(arkode_mem, EStab, estab_data);
  if (wrapper_data == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "arkSetStabilityFn",
                    "ARKUserStabilkty failure");
    return(ARK_MEM_FAIL);
  }

  /* Attach ARKODE function and wrapper structure to heuristics object */
  retval = SUNHeuristicsSetExpStabFn(ark_mem->hconstraints,
                                     ARKControlExpStab, wrapper_data);
  if (retval != SUNHEURISTICS_SUCCESS) {
    arkProcessError(ark_mem, ARK_HEURISTICS_ERR, "ARKODE",
                    "arkSetStabilityFn", "SUNHeuristicsSetExpStabFn failure");
    return(ARK_HEURISTICS_ERR);
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
  long int nst_exp, nst_acc;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkPrintAllStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  retval = SUNHeuristicsGetNumExpSteps(ark_mem->hconstraints, &nst_exp);
  if (retval != SUNHEURISTICS_SUCCESS) {
    arkProcessError(ark_mem, ARK_HEURISTICS_ERR, "ARKODE", "arkPrintAllStats",
                    "Error in SUNHeuristicsGetNumExpSteps");
    return(ARK_HEURISTICS_ERR);
  }
  retval = SUNHeuristicsGetNumAccSteps(ark_mem->hconstraints, &nst_acc);
  if (retval != SUNHEURISTICS_SUCCESS) {
    arkProcessError(ark_mem, ARK_HEURISTICS_ERR, "ARKODE", "arkPrintAllStats",
                    "Error in SUNHeuristicsGetNumAccSteps");
    return(ARK_HEURISTICS_ERR);
  }

  switch(fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    fprintf(outfile, "Current time                 = %"RSYM"\n", ark_mem->tcur);
    fprintf(outfile, "Steps                        = %ld\n", ark_mem->nst);
    fprintf(outfile, "Step attempts                = %ld\n", ark_mem->nst_attempts);
    fprintf(outfile, "Stability limited steps      = %ld\n", nst_exp);
    fprintf(outfile, "Accuracy limited steps       = %ld\n", nst_acc);
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
    fprintf(outfile, ",Stability limited steps,%ld", nst_exp);
    fprintf(outfile, ",Accuracy limited steps,%ld", nst_acc);
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
  case ARK_HEURISTICS_ERR:
    sprintf(name,"ARK_HEURISTICS_ERR");
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
  fprintf(fp, "  Maximum number of error test failures = %i\n",ark_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",ark_mem->maxncf);
  fprintf(fp, "\n");
  (void) SUNHeuristicsWrite(ark_mem->hconstraints, fp);
  (void) SUNControlWrite(ark_mem->hcontroller, fp);
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
