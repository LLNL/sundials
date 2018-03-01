/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * output functions for the ARKODE solver.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
 ARKODE optional input functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeSetDefaults:

 Resets all optional inputs to ARKode default values.  Does not 
 change problem-defining function pointers fe and fi or 
 user_data pointer.  Also leaves alone any data 
 structures/options related to root-finding (those can be reset 
 using ARKodeRootInit).
---------------------------------------------------------------*/
int ARKodeSetDefaults(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->dense_q          = QDENSE_DEF;     /* dense output order */
  ark_mem->fixedstep        = SUNFALSE;       /* default to use adaptive steps */
  ark_mem->reltol           = 1.e-4;          /* relative tolerance */
  ark_mem->itol             = ARK_SS;         /* scalar-scalar solution tolerances */
  ark_mem->ritol            = ARK_SS;         /* scalar-scalar residual tolerances */
  ark_mem->Sabstol          = 1.e-9;          /* solution absolute tolerance */
  ark_mem->SRabstol         = 1.e-9;          /* residual absolute tolerance */
  ark_mem->user_efun        = SUNFALSE;       /* no user-supplied ewt function */
  ark_mem->efun             = arkEwtSet;      /* built-in ewt function */
  ark_mem->e_data           = NULL;           /* ewt function data */
  ark_mem->user_rfun        = SUNFALSE;       /* no user-supplied rwt function */
  ark_mem->rfun             = arkRwtSet;      /* built-in rwt function */
  ark_mem->r_data           = NULL;           /* rwt function data */
  ark_mem->ehfun            = arkErrHandler;  /* default error handler fn */
  ark_mem->eh_data          = ark_mem;        /* error handler data */
  ark_mem->errfp            = stderr;         /* output stream for errors */
  ark_mem->mxstep           = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->mxhnil           = MXHNIL;         /* max warns of t+h==t */
  ark_mem->hin              = ZERO;           /* determine initial step on-the-fly */
  ark_mem->hmin             = ZERO;           /* no minimum step size */
  ark_mem->hmax_inv         = ZERO;           /* no maximum step size */
  ark_mem->tstopset         = SUNFALSE;       /* no stop time set */
  ark_mem->tstop            = ZERO;           /* no fixed stop time */
  ark_mem->diagfp           = NULL;           /* no solver diagnostics file */
  ark_mem->report           = SUNFALSE;       /* don't report solver diagnostics */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSodeSetDenseOrder:

 Specifies the polynomial order for dense output.  Positive 
 values are sent to the interpolation module; negative values 
 imply to use the default.
---------------------------------------------------------------*/
int ARKodeSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetDenseOrder", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided value, or default, depending on argument */
  if (dord < 0) {
    ark_mem->dense_q = QDENSE_DEF;
  } else {
    ark_mem->dense_q = dord;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrHandlerFn:

 Specifies the error handler function
---------------------------------------------------------------*/
int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, 
                          void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetErrHandlerFn", MSGARK_NO_MEM);
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
 ARKodeSetErrFile:

 Specifies the FILE pointer for output (NULL means no messages)
---------------------------------------------------------------*/
int ARKodeSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetErrFile", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->errfp = errfp;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetUserData:

 Specifies the user data pointer for f
---------------------------------------------------------------*/
int ARKodeSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetUserData", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->user_data = user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetDiagnostics:

 Specifies to enable solver diagnostics, and specifies the FILE
 pointer for output (diagfp==NULL disables output)
---------------------------------------------------------------*/
int ARKodeSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetDiagnostics", MSGARK_NO_MEM);
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
 ARKodeSetMaxNumSteps:

 Specifies the maximum number of integration steps
---------------------------------------------------------------*/
int ARKodeSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetMaxNumSteps", MSGARK_NO_MEM);
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
 ARKodeSetMaxHnilWarns:

 Specifies the maximum number of warnings for small h
---------------------------------------------------------------*/
int ARKodeSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetMaxHnilWarns", MSGARK_NO_MEM);
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
 ARKodeSetInitStep:

 Specifies the initial step size
---------------------------------------------------------------*/
int ARKodeSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) {
    ark_mem->hin = ZERO;
  } else {
    ark_mem->hin = hin;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMinStep:

 Specifies the minimum step size
---------------------------------------------------------------*/
int ARKodeSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetMinStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmin <= ZERO) {
    ark_mem->hmin = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * ark_mem->hmax_inv > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKodeSetMinStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmin = hmin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxStep:

 Specifies the maximum step size
---------------------------------------------------------------*/
int ARKodeSetMaxStep(void *arkode_mem, realtype hmax)
{
  realtype hmax_inv;
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetMaxStep", MSGARK_NO_MEM);
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
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKodeSetMaxStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmax_inv = hmax_inv;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetStopTime:

 Specifies the time beyond which the integration is not to proceed.
---------------------------------------------------------------*/
int ARKodeSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetStopTime", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If ARKode was called at least once, test if tstop is legal
     (i.e. if it was not already passed).
     If ARKodeSetStopTime is called before the first call to ARKode,
     tstop will be checked in ARKode. */
  if (ark_mem->nst > 0) {
    if ( (tstop - ark_mem->tcur) * ark_mem->h < ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                      "ARKodeSetStopTime", MSGARK_BAD_TSTOP, 
                      tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  ark_mem->tstop    = tstop;
  ark_mem->tstopset = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetFixedStep:

 Specifies to use a fixed time step size instead of performing 
 any form of temporal adaptivity.  ARKode will use this step size 
 for all steps (unless tstop is set, in which case it may need to 
 modify that last step approaching tstop.  If any (non)linear
 solver failure occurs, ARKode will immediately return with an 
 error message since the time step size cannot be modified.  

 Any nonzero argument will result in the use of that fixed step 
 size; an argument of 0 will re-enable temporal adaptivity.
---------------------------------------------------------------*/
int ARKodeSetFixedStep(void *arkode_mem, realtype hfixed)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetFixedStep", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set ark_mem entry */
  if (hfixed != ZERO) {
    ark_mem->fixedstep = SUNTRUE;
    ark_mem->hin = hfixed;
  } else {
    ark_mem->fixedstep = SUNFALSE;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetRootDirection:

 Specifies the direction of zero-crossings to be monitored.
 The default is to monitor both crossings.
---------------------------------------------------------------*/
int ARKodeSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  int i;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetRootDirection", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetRootDirection", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  if (ark_root_mem->nrtfn == 0) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE::ARKStep", 
                    "ARKodeSetRootDirection", MSGARK_NO_ROOT);
    return(ARK_ILL_INPUT);    
  }

  for(i=0; i<ark_root_mem->nrtfn; i++) 
    ark_root_mem->rootdir[i] = rootdir[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNoInactiveRootWarn:

 Disables issuing a warning if some root function appears
 to be identically zero at the beginning of the integration
---------------------------------------------------------------*/
int ARKodeSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetNoInactiveRootWarn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetNoInactiveRootWarn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  ark_root_mem->mxgnull = 0;
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetPostprocessStepFn:

 Specifies a user-provided step postprocessing function having 
 type ARKPostProcessStepFn.  A NULL input function disables step
 postprocessing.

 IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA, 
 THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND 
 STABILITY ARE LOST.
---------------------------------------------------------------*/
int ARKodeSetPostprocessStepFn(void *arkode_mem, 
                               ARKPostProcessStepFn ProcessStep)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeSetPostprocessStepFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStep = ProcessStep;
  return(ARK_SUCCESS);
}


/*===============================================================
 ARKODE optional output functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeGetNumSteps:

 Returns the current number of integration steps
---------------------------------------------------------------*/
int ARKodeGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->nst;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetActualInitStep:

 Returns the step size used on the first step
---------------------------------------------------------------*/
int ARKodeGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetActualInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *hinused = ark_mem->h0u;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetLastStep:

 Returns the step size used on the last successful step
---------------------------------------------------------------*/
int ARKodeGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetLastStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hlast = ark_mem->hold;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentStep:

 Returns the step size to be attempted on the next step
---------------------------------------------------------------*/
int ARKodeGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetCurrentStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  
  *hcur = ark_mem->next_h;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentTime:

 Returns the current value of the independent variable
---------------------------------------------------------------*/
int ARKodeGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetCurrentTime", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tcur = ark_mem->tcur;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetTolScaleFactor:

 Returns a suggested factor for scaling tolerances
---------------------------------------------------------------*/
int ARKodeGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetTolScaleFactor", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tolsfact = ark_mem->tolsf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetErrWeights:

 This routine returns the current error weight vector.
---------------------------------------------------------------*/
int ARKodeGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetErrWeights", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ewt, eweight);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetResWeights:

 This routine returns the current residual weight vector.
---------------------------------------------------------------*/
int ARKodeGetResWeights(void *arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetResWeights", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->rwt, rweight);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetWorkSpace:

 Returns integrator work space requirements
---------------------------------------------------------------*/
int ARKodeGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetWorkSpace", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *leniw = ark_mem->liw;
  *lenrw = ark_mem->lrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumGEvals:

 Returns the current number of calls to g (for rootfinding)
---------------------------------------------------------------*/
int ARKodeGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetNumGEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetNumGEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  *ngevals = ark_root_mem->nge;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetRootInfo:

 Returns pointer to array rootsfound showing roots found
---------------------------------------------------------------*/
int ARKodeGetRootInfo(void *arkode_mem, int *rootsfound)
{
  int i;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetRootInfo", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetRootInfo", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  for (i=0; i<ark_root_mem->nrtfn; i++) 
    rootsfound[i] = ark_root_mem->iroots[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetStepStats:

 Returns step statistics
---------------------------------------------------------------*/
int ARKodeGetStepStats(void *arkode_mem, long int *nsteps,
                       realtype *hinused, realtype *hlast,
                       realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeGetStepStats", MSGARK_NO_MEM);
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


/*-----------------------------------------------------------------*/

char *ARKodeGetReturnFlagName(long int flag)
{
  char *name;
  name = (char *)malloc(24*sizeof(char));

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
  default:
    sprintf(name,"NONE");
  }

  return(name);
}



/*===============================================================
  ARKODE parameter output
===============================================================*/

/*---------------------------------------------------------------
 ARKodeWriteParameters:

 Outputs all solver parameters to the provided file pointer.
---------------------------------------------------------------*/
int ARKodeWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::ARKStep", 
                    "ARKodeWriteParameters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ARKode solver parameters:\n");
  fprintf(fp, "  Dense output order %i\n",ark_mem->dense_q);
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

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
      EOF
---------------------------------------------------------------*/
