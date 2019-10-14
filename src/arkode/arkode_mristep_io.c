/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the optional input and output functions
 * for the ARKode MRIStep time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  MRIStep Optional input functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetDenseOrder: Specifies the polynomial order for dense
  output.  Positive values are sent to the interpolation module;
  negative values imply to use the default.
  ---------------------------------------------------------------*/
int MRIStepSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetDenseOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDenseOrder(ark_mem, dord));
}

/*---------------------------------------------------------------
  MRIStepSetErrHandlerFn: Specifies the error handler function
  ---------------------------------------------------------------*/
int MRIStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrHandlerFn(ark_mem, ehfun, eh_data));
}

/*---------------------------------------------------------------
  MRIStepSetErrFile: Specifies the FILE pointer for output (NULL
  means no messages)
  ---------------------------------------------------------------*/
int MRIStepSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrFile(ark_mem, errfp));
}

/*---------------------------------------------------------------
  MRIStepSetUserData: Specifies the user data pointer

  NOTE: This function assumes that the inner ARKodeMem user_data
  variable is only passed to the RHS function(s) and that separate
  variables (e.g., e_data, J_data, etc.) are used to pass the user
  data to all other user functions. This ensures that re-attaching
  the outer stepper to the inner user_data variable does not
  clobber the pointers passed to other user functions.
  ---------------------------------------------------------------*/
int MRIStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem         ark_mem;       /* outer ARKode memory  */
  ARKodeMRIStepMem  step_mem;      /* outer stepper memory */
  ARKodeMem         inner_ark_mem; /* inner ARKode memory  */
  int               retval;

  /* Access MRIStep memory */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetUserData",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check for inner memory */
  if (step_mem->inner_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetUserData", "The inner stepper memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* Attach user_data to the outer stepper */
  retval = arkSetUserData(ark_mem, user_data);
  if (retval != ARK_SUCCESS) return(retval);

  /* Attach user_data to the inner stepper and re-attach the outer
     stepper to the inner stepper ARKodeMem user_data pointer. */
  switch (step_mem->inner_stepper_id) {
  case MRISTEP_ARKSTEP:
    retval = ARKStepSetUserData(step_mem->inner_mem, user_data);
    if (retval != ARK_SUCCESS) return(retval);
    inner_ark_mem = (ARKodeMem)(step_mem->inner_mem);
    inner_ark_mem->user_data = arkode_mem;
    break;
  }

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetDiagnostics: Specifies to enable solver diagnostics,
  and specifies the FILE pointer for output (diagfp==NULL
  disables output)
  ---------------------------------------------------------------*/
int MRIStepSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDiagnostics(ark_mem, diagfp));
}

/*---------------------------------------------------------------
  MRIStepSetMaxNumSteps: Specifies the maximum number of
  integration steps
  ---------------------------------------------------------------*/
int MRIStepSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxNumSteps(ark_mem, mxsteps));
}

/*---------------------------------------------------------------
  MRIStepSetMaxHnilWarns: Specifies the maximum number of warnings
  for small h
  ---------------------------------------------------------------*/
int MRIStepSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxHnilWarns(ark_mem, mxhnil));
}

/*---------------------------------------------------------------
  MRIStepSetStopTime: Specifies the time beyond which the
  integration is not to proceed.
  ---------------------------------------------------------------*/
int MRIStepSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetStopTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetStopTime(ark_mem, tstop));
}

/*---------------------------------------------------------------
  MRIStepSetFixedStep: Specifies the fixed time step sizes to use
  with MRIStep. MRIStep will use this step size for all steps
  (unless tstop is set, in which case it may need to modify that
  last step approaching tstop. If any solver failure occurs in the
  timestepping module, MRIStep will typically immediately return
  with an error message indicating that the selected step size
  cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size.
  ---------------------------------------------------------------*/
int MRIStepSetFixedStep(void *arkode_mem, realtype hsfixed)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetFixedStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (hsfixed == ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetFixedStep",
                    "MIRStep does not support adaptive steps at this time.");
    return(ARK_ILL_INPUT);
  }

  /* using an explicit method with fixed step sizes, enforce use of
     arkEwtSmallReal to compute error weight vector */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSetSmallReal;
  ark_mem->e_data    = ark_mem;

  return(arkSetFixedStep(ark_mem, hsfixed));
}

/*---------------------------------------------------------------
  MRIStepSetRootDirection: Specifies the direction of zero-crossings
  to be monitored.  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int MRIStepSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetRootDirection(ark_mem, rootdir));
}

/*---------------------------------------------------------------
  MRIStepSetNoInactiveRootWarn:  Disables issuing a warning if
  some root function appears to be identically zero at the
  beginning of the integration
  ---------------------------------------------------------------*/
int MRIStepSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetNoInactiveRootWarn(ark_mem));
}

/*---------------------------------------------------------------
  MRIStepSetPostprocessStepFn:  Specifies a user-provided step
  postprocessing function having type ARKPostProcessStepFn.  A
  NULL input function disables step postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int MRIStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessStepFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetPostprocessStepFn(ark_mem, ProcessStep));
}


/*===============================================================
  MRIStep Optional output functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepGetNumSteps:  Returns the current number of integration
  steps
  ---------------------------------------------------------------*/
int MRIStepGetNumSteps(void *arkode_mem, long int *nssteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumSteps(ark_mem, nssteps));
}

/*---------------------------------------------------------------
  MRIStepGetLastStep: Returns the step size used on the last
  successful step
  ---------------------------------------------------------------*/
int MRIStepGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetLastStep(ark_mem, hlast));
}

/*---------------------------------------------------------------
  MRIStepGetCurrentTime: Returns the current value of the
  independent variable
  ---------------------------------------------------------------*/
int MRIStepGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentTime(ark_mem, tcur));
}

/*---------------------------------------------------------------
  MRIStepGetWorkSpace: Returns integrator work space requirements
  ---------------------------------------------------------------*/
int MRIStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetWorkSpace(ark_mem, lenrw, leniw));
}

/*---------------------------------------------------------------
  MRIStepGetNumGEvals: Returns the current number of calls to g
  (for rootfinding)
  ---------------------------------------------------------------*/
int MRIStepGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumGEvals(ark_mem, ngevals));
}

/*---------------------------------------------------------------
  MRIStepGetRootInfo: Returns pointer to array rootsfound showing
  roots found
  ---------------------------------------------------------------*/
int MRIStepGetRootInfo(void *arkode_mem, int *rootsfound)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetRootInfo(ark_mem, rootsfound));
}

/*---------------------------------------------------------------
  MRIStepGetReturnFlagName: translates from return flags IDs to
  names
  ---------------------------------------------------------------*/
char *MRIStepGetReturnFlagName(long int flag)
{ return(arkGetReturnFlagName(flag)); }



/*===============================================================
  MRIStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetDefaults:

  Resets all MRIStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int MRIStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set default values for integrator optional inputs */
  step_mem->q         = 3;              /* method order */
  step_mem->p         = 0;              /* embedding order */
  step_mem->stages    = 0;              /* no stages */
  step_mem->B         = NULL;           /* no Butcher table */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.
  ---------------------------------------------------------------*/
int MRIStepSetTable(void *arkode_mem, int q, ARKodeButcherTable B)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for illegal inputs */
  if (B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTables", MSG_ARK_NO_MEM);
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;

  /* set the relevant parameters */
  step_mem->stages = B->stages;
  step_mem->q = B->q;
  step_mem->p = 0; /* assume fixed stepping */

  /* copy the table in step memory */
  step_mem->B = ARKodeButcherTable_Copy(B);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetTableNum:

  Specifies to use pre-existing Butcher tables for the problem,
  based on the integer flags passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int MRIStepSetTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that argument specifies an explicit table (assume explicit) */
  if (itable < MIN_ERK_NUM || itable > MAX_ERK_NUM ) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTableNum",
                    "Illegal MRI table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B);  step_mem->B = NULL;

  /* fill in table based on argument */
  step_mem->B = ARKodeButcherTable_LoadERK(itable);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTableNum",
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  step_mem->stages = step_mem->B->stages;
  step_mem->q = step_mem->B->q;
  step_mem->p = step_mem->B->p;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPreInnerFn:

  Sets the user-supplied function called BEFORE the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPreInnerFn(void *arkode_mem, MRIStepPreInnerFn prefn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->pre_inner_evolve = prefn;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPostInnerFn:

  Sets the user-supplied function called AFTER the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPostInnerFn(void *arkode_mem, MRIStepPostInnerFn postfn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->post_inner_evolve = postfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  MRIStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepGetLastInnerStepFlag:

  Returns the last return value from the inner stepper.
  ---------------------------------------------------------------*/
int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetLastInnerStepFlag",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get the last return value from the inner stepper */
  *flag = step_mem->inner_retval;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNumRhsEvals:

  Returns the current number of calls to fs and ff
  ---------------------------------------------------------------*/
int MRIStepGetNumRhsEvals(void *arkode_mem, long int *nfs_evals)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get number of fs evals from step_mem */
  *nfs_evals = step_mem->nfs;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetCurrentButcherTables:

  Sets pointers to the slow and fast Butcher tables currently in
  use.
  ---------------------------------------------------------------*/
int MRIStepGetCurrentButcherTables(void *arkode_mem, ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetCurrentButcherTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get table from step_mem */
  *B = step_mem->B;

  return(ARK_SUCCESS);
}


/*===============================================================
  MRIStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKode infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteButcher(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteButcher", "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* wrie outer Butcher table */
  fprintf(fp, "\nMRIStep Butcher tables:\n");
  if (step_mem->B != NULL) {
    fprintf(fp, "  Slow Butcher table (stages = %i):\n", step_mem->stages);
    ARKodeButcherTable_Write(step_mem->B, fp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
