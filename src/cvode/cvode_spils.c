/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * -----------------------------------------------------------------
 * This is the implementation file for the CVSPILS linear solver
 * interface.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_spils_impl.h"
#include <sundials/sundials_math.h>

/* Private constants */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)

/* Algorithmic constants */

#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */


/*===============================================================
  CVSPILS Exported functions -- Required
===============================================================*/

/*---------------------------------------------------------------
 CVSpilsSetLinearSolver specifies the iterative linear solver.
---------------------------------------------------------------*/
int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS)
{
  int retval;
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if any input is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
		    "CVSpilsSetLinearSolver", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  if (LS == NULL) {
    cvProcessError(NULL, CVSPILS_ILL_INPUT, "CVSPILS", 
		    "CVSpilsSetLinearSolver", 
                    "LS must be non-NULL");
    return(CVSPILS_ILL_INPUT);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if solver and vector are compatible with SPILS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_ITERATIVE) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", 
                    "CVSpilsSetLinearSolver", 
                    "Non-iterative LS supplied to CVSpils interface");
    return(CVSPILS_ILL_INPUT);
  }
  if ( (cv_mem->cv_tempv->ops->nvlinearsum == NULL) ||
       (cv_mem->cv_tempv->ops->nvconst == NULL) ||
       (cv_mem->cv_tempv->ops->nvdotprod == NULL) ){
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS", 
                    "CVSpilsSetLinearSolver", MSGS_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  /* free any existing system solver attached to CVode */
  if (cv_mem->cv_lfree)  cv_mem->cv_lfree(cv_mem);

  /* Set four main system linear solver function fields in cv_mem */
  cv_mem->cv_linit  = cvSpilsInitialize;
  cv_mem->cv_lsetup = cvSpilsSetup;
  cv_mem->cv_lsolve = cvSpilsSolve;
  cv_mem->cv_lfree  = cvSpilsFree;
  
  /* Get memory for CVSpilsMemRec */
  cvspils_mem = NULL;
  cvspils_mem = (CVSpilsMem) malloc(sizeof(struct CVSpilsMemRec));
  if (cvspils_mem == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPILS", 
                    "CVSpilsSetLinearSolver", MSGS_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  cvspils_mem->LS = LS;
  
  /* Set defaults for Jacobian-related fields */
  cvspils_mem->jtimesDQ = SUNTRUE;
  cvspils_mem->jtsetup = NULL;
  cvspils_mem->jtimes = CVSpilsDQJtimes;
  cvspils_mem->j_data = cv_mem;

  /* Set defaults for preconditioner-related fields */
  cvspils_mem->pset   = NULL;
  cvspils_mem->psolve = NULL;
  cvspils_mem->pfree  = NULL;
  cvspils_mem->P_data = cv_mem->cv_user_data;

  /* Initialize counters */
  cvSpilsInitializeCounters(cvspils_mem);

  /* Set default values for the rest of the SPILS parameters */
  cvspils_mem->jbad = SUNTRUE;
  cvspils_mem->eplifac = CVSPILS_EPLIN;
  cvspils_mem->last_flag = CVSPILS_SUCCESS;

  /* Attach default CVSpils interface routines to iterative LS */
  retval = SUNLinSolSetATimes(LS, cv_mem, CVSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVSPILS_SUNLS_FAIL, "CVSPILS", 
                    "CVSpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetATimes");
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_SUNLS_FAIL);
  }
  retval = SUNLinSolSetPreconditioner(LS, cv_mem, NULL, NULL);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVSPILS_SUNLS_FAIL, "CVSPILS", 
                    "CVSpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetPreconditioner");
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_SUNLS_FAIL);
  }

  /* Allocate memory for ytemp and x */
  cvspils_mem->ytemp = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->ytemp == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPILS", 
                    "CVSpilsSetLinearSolver", MSGS_MEM_FAIL);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }

  cvspils_mem->x = N_VClone(cv_mem->cv_tempv);
  if (cvspils_mem->x == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVSPILS", 
                    "CVSpilsSetLinearSolver", MSGS_MEM_FAIL);
    N_VDestroy(cvspils_mem->ytemp);
    free(cvspils_mem); cvspils_mem = NULL;
    return(CVSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, cvspils_mem->ytemp);
  cvspils_mem->sqrtN = SUNRsqrt( N_VDotProd(cvspils_mem->ytemp, 
                                            cvspils_mem->ytemp) );

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvspils_mem;

  return(CVSPILS_SUCCESS);
}


/*===============================================================
  CVSPILS Exported functions -- Optional input/output
===============================================================*/


int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsSetEpsLin", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsSetEpsLin", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVSPILS",
                   "CVSpilsSetEpsLin", MSGS_BAD_EPLIN);
    return(CVSPILS_ILL_INPUT);
  }

  cvspils_mem->eplifac = (eplifac == ZERO) ? CVSPILS_EPLIN : eplifac;

  return(CVSPILS_SUCCESS);
}


int CVSpilsSetPreconditioner(void *cvode_mem, 
                             CVSpilsPrecSetupFn psetup,
                             CVSpilsPrecSolveFn psolve)
{
  int retval;
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  PSetupFn cvspils_psetup;
  PSolveFn cvspils_psolve;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsSetPreconditioner", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* store function pointers for user-supplied routines in CVSpils interface */
  cvspils_mem->pset   = psetup;
  cvspils_mem->psolve = psolve;

  /* notify iterative linear solver to call CVSpils interface routines */
  cvspils_psetup = (psetup == NULL) ? NULL : CVSpilsPSetup;
  cvspils_psolve = (psolve == NULL) ? NULL : CVSpilsPSolve;
  retval = SUNLinSolSetPreconditioner(cvspils_mem->LS, cv_mem, 
                                      cvspils_psetup, cvspils_psolve);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVSPILS_SUNLS_FAIL, "CVSPILS", 
                   "CVSpilsSetPreconditioner", 
                   "Error in calling SUNLinSolSetPreconditioner");
    return(CVSPILS_SUNLS_FAIL);
  }

  return(CVSPILS_SUCCESS);
}


int CVSpilsSetJacTimes(void *cvode_mem,
                       CVSpilsJacTimesSetupFn jtsetup,
                       CVSpilsJacTimesVecFn jtimes)
{
  int retval;
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsSetJacTimes", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsSetJacTimes", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* store function pointers for user-supplied routines in CVSpils 
     interface (NULL jtimes implies use of DQ default) */
  if (jtimes != NULL) {
    cvspils_mem->jtimesDQ = SUNFALSE;
    cvspils_mem->jtimes   = jtimes;
  } else {
    cvspils_mem->jtimesDQ = SUNTRUE;
  }
  cvspils_mem->jtsetup = jtsetup;

  /* notify iterative linear solver to call CVSpils interface routines */
  retval = SUNLinSolSetATimes(cvspils_mem->LS, cv_mem, CVSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVSPILS_SUNLS_FAIL, "CVSPILS", 
                    "CVSpilsSetJacTimes", 
                    "Error in calling SUNLinSolSetATimes");
    return(CVSPILS_SUNLS_FAIL);
  }

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS,
                        long int *leniwLS)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetWorkSpace", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* start with fixed sizes plus NVectors */
  *lenrwLS = 4;
  *leniwLS = 10;

  /* add NVector sizes */
  if (cv_mem->cv_tempv->ops->nvspace) {
    N_VSpace(cv_mem->cv_tempv, &lrw1, &liw1);
    *lenrwLS += 2*lrw1;
    *leniwLS += 2*liw1;
  }

  /* add LS sizes */
  if (cvspils_mem->LS->ops->space) {
    flag = SUNLinSolSpace(cvspils_mem->LS, &lrw, &liw);
    *lenrwLS += lrw;
    *leniwLS += liw;
  }

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumPrecEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *npevals = cvspils_mem->npe;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumPrecSolves", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *npsolves = cvspils_mem->nps;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumLinIters", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *nliters = cvspils_mem->nli;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumConvFails", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *nlcfails = cvspils_mem->ncfl;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumJTSetupEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumJTSetupEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *njtsetups = cvspils_mem->njtsetup;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumJtimesEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *njvevals = cvspils_mem->njtimes;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetNumRhsEvals", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *nfevalsLS = cvspils_mem->nfes;

  return(CVSPILS_SUCCESS);
}


int CVSpilsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS",
                   "CVSpilsGetLastFlag", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS",
                   "CVSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  *flag = cvspils_mem->last_flag;

  return(CVSPILS_SUCCESS);
}


char *CVSpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVSPILS_SUCCESS:
    sprintf(name,"CVSPILS_SUCCESS");
    break; 
  case CVSPILS_MEM_NULL:
    sprintf(name,"CVSPILS_MEM_NULL");
    break;
  case CVSPILS_LMEM_NULL:
    sprintf(name,"CVSPILS_LMEM_NULL");
    break;
  case CVSPILS_ILL_INPUT:
    sprintf(name,"CVSPILS_ILL_INPUT");
    break;
  case CVSPILS_MEM_FAIL:
    sprintf(name,"CVSPILS_MEM_FAIL");
    break;
  case CVSPILS_PMEM_NULL:
    sprintf(name,"CVSPILS_PMEM_NULL");
    break;
  case CVSPILS_SUNLS_FAIL:
    sprintf(name,"CVSPILS_SUNLS_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*=================================================================
  CVSPILS private functions
  =================================================================*/

/*-----------------------------------------------------------------
  CVSpilsATimes

  This routine generates the matrix-vector product z = Mv, where
  M = I - gamma*J. The product J*v is obtained by calling the jtimes 
  routine. It is then scaled by -gamma and added to v to obtain M*v.
  The return value is the same as the value returned by jtimes --
  0 if successful, nonzero otherwise.
  -----------------------------------------------------------------*/
int CVSpilsATimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem   cv_mem;
  CVSpilsMem cvspils_mem;
  int jtflag;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "CVSpilsATimes", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "CVSpilsATimes", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  jtflag = cvspils_mem->jtimes(v, z, cv_mem->cv_tn,
                               cvspils_mem->ycur,
                               cvspils_mem->fcur,
                               cvspils_mem->j_data,
                               cvspils_mem->ytemp);
  cvspils_mem->njtimes++;
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -cv_mem->cv_gamma, z, z);

  return(0);
}


/*---------------------------------------------------------------
 CVSpilsPSetup:

 This routine interfaces between the generic iterative linear 
 solvers and the user's psetup routine.  It passes to psetup all 
 required state information from cvode_mem.  Its return value 
 is the same as that returned by psetup. Note that the generic
 iterative linear solvers guarantee that CVSpilsPSetup will only
 be called in the case that the user's psetup routine is non-NULL.
---------------------------------------------------------------*/
int CVSpilsPSetup(void *cvode_mem)
{
  int        retval;
  CVodeMem   cv_mem;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "CVSpilsPSetup", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "CVSpilsPSetup", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Call user pset routine to update preconditioner and possibly 
     reset jcur (pass !jbad as update suggestion) */
  retval = cvspils_mem->pset(cv_mem->cv_tn, 
                             cvspils_mem->ycur, 
                             cvspils_mem->fcur, 
                             !(cvspils_mem->jbad),
                             &cv_mem->cv_jcur,
                             cv_mem->cv_gamma, 
                             cvspils_mem->P_data);
  return(retval);     
}


/*-----------------------------------------------------------------
  CVSpilsPSolve

  This routine interfaces between the generic SUNLinSolSolve 
  routine and the user's psolve routine.  It passes to psolve all
  required state information from cvode_mem.  Its return value is 
  the same as that returned by psolve. Note that the generic 
  SUNLinSol solver guarantees that CVSpilsPSolve will not be called 
  in the case in which preconditioning is not done. This is the 
  only case in which the user's psolve routine is allowed to be 
  NULL.
  -----------------------------------------------------------------*/
int CVSpilsPSolve(void *cvode_mem, N_Vector r, N_Vector z,
                  realtype tol, int lr)
{
  CVodeMem   cv_mem;
  CVSpilsMem cvspils_mem;
  int retval;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "CVSpilsPSolve", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "CVSpilsPSolve", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* call the user-supplied psolve routine, and accumulate count */
  retval = cvspils_mem->psolve(cv_mem->cv_tn, cvspils_mem->ycur,
                               cvspils_mem->fcur, r, z,
                               cv_mem->cv_gamma, tol, lr,
                               cvspils_mem->P_data);
  cvspils_mem->nps++;
  return(retval);     
}


/*-----------------------------------------------------------------
  CVSpilsDQJtimes

  This routine generates a difference quotient approximation to
  the Jacobian times vector f_y(t,y) * v. The approximation is 
  Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
  input, i.e. the WRMS norm of v/vnrm is 1.
  -----------------------------------------------------------------*/
int CVSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
                    N_Vector y, N_Vector fy,
                    void *cvode_mem, N_Vector work)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* Return immediately if cvode_mem or cv_mem->cv_lmem are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "CVSpilsDQJtimes", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "CVSpilsDQJtimes", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, cv_mem->cv_ewt);

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = cv_mem->cv_f(t, work, Jv, cv_mem->cv_user_data); 
    cvspils_mem->nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}


/*-----------------------------------------------------------------
  cvSpilsInitialize

  This routine performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)
  -----------------------------------------------------------------*/
int cvSpilsInitialize(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;

  /* Return immediately if cv_mem or cv_mem->cv_lmem are NULL */
  if (cv_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "cvSpilsInitialize", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "cvSpilsInitialize", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  
  cvSpilsInitializeCounters(cvspils_mem);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (cvspils_mem->jtimesDQ) {
    cvspils_mem->jtsetup = NULL;
    cvspils_mem->jtimes = CVSpilsDQJtimes;
    cvspils_mem->j_data = cv_mem;
  } else {
    cvspils_mem->j_data = cv_mem->cv_user_data;
  }

  /* if psetup is not present, then cvSpilsSetup does not need to be 
     called, so set the lsetup function to NULL */
  if (cvspils_mem->pset == NULL)  cv_mem->cv_lsetup = NULL;

  /* Call LS initialize routine */
  cvspils_mem->last_flag = SUNLinSolInitialize(cvspils_mem->LS);
  return(cvspils_mem->last_flag);
}


/*-----------------------------------------------------------------
  cvSpilsSetup

  This routine calls the LS 'setup' routine.
  -----------------------------------------------------------------*/
int cvSpilsSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                 N_Vector fpred, booleantype *jcurPtr,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype dgamma;
  int  retval;
  CVSpilsMem cvspils_mem;

  /* Return immediately if cv_mem or cv_mem->cv_lmem are NULL */
  if (cv_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "cvSpilsSetup", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "cvSpilsSetup", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Set CVSpils N_Vector pointers to current solution and rhs */
  cvspils_mem->ycur = ypred;
  cvspils_mem->fcur = fpred;

  /* Use nst, gamma/gammap, and convfail to set J/P eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  cvspils_mem->jbad = (cv_mem->cv_nst == 0) || 
    (cv_mem->cv_nst > cvspils_mem->nstlpre + CVSPILS_MSBPRE) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVSPILS_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  *jcurPtr = cvspils_mem->jbad;
  
  /* Call LS setup routine -- the LS will call CVSpilsPSetup, who will 
     pass the heuristic suggestions above to the user code(s) */
  retval = SUNLinSolSetup(cvspils_mem->LS, NULL);

  /* If user set jcur to SUNTRUE, increment npe and save nst value */
  if (*jcurPtr) {
    cvspils_mem->npe++;
    cvspils_mem->nstlpre = cv_mem->cv_nst;
  }
  
  /* Update jcur flag if we suggested an update */
  if (cvspils_mem->jbad) *jcurPtr = SUNTRUE;

  return(retval);
}


/*-----------------------------------------------------------------
  cvSpilsSolve

  This routine interfaces between CVode and the generic 
  SUNLinearSolver object LS, by setting the appropriate tolerance 
  and scaling vectors, calling the solver, and accumulating 
  statistics from the solve for use/reporting by CVode.
  -----------------------------------------------------------------*/
int cvSpilsSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                 N_Vector ynow, N_Vector fnow)
{
  realtype bnorm;
  CVSpilsMem cvspils_mem;
  int nli_inc, retval;
  
  /* Return immediately if cv_mem or cv_mem->cv_lmem are NULL */
  if (cv_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", 
                   "cvSpilsSolve", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", 
                   "cvSpilsSolve", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  cvspils_mem->deltar = cvspils_mem->eplifac * cv_mem->cv_tq[4]; 
  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= cvspils_mem->deltar) {
    if (cv_mem->cv_mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve 
     interface routines */
  cvspils_mem->ycur = ynow;
  cvspils_mem->fcur = fnow;

  /* Set input tolerance and initial guess x = 0 to LS */  
  cvspils_mem->delta = cvspils_mem->deltar * cvspils_mem->sqrtN;
  N_VConst(ZERO, cvspils_mem->x);

  /* Set scaling vectors for LS to use */
  retval = SUNLinSolSetScalingVectors(cvspils_mem->LS,
                                      weight,
                                      weight);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVSPILS_SUNLS_FAIL, "CVSPILS", "cvSpilsSolve", 
                    "Error in calling SUNLinSolSetScalingVectors");
    return(CVSPILS_SUNLS_FAIL);
  }

  /* If a user-provided jtsetup routine is supplied, call that here */
  if (cvspils_mem->jtsetup) {
    retval = cvspils_mem->jtsetup(cv_mem->cv_tn, ynow, fnow, 
                                  cvspils_mem->j_data);
    cvspils_mem->njtsetup++;
    if (retval != 0) {
      cvProcessError(cv_mem, retval, "CVSPILS", 
                      "cvSpilsSolve", MSGS_JTSETUP_FAILED);
      return(retval);
    }
  }
  
  /* Call solver, and copy x to b */
  retval = SUNLinSolSolve(cvspils_mem->LS, NULL, cvspils_mem->x,
                          b, cvspils_mem->delta);
  N_VScale(ONE, cvspils_mem->x, b);

  /* Retrieve solver statistics */
  nli_inc = SUNLinSolNumIters(cvspils_mem->LS);
  
  /* Increment counters nli and ncfl */
  cvspils_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) cvspils_mem->ncfl++;

  /* Interpret solver return value  */
  cvspils_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
    /* allow reduction but not solution on first Newton iteration, 
       otherwise return with a recoverable failure */
    if (cv_mem->cv_mnewt == 0) return(0);
    else                       return(1);
    break;
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL:
    return(1);
    break;
  case SUNLS_MEM_NULL:
  case SUNLS_ILL_INPUT:
  case SUNLS_MEM_FAIL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL:
    return(-1);
    break;
  case SUNLS_PACKAGE_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_PACKAGE_FAIL_UNREC, "CVSPILS", 
                   "cvSpilsSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_ATIMES_FAIL_UNREC, "CVSPILS", 
                   "cvSpilsSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_PSOLVE_FAIL_UNREC, "CVSPILS", 
                   "cvSpilsSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }
  
  return(0); 
}


/*-----------------------------------------------------------------
  cvSpilsFree

  This routine frees memory associates with the CVSpils system
  solver interface.
  -----------------------------------------------------------------*/
int cvSpilsFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;

  /* Return immediately if cv_mem or cv_mem->cv_lmem are NULL */
  if (cv_mem == NULL)  return (CVSPILS_SUCCESS);
  if (cv_mem->cv_lmem == NULL)  return(CVSPILS_SUCCESS);
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Free N_Vector memory */
  if (cvspils_mem->ytemp) {
    N_VDestroy(cvspils_mem->ytemp);
    cvspils_mem->ytemp = NULL;
  }
  if (cvspils_mem->x) {
    N_VDestroy(cvspils_mem->x);
    cvspils_mem->x = NULL;
  }

  /* Nullify other N_Vector pointers */
  cvspils_mem->ycur = NULL;
  cvspils_mem->fcur = NULL;

  /* Free preconditioner memory (if applicable) */
  if (cvspils_mem->pfree)  cvspils_mem->pfree(cv_mem);
  
  /* free CVSpils interface structure */
  free(cv_mem->cv_lmem);
  
  return(CVSPILS_SUCCESS);
}


/*-----------------------------------------------------------------
  cvSpilsInitializeCounters

  This routine resets all counters from an CVSpilsMem structure.
  -----------------------------------------------------------------*/
int cvSpilsInitializeCounters(CVSpilsMem cvspils_mem)
{
  cvspils_mem->npe      = 0;
  cvspils_mem->nli      = 0;
  cvspils_mem->nps      = 0;
  cvspils_mem->ncfl     = 0;
  cvspils_mem->nstlpre  = 0;
  cvspils_mem->njtsetup = 0;
  cvspils_mem->njtimes  = 0;
  cvspils_mem->nfes     = 0;
  return(0);
}
