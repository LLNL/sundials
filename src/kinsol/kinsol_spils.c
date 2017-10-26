/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the KINSPILS linear solvers.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include "kinsol_spils_impl.h"

#include <sundials/sundials_math.h>

/*------------------------------------------------------------------
  private constants
  ------------------------------------------------------------------*/

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)


/*==================================================================
  KINSPILS Exported functions -- Required
  ==================================================================*/

/*---------------------------------------------------------------
  KINSpilsSetLinearSolver specifies the iterative linear solver
  ---------------------------------------------------------------*/
int KINSpilsSetLinearSolver(void *kinmem, SUNLinearSolver LS)
{
  int retval;
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if any input is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
		    "KINSpilsSetLinearSolver", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  if (LS == NULL) {
    KINProcessError(NULL, KINSPILS_ILL_INPUT, "KINSPILS",
		    "KINSpilsSetLinearSolver",
                    "LS must be non-NULL");
    return(KINSPILS_ILL_INPUT);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if solver and vector are compatible with SPILS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_ITERATIVE) {
    KINProcessError(kin_mem, KINSPILS_ILL_INPUT, "KINSPILS",
                    "KINSpilsSetLinearSolver",
                    "Non-iterative LS supplied to KINSpils interface");
    return(KINSPILS_ILL_INPUT);
  }

  /*
    check for required vector operations

    Note: do NOT need to check for N_VLinearSum, N_VProd, N_VScale,
    N_VDiv, or N_VWL2Norm because they are required by KINSOL.
  */
  if ( (kin_mem->kin_vtemp1->ops->nvconst == NULL) ||
       (kin_mem->kin_vtemp1->ops->nvdotprod == NULL) ||
       (kin_mem->kin_vtemp1->ops->nvl1norm == NULL) ) {
    KINProcessError(kin_mem, KINSPILS_ILL_INPUT, "KINSPILS",
                    "KINSpilsSetLinearSolver", MSGS_BAD_NVECTOR);
    return(KINSPILS_ILL_INPUT);
  }

  /* free any existing system solver attached to KIN */
  if (kin_mem->kin_lfree) kin_mem->kin_lfree(kin_mem);

  /* This is an iterative linear solver */
  kin_mem->kin_inexact_ls = SUNTRUE;

  /* Set four main system linear solver function fields in kin_mem */
  kin_mem->kin_linit  = kinSpilsInitialize;
  kin_mem->kin_lsetup = kinSpilsSetup;
  kin_mem->kin_lsolve = kinSpilsSolve;
  kin_mem->kin_lfree  = kinSpilsFree;
  
  /* Get memory for KINSpilsMemRec */
  kinspils_mem = NULL;
  kinspils_mem = (KINSpilsMem) malloc(sizeof(struct KINSpilsMemRec));
  if (kinspils_mem == NULL) {
    KINProcessError(kin_mem, KINSPILS_MEM_FAIL, "KINSPILS",
                    "KINSpilsSetLinearSolver", MSGS_MEM_FAIL);
    return(KINSPILS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  kinspils_mem->LS = LS;

  /* Set defaults for Jacobian-related fields */
  kinspils_mem->jtimesDQ = SUNTRUE;
  kinspils_mem->jtimes   = KINSpilsDQJtimes;
  kinspils_mem->jdata    = kin_mem;

  /* Set defaults for preconditioner-related fields */
  kinspils_mem->pset   = NULL;
  kinspils_mem->psolve = NULL;
  kinspils_mem->pfree  = NULL;
  kinspils_mem->pdata  = kin_mem->kin_user_data;

  /* Initialize counters */
  kinSpilsInitializeCounters(kinspils_mem);

  /* Set default values for the rest of the SPILS parameters */
  kinspils_mem->last_flag = KINSPILS_SUCCESS;

  /* Attach default KINSpils interface routines to iterative LS */
  retval = SUNLinSolSetATimes(LS, kin_mem, KINSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    KINProcessError(kin_mem, KINSPILS_SUNLS_FAIL, "KINSPILS",
                    "KINSpilsSetLinearSolver",
                    "Error in calling SUNLinSolSetATimes");
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_SUNLS_FAIL);
  }
  retval = SUNLinSolSetPreconditioner(LS, kin_mem, NULL, NULL);
  if (retval != SUNLS_SUCCESS) {
    KINProcessError(kin_mem, KINSPILS_SUNLS_FAIL, "KINSPILS",
                    "KINSpilsSetLinearSolver",
                    "Error in calling SUNLinSolSetPreconditioner");
    free(kinspils_mem); kinspils_mem = NULL;
    return(KINSPILS_SUNLS_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  kin_mem->kin_lmem = kinspils_mem;

  return(KINSPILS_SUCCESS);
}


/*==================================================================
  KINSPILS Exported functions -- Optional input/output
  ==================================================================*/

/*------------------------------------------------------------------
  KINSpilsSetPreconditioner sets the preconditioner setup and solve
  functions
  ------------------------------------------------------------------*/
int KINSpilsSetPreconditioner(void *kinmem,
			      KINSpilsPrecSetupFn psetup,
                              KINSpilsPrecSolveFn psolve)
{
  int retval;
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  PSetupFn kinspils_psetup;
  PSolveFn kinspils_psolve;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsSetPreconditioner", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  kinspils_mem->pset   = psetup;
  kinspils_mem->psolve = psolve;

  /* notify iterative linear solver to call KINSpils interface routines */
  kinspils_psetup = (psetup == NULL) ? NULL : KINSpilsPSetup;
  kinspils_psolve = (psolve == NULL) ? NULL : KINSpilsPSolve;
  retval = SUNLinSolSetPreconditioner(kinspils_mem->LS, kin_mem,
                                      kinspils_psetup, kinspils_psolve);
  if (retval != SUNLS_SUCCESS) {
    KINProcessError(kin_mem, KINSPILS_SUNLS_FAIL, "KINSPILS", 
                    "KINSpilsSetPreconditioner", 
                    "Error in calling SUNLinSolSetPreconditioner");
    return(KINSPILS_SUNLS_FAIL);
  }

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsSetJacTimesVecFn sets the matrix-vector product function
  ------------------------------------------------------------------*/
int KINSpilsSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jtv)
{
  int retval;
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsSetJacTimesVecFn", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* store function pointers for user-supplied routine in KINSpils 
     interface (NULL jtimes implies use of DQ default) */
  if (jtv != NULL) {
    kinspils_mem->jtimesDQ = SUNFALSE;
    kinspils_mem->jtimes   = jtv;
  } else {
    kinspils_mem->jtimesDQ = SUNTRUE;
  }

  /* notify iterative linear solver to call KINSpils interface routines */
  retval = SUNLinSolSetATimes(kinspils_mem->LS, kin_mem, KINSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    KINProcessError(kin_mem, KINSPILS_SUNLS_FAIL, "KINSPILS", 
                    "KINSpilsSetJacTimes", 
                    "Error in calling SUNLinSolSetATimes");
    return(KINSPILS_SUNLS_FAIL);
  }

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetWorkSpace returns the integer and real workspace size
  ------------------------------------------------------------------*/
int KINSpilsGetWorkSpace(void *kinmem, long int *lenrwLS,
                         long int *leniwLS)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetWorkSpace", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL,
                    "KINSPILS", "KINSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* start with fixed sizes */
  *lenrwLS = 0;
  *leniwLS = 9;

  /* add N_Vector sizes */
  if (kin_mem->kin_vtemp1->ops->nvspace) {
    N_VSpace(kin_mem->kin_vtemp1, &lrw1, &liw1);
    *lenrwLS += lrw1;
    *leniwLS += liw1;
  }

  /* add LS sizes */
  if (kinspils_mem->LS->ops->space) {
    flag = SUNLinSolSpace(kinspils_mem->LS, &lrw, &liw);
    *lenrwLS += lrw;
    *leniwLS += liw;
  }

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumPrecEvals returns the total number of preconditioner
  evaluations
  ------------------------------------------------------------------*/
int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumPrecEvals", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *npevals = kinspils_mem->npe;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumPrecSolves returns the total number of times the
  preconditioner was applied
  ------------------------------------------------------------------*/
int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumPrecSolves", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *npsolves = kinspils_mem->nps;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumLinIters returns the total number of linear 
  iterations
  ------------------------------------------------------------------*/
int KINSpilsGetNumLinIters(void *kinmem, long int *nliters)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumLinIters", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *nliters = kinspils_mem->nli;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumConvFails returns the total numbe of convergence 
  failures
  ------------------------------------------------------------------*/
int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumConvFails", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *nlcfails = kinspils_mem->ncfl;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumJtimesEvals returns the number of times the matrix
  vector product was computed
  ------------------------------------------------------------------*/
int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njvevals)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumJtimesEvals", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *njvevals = kinspils_mem->njtimes;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetNumFuncEvals returns the number of calls to the user's
  F routine by the linear solver module
  ------------------------------------------------------------------*/
int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevals)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetNumFuncEvals", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetNumFuncEvals", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *nfevals = kinspils_mem->nfes;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetLastFlag returns the last flag set in the KINSPILS 
  function
  ------------------------------------------------------------------*/
int KINSpilsGetLastFlag(void *kinmem, long int *flag)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kinmem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
                    "KINSpilsGetLastFlag", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
                    "KINSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  *flag = kinspils_mem->last_flag;

  return(KINSPILS_SUCCESS);
}

/*------------------------------------------------------------------
  KINSpilsGetReturnFlagName
  ------------------------------------------------------------------*/
char *KINSpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINSPILS_SUCCESS:
    sprintf(name, "KINSPILS_SUCCESS");
    break;
  case KINSPILS_MEM_NULL:
    sprintf(name, "KINSPILS_MEM_NULL");
    break;
  case KINSPILS_LMEM_NULL:
    sprintf(name, "KINSPILS_LMEM_NULL");
    break;
  case KINSPILS_ILL_INPUT:
    sprintf(name, "KINSPILS_ILL_INPUT");
    break;
  case KINSPILS_MEM_FAIL:
    sprintf(name, "KINSPILS_MEM_FAIL");
    break;
  case KINSPILS_PMEM_NULL:
    sprintf(name, "KINSPILS_PMEM_NULL");
    break;
  case KINSPILS_SUNLS_FAIL:
    sprintf(name,"KINSPILS_SUNLS_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}


/*==================================================================
  KINSPILS Private functions
  ==================================================================*/

/*------------------------------------------------------------------
  KINSpilsATimes
  
  This routine coordinates the generation of the matrix-vector
  product z = J*v by calling either KINSpilsDQJtimes, which uses
  a difference quotient approximation for J*v, or by calling the
  user-supplied routine KINSpilsJacTimesVecFn if it is non-null.
  ------------------------------------------------------------------*/
int KINSpilsATimes(void *kinmem, N_Vector v, N_Vector z)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", 
                    "KINSpilsATimes", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS", 
                    "KINSpilsATimes", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  retval = kinspils_mem->jtimes(v, z, kin_mem->kin_uu,
                                &(kinspils_mem->new_uu),
                                kinspils_mem->jdata);
  kinspils_mem->njtimes++;

  return(retval);
}


/*---------------------------------------------------------------
  KINSpilsPSetup:

  This routine interfaces between the generic iterative linear 
  solvers and the user's psetup routine. It passes to psetup all 
  required state information from kin_mem. Its return value 
  is the same as that returned by psetup. Note that the generic
  iterative linear solvers guarantee that KINSpilsPSetup will only
  be called in the case that the user's psetup routine is non-NULL.
  ---------------------------------------------------------------*/
int KINSpilsPSetup(void *kinmem)
{
  KINMem      kin_mem;
  KINSpilsMem kinspils_mem;
  int         retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", 
		    "KINSpilsPSetup", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS", 
		    "KINSpilsPSetup", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* Call user pset routine to update preconditioner */
  retval = kinspils_mem->pset(kin_mem->kin_uu, kin_mem->kin_uscale,
                              kin_mem->kin_fval, kin_mem->kin_fscale, 
                              kinspils_mem->pdata);
  kinspils_mem->npe++;
  return(retval);
}


/*------------------------------------------------------------------
  KINSpilsPSolve

  This routine interfaces between the generic iterative linear
  solvers and the user's psolve routine. It passes to psolve all
  required state information from kinsol_mem. Its return value is
  the same as that returned by psolve. Note that the generic
  SUNLinSol solver guarantees that KINSpilsPSolve will not be called
  in the case in which preconditioning is not done. This is the only
  case in which the user's psolve routine is allowed to be NULL.
  ------------------------------------------------------------------*/
int KINSpilsPSolve(void *kinmem, N_Vector r, N_Vector z,
                   realtype toldummy, int lrdummy)
{
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", 
		    "KINSpilsPSolve", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS", 
		    "KINSpilsPSolve", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;
  
  /* copy the rhs into z before the psolve call */   
  /* Note: z returns with the solution */
  N_VScale(ONE, r, z);

  retval = kinspils_mem->psolve(kin_mem->kin_uu, kin_mem->kin_uscale,
                                kin_mem->kin_fval, kin_mem->kin_fscale,
                                z, kinspils_mem->pdata);

  kinspils_mem->nps++;
  return(retval);
}


/*------------------------------------------------------------------
  KINSpilsDQJtimes

  This routine generates the matrix-vector product z = J*v using a
  difference quotient approximation. The approximation is 
  J*v = [func(uu + sigma*v) - func(uu)]/sigma. Here sigma is based
  on the dot products (uscale*uu, uscale*v) and
  (uscale*v, uscale*v), the L1Norm(uscale*v), and on sqrt_relfunc
  (the square root of the relative error in the function). Note
  that v in the argument list has already been both preconditioned
  and unscaled.
 
  NOTE: Unlike the DQ Jacobian functions for direct linear solvers
        (which are called from within the lsetup function), this
        function is called from within the lsolve function and thus
        a recovery may still be possible even if the system function
        fails (recoverably).
  ------------------------------------------------------------------*/
int KINSpilsDQJtimes(N_Vector v, N_Vector Jv, N_Vector u, 
                     booleantype *new_u, void *kinmem)
{
  realtype sigma, sigma_inv, sutsv, sq1norm, sign, vtv;
  KINMem kin_mem;
  KINSpilsMem kinspils_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS", 
		    "KINSpilsDQJtimes", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS", 
		    "KINSpilsDQJtimes", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* scale the vector v and put Du*v into vtemp1 */
  N_VProd(v, kin_mem->kin_uscale, kin_mem->kin_vtemp1);

  /* scale u and put into Jv (used as a temporary storage) */
  N_VProd(u, kin_mem->kin_uscale, Jv);

  /* compute dot product (Du*u).(Du*v) */
  sutsv = N_VDotProd(Jv, kin_mem->kin_vtemp1);

  /* compute dot product (Du*v).(Du*v) */
  vtv = N_VDotProd(kin_mem->kin_vtemp1, kin_mem->kin_vtemp1);

  sq1norm = N_VL1Norm(kin_mem->kin_vtemp1);

  sign = (sutsv >= ZERO) ? ONE : -ONE ;
 
  /* this expression for sigma is from p. 469, Brown and Saad paper */
  sigma = sign*(kin_mem->kin_sqrt_relfunc)*SUNMAX(SUNRabs(sutsv),sq1norm)/vtv;

  sigma_inv = ONE/sigma;

  /* compute the u-prime at which to evaluate the function func */
  N_VLinearSum(ONE, u, sigma, v, kin_mem->kin_vtemp1);
 
  /* call the system function to calculate func(u+sigma*v) */
  retval = kin_mem->kin_func(kin_mem->kin_vtemp1, kin_mem->kin_vtemp2,
                             kin_mem->kin_user_data);    
  kinspils_mem->nfes++;
  if (retval != 0) return(retval);

  /* finish the computation of the difference quotient */
  N_VLinearSum(sigma_inv, kin_mem->kin_vtemp2, -sigma_inv, kin_mem->kin_fval, Jv);

  return(0);
}


/*------------------------------------------------------------------
  kinSpilsInitialize performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)
  ------------------------------------------------------------------*/
int kinSpilsInitialize(KINMem kin_mem)
{
  int retval;
  KINSpilsMem kinspils_mem;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
		    "kinSpilsInitialize", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
		    "kinSpilsInitialize", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;
  
  /* initialize counters */
  kinSpilsInitializeCounters(kinspils_mem);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (kinspils_mem->jtimesDQ) {
    kinspils_mem->jtimes  = KINSpilsDQJtimes;
    kinspils_mem->jdata   = kin_mem;
  } else {
    kinspils_mem->jdata   = kin_mem->kin_user_data;
  }

  /* the Picard iteration is incompatible with a difference-quotient J*v */
  if ( (kin_mem->kin_globalstrategy == KIN_PICARD) && kinspils_mem->jtimesDQ ) {
    KINProcessError(kin_mem, KIN_ILL_INPUT, "KINSOL",
                    "kinSpilsInitialize", MSG_NOL_FAIL);
    return(KIN_ILL_INPUT);
  }

  /* if NOT preconditioning or do NOT need to setup the
     preconditioner, then set the lsetup function to NULL */
  if ((kinspils_mem->psolve == NULL) || 
      (kinspils_mem->pset == NULL)) {
    kin_mem->kin_lsetup = NULL;
  }

  /* Set scaling vectors assuming RIGHT preconditioning */
  /* NOTE: retval is non-zero only if LS == NULL        */
  retval = SUNLinSolSetScalingVectors(kinspils_mem->LS,
                                      kin_mem->kin_fscale,
                                      kin_mem->kin_fscale);
  if (retval != SUNLS_SUCCESS) {
    KINProcessError(kin_mem, KINSPILS_SUNLS_FAIL, "KINSPILS", "kinSpilsInitialize", 
		    "Error in calling SUNLinSolSetScalingVectors");
    return(KINSPILS_SUNLS_FAIL);
  }

  /* Call LS initialize routine */
  kinspils_mem->last_flag = SUNLinSolInitialize(kinspils_mem->LS);
  return(kinspils_mem->last_flag);
}


/*------------------------------------------------------------------
  kinSpilsSetup call the LS setup routine
  ------------------------------------------------------------------*/
int kinSpilsSetup(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;
  int retval;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
		    "kinSpilsSetup", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
		    "kinSpilsSetup", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* Call LS setup routine -- the LS will call KINSpilsPSetup
     if applicable */
  retval = SUNLinSolSetup(kinspils_mem->LS, NULL);

  /* save nni value from most recent lsetup call */
  kin_mem->kin_nnilset = kin_mem->kin_nni;

  return(retval);
}


/*------------------------------------------------------------------
  kinSpilsSolve interfaces between KINSOL and the generic
  SUNLinearSolver object
  ------------------------------------------------------------------*/
int kinSpilsSolve(KINMem kin_mem, N_Vector xx, N_Vector bb,
                  realtype *sJpnorm, realtype *sFdotJp)
{
  KINSpilsMem kinspils_mem;
  int nli_inc, retval;
  realtype res_norm;
  
  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) {
    KINProcessError(NULL, KINSPILS_MEM_NULL, "KINSPILS",
		    "kinSpilsSolve", MSGS_KINMEM_NULL);
    return(KINSPILS_MEM_NULL);
  }
  if (kin_mem->kin_lmem == NULL) {
    KINProcessError(kin_mem, KINSPILS_LMEM_NULL, "KINSPILS",
		    "kinSpilsSolve", MSGS_LMEM_NULL);
    return(KINSPILS_LMEM_NULL);
  }
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* Set initial guess x = 0 to LS */  
  N_VConst(ZERO, xx);
  
  /* set flag required for user Jacobian routine */
  kinspils_mem->new_uu = SUNTRUE;

  /* Call solver */
  retval = SUNLinSolSolve(kinspils_mem->LS, NULL, xx, bb, kin_mem->kin_eps);

  /* Retrieve solver statistics */
  res_norm = SUNLinSolResNorm(kinspils_mem->LS);
  nli_inc  = SUNLinSolNumIters(kinspils_mem->LS);

  if (kin_mem->kin_printfl > 2) 
    KINPrintInfo(kin_mem, PRNT_NLI, "KINSPILS", "kinSpilsSolve",
                 INFO_NLI, nli_inc);

  /* Increment counters nli and ncfl */
  kinspils_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) kinspils_mem->ncfl++;

  /* Interpret solver return value */
  kinspils_mem->last_flag = retval;

  if ( (retval != 0) && (retval != SUNLS_RES_REDUCED) ) {
    
    switch(retval) {
    case SUNLS_ATIMES_FAIL_REC:
    case SUNLS_PSOLVE_FAIL_REC:
      return(1);
      break;
    case SUNLS_MEM_NULL:
    case SUNLS_ILL_INPUT:
    case SUNLS_MEM_FAIL:
    case SUNLS_GS_FAIL:
    case SUNLS_CONV_FAIL:
    case SUNLS_QRFACT_FAIL:
    case SUNLS_LUFACT_FAIL:
    case SUNLS_QRSOL_FAIL:
      return(-1);
      break;
    case SUNLS_PACKAGE_FAIL_REC:
      KINProcessError(kin_mem, SUNLS_PACKAGE_FAIL_REC, "KINSPILS",
                      "kinSpilsSolve",
                      "Failure in SUNLinSol external package");
      return(-1);
      break;
    case SUNLS_PACKAGE_FAIL_UNREC:
      KINProcessError(kin_mem, SUNLS_PACKAGE_FAIL_UNREC, "KINSPILS", 
                      "kinSpilsSolve",
                      "Failure in SUNLinSol external package");
      return(-1);
      break;
    case SUNLS_ATIMES_FAIL_UNREC:
      KINProcessError(kin_mem, SUNLS_ATIMES_FAIL_UNREC, "KINSPILS", 
                      "kinSpilsSolve", MSGS_JTIMES_FAILED);    
      return(-1);
      break;
    case SUNLS_PSOLVE_FAIL_UNREC:
      KINProcessError(kin_mem, SUNLS_PSOLVE_FAIL_UNREC, "KINSPILS", 
                      "kinSpilsSolve", MSGS_PSOLVE_FAILED);
      return(-1);
      break;
    }
  }

  /*
    SUNLinSolSolve returned SUNLS_SUCCESS or SUNLS_RES_REDUCED
 
    Compute the terms sJpnorm and sFdotJp for use in the linesearch
    routine and in KINForcingTerm.  Both of these terms are subsequently
    corrected if the step is reduced by constraints or the linesearch.

    sJpnorm is the norm of the scaled product (scaled by fscale) of the
    current Jacobian matrix J and the step vector p (= solution vector xx).

    sFdotJp is the dot product of the scaled f vector and the scaled
    vector J*p, where the scaling uses fscale.
  */

  retval = KINSpilsATimes(kin_mem, xx, bb);
  if (retval > 0) {
    kinspils_mem->last_flag = SUNLS_ATIMES_FAIL_REC;
    return(1);
  }      
  else if (retval < 0) {
    kinspils_mem->last_flag = SUNLS_ATIMES_FAIL_UNREC;
    return(-1);
  }

  *sJpnorm = N_VWL2Norm(bb, kin_mem->kin_fscale);
  N_VProd(bb, kin_mem->kin_fscale, bb);
  N_VProd(bb, kin_mem->kin_fscale, bb);
  *sFdotJp = N_VDotProd(kin_mem->kin_fval, bb);

  if (kin_mem->kin_printfl > 2)
    KINPrintInfo(kin_mem, PRNT_EPS, "KINSPILS", "kinSpilsSolve",
                 INFO_EPS, res_norm, kin_mem->kin_eps);

  return(0);
}


/*------------------------------------------------------------------
  kinSpilsFree frees memory associated with the KINSpils system
  solver interface
  ------------------------------------------------------------------*/
int kinSpilsFree(KINMem kin_mem)
{
  KINSpilsMem kinspils_mem;

  /* Return immediately if kin_mem or kin_mem->kin_lmem are NULL */
  if (kin_mem == NULL) return (KINSPILS_SUCCESS);
  if (kin_mem->kin_lmem == NULL) return(KINSPILS_SUCCESS);
  kinspils_mem = (KINSpilsMem) kin_mem->kin_lmem;

  /* Free preconditioner memory (if applicable) */
  if (kinspils_mem->pfree) kinspils_mem->pfree(kin_mem);

  /* free KINSpils interface structure */
  free(kin_mem->kin_lmem);

  return(KINSPILS_SUCCESS);
}


/*------------------------------------------------------------------
  kinSpilsInitializeCounters resets counters for the SPILS interface 
  ------------------------------------------------------------------*/
int kinSpilsInitializeCounters(KINSpilsMem kinspils_mem)
{
  kinspils_mem->npe     = 0;
  kinspils_mem->nli     = 0;
  kinspils_mem->nps     = 0;
  kinspils_mem->ncfl    = 0;
  kinspils_mem->njtimes = 0;
  kinspils_mem->nfes    = 0;
  
  return(0);
}
