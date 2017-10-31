/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------
 * This is the implementation file for the IDASPILS linear solver
 * interface
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_spils_impl.h"
#include <sundials/sundials_math.h>

/* Private constants */
#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT05   RCONST(0.05)
#define PT9    RCONST(0.9)
#define ONE    RCONST(1.0)

/* Algorithmic constants */
#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */


/*===============================================================
  IDASPILS Exported functions -- Required
===============================================================*/

/*---------------------------------------------------------------
  IDASpilsSetLinearSolver specifies the iterative linear solver
  ---------------------------------------------------------------*/
int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS)
{
  int retval;
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if any input is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "IDASpilsSetLinearSolver", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  if (LS == NULL) {
    IDAProcessError(NULL, IDASPILS_ILL_INPUT, "IDASPILS", 
		    "IDASpilsSetLinearSolver", 
                    "LS must be non-NULL");
    return(IDASPILS_ILL_INPUT);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if solver and vector are compatible with SPILS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_ITERATIVE) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASPILS", 
                    "IDASpilsSetLinearSolver", 
                    "Non-iterative LS supplied to IDASpils interface");
    return(IDASPILS_ILL_INPUT);
  }
  if ( (IDA_mem->ida_tempv1->ops->nvdotprod == NULL) ||
       (IDA_mem->ida_tempv1->ops->nvconst == NULL) ||
       (IDA_mem->ida_tempv1->ops->nvscale == NULL) ||
       (IDA_mem->ida_tempv1->ops->nvlinearsum == NULL) ) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASPILS", 
                    "IDASpilsSetLinearSolver", MSGS_BAD_NVECTOR);
    return(IDASPILS_ILL_INPUT);
  }

  /* free any existing system solver attached to IDA */
  if (IDA_mem->ida_lfree)  IDA_mem->ida_lfree(IDA_mem);

  /* Set four main system linear solver function fields in IDA_mem */
  IDA_mem->ida_linit  = idaSpilsInitialize;
  IDA_mem->ida_lsetup = idaSpilsSetup;
  IDA_mem->ida_lsolve = idaSpilsSolve;
  IDA_mem->ida_lperf  = idaSpilsPerf;
  IDA_mem->ida_lfree  = idaSpilsFree;
  
  /* Get memory for IDASpilsMemRec */
  idaspils_mem = NULL;
  idaspils_mem = (IDASpilsMem) malloc(sizeof(struct IDASpilsMemRec));
  if (idaspils_mem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_MEM_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", MSGS_MEM_FAIL);
    return(IDASPILS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  idaspils_mem->LS = LS;
  
  /* Set defaults for Jacobian-related fields */
  idaspils_mem->jtimesDQ = SUNTRUE;
  idaspils_mem->jtsetup  = NULL;
  idaspils_mem->jtimes   = IDASpilsDQJtimes;
  idaspils_mem->jdata    = IDA_mem;

  /* Set defaults for preconditioner-related fields */
  idaspils_mem->pset   = NULL;
  idaspils_mem->psolve = NULL;
  idaspils_mem->pfree  = NULL;
  idaspils_mem->pdata  = IDA_mem->ida_user_data;

  /* Set default values for the rest of the Spils parameters */
  idaspils_mem->eplifac  = PT05;
  idaspils_mem->dqincfac = ONE;

  /* Initialize counters */
  idaSpilsInitializeCounters(idaspils_mem);

  /* Set default values for the rest of the SPILS parameters */
  idaspils_mem->last_flag = IDASPILS_SUCCESS;

  /* Attach default IDASpils interface routines to iterative LS */
  retval = SUNLinSolSetATimes(LS, IDA_mem, IDASpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    IDAProcessError(IDA_mem, IDASPILS_SUNLS_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetATimes");
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_SUNLS_FAIL);
  }
  retval = SUNLinSolSetPreconditioner(LS, IDA_mem, NULL, NULL);
  if (retval != SUNLS_SUCCESS) {
    IDAProcessError(IDA_mem, IDASPILS_SUNLS_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetPreconditioner");
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_SUNLS_FAIL);
  }

  /* Allocate memory for ytemp, yptemp and x */
  idaspils_mem->ytemp = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->ytemp == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_MEM_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", MSGS_MEM_FAIL);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  idaspils_mem->yptemp = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->yptemp == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_MEM_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", MSGS_MEM_FAIL);
    N_VDestroy(idaspils_mem->ytemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  idaspils_mem->x = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->x == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_MEM_FAIL, "IDASPILS", 
                    "IDASpilsSetLinearSolver", MSGS_MEM_FAIL);
    N_VDestroy(idaspils_mem->ytemp);
    N_VDestroy(idaspils_mem->yptemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, idaspils_mem->ytemp);
  idaspils_mem->sqrtN = SUNRsqrt( N_VDotProd(idaspils_mem->ytemp, 
                                             idaspils_mem->ytemp) );

  /* Attach linear solver memory to integrator memory */
  IDA_mem->ida_lmem = idaspils_mem;

  return(IDASPILS_SUCCESS);
}


/*===============================================================
  IDASPILS Exported functions -- Optional input/output
  ===============================================================*/


/*---------------------------------------------------------------*/
int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsSetEpsLin", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsSetEpsLin", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Check for legal eplifac */
  if (eplifac < ZERO) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASPILS",
                    "IDASpilsSetEpsLin", MSGS_NEG_EPLIFAC);
    return(IDASPILS_ILL_INPUT);
  }

  if (eplifac == ZERO)
    idaspils_mem->eplifac = PT05;
  else
    idaspils_mem->eplifac = eplifac;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsSetIncrementFactor(void *ida_mem, realtype dqincfac)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsSetIncrementFactor", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsSetIncrementFactor", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Check for legal maxrs */
  if (dqincfac <= ZERO) {
    IDAProcessError(IDA_mem, IDASPILS_ILL_INPUT, "IDASPILS",
                    "IDASpilsSetIncrementFactor", MSGS_NEG_DQINCFAC);
    return(IDASPILS_ILL_INPUT);
  }

  idaspils_mem->dqincfac = dqincfac;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsSetPreconditioner(void *ida_mem,
                              IDASpilsPrecSetupFn psetup,
                              IDASpilsPrecSolveFn psolve)
{
  int retval;
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  PSetupFn idaspils_psetup;
  PSolveFn idaspils_psolve;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsSetPreconditioner", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  idaspils_mem->pset   = psetup;
  idaspils_mem->psolve = psolve;

  /* notify iterative linear solver to call IDASpils interface routines */
  idaspils_psetup = (psetup == NULL) ? NULL : IDASpilsPSetup;
  idaspils_psolve = (psolve == NULL) ? NULL : IDASpilsPSolve;
  retval = SUNLinSolSetPreconditioner(idaspils_mem->LS, IDA_mem, 
                                      idaspils_psetup, idaspils_psolve);
  if (retval != SUNLS_SUCCESS) {
    IDAProcessError(IDA_mem, IDASPILS_SUNLS_FAIL, "IDASPILS", 
                    "IDASpilsSetPreconditioner", 
                    "Error in calling SUNLinSolSetPreconditioner");
    return(IDASPILS_SUNLS_FAIL);
  }

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsSetJacTimes(void *ida_mem,
                        IDASpilsJacTimesSetupFn jtsetup,
                        IDASpilsJacTimesVecFn jtimes)
{
  int retval;
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsSetJacTimes", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsSetJacTimes", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* store function pointers for user-supplied routines in IDASpils 
     interface (NULL jtimes implies use of DQ default) */
  if (jtimes != NULL) {
    idaspils_mem->jtimesDQ = SUNFALSE;
    idaspils_mem->jtimes   = jtimes;
  } else {
    idaspils_mem->jtimesDQ = SUNTRUE;
  }
  idaspils_mem->jtsetup = jtsetup;

  /* notify iterative linear solver to call IDASpils interface routine */
  retval = SUNLinSolSetATimes(idaspils_mem->LS, IDA_mem, IDASpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    IDAProcessError(IDA_mem, IDASPILS_SUNLS_FAIL, "IDASPILS", 
                    "IDASpilsSetJacTimes", 
                    "Error in calling SUNLinSolSetATimes");
    return(IDASPILS_SUNLS_FAIL);
  }

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetWorkSpace(void *ida_mem, long int *lenrwLS,
                         long int *leniwLS)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetWorkSpace", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* start with fixed sizes */
  *lenrwLS = 4;
  *leniwLS = 14;

  /* add N_Vector sizes */
  if (IDA_mem->ida_tempv1->ops->nvspace) {
    N_VSpace(IDA_mem->ida_tempv1, &lrw1, &liw1);
    *lenrwLS += 3*lrw1;
    *leniwLS += 3*liw1;
  }

  /* add LS sizes */
  if (idaspils_mem->LS->ops->space) {
    flag = SUNLinSolSpace(idaspils_mem->LS, &lrw, &liw);
    *lenrwLS += lrw;
    *leniwLS += liw;
  }

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumPrecEvals(void *ida_mem, long int *npevals)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumPrecEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *npevals = idaspils_mem->npe;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumPrecSolves(void *ida_mem, long int *npsolves)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumPrecSolves", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *npsolves = idaspils_mem->nps;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumLinIters(void *ida_mem, long int *nliters)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumLinIters", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *nliters = idaspils_mem->nli;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumConvFails(void *ida_mem, long int *nlcfails)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumConvFails", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *nlcfails = idaspils_mem->ncfl;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumJTSetupEvals(void *ida_mem, long int *njtsetups)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumJTSetupEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumJTSetupEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *njtsetups = idaspils_mem->njtsetup;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumJtimesEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *njvevals = idaspils_mem->njtimes;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetNumResEvals", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetNumResEvals", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *nrevalsLS = idaspils_mem->nres;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
int IDASpilsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS",
                    "IDASpilsGetLastFlag", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS",
                    "IDASpilsGetLastFlag", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  *flag = idaspils_mem->last_flag;

  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------*/
char *IDASpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDASPILS_SUCCESS:
    sprintf(name,"IDASPILS_SUCCESS");
    break; 
  case IDASPILS_MEM_NULL:
    sprintf(name,"IDASPILS_MEM_NULL");
    break;
  case IDASPILS_LMEM_NULL:
    sprintf(name,"IDASPILS_LMEM_NULL");
    break;
  case IDASPILS_ILL_INPUT:
    sprintf(name,"IDASPILS_ILL_INPUT");
    break;
  case IDASPILS_MEM_FAIL:
    sprintf(name,"IDASPILS_MEM_FAIL");
    break;
  case IDASPILS_PMEM_NULL:
    sprintf(name,"IDASPILS_PMEM_NULL");
    break;
  case IDASPILS_SUNLS_FAIL:
    sprintf(name,"IDASPILS_SUNLS_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*===============================================================
  IDASPILS Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  IDASpilsATimes:

  This routine generates the matrix-vector product z = Jv, where
  J is the system Jacobian, by calling either the user provided
  routine or the internal DQ routine.  The return value is 
  the same as the value returned by jtimes -- 
  0 if successful, nonzero otherwise.
  ---------------------------------------------------------------*/
int IDASpilsATimes(void *ida_mem, N_Vector v, N_Vector z)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  int jtflag;

  /* Return immediately if ida_mem or ida_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
                    "IDASpilsATimes", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
                    "IDASpilsATimes", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  jtflag = idaspils_mem->jtimes(IDA_mem->ida_tn,
                                idaspils_mem->ycur,
                                idaspils_mem->ypcur,
                                idaspils_mem->rcur, v, z,
                                IDA_mem->ida_cj,
                                idaspils_mem->jdata,
                                idaspils_mem->ytemp,
                                idaspils_mem->yptemp);
  idaspils_mem->njtimes++;
  return(jtflag);
}


/*---------------------------------------------------------------
  IDASpilsPSetup:

  This routine interfaces between the generic iterative linear 
  solvers and the user's psetup routine.  It passes to psetup all 
  required state information from ida_mem.  Its return value 
  is the same as that returned by psetup. Note that the generic
  iterative linear solvers guarantee that IDASpilsPSetup will only
  be called in the case that the user's psetup routine is non-NULL.
  ---------------------------------------------------------------*/
int IDASpilsPSetup(void *ida_mem)
{
  int         retval;
  IDAMem      IDA_mem;
  IDASpilsMem idaspils_mem;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "IDASpilsPSetup", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "IDASpilsPSetup", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Call user pset routine to update preconditioner and possibly 
     reset jcur (pass !jbad as update suggestion) */
  retval = idaspils_mem->pset(IDA_mem->ida_tn, 
                              idaspils_mem->ycur, 
                              idaspils_mem->ypcur, 
                              idaspils_mem->rcur, 
                              IDA_mem->ida_cj, 
                              idaspils_mem->pdata);
  idaspils_mem->npe++;
  return(retval);
}


/*---------------------------------------------------------------
  IDASpilsPSolve:

  This routine interfaces between the generic SUNLinSolSolve 
  routine and the user's psolve routine.  It passes to psolve all 
  required state information from ida_mem.  Its return value is
  the same as that returned by psolve.  Note that the generic 
  SUNLinSol solver guarantees that IDASilsPSolve will not be 
  called in the case in which preconditioning is not done. This 
  is the only case in which the user's psolve routine is allowed 
  to be NULL.
  ---------------------------------------------------------------*/
int IDASpilsPSolve(void *ida_mem, N_Vector r, N_Vector z,
                   realtype tol, int lr)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  int retval;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "IDASpilsPSolve", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "IDASpilsPSolve", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  retval = idaspils_mem->psolve(IDA_mem->ida_tn,
                                idaspils_mem->ycur,
                                idaspils_mem->ypcur,
                                idaspils_mem->rcur, r, z,
                                IDA_mem->ida_cj, tol,
                                idaspils_mem->pdata);
  idaspils_mem->nps++;
  return(retval);
}


/*---------------------------------------------------------------
  IDASpilsDQJtimes:

  This routine generates a difference quotient approximation to 
  the matrix-vector product z = Jv, where J is the system 
  Jacobian. The approximation is 
       Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,  
  where
       y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
       sigma = sqrt(Neq)*dqincfac.
  The return value from the call to res is saved in order to set 
  the return flag from IDASpilsSolve.
  ---------------------------------------------------------------*/
int IDASpilsDQJtimes(realtype tt, N_Vector yy, N_Vector yp,
                     N_Vector rr, N_Vector v, N_Vector Jv, 
                     realtype c_j, void *ida_mem, N_Vector work1, 
                     N_Vector work2)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  N_Vector y_tmp, yp_tmp;
  realtype sig=ZERO, siginv;
  int iter, retval;

  /* Return immediately if ida_mem or IDA_mem->ida_lmem are NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "IDASpilsDQJtimes", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "IDASpilsDQJtimes", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  sig = idaspils_mem->sqrtN * idaspils_mem->dqincfac;  /* GMRES */
  /*sig = idaspils_mem->dqincfac / N_VWrmsNorm(v, IDA_mem->ida_ewt);*/  /* BiCGStab/TFQMR */ 

  /* Rename work1 and work2 for readibility */
  y_tmp  = work1;
  yp_tmp = work2;

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
    N_VLinearSum(sig, v, ONE, yy, y_tmp);
    N_VLinearSum(c_j*sig, v, ONE, yp, yp_tmp);
    
    /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
    retval = IDA_mem->ida_res(tt, y_tmp, yp_tmp, Jv, IDA_mem->ida_user_data); 
    idaspils_mem->nres++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Set Jv to [Jv - rr]/sig and return. */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, rr, Jv);

  return(0);
}


/*---------------------------------------------------------------
 idaSpilsInitialize performs remaining initializations specific
 to the iterative linear solver interface (and solver itself)
---------------------------------------------------------------*/
int idaSpilsInitialize(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "idaSpilsInitialize", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "idaSpilsInitialize", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;
  
  idaSpilsInitializeCounters(idaspils_mem);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (idaspils_mem->jtimesDQ) {
    idaspils_mem->jtsetup = NULL;
    idaspils_mem->jtimes  = IDASpilsDQJtimes;
    idaspils_mem->jdata   = IDA_mem;
  } else {
    idaspils_mem->jdata   = IDA_mem->ida_user_data;
  }

  /* if psetup is not present, then idaSpilsSetup does not need to be 
     called, so set the lsetup function to NULL */
  if (idaspils_mem->pset == NULL)  IDA_mem->ida_lsetup = NULL;

  /* Call LS initialize routine */
  idaspils_mem->last_flag = SUNLinSolInitialize(idaspils_mem->LS);
  return(idaspils_mem->last_flag);
}


/*---------------------------------------------------------------
 idaSpilsSetup calls the LS 'setup' routine.
---------------------------------------------------------------*/
int idaSpilsSetup(IDAMem IDA_mem, N_Vector y, N_Vector yp, N_Vector r, 
                  N_Vector vt1, N_Vector vt2, N_Vector vt3)
{
  int  retval;
  IDASpilsMem idaspils_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "idaSpilsSetup", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "idaSpilsSetup", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Set IDASpils N_Vector pointers to inputs */
  idaspils_mem->ycur  = y;
  idaspils_mem->ypcur = yp;
  idaspils_mem->rcur  = r;

  /* Call LS setup routine -- the LS will call IDASpilsPSetup 
     if applicable */
  retval = SUNLinSolSetup(idaspils_mem->LS, NULL);
  return(retval);
}


/*---------------------------------------------------------------
 idaSpilsSolve: interfaces between IDA and the generic 
 SUNLinearSolver object LS, by setting the appropriate tolerance 
 and scaling vectors, calling the solver, and accumulating 
 statistics from the solve for use/reporting by IDA.
---------------------------------------------------------------*/
int idaSpilsSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                  N_Vector ycur, N_Vector ypcur, N_Vector rescur)
{
  IDASpilsMem idaspils_mem;
  int nli_inc, retval;
  
  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "idaSpilsSolve", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "idaSpilsSolve", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Set convergence test constant epslin, in terms of the
     Newton convergence test constant epsNewt and safety factors. The factor
     sqrt(Neq) assures that the convergence test is applied to the WRMS norm 
     of the residual vector, rather than the weighted L2 norm. */
  idaspils_mem->epslin = idaspils_mem->sqrtN *
    idaspils_mem->eplifac * IDA_mem->ida_epsNewt;

  /* Set vectors ycur, ypcur and rcur for use by the Atimes and 
     Psolve interface routines */
  idaspils_mem->ycur  = ycur;
  idaspils_mem->ypcur = ypcur;
  idaspils_mem->rcur  = rescur;

  /* Set initial guess x = 0 to LS */  
  N_VConst(ZERO, idaspils_mem->x);

  /* Set scaling vectors for LS to use */
  retval = SUNLinSolSetScalingVectors(idaspils_mem->LS, weight, weight);
  if (retval != SUNLS_SUCCESS) {
    IDAProcessError(IDA_mem, IDASPILS_SUNLS_FAIL, "IDASPILS", "idaSpilsSolve", 
                    "Error in calling SUNLinSolSetScalingVectors");
    return(IDASPILS_SUNLS_FAIL);
  }

  /* If a user-provided jtsetup routine is supplied, call that here */
  if (idaspils_mem->jtsetup) {
    retval = idaspils_mem->jtsetup(IDA_mem->ida_tn, ycur, ypcur, rescur,
                                   IDA_mem->ida_cj, idaspils_mem->jdata);
    idaspils_mem->njtsetup++;
    if (retval != 0) {
      IDAProcessError(IDA_mem, retval, "IDASPILS", 
                      "idaSpilsSolve", MSGS_JTSETUP_FAILED);
      return(retval);
    }
  }
  
  /* Call solver */
  retval = SUNLinSolSolve(idaspils_mem->LS, NULL, idaspils_mem->x,
                          b, idaspils_mem->epslin);

  /* Retrieve solver statistics */
  nli_inc = SUNLinSolNumIters(idaspils_mem->LS);
  
  /* Copy x (or preconditioned residual vector if no iterations required) to b */
  if (nli_inc == 0) N_VScale(ONE, SUNLinSolResid(idaspils_mem->LS), b);
  else N_VScale(ONE, idaspils_mem->x, b);

  /* Increment counters nli and ncfl */
  idaspils_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) idaspils_mem->ncfl++;

  /* Interpret solver return value  */
  idaspils_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
  case SUNLS_CONV_FAIL:
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
    IDAProcessError(IDA_mem, SUNLS_PACKAGE_FAIL_UNREC, "IDASPILS", 
                    "idaSpilsSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    IDAProcessError(IDA_mem, SUNLS_PSOLVE_FAIL_UNREC, "IDASPILS", 
                    "idaSpilsSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }
  
  return(0); 
}


/*---------------------------------------------------------------
 idaSpilsPerf: accumulates performance statistics information 
 for IDA
---------------------------------------------------------------*/
int idaSpilsPerf(IDAMem IDA_mem, int perftask)
{
  IDASpilsMem idaspils_mem;
  realtype rcfn, rcfl;
  long int nstd, nnid;
  booleantype lcfn, lcfl;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPILS", 
		    "idaSpilsPerf", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASPILS_LMEM_NULL, "IDASPILS", 
		    "idaSpilsPerf", MSGS_LMEM_NULL);
    return(IDASPILS_LMEM_NULL);
  }
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* when perftask == 0, store current performance statistics */
  if (perftask == 0) {
    idaspils_mem->nst0  = IDA_mem->ida_nst;
    idaspils_mem->nni0  = IDA_mem->ida_nni;
    idaspils_mem->ncfn0 = IDA_mem->ida_ncfn;
    idaspils_mem->ncfl0 = idaspils_mem->ncfl;  
    idaspils_mem->nwarn = 0;
    return(0);
  }

  /* Compute statistics since last call

     Note: the performance monitor that checked whether the average 
       number of linear iterations was too close to maxl has been 
       removed, since the 'maxl' value is no longer owned by the 
       IDASpils interface.
   */
  nstd = IDA_mem->ida_nst - idaspils_mem->nst0;
  nnid = IDA_mem->ida_nni - idaspils_mem->nni0;
  if (nstd == 0 || nnid == 0) return(0);
  
  rcfn = (realtype) ( (IDA_mem->ida_ncfn - idaspils_mem->ncfn0) /
                      ((realtype) nstd) );
  rcfl = (realtype) ( (idaspils_mem->ncfl - idaspils_mem->ncfl0) /
                      ((realtype) nnid) );
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lcfn || lcfl)) return(0);
  idaspils_mem->nwarn++;
  if (idaspils_mem->nwarn > 10) return(1);
  if (lcfn) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPILS", "idaSpilsPerf",
                    MSGS_CFN_WARN, IDA_mem->ida_tn, rcfn);
  if (lcfl) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPILS", "idaSpilsPerf",
                    MSGS_CFL_WARN, IDA_mem->ida_tn, rcfl);
  return(0);
}


/*---------------------------------------------------------------
 idaSpilsFree frees memory associates with the IDASpils system
 solver interface.
---------------------------------------------------------------*/
int idaSpilsFree(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;

  /* Return immediately if IDA_mem or IDA_mem->ida_lmem are NULL */
  if (IDA_mem == NULL)  return (IDASPILS_SUCCESS);
  if (IDA_mem->ida_lmem == NULL)  return(IDASPILS_SUCCESS);
  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Free N_Vector memory */
  if (idaspils_mem->ytemp) {
    N_VDestroy(idaspils_mem->ytemp);
    idaspils_mem->ytemp = NULL;
  }
  if (idaspils_mem->yptemp) {
    N_VDestroy(idaspils_mem->yptemp);
    idaspils_mem->yptemp = NULL;
  }
  if (idaspils_mem->x) {
    N_VDestroy(idaspils_mem->x);
    idaspils_mem->x = NULL;
  }

  /* Nullify other N_Vector pointers */
  idaspils_mem->ycur  = NULL;
  idaspils_mem->ypcur = NULL;
  idaspils_mem->rcur  = NULL;

  /* Free preconditioner memory (if applicable) */
  if (idaspils_mem->pfree)  idaspils_mem->pfree(IDA_mem);
  
  /* free IDASpils interface structure */
  free(IDA_mem->ida_lmem);
  
  return(IDASPILS_SUCCESS);
}


/*---------------------------------------------------------------
 idaSpilsInitializeCounters resets all counters from an 
 IDASpilsMem structure.
---------------------------------------------------------------*/
int idaSpilsInitializeCounters(IDASpilsMem idaspils_mem)
{
  idaspils_mem->npe      = 0;
  idaspils_mem->nli      = 0;
  idaspils_mem->nps      = 0;
  idaspils_mem->ncfl     = 0;
  idaspils_mem->nres     = 0;
  idaspils_mem->njtsetup = 0;
  idaspils_mem->njtimes  = 0;

  return(0);
}
