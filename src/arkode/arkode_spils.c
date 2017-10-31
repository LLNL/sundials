/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This is the implementation file for the ARKSPILS linear solver
 * interface.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_spils_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define MAX_DQITERS  3  /* max. # of attempts to recover in DQ J*v */
#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)


/*===============================================================
  ARKSPILS Exported functions -- Required
===============================================================*/

/*---------------------------------------------------------------
 ARKSpilsSetLinearSolver specifies the iterative linear solver.
---------------------------------------------------------------*/
int ARKSpilsSetLinearSolver(void *arkode_mem, SUNLinearSolver LS)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if any input is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetLinearSolver", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (LS == NULL) {
    arkProcessError(NULL, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetLinearSolver", 
                    "LS must be non-NULL");
    return(ARKSPILS_ILL_INPUT);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if solver and vector are compatible with SPILS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_ITERATIVE) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", 
                    "Non-iterative LS supplied to ARKSpils interface");
    return(ARKSPILS_ILL_INPUT);
  }
  if ( (ark_mem->ark_tempv->ops->nvlinearsum == NULL) ||
       (ark_mem->ark_tempv->ops->nvconst == NULL) ||
       (ark_mem->ark_tempv->ops->nvdotprod == NULL) ){
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  /* free any existing system solver attached to ARKode */
  if (ark_mem->ark_lfree)  ark_mem->ark_lfree(ark_mem);

  /* Set four main system linear solver function fields in ark_mem */
  ark_mem->ark_linit  = arkSpilsInitialize;
  ark_mem->ark_lsetup = arkSpilsSetup;
  ark_mem->ark_lsolve = arkSpilsSolve;
  ark_mem->ark_lfree  = arkSpilsFree;
  ark_mem->ark_lsolve_type = 0;
  
  /* Get memory for ARKSpilsMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMem) malloc(sizeof(struct ARKSpilsMemRec));
  if (arkspils_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  arkspils_mem->LS = LS;
  
  /* Set defaults for Jacobian-related fields */
  arkspils_mem->jtimesDQ = SUNTRUE;
  arkspils_mem->jtsetup = NULL;
  arkspils_mem->jtimes = ARKSpilsDQJtimes;
  arkspils_mem->j_data = ark_mem;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->pset   = NULL;
  arkspils_mem->psolve = NULL;
  arkspils_mem->pfree  = NULL;
  arkspils_mem->P_data = ark_mem->ark_user_data;

  /* Initialize counters */
  arkSpilsInitializeCounters(arkspils_mem);

  /* Set default values for the rest of the SPILS parameters */
  arkspils_mem->jbad = SUNTRUE;
  arkspils_mem->eplifac = ARKSPILS_EPLIN;
  arkspils_mem->last_flag = ARKSPILS_SUCCESS;

  /* Attach default ARKSpils interface routines to iterative LS */
  retval = SUNLinSolSetATimes(LS, ark_mem, ARKSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetATimes");
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_SUNLS_FAIL);
  }
  retval = SUNLinSolSetPreconditioner(LS, ark_mem, NULL, NULL);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", 
                    "Error in calling SUNLinSolSetPreconditioner");
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_SUNLS_FAIL);
  }

  /* Allocate memory for ytemp and x */
  arkspils_mem->ytemp = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->ytemp == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  arkspils_mem->x = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->x == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPILS", 
                    "ARKSpilsSetLinearSolver", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->ytemp);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, arkspils_mem->ytemp);
  arkspils_mem->sqrtN = SUNRsqrt( N_VDotProd(arkspils_mem->ytemp, 
                                             arkspils_mem->ytemp) );

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassLinearSolver specifies the iterative mass-matrix
 linear solver and user-supplied routine to perform the 
 mass-matrix-vector product.
---------------------------------------------------------------*/
int ARKSpilsSetMassLinearSolver(void *arkode_mem,
                                SUNLinearSolver LS,
                                booleantype time_dep)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem, LS or mtimes are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassLinearSolver", 
                    MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (LS == NULL) {
    arkProcessError(NULL, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetLinearSolver", 
                    "LS must be non-NULL");
    return(ARKSPILS_ILL_INPUT);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if solver and vector are compatible with SPILS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_ITERATIVE) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", 
                    "Non-iterative LS supplied to ARKSpils interface");
    return(ARKSPILS_ILL_INPUT);
  }
  if ( (ark_mem->ark_tempv->ops->nvconst == NULL) ||
       (ark_mem->ark_tempv->ops->nvdotprod == NULL) ){
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  /* free any existing system solver attached to ARKode */
  if (ark_mem->ark_mfree)  ark_mem->ark_mfree(ark_mem);

  /* Set four main system linear solver function fields in ark_mem */
  ark_mem->ark_minit  = arkSpilsMassInitialize;
  ark_mem->ark_msetup = arkSpilsMassSetup;
  ark_mem->ark_mmult  = arkSpilsMassMult;
  ark_mem->ark_msolve = arkSpilsMassSolve;
  ark_mem->ark_mfree  = arkSpilsMassFree;
  ark_mem->ark_msolve_type = 0;
  
  /* notify arkode of non-identity mass matrix */
  ark_mem->ark_mass_matrix = SUNTRUE;

  /* Get memory for ARKSpilsMassMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMassMem) malloc(sizeof(struct ARKSpilsMassMemRec));
  if (arkspils_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer; flag indicating time-dependence */
  arkspils_mem->LS = LS;
  arkspils_mem->time_dependent = time_dep;

  /* Set mass-matrix-vector product routines to NULL */
  arkspils_mem->mtsetup = NULL;
  arkspils_mem->mtimes  = NULL;
  arkspils_mem->mt_data = NULL;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->pset   = NULL;
  arkspils_mem->psolve = NULL;
  arkspils_mem->pfree  = NULL;
  arkspils_mem->P_data = ark_mem->ark_user_data;

  /* Initialize counters */
  arkSpilsInitializeMassCounters(arkspils_mem);

  /* Set default values for the rest of the SPILS parameters */
  arkspils_mem->eplifac   = ARKSPILS_EPLIN;
  arkspils_mem->last_flag = ARKSPILS_SUCCESS;

  /* Attach default ARKSpils interface routines to iterative LS */
  retval = SUNLinSolSetATimes(LS, ark_mem, NULL);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", 
                    "Error in calling SUNLinSolSetATimes");
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_SUNLS_FAIL);
  }
  retval = SUNLinSolSetPreconditioner(LS, ark_mem, NULL, NULL);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", 
                    "Error in calling SUNLinSolSetPreconditioner");
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_SUNLS_FAIL);
  }

  /* Allocate memory for x */
  arkspils_mem->x = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->x == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassLinearSolver", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, arkspils_mem->x);
  arkspils_mem->sqrtN = SUNRsqrt( N_VDotProd(arkspils_mem->x, 
                                             arkspils_mem->x) );

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}


/*===============================================================
  ARKSPILS Exported functions -- Optional input/output
===============================================================*/

/*---------------------------------------------------------------
 ARKSpilsSetEpsLin
---------------------------------------------------------------*/
int ARKSpilsSetEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetEpsLin", MSGS_BAD_EPLIN);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->eplifac = (eplifac == ZERO) ? ARKSPILS_EPLIN : eplifac;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetPreconditioner
---------------------------------------------------------------*/
int ARKSpilsSetPreconditioner(void *arkode_mem, 
			      ARKSpilsPrecSetupFn psetup, 
			      ARKSpilsPrecSolveFn psolve)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  PSetupFn arkspils_psetup;
  PSolveFn arkspils_psolve;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPreconditioner", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* store function pointers for user-supplied routines in ARKSpils interface */
  arkspils_mem->pset   = psetup;
  arkspils_mem->psolve = psolve;

  /* notify iterative linear solver to call ARKSpils interface routines */
  arkspils_psetup = (psetup == NULL) ? NULL : ARKSpilsPSetup;
  arkspils_psolve = (psolve == NULL) ? NULL : ARKSpilsPSolve;
  retval = SUNLinSolSetPreconditioner(arkspils_mem->LS, ark_mem, 
                                      arkspils_psetup, arkspils_psolve);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetPreconditioner", 
                    "Error in calling SUNLinSolSetPreconditioner");
    return(ARKSPILS_SUNLS_FAIL);
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetJacTimes
---------------------------------------------------------------*/
int ARKSpilsSetJacTimes(void *arkode_mem, 
                        ARKSpilsJacTimesSetupFn jtsetup,
                        ARKSpilsJacTimesVecFn jtimes)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsSetJacTimes", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsSetJacTimes", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* store function pointers for user-supplied routines in ARKSpils 
     interface (NULL jtimes implies use of DQ default) */
  if (jtimes != NULL) {
    arkspils_mem->jtimesDQ = SUNFALSE;
    arkspils_mem->jtimes   = jtimes;
  } else {
    arkspils_mem->jtimesDQ = SUNTRUE;
  }
  arkspils_mem->jtsetup = jtsetup;

  /* notify iterative linear solver to call ARKSpils interface routines */
  retval = SUNLinSolSetATimes(arkspils_mem->LS, ark_mem, ARKSpilsATimes);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetJacTimes", 
                    "Error in calling SUNLinSolSetATimes");
    return(ARKSPILS_SUNLS_FAIL);
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetWorkSpace
---------------------------------------------------------------*/
int ARKSpilsGetWorkSpace(void *arkode_mem, long int *lenrw, 
			 long int *leniw)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetWorkSpace", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* start with fixed sizes plus NVectors */
  *lenrw = 4;
  *leniw = 10;

  /* add NVector sizes */
  if (ark_mem->ark_tempv->ops->nvspace) {
    N_VSpace(ark_mem->ark_tempv, &lrw1, &liw1);
    *lenrw += 2*lrw1;
    *leniw += 2*liw1;
  }

  /* add LS sizes */
  if (arkspils_mem->LS->ops->space) {
    flag = SUNLinSolSpace(arkspils_mem->LS, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumPrecEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumPrecEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *npevals = arkspils_mem->npe;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumPrecSolves
---------------------------------------------------------------*/
int ARKSpilsGetNumPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumPrecSolves", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *npsolves = arkspils_mem->nps;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumLinIters
---------------------------------------------------------------*/
int ARKSpilsGetNumLinIters(void *arkode_mem, long int *nliters)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumLinIters", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nliters = arkspils_mem->nli;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumConvFails
---------------------------------------------------------------*/
int ARKSpilsGetNumConvFails(void *arkode_mem, long int *nlcfails)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumConvFails", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nlcfails = arkspils_mem->ncfl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumJTSetupEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumJTSetupEvals(void *arkode_mem, long int *njtsetups)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumJTSetupEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumJTSetupEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *njtsetups = arkspils_mem->njtsetup;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumJtimesEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumJtimesEvals(void *arkode_mem, long int *njvevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumJtimesEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *njvevals = arkspils_mem->njtimes;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumRhsEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumRhsEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *nfevalsLS = arkspils_mem->nfes;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetLastFlag
---------------------------------------------------------------*/
int ARKSpilsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetLastFlag", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  *flag = arkspils_mem->last_flag;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetReturnFlagName -- ADD SIMILAR ROUTINE TO GENERIC SUNLinearSolver MODULE
---------------------------------------------------------------*/
char *ARKSpilsGetReturnFlagName(long int flag)
{
  char *name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKSPILS_SUCCESS:
    sprintf(name,"ARKSPILS_SUCCESS");
    break; 
  case ARKSPILS_MEM_NULL:
    sprintf(name,"ARKSPILS_MEM_NULL");
    break;
  case ARKSPILS_LMEM_NULL:
    sprintf(name,"ARKSPILS_LMEM_NULL");
    break;
  case ARKSPILS_ILL_INPUT:
    sprintf(name,"ARKSPILS_ILL_INPUT");
    break;
  case ARKSPILS_MEM_FAIL:
    sprintf(name,"ARKSPILS_MEM_FAIL");
    break;
  case ARKSPILS_PMEM_NULL:
    sprintf(name,"ARKSPILS_PMEM_NULL");
    break;
  case ARKSPILS_MASSMEM_NULL:
    sprintf(name,"ARKSPILS_MASSMEM_NULL");
    break;
  case ARKSPILS_SUNLS_FAIL:
    sprintf(name,"ARKSPILS_SUNLS_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassEpsLin
---------------------------------------------------------------*/
int ARKSpilsSetMassEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassEpsLin", MSGS_BAD_EPLIN);
    return(ARKSPILS_ILL_INPUT);
  }
  arkspils_mem->eplifac = (eplifac == ZERO) ? ARKSPILS_EPLIN : eplifac;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassPreconditioner
---------------------------------------------------------------*/
int ARKSpilsSetMassPreconditioner(void *arkode_mem, 
				  ARKSpilsMassPrecSetupFn psetup, 
				  ARKSpilsMassPrecSolveFn psolve)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;
  PSetupFn arkspils_mpsetup;
  PSolveFn arkspils_mpsolve;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPreconditioner", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassPreconditioner", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* store function pointers for user-supplied routines in ARKSpils interface */
  arkspils_mem->pset   = psetup;
  arkspils_mem->psolve = psolve;

  /* notify iterative linear solver to call ARKSpils interface routines */
  arkspils_mpsetup = (psetup == NULL) ? NULL : ARKSpilsMPSetup;
  arkspils_mpsolve = (psolve == NULL) ? NULL : ARKSpilsMPSolve;
  retval = SUNLinSolSetPreconditioner(arkspils_mem->LS, ark_mem, 
                                      arkspils_mpsetup, arkspils_mpsolve);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassPreconditioner", 
                    "Error in calling SUNLinSolSetPreconditioner");
    return(ARKSPILS_SUNLS_FAIL);
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsSetMassTimes
---------------------------------------------------------------*/
int ARKSpilsSetMassTimes(void *arkode_mem, 
                         ARKSpilsMassTimesSetupFn mtsetup,
                         ARKSpilsMassTimesVecFn mtimes,
                         void *mtimes_data)
{
  int retval;
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem, ark_mem->ark_mass_mem or mtimes are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassTimes", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsSetMassTimes", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;
  if (mtimes == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		    "ARKSpilsSetMassTimes", 
                    "non-NULL mtimes function must be supplied");
    return(ARKSPILS_ILL_INPUT);
  }

  /* store pointers for user-supplied routines and data structure
     in ARKSpils interface */
  arkspils_mem->mtsetup = mtsetup;
  arkspils_mem->mtimes  = mtimes;
  arkspils_mem->mt_data = mtimes_data;

  /* notify iterative linear solver to call ARKSpils interface routines */
  retval = SUNLinSolSetATimes(arkspils_mem->LS, ark_mem, ARKSpilsMTimes);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", 
                    "ARKSpilsSetMassTimes", 
                    "Error in calling SUNLinSolSetATimes");
    return(ARKSPILS_SUNLS_FAIL);
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetMassWorkSpace
---------------------------------------------------------------*/
int ARKSpilsGetMassWorkSpace(void *arkode_mem, long int *lenrw, 
			     long int *leniw)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetMassWorkSpace", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetMassWorkSpace", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* start with fixed sizes */
  *lenrw = 4;
  *leniw = 8;

  /* add NVector sizes */
  if (ark_mem->ark_tempv->ops->nvspace) {
    N_VSpace(ark_mem->ark_tempv, &lrw1, &liw1);
    *lenrw += lrw1;
    *leniw += liw1;
  }
  
  /* add LS sizes */
  if (arkspils_mem->LS->ops->space) {
    flag = SUNLinSolSpace(arkspils_mem->LS, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassPrecEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumMassPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecEvals", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *npevals = arkspils_mem->npe;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassPrecSolves
---------------------------------------------------------------*/
int ARKSpilsGetNumMassPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecSolves", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassPrecSolves", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *npsolves = arkspils_mem->nps;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassIters
---------------------------------------------------------------*/
int ARKSpilsGetNumMassIters(void *arkode_mem, long int *nmiters)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassLinIters", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassLinIters", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmiters = arkspils_mem->nli;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMassConvFails
---------------------------------------------------------------*/
int ARKSpilsGetNumMassConvFails(void *arkode_mem, long int *nmcfails)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassConvFails", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMassConvFails", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmcfails = arkspils_mem->ncfl;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMTSetups
---------------------------------------------------------------*/
int ARKSpilsGetNumMTSetups(void *arkode_mem, long int *nmtsetups)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMTSetups", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMTSetups", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmtsetups = arkspils_mem->nmtsetup;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetNumMtimesEvals
---------------------------------------------------------------*/
int ARKSpilsGetNumMtimesEvals(void *arkode_mem, long int *nmvevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMtimesEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetNumMtimesEvals", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *nmvevals = arkspils_mem->nmtimes;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSpilsGetLastMassFlag
---------------------------------------------------------------*/
int ARKSpilsGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastMassFlag", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsGetLastMassFlag", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  *flag = arkspils_mem->last_flag;

  return(ARKSPILS_SUCCESS);
}


/*===============================================================
  ARKSPILS Private functions
===============================================================*/

/*---------------------------------------------------------------
 ARKSpilsATimes:

 This routine generates the matrix-vector product z = Av, where
 A = M - gamma*J. The product M*v is obtained either by calling 
 the mtimes routine or by just using v (if M=I).  The product 
 J*v is obtained by calling the jtimes routine. It is then scaled 
 by -gamma and added to M*v to obtain A*v. The return value is 
 the same as the values returned by jtimes and mtimes -- 
 0 if successful, nonzero otherwise.
---------------------------------------------------------------*/
int ARKSpilsATimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int jtflag, mtflag;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
                    "ARKSpilsATimes", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
                    "ARKSpilsATimes", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  jtflag = arkspils_mem->jtimes(v, z, ark_mem->ark_tn, 
                                arkspils_mem->ycur, 
                                arkspils_mem->fcur, 
                                arkspils_mem->j_data, 
                                arkspils_mem->ytemp);
  arkspils_mem->njtimes++;
  if (jtflag != 0) return(jtflag);

  /* Compute mass matrix vector product and add to result */
  if (ark_mem->ark_mass_matrix) {
    mtflag = ARKSpilsMTimes(arkode_mem, v, arkspils_mem->ytemp);
    if (mtflag != 0) return(mtflag);
    N_VLinearSum(ONE, arkspils_mem->ytemp, -ark_mem->ark_gamma, z, z);
  } else {
    N_VLinearSum(ONE, v, -ark_mem->ark_gamma, z, z);
  }

  return(0);
}

/*---------------------------------------------------------------
 ARKSpilsPSetup:

 This routine interfaces between the generic iterative linear 
 solvers and the user's psetup routine.  It passes to psetup all 
 required state information from arkode_mem.  Its return value 
 is the same as that returned by psetup. Note that the generic
 iterative linear solvers guarantee that ARKSpilsPSetup will only
 be called in the case that the user's psetup routine is non-NULL.
---------------------------------------------------------------*/
int ARKSpilsPSetup(void *arkode_mem)
{
  int         retval;
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsPSetup", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsPSetup", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Call user pset routine to update preconditioner and possibly 
     reset jcur (pass !jbad as update suggestion) */
  retval = arkspils_mem->pset(ark_mem->ark_tn, 
                              arkspils_mem->ycur, 
                              arkspils_mem->fcur, 
                              !(arkspils_mem->jbad),
                              &ark_mem->ark_jcur,
                              ark_mem->ark_gamma, 
                              arkspils_mem->P_data);
  return(retval);     
}

/*---------------------------------------------------------------
 ARKSpilsPSolve:

 This routine interfaces between the generic SUNLinSolSolve 
 routine and the user's psolve routine.  It passes to psolve all 
 required state information from arkode_mem.  Its return value 
 is the same as that returned by psolve. Note that the generic 
 SUNLinSol solver guarantees that ARKSpilsPSolve will not be 
 called in the case in which preconditioning is not done. This 
 is the only case in which the user's psolve routine is allowed 
 to be NULL.
---------------------------------------------------------------*/
int ARKSpilsPSolve(void *arkode_mem, N_Vector r, N_Vector z,
                   realtype tol, int lr)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int retval;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsPSolve", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsPSolve", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* call the user-supplied psolve routine, and accumulate count */
  retval = arkspils_mem->psolve(ark_mem->ark_tn, 
                                arkspils_mem->ycur, 
                                arkspils_mem->fcur, r, z, 
                                ark_mem->ark_gamma, 
                                tol, lr, 
                                arkspils_mem->P_data);
  arkspils_mem->nps++;

  return(retval);     
}

/*---------------------------------------------------------------
 ARKSpilsMTimes:

 This routine generates the matrix-vector product z = Mv, where
 M is the system mass matrix, by calling the user-supplied mtimes
 routine.. The return value is the same as the value returned 
 by mtimes -- 0 if successful, nonzero otherwise.
---------------------------------------------------------------*/
int ARKSpilsMTimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem       ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int retval;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsMTimes", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsMTimes", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* call user-supplied mtimes routine and increment counter */
  retval = arkspils_mem->mtimes(v, z, ark_mem->ark_tn, 
                                arkspils_mem->mt_data);
  arkspils_mem->nmtimes++;
  return(retval);
}


/*---------------------------------------------------------------
 ARKSpilsMPSetup:

 This routine interfaces between the generic iterative linear 
 solver and the user's mass matrix psetup routine.  It passes to
 psetup all required state information from arkode_mem.  Its 
 return value is the same as that returned by psetup. Note that 
 the generic iterative linear solvers guarantee that 
 ARKSpilsMPSetup will only be called if the user's psetup 
 routine is non-NULL.
---------------------------------------------------------------*/
int ARKSpilsMPSetup(void *arkode_mem)
{
  ARKodeMem       ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int retval;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsMPSetup", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsMPSetup", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* only proceed if the mass matrix is time-dependent or if 
     mpsetup has not been called previously */
  if (!arkspils_mem->time_dependent && arkspils_mem->npe) 
    return(0);
  
  /* call user-supplied pset routine and increment counter */
  retval = arkspils_mem->pset(ark_mem->ark_tn, 
                              arkspils_mem->P_data);  
  arkspils_mem->npe++;
  return(retval);     
}


/*---------------------------------------------------------------
 ARKSpilsMPSolve:

 This routine interfaces between the generic Sp***Solve routine
 (within the SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG solver) and the 
 user's mass matrix psolve routine.  It passes to psolve all 
 required state information from arkode_mem.  Its return value is
 the same as that returned by psolve. Note that the generic SP*** 
 solver guarantees that ARKSpilsMPSolve will not be called in the 
 case in which preconditioning is not done. This is the only case 
 in which the user's psolve routine is allowed to be NULL.
---------------------------------------------------------------*/
int ARKSpilsMPSolve(void *arkode_mem, N_Vector r, N_Vector z,
                    realtype tol, int lr)
{
  ARKodeMem       ark_mem;
  ARKSpilsMassMem arkspils_mem;
  int retval;

  /* Return immediately if arkode_mem or ark_mem->ark_mass_mem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsMPSolve", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "ARKSpilsMPSolve", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* call the user-supplied psolve routine, and accumulate count */
  retval = arkspils_mem->psolve(ark_mem->ark_tn, r, z, tol, lr, 
                                arkspils_mem->P_data);
  arkspils_mem->nps++;
  return(retval);     
}


/*---------------------------------------------------------------
 ARKSpilsDQJtimes:

 This routine generates a difference quotient approximation to
 the Jacobian times vector f_y(t,y) * v. The approximation is 
 Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) 
 is input, i.e. the WRMS norm of v/vnrm is 1.
---------------------------------------------------------------*/
int ARKSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
		     N_Vector y, N_Vector fy,
		     void *arkode_mem, N_Vector work)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* Return immediately if arkode_mem or ark_mem->ark_lmem are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "ARKSpilsDQJtimes", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "ARKSpilsDQJtimes", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ark_mem->ark_ewt);

  for (iter=0; iter<MAX_DQITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = ark_mem->ark_fi(t, work, Jv, ark_mem->ark_user_data); 
    arkspils_mem->nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    /* If fi failed recoverably, shrink sig and retry */
    sig *= PT25;

  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}


/*---------------------------------------------------------------
 arkSpilsInitialize performs remaining initializations specific
 to the iterative linear solver interface (and solver itself)
---------------------------------------------------------------*/
int arkSpilsInitialize(ARKodeMem ark_mem)
{
  int retval;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsInitialize", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "arkSpilsInitialize", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  
  arkSpilsInitializeCounters(arkspils_mem);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (arkspils_mem->jtimesDQ) {
    arkspils_mem->jtsetup = NULL;
    arkspils_mem->jtimes = ARKSpilsDQJtimes;
    arkspils_mem->j_data = ark_mem;
  } else {
    arkspils_mem->j_data = ark_mem->ark_user_data;
  }

  /* if psetup is not present, then arkSpilsSetup does not need to be 
     called, so set the lsetup function to NULL */
  if (arkspils_mem->pset == NULL)  ark_mem->ark_lsetup = NULL;

  /* Ensure that if a mass matrix / solver are used, that system 
     and mass matrix solver types match (and both are set up) */
  if (ark_mem->ark_mass_matrix) {

    /* access mass matrix solver interface */
    if (ark_mem->ark_mass_mem == NULL) {
      arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
                      "arkSpilsInitialize", MSGS_MASSMEM_NULL);
      return(ARKSPILS_MASSMEM_NULL);
    }

    /* check that ark_msolve_type is compatible */
    if (ark_mem->ark_msolve_type != 0) {
      arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
		      "arkSpilsInitialize",
                      "Spils and Spils solvers cannot be combined");
      arkspils_mem->last_flag = ARKSPILS_ILL_INPUT;
      return(-1);
    }

    /* initialize mass matrix linear solver */
    retval = arkSpilsMassInitialize(ark_mem);
    if (retval != ARKSPILS_SUCCESS) {
      arkspils_mem->last_flag = retval;
      return(arkspils_mem->last_flag);
    }
  }

  /* Call LS initialize routine */
  arkspils_mem->last_flag = SUNLinSolInitialize(arkspils_mem->LS);
  return(arkspils_mem->last_flag);
}


/*---------------------------------------------------------------
 arkSpilsSetup calls the LS 'setup' routine.
---------------------------------------------------------------*/
int arkSpilsSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                  N_Vector fpred, booleantype *jcurPtr, 
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype dgamma;
  int  retval;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsSetup", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "arkSpilsSetup", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Set ARKSpils N_Vector pointers to current solution and rhs */
  arkspils_mem->ycur = ypred;
  arkspils_mem->fcur = fpred;

  /* Use nst, gamma/gammap, and convfail to set J/P eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  arkspils_mem->jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkspils_mem->nstlpre + ARKSPILS_MSBPRE) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKSPILS_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  *jcurPtr = arkspils_mem->jbad;
  
  /* Call LS setup routine -- the LS will call ARKSpilsPSetup, who will 
     pass the heuristic suggestions above to the user code(s) */
  retval = SUNLinSolSetup(arkspils_mem->LS, NULL);

  /* If user set jcur to SUNTRUE, increment npe and save nst value */
  if (*jcurPtr) {
    arkspils_mem->npe++;
    arkspils_mem->nstlpre = ark_mem->ark_nst;
  }

  /* Update jcurPtr flag if we suggested an update */
  if (arkspils_mem->jbad) *jcurPtr = SUNTRUE;

  return(retval);
}

/*---------------------------------------------------------------
 arkSpilsSolve: interfaces between ARKode and the generic 
 SUNLinearSolver object LS, by setting the appropriate tolerance 
 and scaling vectors, calling the solver, and accumulating 
 statistics from the solve for use/reporting by ARKode.
---------------------------------------------------------------*/
int arkSpilsSolve(ARKodeMem ark_mem, N_Vector b, N_Vector ynow, 
                  N_Vector fnow)
{
  realtype bnorm, res_norm;
  ARKSpilsMem arkspils_mem;
  int nli_inc, nps_inc, retval;
  
  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsSolve", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", 
		    "arkSpilsSolve", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  arkspils_mem->deltar = arkspils_mem->eplifac * ark_mem->ark_eRNrm; 
  bnorm = N_VWrmsNorm(b, ark_mem->ark_rwt);
  if (bnorm <= arkspils_mem->deltar) {
    if (ark_mem->ark_mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve 
     interface routines */
  arkspils_mem->ycur = ynow;
  arkspils_mem->fcur = fnow;

  /* Set input tolerance and initial guess x = 0 to LS */  
  arkspils_mem->delta = arkspils_mem->deltar * arkspils_mem->sqrtN;
  N_VConst(ZERO, arkspils_mem->x);

  /* Set scaling vectors for LS to use */
  retval = SUNLinSolSetScalingVectors(arkspils_mem->LS,
                                      ark_mem->ark_ewt,
                                      ark_mem->ark_rwt);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKSPILS_SUNLS_FAIL, "ARKSPILS", "arkSpilsSolve", 
                    "Error in calling SUNLinSolSetScalingVectors");
    return(ARKSPILS_SUNLS_FAIL);
  }

  /* Store previous nps value in nps_inc */
  nps_inc = arkspils_mem->nps;

  /* If a user-provided jtsetup routine is supplied, call that here */
  if (arkspils_mem->jtsetup) {
    retval = arkspils_mem->jtsetup(ark_mem->ark_tn, ynow, fnow, 
                                 arkspils_mem->j_data);
    arkspils_mem->njtsetup++;
    if (retval != 0) {
      arkProcessError(ark_mem, retval, "ARKSPILS", 
                      "arkSpilsSolve", MSGS_JTSETUP_FAILED);
      return(retval);
    }
  }
  
  /* Call solver, and copy x to b */
  retval = SUNLinSolSolve(arkspils_mem->LS, NULL, arkspils_mem->x,
                          b, arkspils_mem->delta);
  N_VScale(ONE, arkspils_mem->x, b);

  /* Retrieve solver statistics */
  res_norm = SUNLinSolResNorm(arkspils_mem->LS);
  nli_inc  = SUNLinSolNumIters(arkspils_mem->LS);
  nps_inc  = arkspils_mem->nps - nps_inc;
  
  /* Increment counters nli and ncfl */
  arkspils_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) arkspils_mem->ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "      kry  %"RSYM"  %"RSYM"  %i  %i\n", 
            bnorm, res_norm, nli_inc, nps_inc);
  
  /* Interpret solver return value  */
  arkspils_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
    /* allow reduction but not solution on first Newton iteration, 
       otherwise return with a recoverable failure */
    if (ark_mem->ark_mnewt == 0) return(0);
    else                         return(1);
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
    arkProcessError(ark_mem, SUNLS_PACKAGE_FAIL_UNREC, "ARKSPILS", 
                    "arkSpilsSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_ATIMES_FAIL_UNREC, "ARKSPILS", 
                    "arkSpilsSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PSOLVE_FAIL_UNREC, "ARKSPILS", 
                    "arkSpilsSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }
  
  return(0); 
}


/*---------------------------------------------------------------
 arkSpilsFree frees memory associates with the ARKSpils system
 solver interface.
---------------------------------------------------------------*/
int arkSpilsFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL)  return (ARKSPILS_SUCCESS);
  if (ark_mem->ark_lmem == NULL)  return(ARKSPILS_SUCCESS);
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Free N_Vector memory */
  if (arkspils_mem->ytemp) {
    N_VDestroy(arkspils_mem->ytemp);
    arkspils_mem->ytemp = NULL;
  }
  if (arkspils_mem->x) {
    N_VDestroy(arkspils_mem->x);
    arkspils_mem->x = NULL;
  }

  /* Nullify other N_Vector pointers */
  arkspils_mem->ycur = NULL;
  arkspils_mem->fcur = NULL;

  /* Free preconditioner memory (if applicable) */
  if (arkspils_mem->pfree)  arkspils_mem->pfree(ark_mem);
  
  /* free ARKSpils interface structure */
  free(ark_mem->ark_lmem);
  
  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSpilsMassInitialize performs remaining initializations specific
 to the iterative linear solver interface (and solver itself)
---------------------------------------------------------------*/
int arkSpilsMassInitialize(ARKodeMem ark_mem)
{
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsMassInitialize", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "arkSpilsMassInitialize", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;
  
  arkSpilsInitializeMassCounters(arkspils_mem);

  /* ensure that mtimes routine and mass matrix solver exist */
  if (arkspils_mem->mtimes == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "arkSpilsMassInitialize",
                    "SpilsMass solver cannot run without user-provided mass matrix-vector product routine");
    arkspils_mem->last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }
  if (arkspils_mem->LS == NULL) {
    arkProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", 
                    "arkSpilsMassInitialize",
                    "SpilsMass solver cannot run without SUNLinearSolver object");
    arkspils_mem->last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }

  /* if psetup is not present, then arkSpilsMassSetup does not need 
     to be called, so set the msetup function to NULL */
  if (arkspils_mem->pset == NULL)  ark_mem->ark_msetup = NULL;

  /* Call LS initialize routine */
  arkspils_mem->last_flag = SUNLinSolInitialize(arkspils_mem->LS);
  return(arkspils_mem->last_flag);
}


/*---------------------------------------------------------------
 arkSpilsMassSetup calls the LS 'setup' routine.
---------------------------------------------------------------*/
int arkSpilsMassSetup(ARKodeMem ark_mem, N_Vector vtemp1,
                      N_Vector vtemp2, N_Vector vtemp3)
{
  int retval;
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsMassSetup", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "arkSpilsMassSetup", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Call LS setup routine */
  retval = SUNLinSolSetup(arkspils_mem->LS, NULL);
  return(retval);
}


/*---------------------------------------------------------------
 arkSpilsMassMult just calls ARKSpilsMTimes routine.
---------------------------------------------------------------*/
int arkSpilsMassMult(ARKodeMem ark_mem, N_Vector v, N_Vector Mv)
{
  int retval;
  /* Return immediately if ark_mem is  NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsMassMult", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  retval = ARKSpilsMTimes(ark_mem, v, Mv);
  return(retval);
}


/*---------------------------------------------------------------
 arkSpilsMassSolve: interfaces between ARKode and the generic 
 SUNLinearSolver object LS, by setting the appropriate tolerance 
 and scaling vectors, calling the solver, and accumulating 
 statistics from the solve for use/reporting by ARKode.
---------------------------------------------------------------*/
int arkSpilsMassSolve(ARKodeMem ark_mem, N_Vector b)
{
  realtype res_norm;
  ARKSpilsMassMem arkspils_mem;
  int nli_inc, nps_inc, retval;
  
  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", 
		    "arkSpilsMassSolve", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSPILS_MASSMEM_NULL, "ARKSPILS", 
		    "arkSpilsMassSolve", MSGS_MASSMEM_NULL);
    return(ARKSPILS_MASSMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* Set input tolerance and initial guess x = 0 to LS */  
  arkspils_mem->deltar = arkspils_mem->eplifac * ark_mem->ark_nlscoef; 
  arkspils_mem->delta  = arkspils_mem->deltar * arkspils_mem->sqrtN;
  N_VConst(ZERO, arkspils_mem->x);

  /* Set scaling vectors for LS to use */
  retval = SUNLinSolSetScalingVectors(arkspils_mem->LS,
                                      ark_mem->ark_ewt,
                                      ark_mem->ark_rwt);
  if (retval != SUNLS_SUCCESS) {
    arkspils_mem->last_flag = retval;
    return(retval);
  }

  /* Store previous nps value in nps_inc */
  nps_inc = arkspils_mem->nps;

  /* If a user-provided mtsetup routine is supplied, the mass 
     matrix is time-dependent, call that here */
  if (!arkspils_mem->time_dependent && arkspils_mem->nmtsetup) {
    retval = arkspils_mem->mtsetup(ark_mem->ark_tn, 
                                   arkspils_mem->mt_data);
    arkspils_mem->nmtsetup++;
    if (retval != 0) {
      arkProcessError(ark_mem, retval, "ARKSPILS", 
                      "arkSpilsMassSolve", MSGS_MTSETUP_FAILED);
      return(retval);
    }
  }
  
  /* Call solver, and copy x to b */
  retval = SUNLinSolSolve(arkspils_mem->LS, NULL, arkspils_mem->x,
                          b, arkspils_mem->delta);
  if (retval != SUNLS_SUCCESS) {
    arkspils_mem->last_flag = retval;
    return(retval);
  }
  N_VScale(ONE, arkspils_mem->x, b);

  /* Retrieve solver statistics */
  res_norm = SUNLinSolResNorm(arkspils_mem->LS);
  nli_inc  = SUNLinSolNumIters(arkspils_mem->LS);
  nps_inc  = arkspils_mem->nps - nps_inc;
  
  /* Increment counters nli and ncfl */
  arkspils_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) arkspils_mem->ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->ark_report) 
    fprintf(ark_mem->ark_diagfp, "      mass  %"RSYM"  %i  %i\n", 
	    res_norm, nli_inc, nps_inc);
  
  /* Interpret solver return value  */
  arkspils_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
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
    arkProcessError(ark_mem, SUNLS_PACKAGE_FAIL_UNREC, "ARKSPILS", 
		    "arkSpilsMassSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_ATIMES_FAIL_UNREC, "ARKSPILS", 
		    "arkSpilsMassSolve", MSGS_MTIMES_FAILED);    
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PSOLVE_FAIL_UNREC, "ARKSPILS", 
		    "arkSpilsMassSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }
  
  return(0); 
}


/*---------------------------------------------------------------
 arkSpilsMassFree frees memory associates with the ARKSpils mass
 matrix solver interface.
---------------------------------------------------------------*/
int arkSpilsMassFree(ARKodeMem ark_mem)
{
  ARKSpilsMassMem arkspils_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL)  return (ARKSPILS_SUCCESS);
  if (ark_mem->ark_mass_mem == NULL)  return(ARKSPILS_SUCCESS);
  arkspils_mem = (ARKSpilsMassMem) ark_mem->ark_mass_mem;

  /* detach ARKSpils interface routines from LS object (ignore return values) */
  SUNLinSolSetATimes(arkspils_mem->LS, NULL, NULL);
  SUNLinSolSetPreconditioner(arkspils_mem->LS, NULL, NULL, NULL);

  /* Free N_Vector memory */
  if (arkspils_mem->x) {
    N_VDestroy(arkspils_mem->x);
    arkspils_mem->x = NULL;
  }

  /* Nullify other N_Vector pointers */
  arkspils_mem->ycur = NULL;

  /* free ARKSpils interface structure */
  free(ark_mem->ark_mass_mem);
  
  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 arkSpilsInitializeCounters and arkSpilsInitializeMassCounters:

 These routines reset all counters from an ARKSpilsMem or 
 ARKSpilsMassMem structure.
---------------------------------------------------------------*/
int arkSpilsInitializeCounters(ARKSpilsMem arkspils_mem)
{
  arkspils_mem->npe      = 0;
  arkspils_mem->nli      = 0;
  arkspils_mem->nps      = 0;
  arkspils_mem->ncfl     = 0;
  arkspils_mem->nstlpre  = 0;
  arkspils_mem->njtsetup = 0;
  arkspils_mem->njtimes  = 0;
  arkspils_mem->nfes     = 0;
  return(0); 
}

int arkSpilsInitializeMassCounters(ARKSpilsMassMem arkspils_mem)
{
  arkspils_mem->npe      = 0;
  arkspils_mem->nli      = 0;
  arkspils_mem->nps      = 0;
  arkspils_mem->ncfl     = 0;
  arkspils_mem->nmtsetup = 0;
  arkspils_mem->nmtimes  = 0;
  return(0); 
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
