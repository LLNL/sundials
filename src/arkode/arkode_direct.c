/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds, Ashley L. Crawford @ SMU
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
 * This is the implementation file for the ARKDLS linear solver 
 * interface
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_direct_impl.h"
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

/*===============================================================
 FUNCTION SPECIFIC CONSTANTS
===============================================================*/

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/*===============================================================
  ARKDLS Exported functions -- Required
===============================================================*/

/*---------------------------------------------------------------
 ARKDlsSetLinearSolver specifies the direct linear solver.
---------------------------------------------------------------*/
int ARKDlsSetLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                          SUNMatrix A)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if any input is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
                    "ARKDlsSetLinearSolver", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if ( (LS == NULL)  || (A == NULL) ) {
    arkProcessError(NULL, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "ARKDlsSetLinearSolver",
                    "Both LS and A must be non-NULL");
    return(ARKDLS_ILL_INPUT);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if solver and vector are compatible with DLS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_DIRECT) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "ARKDlsSetLinearSolver", 
                    "Non-direct LS supplied to ARKDls interface");
    return(ARKDLS_ILL_INPUT);
  }
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL ||
      ark_mem->ark_tempv->ops->nvsetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "ARKDlsSetLinearSolver", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  /* free any existing system solver attached to ARKode */
  if (ark_mem->ark_lfree)  ark_mem->ark_lfree(ark_mem);

  /* Set four main system linear solver function fields in ark_mem */
  ark_mem->ark_linit  = arkDlsInitialize;
  ark_mem->ark_lsetup = arkDlsSetup;
  ark_mem->ark_lsolve = arkDlsSolve;
  ark_mem->ark_lfree  = arkDlsFree;
  ark_mem->ark_lsolve_type = 1;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetLinearSolver", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer */
  arkdls_mem->LS = LS;
  
  /* Initialize Jacobian-related data */
  arkdls_mem->jacDQ = SUNTRUE;
  arkdls_mem->jac = arkDlsDQJac;
  arkdls_mem->J_data = ark_mem;
  arkdls_mem->last_flag = ARKDLS_SUCCESS;

  /* Initialize counters */
  arkDlsInitializeCounters(arkdls_mem);

  /* Store pointer to A and create saved_J */
  arkdls_mem->A = A;
  arkdls_mem->savedJ = SUNMatClone(A);
  if (arkdls_mem->savedJ == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetLinearSolver", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Allocate memory for x */
  arkdls_mem->x = N_VClone(ark_mem->ark_tempv);
  if (arkdls_mem->x == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetLinearSolver", MSGD_MEM_FAIL);
    SUNMatDestroy(arkdls_mem->savedJ);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsSetMassLinearSolver specifies the direct mass-matrix
 linear solver and user-supplied routine to fill the mass matrix
---------------------------------------------------------------*/
int ARKDlsSetMassLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                              SUNMatrix M, booleantype time_dep)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if any input is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetMassLinearSolver", 
                    MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if ( (LS == NULL) || (M == NULL) ) {
    arkProcessError(NULL, ARKDLS_ILL_INPUT, "ARKDLS", 
		    "ARKDlsSetMassLinearSolver",
                    "Both LS and M must be non-NULL");
    return(ARKDLS_ILL_INPUT);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if solver and vector are compatible with DLS */
  if (SUNLinSolGetType(LS) != SUNLINEARSOLVER_DIRECT) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "ARKDlsSetMassLinearSolver", 
                    "Non-direct LS supplied to ARKDls interface");
    return(ARKDLS_ILL_INPUT);
  }
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL ||
      ark_mem->ark_tempv->ops->nvsetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "ARKDlsSetMassLinearSolver", 
                    MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  /* free any existing system solver attached to ARKode */
  if (ark_mem->ark_mfree)  ark_mem->ark_mfree(ark_mem);

  /* Set four main system linear solver function fields in ark_mem */
  ark_mem->ark_minit  = arkDlsMassInitialize;
  ark_mem->ark_msetup = arkDlsMassSetup;
  ark_mem->ark_mmult  = arkDlsMassMult;
  ark_mem->ark_msolve = arkDlsMassSolve;
  ark_mem->ark_mfree  = arkDlsMassFree;
  ark_mem->ark_msolve_type = 1;
  
  /* notify arkode of non-identity mass matrix */
  ark_mem->ark_mass_matrix = SUNTRUE;

  /* Get memory for ARKDlsMassMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMassMem) malloc(sizeof(struct ARKDlsMassMemRec));
  if (arkdls_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetMassLinearSolver", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* set SUNLinearSolver pointer; flag indicating time-dependence */
  arkdls_mem->LS = LS;
  arkdls_mem->time_dependent = time_dep;

  /* Initialize mass-matrix-related data */
  arkdls_mem->mass = NULL;
  arkdls_mem->last_flag = ARKDLS_SUCCESS;

  /* Initialize counters */
  arkDlsInitializeMassCounters(arkdls_mem);

  /* Store pointer to M and create M_lu */
  arkdls_mem->M = M;
  arkdls_mem->M_lu = SUNMatClone(M);
  if (arkdls_mem->M_lu == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetMassLinearSolver", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Allocate memory for x */
  arkdls_mem->x = N_VClone(ark_mem->ark_tempv);
  if (arkdls_mem->x == NULL) {
    arkProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDLS", 
                    "ARKDlsSetMassLinearSolver", MSGD_MEM_FAIL);
    SUNMatDestroy(arkdls_mem->M_lu);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_mass_mem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*===============================================================
  ARKDLS Exported functions -- Optional input/output
===============================================================*/

/*---------------------------------------------------------------
 ARKDlsSetJacFn specifies the Jacobian function.
---------------------------------------------------------------*/
int ARKDlsSetJacFn(void *arkode_mem, ARKDlsJacFn jac)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetJacFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsSetJacFn", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  if (jac != NULL) {
    arkdls_mem->jacDQ  = SUNFALSE;
    arkdls_mem->jac    = jac;
    arkdls_mem->J_data = ark_mem->ark_user_data;
  } else {
    arkdls_mem->jacDQ  = SUNTRUE;
    arkdls_mem->jac    = arkDlsDQJac;
    arkdls_mem->J_data = ark_mem;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsSetMassFn specifies the mass matrix function.
---------------------------------------------------------------*/
int ARKDlsSetMassFn(void *arkode_mem, ARKDlsMassFn mass)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  if (mass != NULL) {
    arkdls_mem->mass = mass;
  } else {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
		    "ARKDlsSetDenseMassFn", "MassFn must be non-NULL");
    return(ARKDLS_ILL_INPUT);
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetWorkSpace returns the length of workspace allocated for 
 the ARKDLS linear solver interface.
---------------------------------------------------------------*/
int ARKDlsGetWorkSpace(void *arkode_mem, long int *lenrw, 
                       long int *leniw)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetWorkSpace", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetWorkSpace", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* initialize outputs with requirements from ARKDlsMem structure */
  *lenrw = 0;
  *leniw = 4;

  /* add NVector sizes */
  if (arkdls_mem->x->ops->nvspace) {
    N_VSpace(arkdls_mem->x, &lrw1, &liw1);
    *lenrw += lrw1;
    *leniw += liw1;
  }
  
  /* add SUNMatrix size (only account for the one owned by Dls interface) */
  if (arkdls_mem->savedJ->ops->space) {
    flag = SUNMatSpace(arkdls_mem->savedJ, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }

  /* add LS sizes */
  if (arkdls_mem->LS->ops->space) {
    flag = SUNLinSolSpace(arkdls_mem->LS, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumJacEvals returns the number of Jacobian evaluations.
---------------------------------------------------------------*/
int ARKDlsGetNumJacEvals(void *arkode_mem, long int *njevals)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumJacEvals", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *njevals = arkdls_mem->nje;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumRhsEvals returns the number of calls to the ODE function
 needed for the DQ Jacobian approximation.
---------------------------------------------------------------*/
int ARKDlsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumRhsEvals", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *nfevalsLS = arkdls_mem->nfeDQ;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetReturnFlagName returns the name associated with a ARKDLS
 return value.
---------------------------------------------------------------*/
char *ARKDlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKDLS_SUCCESS:
    sprintf(name,"ARKDLS_SUCCESS");
    break;   
  case ARKDLS_MEM_NULL:
    sprintf(name,"ARKDLS_MEM_NULL");
    break;
  case ARKDLS_LMEM_NULL:
    sprintf(name,"ARKDLS_LMEM_NULL");
    break;
  case ARKDLS_ILL_INPUT:
    sprintf(name,"ARKDLS_ILL_INPUT");
    break;
  case ARKDLS_MEM_FAIL:
    sprintf(name,"ARKDLS_MEM_FAIL");
    break;
  case ARKDLS_MASSMEM_NULL:
    sprintf(name,"ARKDLS_MASSMEM_NULL");
    break;
  case ARKDLS_JACFUNC_UNRECVR:
    sprintf(name,"ARKDLS_JACFUNC_UNRECVR");
    break;
  case ARKDLS_JACFUNC_RECVR:
    sprintf(name,"ARKDLS_JACFUNC_RECVR");
    break;
  case ARKDLS_MASSFUNC_UNRECVR:
    sprintf(name,"ARKDLS_MASSFUNC_UNRECVR");
    break;
  case ARKDLS_MASSFUNC_RECVR:
    sprintf(name,"ARKDLS_MASSFUNC_RECVR");
    break;
  case ARKDLS_SUNMAT_FAIL:
    sprintf(name,"ARKDLS_SUNMAT_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
 ARKDlsGetLastFlag returns the last flag set in a ARKDLS function.
---------------------------------------------------------------*/
int ARKDlsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastFlag", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastFlag", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  *flag = arkdls_mem->last_flag;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetMassWorkSpace returns the length of workspace allocated 
 for the ARKDLS mass matrix linear solver.
---------------------------------------------------------------*/
int ARKDlsGetMassWorkSpace(void *arkode_mem, long int *lenrw, 
                           long int *leniw)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;
  sunindextype lrw1, liw1;
  long int lrw, liw;
  int flag;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetMassWorkSpace", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetMassWorkSpace", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* initialize outputs with requirements from ARKDlsMem structure */
  *lenrw = 0;
  *leniw = 6;

  /* add NVector sizes */
  if (arkdls_mem->x->ops->nvspace) {
    N_VSpace(arkdls_mem->x, &lrw1, &liw1);
    *lenrw += lrw1;
    *leniw += liw1;
  }
  
  /* add SUNMatrix size (only account for the one owned by Dls interface) */
  if (arkdls_mem->M_lu->ops->space) {
    flag = SUNMatSpace(arkdls_mem->M_lu, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }

  /* add LS sizes */
  if (arkdls_mem->LS->ops->space) {
    flag = SUNLinSolSpace(arkdls_mem->LS, &lrw, &liw);
    *lenrw += lrw;
    *leniw += liw;
  }
  
  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumMassSetups returns the number of mass matrix solver
 'setup' calls
---------------------------------------------------------------*/
int ARKDlsGetNumMassSetups(void *arkode_mem, long int *nmsetups)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassSetups", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassSetups", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *nmsetups = arkdls_mem->mass_setups;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumMassSolves returns the number of mass matrix solver
 'solve' calls
---------------------------------------------------------------*/
int ARKDlsGetNumMassSolves(void *arkode_mem, long int *nmsolves)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassSolves", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassSolves", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *nmsolves = arkdls_mem->mass_solves;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetNumMassMult returns the number of mass matrix-vector
 products
---------------------------------------------------------------*/
int ARKDlsGetNumMassMult(void *arkode_mem, long int *nmmults)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassMult", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetNumMassMult", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *nmmults = arkdls_mem->mass_mults;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDlsGetLastMassFlag returns the last flag set in a ARKDLS mass
  matrix function.
---------------------------------------------------------------*/
int ARKDlsGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastMassFlag", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Access the ARKDlsMem structure */
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "ARKDlsGetLastMassFlag", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  *flag = arkdls_mem->last_flag;

  return(ARKDLS_SUCCESS);
}


/*===============================================================
  ARKDLS Private functions
===============================================================*/

/*---------------------------------------------------------------
 arkDlsDQJac:

 This routine is a wrapper for the Dense and Band
 implementations of the difference quotient Jacobian 
 approximation routines.
---------------------------------------------------------------*/
int arkDlsDQJac(realtype t, N_Vector y, N_Vector fy, 
                SUNMatrix Jac, void *arkode_mem, N_Vector tmp1, 
                N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  ARKodeMem ark_mem;
  ark_mem = (ARKodeMem) arkode_mem;

  /* verify that Jac is non-NULL */
  if (Jac == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "arkDlsDQJac", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }

  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
    retval = arkDlsDenseDQJac(t, y, fy, Jac, ark_mem, tmp1);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
    retval = arkDlsBandDQJac(t, y, fy, Jac, ark_mem, tmp1, tmp2);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_SPARSE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKDLS", 
		    "arkDlsDQJac", 
                    "arkDlsDQJac not implemented for SUNMATRIX_SPARSE");
    retval = ARK_ILL_INPUT;
  } else {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKDLS", 
		    "arkDlsDQJac", 
                    "unrecognized matrix type for arkDlsDQJac");
    retval = ARK_ILL_INPUT;
  }
  return(retval);
}

/*---------------------------------------------------------------
 arkDlsDenseDQJac:

 This routine generates a dense difference quotient approximation 
 to the Jacobian of f(t,y). It assumes a dense SUNMatrix input
 (stored column-wise, and that elements within each column are
 contiguous). The address of the jth column of J is obtained via 
 the function SUNDenseMatrix_Column() and this pointer is 
 associated with an N_Vector using the 
 N_VGetArrayPointer/N_VSetArrayPointer functions.  Finally, the
 actual computation of the jth column of the Jacobian is done 
 with a call to N_VLinearSum.
---------------------------------------------------------------*/
int arkDlsDenseDQJac(realtype t, N_Vector y, N_Vector fy, 
                     SUNMatrix Jac, ARKodeMem ark_mem, N_Vector tmp1)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  sunindextype j, N;
  int retval = 0;
  ARKDlsMem arkdls_mem;

  /* access DlsMem interface structure */
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* access matrix dimension */
  N = SUNDenseMatrix_Rows(Jac);

  /* Rename work vector for readibility */
  ftemp = tmp1;

  /* Create an empty vector for matrix column calculations */
  jthCol = N_VCloneEmpty(tmp1);

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ark_mem->ark_ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->ark_uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->ark_rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) * ark_mem->ark_uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = SUNMAX(srur*SUNRabs(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = ark_mem->ark_fi(t, y, ftemp, ark_mem->ark_user_data);
    arkdls_mem->nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    /* DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol); */  /* UNNECESSARY?? */
  }

  /* Destroy jthCol vector */
  N_VSetArrayPointer(NULL, jthCol);  /* SHOULDN'T BE NEEDED */
  N_VDestroy(jthCol);

  return(retval);
}


/*---------------------------------------------------------------
 arkDlsBandDQJac:

 This routine generates a banded difference quotient approximation 
 to the Jacobian of f(t,y).  It assumes a band SUNMatrix input 
 (stored column-wise, and that elements within each column are 
 contiguous). This makes it possible to get the address 
 of a column of J via the function SUNBandMatrix_Column() and to
 write a simple for loop to set each of the elements of a column 
 in succession.
---------------------------------------------------------------*/
int arkDlsBandDQJac(realtype t, N_Vector y, N_Vector fy, 
                    SUNMatrix Jac, ARKodeMem ark_mem,
                    N_Vector tmp1, N_Vector tmp2)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  sunindextype group, i, j, width, ngroups, i1, i2;
  int retval = 0;
  sunindextype N, mupper, mlower;
  ARKDlsMem arkdls_mem;

  /* access DlsMem interface structure */
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* access matrix dimensions */
  N = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ark_mem->ark_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->ark_uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->ark_rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->ark_h) * ark_mem->ark_uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ytemp, ftemp, 
                             ark_mem->ark_user_data);
    arkdls_mem->nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(Jac, j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
  
  return(retval);
}


/*---------------------------------------------------------------
 arkDlsInitialize performs remaining initializations specific
 to the direct linear solver interface (and solver itself)
---------------------------------------------------------------*/
int arkDlsInitialize(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;
  ARKDlsMassMem arkdls_mass_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
                    "arkDlsInitialize", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "arkDlsInitialize", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  arkDlsInitializeCounters(arkdls_mem);

  /* Set Jacobian function and data, depending on jacDQ (in case 
     it has changed based on user input) */
  if (arkdls_mem->jacDQ) {
    arkdls_mem->jac    = arkDlsDQJac;
    arkdls_mem->J_data = ark_mem;
  } else {
    arkdls_mem->J_data = ark_mem->ark_user_data;
  }

  /* Ensure that if a mass matrix / solver are used, that system 
     and mass matrix solvers and types match */
  if (ark_mem->ark_mass_matrix) {

    /* access mass matrix solver interface */
    if (ark_mem->ark_mass_mem == NULL) {
      arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
                      "arkDlsInitialize", MSGD_MASSMEM_NULL);
      return(ARKDLS_MASSMEM_NULL);
    }
    arkdls_mass_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

    /* check that ark_msolve_type is compatible */
    if (ark_mem->ark_msolve_type != 1) {
      arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
		      "arkDlsInitialize",
                      "Dls and Spils solvers cannot be combined");
      arkdls_mem->last_flag = ARKDLS_ILL_INPUT;
      return(-1);
    }

    /* check that system and mass matrix types are compatible */
    if (SUNMatGetID(arkdls_mem->A) != SUNMatGetID(arkdls_mass_mem->M)) {
      arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
		      "arkDlsInitialize",
                      "System and mass matrices must have the same type");
      arkdls_mem->last_flag = ARKDLS_ILL_INPUT;
      return(-1);
    }
  }

  /* Call LS initialize routine */
  arkdls_mem->last_flag = SUNLinSolInitialize(arkdls_mem->LS);
  return(arkdls_mem->last_flag);
}


/*---------------------------------------------------------------
 arkDlsSetup determines whether to update a Jacobian matrix (or
 use a stored version), based on heuristics regarding previous 
 convergence issues, the number of time steps since it was last
 updated, etc.; it then creates the system matrix from this, the
 'gamma' factor and the mass/identity matrix, 
    A = M-gamma*J.
 This routine then calls the LS 'setup' routine with A.
---------------------------------------------------------------*/
int arkDlsSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                N_Vector fpred, booleantype *jcurPtr, 
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  ARKDlsMem arkdls_mem;
  ARKDlsMassMem arkdls_mass_mem;
  int retval;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
                    "arkDlsSetup", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
                    "arkDlsSetup", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = SUNRabs((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkdls_mem->nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
 
  /* If jok = SUNTRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = SUNFALSE;
    retval = SUNMatCopy(arkdls_mem->savedJ, arkdls_mem->A);
    if (retval) {
      arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                      "arkDlsSetup",  MSGD_MATCOPY_FAILED);
      arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
      return(-1);
    }

  /* If jok = SUNFALSE, call jac routine for new J value */
  } else {
    arkdls_mem->nje++;
    arkdls_mem->nstlj = ark_mem->ark_nst;
    *jcurPtr = SUNTRUE;
    retval = SUNMatZero(arkdls_mem->A);
    if (retval) {
      arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                      "arkDlsSetup",  MSGD_MATZERO_FAILED);
      arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
      return(-1);
    }

    retval = arkdls_mem->jac(ark_mem->ark_tn, ypred, fpred, arkdls_mem->A, 
                             arkdls_mem->J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKDLS", 
                      "arkDlsSetup",  MSGD_JACFUNC_FAILED);
      arkdls_mem->last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    retval = SUNMatCopy(arkdls_mem->A, arkdls_mem->savedJ);
    if (retval) {
      arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                      "arkDlsSetup",  MSGD_MATCOPY_FAILED);
      arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
      return(-1);
    }

  }
  
  /* Scale and add mass matrix to get A = M-gamma*J*/
  if (ark_mem->ark_mass_matrix) {

    /* Access mass matrix solver interface */
    if (ark_mem->ark_mass_mem == NULL) {
      arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
                      "arkDlsSetup", MSGD_MASSMEM_NULL);
      return(ARKDLS_MASSMEM_NULL);
    }
    arkdls_mass_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

    /* Compute mass matrix */
    retval = arkDlsMassSetup(ark_mem, vtemp1, vtemp2, vtemp3);

    /* Perform linear combination A = M-gamma*A */
    retval = SUNMatScaleAdd(-ark_mem->ark_gamma, arkdls_mem->A, 
                            arkdls_mass_mem->M);
    if (retval) {
      arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                      "arkDlsSetup",  MSGD_MATSCALEADD_FAILED);
      arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
      return(-1);
    }

  } else {
    retval = SUNMatScaleAddI(-ark_mem->ark_gamma, arkdls_mem->A);
    if (retval) {
      arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                      "arkDlsSetup",  MSGD_MATSCALEADDI_FAILED);
      arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
      return(-1);
    }
  }

  /* Call generic linear solver 'setup' with this system matrix, and
     return success/failure flag */
  arkdls_mem->last_flag = SUNLinSolSetup(arkdls_mem->LS, arkdls_mem->A);
  return(arkdls_mem->last_flag);
}


/*---------------------------------------------------------------
 arkDlsSolve interfaces between ARKode and the generic 
 SUNLinearSolver object LS, by calling the solver and scaling 
 the solution appropriately when gamrat != 1.
---------------------------------------------------------------*/
int arkDlsSolve(ARKodeMem ark_mem, N_Vector b, N_Vector ycur, N_Vector fcur)
{
  int retval;
  ARKDlsMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "arkDlsSolve", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKDLS_LMEM_NULL, "ARKDLS", 
		    "arkDlsSolve", MSGD_LMEM_NULL);
    return(ARKDLS_LMEM_NULL);
  }
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* call the generic linear system solver, and copy b to x */
  retval = SUNLinSolSolve(arkdls_mem->LS, arkdls_mem->A, arkdls_mem->x, b, ZERO);
  N_VScale(ONE, arkdls_mem->x, b);
  
  /* scale the correction to account for change in gamma */
  if (ark_mem->ark_gamrat != ONE) 
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);
  
  /* store solver return value and return */
  arkdls_mem->last_flag = retval;
  return(retval);
}


/*---------------------------------------------------------------
 arkDlsFree frees memory associates with the ARKDls system
 solver interface.
---------------------------------------------------------------*/
int arkDlsFree(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_lmem are NULL */
  if (ark_mem == NULL)  return (ARKDLS_SUCCESS);
  if (ark_mem->ark_lmem == NULL)  return(ARKDLS_SUCCESS);
  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Free x vector */
  if (arkdls_mem->x) {
    N_VDestroy(arkdls_mem->x);
    arkdls_mem->x = NULL;
  }

  /* Free savedJ memory */
  if (arkdls_mem->savedJ) {
    SUNMatDestroy(arkdls_mem->savedJ);
    arkdls_mem->savedJ = NULL;
  }

  /* Nullify other SUNMatrix pointer */
  arkdls_mem->A = NULL;

  /* free ARKDls interface structure */
  free(ark_mem->ark_lmem);
  
  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkDlsMassInitialize performs remaining initializations specific
 to the direct mass linear solver interface (and solver itself)
---------------------------------------------------------------*/
int arkDlsMassInitialize(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "arkDlsMassInitialize", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "arkDlsMassInitialize", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;
  
  arkDlsInitializeMassCounters(arkdls_mem);

  /* Ensure that mass matrix routine, matrix and solver exist */
  if (arkdls_mem->mass == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "arkDlsMassInitialize",
                    "DlsMass solver cannot run without user-provided mass-matrix routine");
    arkdls_mem->last_flag = ARKDLS_ILL_INPUT;
    return(-1);
  }
  if (arkdls_mem->M == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "arkDlsMassInitialize",
                    "DlsMass solver cannot run without SUNMatrix object");
    arkdls_mem->last_flag = ARKDLS_ILL_INPUT;
    return(-1);
  }
  if (arkdls_mem->LS == NULL) {
    arkProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDLS", 
                    "arkDlsMassInitialize",
                    "DlsMass solver cannot run without SUNLinearSolver object");
    arkdls_mem->last_flag = ARKDLS_ILL_INPUT;
    return(-1);
  }

  /* Call LS initialize routine */
  arkdls_mem->last_flag = SUNLinSolInitialize(arkdls_mem->LS);
  return(arkdls_mem->last_flag);
}


/*---------------------------------------------------------------
 arkDlsMassSetup updates the system mass matrix and calls the LS
 'setup' routine with M.
---------------------------------------------------------------*/
int arkDlsMassSetup(ARKodeMem ark_mem, N_Vector vtemp1,
                    N_Vector vtemp2, N_Vector vtemp3)
{
  ARKDlsMassMem arkdls_mem;
  int retval;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "arkDlsMassSetup", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "arkDlsMassSetup", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* If mass matrix is not time dependent, and if it has been set up 
   previously, just reuse existing M and M_lu */
  if (!arkdls_mem->time_dependent && arkdls_mem->mass_setups) {
    arkdls_mem->last_flag = ARKDLS_SUCCESS;
    return(arkdls_mem->last_flag);
  }
  
  /* Compute mass matrix */
  retval = SUNMatZero(arkdls_mem->M);
  if (retval) {
    arkProcessError(ark_mem, ARKDLS_SUNMAT_FAIL, "ARKDLS", 
                    "arkDlsMassSetup",  MSGD_MATZERO_FAILED);
    arkdls_mem->last_flag = ARKDLS_SUNMAT_FAIL;
    return(-1);
  }

  retval = arkdls_mem->mass(ark_mem->ark_tn, 
                            arkdls_mem->M, 
                            ark_mem->ark_user_data, 
                            vtemp1, vtemp2, vtemp3);
  if (retval < 0) {
    arkProcessError(ark_mem, ARKDLS_MASSFUNC_UNRECVR, "ARKDLS", 
                    "arkDlsMassSetup",  MSGD_MASSFUNC_FAILED);
    arkdls_mem->last_flag = ARKDLS_MASSFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    arkdls_mem->last_flag = ARKDLS_MASSFUNC_RECVR;
    return(1);
  }

  /* Copy M into M_lu for factorization */
  retval = SUNMatCopy(arkdls_mem->M, arkdls_mem->M_lu);
  
  /* Call generic linear solver 'setup' with this system matrix, and
     return success/failure flag */
  arkdls_mem->last_flag = SUNLinSolSetup(arkdls_mem->LS, arkdls_mem->M_lu);
  arkdls_mem->mass_setups++;
  return(arkdls_mem->last_flag);
}


/*---------------------------------------------------------------
 arkDlsMassMult performs a mass-matrix-vector product using 
 the current mass matrix stored in the DLS solver object.  This
 is only used when updating the residual weight vector.
---------------------------------------------------------------*/
int arkDlsMassMult(ARKodeMem ark_mem, N_Vector v, N_Vector Mv)
{
  int retval;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "arkDlsMassMult", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "arkDlsMassMult", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* Use SUNMatrix routine to perform product */
  retval = SUNMatMatvec(arkdls_mem->M, v, Mv);
  arkdls_mem->mass_mults++;
  return(retval);
}


/*---------------------------------------------------------------
 arkDlsMassSolve interfaces between ARKode and the generic 
 SUNLinearSolver object LS.
---------------------------------------------------------------*/
int arkDlsMassSolve(ARKodeMem ark_mem, N_Vector b)
{
  int retval;
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARKDLS_MEM_NULL, "ARKDLS", 
		    "arkDlsMassSolve", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKDLS_MASSMEM_NULL, "ARKDLS", 
		    "arkDlsMassSolve", MSGD_MASSMEM_NULL);
    return(ARKDLS_MASSMEM_NULL);
  }
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* call the generic linear system solver, and copy b to x */
  retval = SUNLinSolSolve(arkdls_mem->LS, arkdls_mem->M_lu, arkdls_mem->x, b, ZERO);
  N_VScale(ONE, arkdls_mem->x, b);
  arkdls_mem->mass_solves++;
  
  /* store solver return value and return */
  arkdls_mem->last_flag = retval;
  return(retval);
}


/*---------------------------------------------------------------
 arkDlsMassFree frees memory associates with the ARKDls mass
 matrix solver interface.
---------------------------------------------------------------*/
int arkDlsMassFree(ARKodeMem ark_mem)
{
  ARKDlsMassMem arkdls_mem;

  /* Return immediately if ark_mem or ark_mem->ark_mass_mem are NULL */
  if (ark_mem == NULL)  return (ARKDLS_SUCCESS);
  if (ark_mem->ark_mass_mem == NULL)  return(ARKDLS_SUCCESS);
  arkdls_mem = (ARKDlsMassMem) ark_mem->ark_mass_mem;

  /* Free x vector */
  if (arkdls_mem->x) {
    N_VDestroy(arkdls_mem->x);
    arkdls_mem->x = NULL;
  }

  /* Free M_lu memory */
  if (arkdls_mem->M_lu) {
    SUNMatDestroy(arkdls_mem->M_lu);
    arkdls_mem->M_lu = NULL;
  }

  /* Nullify other SUNMatrix pointer */
  arkdls_mem->M = NULL;

  /* free ARKDls interface structure */
  free(ark_mem->ark_mass_mem);
  
  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkDlsInitializeCounters and arkDlsInitializeMassCounters reset
 the counters inside the ARKDlsMem and ARKDlsMassMem objects.
---------------------------------------------------------------*/
int arkDlsInitializeCounters(ARKDlsMem arkdls_mem)
{
  arkdls_mem->nje   = 0;
  arkdls_mem->nfeDQ = 0;
  arkdls_mem->nstlj = 0;
  return(0);
}


int arkDlsInitializeMassCounters(ARKDlsMassMem arkdls_mem)
{
  arkdls_mem->mass_setups = 0;
  arkdls_mem->mass_solves = 0;
  arkdls_mem->mass_mults  = 0;
  return(0);
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
