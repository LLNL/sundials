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
 * Implementation file for the generic ARKSLS linear solver module.
 *---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_sparse_impl.h"
#include <sundials/sundials_math.h>

/*===============================================================
 FUNCTION SPECIFIC CONSTANTS (none)
===============================================================*/

/*===============================================================
 EXPORTED FUNCTIONS
===============================================================*/
              
/*---------------------------------------------------------------
 ARKSlsSetSparseJacFn specifies the sparse Jacobian function.
---------------------------------------------------------------*/
int ARKSlsSetSparseJacFn(void *arkode_mem, ARKSlsSparseJacFn jac)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		   "ARKSlsSetSparseJacFn", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSLS_LMEM_NULL, "ARKSLS", 
		    "ARKSlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(ARKSLS_LMEM_NULL);
  }
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;

  arksls_mem->s_Jeval = jac;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSlsSetSparseMassFn specifies the sparse mass matrix function.
---------------------------------------------------------------*/
int ARKSlsSetSparseMassFn(void *arkode_mem, ARKSlsSparseMassFn mass)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		   "ARKSlsSetSparseMassFn", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MASSMEM_NULL, "ARKSLS", 
		    "ARKSlsSetSparseMassFn", MSGSP_MASSMEM_NULL);
    return(ARKSLS_MASSMEM_NULL);
  }
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;

  arksls_mem->s_Meval = mass;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSlsGetNumJacEvals returns the number of Jacobian evaluations.
---------------------------------------------------------------*/
int ARKSlsGetNumJacEvals(void *arkode_mem, long int *njevals)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		   "ARKSlsGetNumJacEvals", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSLS_LMEM_NULL, "ARKSLS", 
		    "ARKSlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(ARKSLS_LMEM_NULL);
  }
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;

  *njevals = arksls_mem->s_nje;

  return(ARKSLS_SUCCESS);
}



/*---------------------------------------------------------------
 ARKSlsGetNumMassEvals returns the number of mass matrix evaluations.
---------------------------------------------------------------*/
int ARKSlsGetNumMassEvals(void *arkode_mem, long int *nmevals)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		   "ARKSlsGetNumMassEvals", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MASSMEM_NULL, "ARKSLS", 
		    "ARKSlsGetNumMassEvals", MSGSP_MASSMEM_NULL);
    return(ARKSLS_MASSMEM_NULL);
  }
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;

  *nmevals = arksls_mem->s_nme;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSlsGetReturnFlagName returns the name associated with a 
 ARKSLS return value.
---------------------------------------------------------------*/
char *ARKSlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKSLS_SUCCESS:
    sprintf(name,"ARKSLS_SUCCESS");
    break;   
  case ARKSLS_MEM_NULL:
    sprintf(name,"ARKSLS_MEM_NULL");
    break;
  case ARKSLS_LMEM_NULL:
    sprintf(name,"ARKSLS_LMEM_NULL");
    break;
  case ARKSLS_ILL_INPUT:
    sprintf(name,"ARKSLS_ILL_INPUT");
    break;
  case ARKSLS_MEM_FAIL:
    sprintf(name,"ARKSLS_MEM_FAIL");
    break;
  case ARKSLS_JAC_NOSET:
    sprintf(name,"ARKSLS_JAC_NOSET");
    break;
  case ARKSLS_PACKAGE_FAIL:
    sprintf(name,"ARKSLS_PACKAGE_FAIL");
    break;
  case ARKSLS_MASSMEM_NULL:
    sprintf(name,"ARKSLS_MASSMEM_NULL");
    break;
  case ARKSLS_JACFUNC_UNRECVR:
    sprintf(name,"ARKSLS_JACFUNC_UNRECVR");
    break;
  case ARKSLS_JACFUNC_RECVR:
    sprintf(name,"ARKSLS_JACFUNC_RECVR");
    break;
  case ARKSLS_MASSFUNC_UNRECVR:
    sprintf(name,"ARKSLS_MASSFUNC_UNRECVR");
    break;
  case ARKSLS_MASSFUNC_RECVR:
    sprintf(name,"ARKSLS_MASSFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
 ARKSlsGetLastFlag returns the last flag set in a ARKSLS function.
---------------------------------------------------------------*/
int ARKSlsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSlsMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKSlsGetLastFlag", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    arkProcessError(ark_mem, ARKSLS_LMEM_NULL, "ARKSLS", 
		    "ARKSlsGetLastFlag", MSGSP_LMEM_NULL);
    return(ARKSLS_LMEM_NULL);
  }
  arksls_mem = (ARKSlsMem) ark_mem->ark_lmem;

  *flag = arksls_mem->s_last_flag;

  return(ARKSLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKSlsGetLastMassFlag returns the last flag set in a ARKSLS mass
 matrix function.
---------------------------------------------------------------*/
int ARKSlsGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSlsMassMem arksls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKSLS_MEM_NULL, "ARKSLS", 
		    "ARKSlsGetLastMassFlag", MSGSP_ARKMEM_NULL);
    return(ARKSLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_mass_mem == NULL) {
    arkProcessError(ark_mem, ARKSLS_MASSMEM_NULL, "ARKSLS", 
		    "ARKSlsGetLastMassFlag", MSGSP_MASSMEM_NULL);
    return(ARKSLS_MASSMEM_NULL);
  }
  arksls_mem = (ARKSlsMassMem) ark_mem->ark_mass_mem;

  *flag = arksls_mem->s_last_flag;

  return(ARKSLS_SUCCESS);
}

/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
