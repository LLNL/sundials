/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * This is the implementation file for an IDASLS linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_sparse_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * IDASlsSetSparseJacFn specifies the sparse Jacobian function.
 */
int IDASlsSetSparseJacFn(void *ida_mem, IDASlsSparseJacFn jac)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDASlsSetSparseJacFn", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASLS", 
		    "IDASlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  idasls_mem->s_jaceval = jac;

  return(IDASLS_SUCCESS);
}

/*
 * IDASlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int IDASlsGetNumJacEvals(void *ida_mem, long int *njevals)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDASlsGetNumJacEvals", MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASLS", 
		    "IDASlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  *njevals = idasls_mem->s_nje;

  return(IDASLS_SUCCESS);
}

/*
 * IDASlsGetReturnFlagName returns the name associated with a IDASLS
 * return value.
 */
char *IDASlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case IDASLS_SUCCESS:
    sprintf(name,"IDASLS_SUCCESS");
    break;   
  case IDASLS_MEM_NULL:
    sprintf(name,"IDASLS_MEM_NULL");
    break;
  case IDASLS_LMEM_NULL:
    sprintf(name,"IDASLS_LMEM_NULL");
    break;
  case IDASLS_ILL_INPUT:
    sprintf(name,"IDASLS_ILL_INPUT");
    break;
  case IDASLS_MEM_FAIL:
    sprintf(name,"IDASLS_MEM_FAIL");
    break;
  case IDASLS_JAC_NOSET:
    sprintf(name,"IDASLS_JAC_NOSET");
    break;
  case IDASLS_JACFUNC_UNRECVR:
    sprintf(name,"IDASLS_JACFUNC_UNRECVR");
    break;
  case IDASLS_JACFUNC_RECVR:
    sprintf(name,"IDASLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * IDASlsGetLastFlag returns the last flag set in a IDASLS function.
 */
int IDASlsGetLastFlag(void *ida_mem, long int *flag)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDASlsGetLastFlag", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_lmem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_LMEM_NULL, "IDASLS", 
		    "IDASlsGetLastFlag", MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;

  *flag = idasls_mem->s_last_flag;

  return(IDASLS_SUCCESS);
}

