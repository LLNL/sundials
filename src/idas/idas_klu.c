/*
 * -----------------------------------------------------------------
 * $Revision: 4951 $
 * $Date: 2016-09-22 10:21:00 -0700 (Thu, 22 Sep 2016) $
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
 * This is the implementation file for the IDASKLU linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "idas/idas_klu.h"
#include "idas_impl.h"
#include "idas_sparse_impl.h"
#include "sundials/sundials_klu_impl.h"
#include "sundials/sundials_math.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)
#define TWOTHIRDS    RCONST(0.6666666666666667)

/* IDASKLU linit, lsetup, lsolve, and lfree routines */
 
static int IDAKLUInit(IDAMem IDA_mem);

static int IDAKLUSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
		       N_Vector rrp, N_Vector tmp1,
		       N_Vector tmp2, N_Vector tmp3);

static int IDAKLUSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
			     N_Vector ycur, N_Vector ypcur, N_Vector rrcur);

static int IDAKLUFree(IDAMem IDA_mem);

/* IDAKLU lfreeB function */

static int IDAKLUFreeB(IDABMem IDAB_mem);


/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDAKLU
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDA / KLU linear solver module.  
 * IDAKLU first calls the existing lfree routine if this is not NULL.
 * Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDAKLUInit, IDAKLUSetup,
 * IDAKLUSolve, NULL, and IDAKLUFree, respectively.
 * It allocates memory for a structure of type IDAkluMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem) to TRUE.
 * Finally, it allocates memory for KLU.
 * The return value is IDASLS_SUCCESS = 0, IDASLS_LMEM_FAIL = -1,
 * or IDASLS_ILL_INPUT = -2.
 *
 * NOTE: The KLU linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDAKLU will first 
 *       test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int IDAKLU(void *ida_mem, int n, int nnz, int sparsetype)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDAKLU", 
                    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if (IDA_mem->ida_tempv1->ops->nvgetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDAKLU", 
                    MSGSP_BAD_NVECTOR);
    return(IDASLS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) flag = IDA_mem->ida_lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  IDA_mem->ida_linit  = IDAKLUInit;
  IDA_mem->ida_lsetup = IDAKLUSetup;
  IDA_mem->ida_lsolve = IDAKLUSolve;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = IDAKLUFree;

  /* Get memory for IDASlsMemRec. */
  idasls_mem = (IDASlsMem) malloc(sizeof(struct IDASlsMemRec));
  if (idasls_mem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
                    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData)malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
                    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  IDA_mem->ida_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  idasls_mem->s_jaceval = NULL;
  idasls_mem->s_jacdata = NULL;
  idasls_mem->sparsetype = sparsetype;

  /* Allocate memory for the sparse Jacobian */
  idasls_mem->s_JacMat = SparseNewMat(n, n, nnz, sparsetype);
  if (idasls_mem->s_JacMat == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
                    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Initialize KLU structures */
  switch (sparsetype) {
    case CSC_MAT:
      klu_data->sun_klu_solve = &klu_solve;
      break;
    case CSR_MAT:
      klu_data->sun_klu_solve = &klu_tsolve;
      break;
    default:
      SparseDestroyMat(idasls_mem->s_JacMat);
      free(klu_data);
      free(idasls_mem);
      return(IDASLS_ILL_INPUT);
  }
  klu_data->s_Symbolic = NULL;
  klu_data->s_Numeric = NULL;

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLU", 
                    MSGSP_PACKAGE_FAIL);
    return(IDASLS_PACKAGE_FAIL);
  }

  /* Set ordering to COLAMD as the idas default use.
     Users can set a different value with IDAKLUSetOrdering,
     and the user-set value is loaded before any call to klu_analyze in
     IDAKLUSetup.  */
  klu_data->s_ordering = 1;
  klu_data->s_Common.ordering = klu_data->s_ordering;

  /* Attach linear solver memory to the integrator memory */
  idasls_mem->s_solver_data = (void *) klu_data;
  IDA_mem->ida_lmem = idasls_mem;

  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  return(IDASLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDAKLUReInit
 * -----------------------------------------------------------------
 * This routine reinitializes memory and flags for a new factorization 
 * (symbolic and numeric) to be conducted at the next solver setup
 * call.  This routine is useful in the cases where the number of nonzeroes 
 * has changed or if the structure of the linear system has changed
 * which would require a new symbolic (and numeric factorization).
 *
 * The reinit_type argumenmt governs the level of reinitialization:
 *
 * reinit_type = 1: The Jacobian matrix will be destroyed and 
 *                  a new one will be allocated based on the nnz
 *                  value passed to this call. New symbolic and
 *                  numeric factorizations will be completed at the next
 *                  solver setup.
 *
 * reinit_type = 2: Only symbolic and numeric factorizations will be 
 *                  completed.  It is assumed that the Jacobian size
 *                  has not exceeded the size of nnz given in the prior
 *                  call to IDAKLU.
 *
 * This routine assumes no other changes to solver use are necessary.
 *
 * The return value is IDASLS_SUCCESS = 0, IDASLS_MEM_NULL = -1, 
 * IDASLS_LMEM_NULL = -2, IDASLS_ILL_INPUT = -3, or IDASLS_MEM_FAIL = -4.
 *
 * -----------------------------------------------------------------
 */

int IDAKLUReInit(void *ida_mem_v, int n, int nnz, int reinit_type)
{
  IDAMem ida_mem;
  IDASlsMem idasls_mem;
  KLUData klu_data;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem_v == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDAKLUReInit", 
                    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  ida_mem = (IDAMem) ida_mem_v;

  /* Return immediately if ida_lmem is NULL. */
  if (ida_mem->ida_lmem == NULL) {
    IDAProcessError(NULL, IDASLS_LMEM_NULL, "IDASLS", "IDAKLUReInit", 
                    MSGSP_LMEM_NULL);
    return(IDASLS_LMEM_NULL);
  }

  idasls_mem = (IDASlsMem) (ida_mem->ida_lmem);
  klu_data = (KLUData) idasls_mem->s_solver_data;

  /* Return if reinit_type is not valid */
  if ((reinit_type != 1) && (reinit_type != 2)) {
    IDAProcessError(NULL, IDASLS_ILL_INPUT, "IDASLS", "IDAKLUReInit", 
                    MSGSP_ILL_INPUT);
    return(IDASLS_ILL_INPUT);
  }

  if (reinit_type == 1) {

    /* Destroy previous Jacobian information */
    if (idasls_mem->s_JacMat) {
      SparseDestroyMat(idasls_mem->s_JacMat);
    }

    /* Allocate memory for the sparse Jacobian */
    idasls_mem->s_JacMat = SparseNewMat(n, n, nnz, idasls_mem->sparsetype);
    if (idasls_mem->s_JacMat == NULL) {
      IDAProcessError(ida_mem, IDASLS_MEM_FAIL, "IDASLS", "IDAKLU", 
                      MSGSP_MEM_FAIL);
      return(IDASLS_MEM_FAIL);
    }
  }

  /* Free the prior factorazation and reset for first factorization */
  if( klu_data->s_Symbolic != NULL)
    klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  if( klu_data->s_Numeric != NULL)
    klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  idasls_mem->s_first_factorize = 1;

  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * IDAKLU interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the IDAKLU
  linear solver module.  
  It returns 0 if successful.
*/

static int IDAKLUInit(IDAMem IDA_mem)
{
  IDASlsMem idasls_mem;

  idasls_mem = (IDASlsMem)IDA_mem->ida_lmem;

  idasls_mem->s_jacdata = IDA_mem->ida_user_data;

  idasls_mem->s_nje = 0;
  /* Force factorization every call to IDAS */
  idasls_mem->s_first_factorize = 1;

  idasls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the IDAKLU linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     IDASLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int IDAKLUSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
		       N_Vector rrp, N_Vector tmp1, N_Vector tmp2,
		       N_Vector tmp3)
{
  int retval;
  realtype tn, cj;
  IDASlsMem idasls_mem;
  IDASlsSparseJacFn jaceval;
  KLUData klu_data;
  SlsMat JacMat;
  void *jacdata;
  
  realtype uround_twothirds;

  uround_twothirds = SUNRpowerR(IDA_mem->ida_uround,TWOTHIRDS);

  idasls_mem = (IDASlsMem) (IDA_mem->ida_lmem);
  tn = IDA_mem->ida_tn; 
  cj = IDA_mem->ida_cj;

  klu_data = (KLUData) idasls_mem->s_solver_data;

  jaceval = idasls_mem->s_jaceval;
  jacdata = idasls_mem->s_jacdata;
  JacMat = idasls_mem->s_JacMat;

  /* Check that Jacobian eval routine is set */
  if (jaceval == NULL) {
    IDAProcessError(IDA_mem, IDASLS_JAC_NOSET, "IDASSLS", "IDAKLUSetup", 
		    MSGSP_JAC_NOSET);
    free(idasls_mem); idasls_mem = NULL;
    return(IDASLS_JAC_NOSET);
  }

  /* Increment nje counter and call Jacobian eval routine. */
  idasls_mem->s_nje++;
  retval = jaceval(tn, cj, yyp, ypp, rrp, JacMat, jacdata, 
		   tmp1, tmp2, tmp3);

  if (retval < 0) {
    IDAProcessError(IDA_mem, IDASLS_JACFUNC_UNRECVR, "IDASSLS", 
		    "IDAKLUSetup", MSGSP_JACFUNC_FAILED);
    idasls_mem->s_last_flag = IDASLS_JACFUNC_UNRECVR;
    return(IDASLS_JACFUNC_UNRECVR);
  }
  if (retval > 0) {
    idasls_mem->s_last_flag = IDASLS_JACFUNC_RECVR;
    return(+1);
  }

  if (idasls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get the symbolic factorization
       ------------------------------------------------------------*/ 
    /* Update the ordering option with any user-updated values from 
       calls to IDAKLUSetOrdering */
    klu_data->s_Common.ordering = klu_data->s_ordering;

    if (klu_data->s_Symbolic != NULL) {
       klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
    }
    klu_data->s_Symbolic = klu_analyze(JacMat->NP, JacMat->indexptrs, 
				       JacMat->indexvals, &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(IDASLS_PACKAGE_FAIL);
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of  the Jacobian.
       ------------------------------------------------------------*/
    if( klu_data->s_Numeric != NULL) {
       klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
    }
    klu_data->s_Numeric = klu_factor(JacMat->indexptrs, JacMat->indexvals, JacMat->data, 
				     klu_data->s_Symbolic, &(klu_data->s_Common));

    if (klu_data->s_Numeric == NULL) {
      IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASLS", "IDAKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(IDASLS_PACKAGE_FAIL);
    }

    idasls_mem->s_first_factorize = 0;
  }
  else {

    retval = klu_refactor(JacMat->indexptrs, JacMat->indexvals, JacMat->data, 
			  klu_data->s_Symbolic, klu_data->s_Numeric,
			  &(klu_data->s_Common));
    if (retval == 0) {
      IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASLS", "idaKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(IDASLS_PACKAGE_FAIL);
    }
    
    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition 
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/
    
    retval = klu_rcond(klu_data->s_Symbolic, klu_data->s_Numeric,
		       &(klu_data->s_Common));
    if (retval == 0) {
      IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASLS", "idaKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(IDASLS_PACKAGE_FAIL);
    }

    if ( (klu_data->s_Common.rcond)  < uround_twothirds ) {
      
      /* Condition number may be getting large.  
	 Compute more accurate estimate */
      retval = klu_condest(JacMat->indexptrs, JacMat->data, 
			   klu_data->s_Symbolic, klu_data->s_Numeric,
			   &(klu_data->s_Common));
      if (retval == 0) {
	IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASLS", "idaKLUSetup", 
			MSGSP_PACKAGE_FAIL);
	return(IDASLS_PACKAGE_FAIL);
      }
      
      if ( (klu_data->s_Common.condest) > 
	   (1.0/uround_twothirds) ) {

	/* More accurate estimate also says condition number is 
	   large, so recompute the numeric factorization */

	klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
	
	klu_data->s_Numeric = klu_factor(JacMat->indexptrs, JacMat->indexvals, 
					 JacMat->data, klu_data->s_Symbolic, 
					 &(klu_data->s_Common));

	if (klu_data->s_Numeric == NULL) {
	  IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASLS", 
			  "IDAKLUSetup", MSGSP_PACKAGE_FAIL);
	  return(IDASLS_PACKAGE_FAIL);
	}
      }
    }
  }

  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the IDAKLU linear
  solver module.  It calls the KLU solve routine, scales the
  solution vector according to cjratio, then returns IDASLS_SUCCESS = 0.
*/

static int IDAKLUSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
		       N_Vector ycur, N_Vector ypcur, N_Vector rrcur)
{
  int flag;
  realtype cjratio;
  IDASlsMem idasls_mem;
  KLUData klu_data;
  SlsMat JacMat;
  realtype *bd;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;
  JacMat = idasls_mem->s_JacMat;
  cjratio = IDA_mem->ida_cjratio;
  klu_data = (KLUData) idasls_mem->s_solver_data;
  bd = N_VGetArrayPointer(b);

  /* Call KLU to solve the linear system */
  flag = klu_data->sun_klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, JacMat->NP, 1, bd, 
                                 &(klu_data->s_Common));
  if (flag == 0) {
    IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLUSolve", 
		    MSGSP_PACKAGE_FAIL);
    return(IDASLS_PACKAGE_FAIL);
  }

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  idasls_mem->s_last_flag = IDASLS_SUCCESS;
  return(IDASLS_SUCCESS);
}

/*
  This routine frees memory specific to the IDAKLU linear solver.
*/

static int IDAKLUFree(IDAMem IDA_mem)
{
  IDASlsMem idasls_mem;
  KLUData klu_data;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;
  klu_data = (KLUData) idasls_mem->s_solver_data;

  if( klu_data->s_Numeric != NULL)
  {
     klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  }
  if( klu_data->s_Symbolic != NULL)
  {
     klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));
  }

  if (idasls_mem->s_JacMat) {
    SparseDestroyMat(idasls_mem->s_JacMat);
    idasls_mem->s_JacMat = NULL;
  }

  free(klu_data); 
  free(IDA_mem->ida_lmem); 

  return(IDASLS_SUCCESS);
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * IDAKLUB is a wrapper around IDAKLU.
 */

int IDAKLUB(void *ida_mem, int which, int n, int nnz, int sparsetype)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  IDASlsMemB idaslsB_mem;
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDAKLUB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDAKLUB",  
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDAKLUB", 
		    MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Alloc memory for IDASlsMemRecB */
  idaslsB_mem = (IDASlsMemB) malloc(sizeof(struct IDASlsMemRecB));
  if (idaslsB_mem == NULL) {
    IDAProcessError(IDAB_mem->IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", 
		    "IDAKLUB", MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  
  }

  /* set matrix type and initialize Jacob function. */
  idaslsB_mem->s_djacB = NULL;

  /* Attach lmemB data and lfreeB function. */
  IDAB_mem->ida_lmem  = idaslsB_mem;
  IDAB_mem->ida_lfree = IDAKLUFreeB;

  /* Call IDAKLU to the IDAS data of the backward problem. */
  ida_memB = (void *)IDAB_mem->IDA_mem;
  flag = IDAKLU(ida_memB, n, nnz, sparsetype);

  if (flag != IDASLS_SUCCESS) {
    free(idaslsB_mem);
    idaslsB_mem = NULL;
  }

  return(flag);
}


/*
 * IDAKLUReInitB is a wrapper around IDAKLUReInit.
 */

int IDAKLUReInitB(void *ida_mem, int which, int n, int nnz, int reinit_type)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  /* IDASlsMemB idaslsB_mem; */
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDAKLUReInitB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDAKLUReInitB",  
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDAKLUReInitB", 
		    MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  ida_memB = (void *)(IDAB_mem->IDA_mem);
  
  flag = IDAKLUReInit(ida_memB, n, nnz, reinit_type);

  return(flag);
}


/*
 * IDAKLUSetOrderingB is a wrapper around IDAKLUSetOrdering.
 */

int IDAKLUSetOrderingB(void *ida_mem, int which, int ordering_choiceB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  /* IDASlsMemB idaslsB_mem; */
  void *ida_memB;
  int flag;
  
  /* Is ida_mem allright? */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDAKLUSetOrderingB", 
		    MSGSP_CAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDASLS_NO_ADJ, "IDASSLS", "IDAKLUSetOrderingB",  
		    MSGSP_NO_ADJ);
    return(IDASLS_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if ( which >= IDAADJ_mem->ia_nbckpbs ) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDAKLUSetOrderingB", 
		    MSGSP_BAD_WHICH);
    return(IDASLS_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL) {
    if( which == IDAB_mem->ida_index ) break;
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  ida_memB = (void *)(IDAB_mem->IDA_mem);
  
  flag = IDAKLUSetOrdering(ida_memB, ordering_choiceB);

  return(flag);
}


/*
 * IDAKLUFreeB frees the linear solver's memory for that backward problem passed 
 * as argument. 
 */

static int IDAKLUFreeB(IDABMem IDAB_mem)
{
  IDASlsMemB idaslsB_mem;

  idaslsB_mem = (IDASlsMemB) IDAB_mem->ida_lmem;

  free(idaslsB_mem);

  return(0);
}


/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * IDAKLUSetOrdering sets the ordering used by KLU for reducing fill.
 * Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
 * The default used in IDA is 1 for COLAMD.
 * -----------------------------------------------------------------
 */

int IDAKLUSetOrdering(void *ida_mem_v, int ordering_choice)
{
  IDAMem ida_mem;
  IDASlsMem idasls_mem;
  KLUData klu_data;

 /* Return immediately if ida_mem_v is NULL */
  if (ida_mem_v == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDAKLUSetOrdering",
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  ida_mem = (IDAMem) ida_mem_v;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) ) {
    IDAProcessError(NULL, IDASLS_ILL_INPUT, "IDASLS", "IDAKLUSetOrdering",
		    MSGSP_ILL_INPUT);
    return(IDASLS_ILL_INPUT);
  }

  idasls_mem = (IDASlsMem) ida_mem->ida_lmem;
  klu_data = (KLUData) idasls_mem->s_solver_data;

  klu_data->s_ordering = ordering_choice;

  return(IDASLS_SUCCESS);
}
