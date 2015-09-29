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
 * Header file for the ARKKLU linear solver module.
 *---------------------------------------------------------------*/

#ifndef _ARKKLU_H
#define _ARKKLU_H

#include "arkode/arkode_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 Function: ARKKLU
-----------------------------------------------------------------
 A call to the ARKKLU function links the main integrator with 
 the ARKKLU linear solver module.
                                                                
 arkode_mem is the pointer to integrator memory returned by        
            ARKodeCreate.             

 n is the size of the linear system (nrows = ncols = n)

 nnz is the maximum number of nonzeros in the sparse matrix
     A = M + gamma*J
                                                                
 The return value of ARKKLU is one of:
     ARKSLS_SUCCESS   if successful
     ARKSLS_MEM_NULL  if the ARKode memory was NULL
     ARKSLS_MEM_FAIL  if there was a memory allocation failure   
     ARKSLS_ILL_INPUT if a required vector operation is missing
                                                                
 NOTE: The KLU linear solver assumes a serial implementation  
       of the NVECTOR package. Therefore, ARKKLU will first
       test for a compatible N_Vector internal representation
       by checking that the function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKKLU(void *arkode_mem, int n, int nnz, int sparsetype); 


/*---------------------------------------------------------------
 Function: ARKMassKLU
-----------------------------------------------------------------
 A call to the ARKMassKLU function links the mass matrix solve
 with the ARKKLU linear solver module.
                                                                
 arkode_mem is the pointer to integrator memory returned by        
            ARKodeCreate.             

 n is the size of the linear system (nrows = ncols = n)

 nnz is the maximum number of nonzeros in the sparse mass matrix
                                                                
 smass is the user-supplied sparse mass matrix setup function.

 The return value of ARKMassKLU is one of:
     ARKSLS_SUCCESS   if successful
     ARKSLS_MEM_NULL  if the ARKode memory was NULL
     ARKSLS_MEM_FAIL  if there was a memory allocation failure   
     ARKSLS_ILL_INPUT if a required vector operation is missing
                                                                
 NOTE: The KLU linear solver assumes a serial implementation  
       of the NVECTOR package. Therefore, ARKKLU will first
       test for a compatible N_Vector internal representation
       by checking that the function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassKLU(void *arkode_mem, int n, int nnz, int sparsetype,
			       ARKSlsSparseMassFn smass); 


/*
 * -----------------------------------------------------------------
 * ARKKLUReInit
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
 *                  call to CVKLU.
 *
 * This routine assumes no other changes to solver use are necessary.
 *
 * The return value is ARKSLS_SUCCESS = 0, ARKSLS_MEM_NULL = -1, 
 * ARKSLS_LMEM_NULL = -2, ARKSLS_ILL_INPUT = -3, or ARKSLS_MEM_FAIL = -4.
 *
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int ARKKLUReInit(void *arkode_mem, int n, int nnz, 
				   int reinit_type);


/*
 * -----------------------------------------------------------------
 * ARKMassKLUReInit
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
 *                  call to CVKLU.
 *
 * This routine assumes no other changes to solver use are necessary.
 *
 * The return value is ARKSLS_SUCCESS = 0, ARKSLS_MEM_NULL = -1, 
 * ARKSLS_MASSMEM_NULL = -8, ARKSLS_ILL_INPUT = -3, or ARKSLS_MEM_FAIL = -4.
 *
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int ARKMassKLUReInit(void *arkode_mem, int n, 
				       int nnz, int reinit_type);

/*================================================================
 Optional Input Specification Functions
=================================================================*/

/*---------------------------------------------------------------
 Function: ARKKLUSetOrdering 
-----------------------------------------------------------------
 This routine sets the ordering used by KLU for reducing fill in 
 the system matrix solve.  Options for ordering_choice are: 
        0 for AMD, 
	1 for COLAMD, and 
	2 for the natural ordering.
 The default used in ARKODE is 1 for COLAMD.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKKLUSetOrdering(void *arkode_mem, 
				      int ordering_choice); 


/*---------------------------------------------------------------
 Function: ARKMassKLUSetOrdering 
-----------------------------------------------------------------
 This routine sets the ordering used by KLU for reducing fill in 
 the mass matrix solve.  Options for ordering_choice are: 
        0 for AMD, 
	1 for COLAMD, and 
	2 for the natural ordering.
 The default used in ARKODE is 1 for COLAMD.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassKLUSetOrdering(void *arkode_mem, 
					  int ordering_choice); 


#ifdef __cplusplus
}
#endif

#endif
