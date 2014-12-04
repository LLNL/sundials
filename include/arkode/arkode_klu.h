/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2014, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 Header file for the ARKKLU linear solver module.
 ---------------------------------------------------------------*/

#ifndef _ARKKLU_H
#define _ARKKLU_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "arkode/arkode_sparse.h"
#include "sundials/sundials_sparse.h"


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
SUNDIALS_EXPORT int ARKKLU(void *arkode_mem, int n, int nnz); 


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
SUNDIALS_EXPORT int ARKMassKLU(void *arkode_mem, int n, int nnz,
			       ARKSlsSparseMassFn smass); 


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
