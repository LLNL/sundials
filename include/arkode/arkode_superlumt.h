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
 * Header file for the ARKSUPERLUMT linear solver module.
 *---------------------------------------------------------------*/

#ifndef _ARKSUPERLUMT_H
#define _ARKSUPERLUMT_H

#include "arkode/arkode_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 Function: ARKSuperLUMT
-----------------------------------------------------------------
 A call to the ARKSuperLUMT function links the main integrator
 with the ARKSUPERLUMT linear solver module.

 arkode_mem is the pointer to integrator memory returned by
            ARKodeCreate.

 num_threads is the number of hardware threads to use

 n is the number of rows in the linear system

 nnz is the maximum number of nonzeros in the sparse matrix
     A = M + gamma*J

 The return value of ARKSuperLUMT is one of:
     ARKSLS_SUCCESS   if successful
     ARKSLS_MEM_NULL  if the ARKode memory was NULL
     ARKSLS_MEM_FAIL  if there was a memory allocation failure
     ARKSLS_ILL_INPUT if a required vector operation is missing

 NOTE: The ARKSUPERLUMT linear solver assumes a serial 
       implementation of the NVECTOR package. Therefore, 
       ARKSuperLUMT will first test for a compatible N_Vector
       internal representation by checking that the function
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSuperLUMT(void *arkode_mem, 
				 int num_threads,
				 int n, int nnz); 


/*---------------------------------------------------------------
 Function: ARKMassSuperLUMT
-----------------------------------------------------------------
 A call to the ARKMassSuperLUMT function links the mass matrix
 solve with the ARKSUPERLUMT linear solver module.

 arkode_mem is the pointer to integrator memory returned by
            ARKodeCreate.

 num_threads is the number of hardware threads to use

 n is the number of rows in the sparse mass matrix

 nnz is the maximum number of nonzeros in the sparse mass matrix

 smass is the user-supplied sparse mass matrix setup function.

 The return value of ARKMassSuperLUMT is one of:
     ARKSLS_SUCCESS   if successful
     ARKSLS_MEM_NULL  if the ARKode memory was NULL
     ARKSLS_MEM_FAIL  if there was a memory allocation failure
     ARKSLS_ILL_INPUT if a required vector operation is missing

 NOTE: The ARKSUPERLUMT linear solver assumes a serial 
       implementation of the NVECTOR package. Therefore, 
       ARKMassSuperLUMT will first test for a compatible 
       N_Vector internal representation by checking that the 
       function N_VGetArrayPointer exists.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassSuperLUMT(void *arkode_mem, 
				     int num_threads,
				     int n, int nnz,
				     ARKSlsSparseMassFn smass); 


/*================================================================
 Optional Input Specification Functions
=================================================================*/

/*---------------------------------------------------------------
 Function: ARKSuperLUMTSetOrdering 
-----------------------------------------------------------------
 Sets the ordering used by ARKSuperLUMT for reducing fill in the 
 system matrix solve.  Options for ordering_choice are:
     0 for natural ordering
     1 for minimal degree ordering on A'*A
     2 for minimal degree ordering on A'+A
     3 for AMD ordering for unsymmetric matrices
 The default used in ARKODE is 3 for COLAMD.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKSuperLUMTSetOrdering(void *ark_mem, 
					    int ordering_choice); 


/*---------------------------------------------------------------
 Function: ARKMassSuperLUMTSetOrdering 
-----------------------------------------------------------------
 Sets the ordering used by ARKMassSuperLUMT for reducing fill in 
 the mass matrix solve.  Options for ordering_choice are:
     0 for natural ordering
     1 for minimal degree ordering on A'*A
     2 for minimal degree ordering on A'+A
     3 for AMD ordering for unsymmetric matrices
 The default used in ARKODE is 3 for COLAMD.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKMassSuperLUMTSetOrdering(void *ark_mem, 
						int ordering_choice); 


#ifdef __cplusplus
}
#endif

#endif
