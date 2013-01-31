/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, Lawrence Livermore National Security
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the IDASuperLUMT linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _IDASSUPERLUMT_H
#define _IDASSUPERLUMT_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <idas/idas_sparse.h>
#include <sundials/sundials_sparse.h>

/*
 * -----------------------------------------------------------------
 * Function : IDASuperLUMT
 * -----------------------------------------------------------------
 * A call to the IDASuperLUMT function links the main integrator      
 * with the IDASuperLUMT linear solver module.                        
 *                                                                
 * ida_mem is the pointer to integrator memory returned by        
 *     IDACreate.             
 *
 * num_threads is the number of threads that SuperLUMT should invoke     
 *                                                                
 * IDASuperLUMT returns:                                              
 *     IDASLU_SUCCESS   = 0  if successful                              
 *     IDASLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     IDASLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The SuperLUMT linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, IDASuperLUMT will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDASuperLUMT(void *ida_mem, int num_threads, 
				   int m, int n, int nnz); 

/*
 * -----------------------------------------------------------------
 * Function: IDASuperLUMTB
 * -----------------------------------------------------------------
 * IDASuperLUMTB links the main IDAS integrator with the IDASuperLUMT
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDASuperLUMTB(void *ida_mem, int num_threads, 
				    int which, int mB, int nB, int nnzB);



  
#ifdef __cplusplus
}
#endif

#endif
