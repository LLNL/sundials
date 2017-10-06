/*
 * -----------------------------------------------------------------
 * $Revision: 4181 $
 * $Date: 2014-07-23 12:43:06 -0700 (Wed, 23 Jul 2014) $
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
 * This is the header file for the KINSuperLUMT linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _KINSUPERLUMT_H
#define _KINSUPERLUMT_H

#include "kinsol/kinsol_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : KINSuperLUMT
 * -----------------------------------------------------------------
 * A call to the KINSuperLUMT function links the main integrator      
 * with the KINSuperLUMT linear solver module.                        
 *                                                                
 * kin_mem is the pointer to integrator memory returned by        
 *     KINCreate.       
 *
 * num_threads is the number of threads that SuperLUMT should invoke         
 *
 *                                                                
 * KINSuperLUMT returns:                                              
 *     KINSLU_SUCCESS   = 0  if successful                              
 *     KINSLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     KINSLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The SuperLUMT linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, KINSuperLUMT will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int KINSuperLUMT(void *kin_mem, int num_threads,
				   int n, int nnz); 


/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * KINSuperLUMTSetOrdering sets the ordering used by SuperLUMT for 
 * reducing fill.
 * Options are: 
 * 0 for natural ordering
 * 1 for minimal degree ordering on A'*A
 * 2 for minimal degree ordering on A'+A
 * 3 for approximate minimal degree ordering for unsymmetric matrices
 * The default used in SUNDIALS is 3 for COLAMD.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int KINSuperLUMTSetOrdering(void *kin_mem, 
					      int ordering_choice); 



#ifdef __cplusplus
}
#endif

#endif
