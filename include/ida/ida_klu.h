/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the IDAKLU linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _IDAKLU_H
#define _IDAKLU_H

#include "ida/ida_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : IDAKLU
 * -----------------------------------------------------------------
 * A call to the IDAKLU function links the main integrator      
 * with the IDAKLU linear solver module.                        
 *                                                                
 * ida_mem is the pointer to integrator memory returned by        
 *     IDACreate.             
 *
 *                                                                
 * IDAKLU returns:                                              
 *     IDASLU_SUCCESS   = 0  if successful                              
 *     IDASLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     IDASLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The KLU linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, IDAKLU will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDAKLU(void *ida_mem, int n, int nnz); 

/*
 * -----------------------------------------------------------------
 * Function: IDAKLUB
 * -----------------------------------------------------------------
 * IDAKLUB links the main IDAS integrator with the IDAKLU
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDAKLUB(void *ida_mem, int which, int nB, int nnzB);

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * IDAKLUSetOrdering sets the ordering used by KLU for reducing fill.
 * Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
 * The default used in KINSOL is 1 for COLAMD.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDAKLUSetOrdering(void *ida_mem, int ordering_choice); 


  
#ifdef __cplusplus
}
#endif

#endif
