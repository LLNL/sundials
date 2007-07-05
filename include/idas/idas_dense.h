/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2007-07-05 19:10:36 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the IDADENSE linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _IDASDENSE_H
#define _IDASDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <idas/idas_direct.h>
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function : IDADense
 * -----------------------------------------------------------------
 * A call to the IDADense function links the main integrator      
 * with the IDADENSE linear solver module.                        
 *                                                                
 * ida_mem is the pointer to integrator memory returned by        
 *     IDACreate.                                                 
 *                                                                
 * Neq  is the problem size                                       
 *                                                                
 * IDADense returns:                                              
 *     IDADIRECT_SUCCESS   = 0  if successful                              
 *     IDADIRECT_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     IDADIRECT_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The dense linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, IDADense will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDADense(void *ida_mem, int Neq); 

/*
 * -----------------------------------------------------------------
 * Function: IDADenseB
 * -----------------------------------------------------------------
 * IDADenseB links the main IDAS integrator with the IDADENSE
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDADenseB(void *ida_mem, int which, int NeqB);
  
#ifdef __cplusplus
}
#endif

#endif
