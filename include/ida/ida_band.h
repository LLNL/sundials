/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, and Radu Serban @ LLNL
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
 * This is the header file for the IDABAND linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _IDABAND_H
#define _IDABAND_H

#include <ida/ida_direct.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : IDABand
 * -----------------------------------------------------------------
 * A call to the IDABand function links the main integrator       
 * with the IDABAND linear solver module.                         
 *                                                                
 * ida_mem is the pointer to the integrator memory returned by    
 *         IDACreate.                                                   
 *                                                                
 * mupper is the upper bandwidth of the banded Jacobian matrix.   
 *                                                                
 * mlower is the lower bandwidth of the banded Jacobian matrix.   
 *                                                                
 * The return values of IDABand are:                              
 *     IDADLS_SUCCESS   = 0  if successful                            
 *     IDADLS_LMEM_FAIL = -1 if there was a memory allocation failure 
 *     IDADLS_ILL_INPUT = -2 if the input was illegal or NVECTOR bad. 
 *                                                                
 * NOTE: The band linear solver assumes a serial implementation   
 *       of the NVECTOR package. Therefore, IDABand will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the N_VGetArrayPointer function exists.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDABand(void *ida_mem, long int Neq, long int mupper, long int mlower);

#ifdef __cplusplus
}
#endif

#endif
