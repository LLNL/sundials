/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2004-07-22 21:14:47 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES diagonal linear
 * solver, CVDIAG. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvdiag_h
#define _cvdiag_h

#include <stdio.h>
#include "nvector.h"
#include "sundialstypes.h"
 
/*
 * -----------------------------------------------------------------
 *                                                                
 * Function : CVDiag                                              
 * -----------------------------------------------------------------
 * A call to the CVDiag function links the main integrator with   
 * the CVDIAG linear solver.                                      
 *                                                                
 * cvode_mem is the pointer to the integrator memory returned by  
 *              CVodeCreate.                                      
 *                                                                
 * The return values of CVDiag are:                               
 *    SUCCESS       =  0 if successful                                
 *    LMEM_FAIL     = -1 if there was a memory allocation failure     
 *    LIN_ILL_INPUT = -2 if a required vector operation is missing
 * -----------------------------------------------------------------
 */                                                                

int CVDiag(void *cvode_mem);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVDIAG linear solver                 
 * -----------------------------------------------------------------
 *                                                                
 * CVDiagGetIntWorkSpace returns the integer workspace used by    
 *     CVDIAG.                                                    
 * CVDiagGetRealWorkSpace returns the real workspace used by      
 *     CVDIAG.                                                    
 * CVDiagGetNumRhsEvals returns the number of calls to the user   
 *     f routine due to finite difference Jacobian evaluation.    
 * Note: the number of diagonal approximate Jacobians formed is   
 * equal to the number of CVDiagSetup calls.                      
 * This number is available through CVodeGetNumLinSolvSetups.     
 * -----------------------------------------------------------------
 */                                                                

int CVDiagGetIntWorkSpace(void *cvode_mem, long int *leniwDI);
int CVDiagGetRealWorkSpace(void *cvode_mem, long int *lenrwDI);
int CVDiagGetNumRhsEvals(void *cvode_mem, long int *nfevalsDI);

#endif

#ifdef __cplusplus
}
#endif
