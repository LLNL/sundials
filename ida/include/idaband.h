/*
 * -----------------------------------------------------------------
 * $Revision: 1.19.2.2 $
 * $Date: 2005-04-06 23:39:54 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the IDA/IDAS band linear solver
 * module, IDABAND. It interfaces between the band module and the
 * integrator when a banded linear solver is appropriate.
 * -----------------------------------------------------------------
 */

#ifndef _IDABAND_H
#define _IDABAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"
 
/*
 * -----------------------------------------------------------------
 * Type : IDABandJacFn
 * -----------------------------------------------------------------
 * A banded Jacobian approximation function bjac must have the    
 * prototype given below. Its parameters are:                     
 *                                                                
 * Neq is the problem size, and length of all vector arguments.   
 *                                                                
 * mupper is the upper bandwidth of the banded Jacobian matrix.   
 *                                                                
 * mlower is the lower bandwidth of the banded Jacobian matrix.   
 *                                                                
 * tt is the current value of the independent variable t.        
 *                                                                
 * yy is the current value of the dependent variable vector,     
 *    namely the predicted value of y(t).                     
 *                                                                
 * yp is the current value of the derivative vector y',          
 *    namely the predicted value of y'(t).                    
 *                                                                
 * rr is the residual vector F(tt,yy,yp).                     
 *                                                                
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.
 *                                                                
 * jac_data  is a pointer to user Jacobian data - the same as the    
 *    jdata parameter passed to IDABand.                      
 *                                                                
 * Jac is the band matrix (of type BandMat) to be loaded by    
 *     an IDABandJacFn routine with an approximation to the    
 *     system Jacobian matrix                                  
 *            J = dF/dy + cj*dF/dy'                             
 *     at the given point (t,y,y'), where the DAE system is    
 *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
 *     the nonzero elements need to be loaded.  See note below.
 *                                                                
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          
 *     N_Vectors which can be used by an IDABandJacFn routine  
 *     as temporary storage or work space.                     
 *                                                                
 * NOTE: The following are two efficient ways to load Jac:
 *                                                                
 * (1) (with macros - no explicit data structure references)      
 *    for (j=0; j < Neq; j++) {                                   
 *       col_j = BAND_COL(Jac,j);                                  
 *       for (i=j-mupper; i <= j+mlower; i++) {                   
 *         generate J_ij = the (i,j)th Jacobian element           
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;                       
 *       }                                                        
 *     }                                                          
 *                                                                
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)     
 *    for (j=0; j < Neq; j++) {                                   
 *       col_j = BAND_COL(Jac,j);                                  
 *       for (k=-mupper; k <= mlower; k++) {                      
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k    
 *         col_j[k] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *       they are accessible as follows: hcur (the current stepsize)
 *       and ewt (the error weight vector) are accessible through   
 *       IDAGetCurrentStep and IDAGetErrWeights, respectively (see  
 *       ida.h). The unit roundoff is available as                  
 *       UNIT_ROUNDOFF defined in sundialstypes.h                   
 *                                                                
 * A third way, using the BAND_ELEM(A,i,j) macro, is much less    
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * The IDABandJacFn should return                                 
 *     0 if successful,                                           
 *     a positive int if a recoverable error occurred, or         
 *     a negative int if a nonrecoverable error occurred.         
 * In the case of a recoverable error return, the integrator will 
 * attempt to recover by reducing the stepsize (which changes cj).
 * -----------------------------------------------------------------
 */
  
typedef int (*IDABandJacFn)(long int Neq, long int mupper, 
                            long int mlower, realtype tt, 
                            N_Vector yy, N_Vector yp, N_Vector rr, 
                            realtype c_j, void *jac_data, BandMat Jac, 
                            N_Vector tmp1, N_Vector tmp2, 
                            N_Vector tmp3);
 
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
 *     IDABAND_SUCCESS   = 0  if successful                            
 *     IDABAND_LMEM_FAIL = -1 if there was a memory allocation failure 
 *     IDABAND_ILL_INPUT = -2 if the input was illegal or NVECTOR bad. 
 *                                                                
 * NOTE: The band linear solver assumes a serial implementation   
 *       of the NVECTOR package. Therefore, IDABand will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the N_VGetArrayPointer function exists.
 * -----------------------------------------------------------------
 */

int IDABand(void *ida_mem, long int Neq, long int mupper, long int mlower);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDABAND linear solver
 * -----------------------------------------------------------------
 *                                                                
 * IDABandSetJacFn specifies the dense Jacobian approximation     
 *         routine to be used. A user-supplied djac routine must  
 *         be of type IDABandJacFn.                               
 *         By default, a difference quotient routine IDABandDQJac,
 *         supplied with this solver is used.                     
 *         It also specifies a pointer to user data which is    
 *         passed to the bjac routine every time it is called.    
 *                                                                
 * The return value of IDABandSet* is one of:
 *    IDABAND_SUCCESS   if successful
 *    IDABAND_MEM_NULL  if the ida memory was NULL
 *    IDABAND_LMEM_NULL if the idaband memory was NULL
 * -----------------------------------------------------------------
 */

int IDABandSetJacFn(void *ida_mem, IDABandJacFn bjac, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDABAND linear solver
 * -----------------------------------------------------------------
 *                                                                
 * IDABandGetWorkSpace returns the real and integer workspace used 
 *     by IDABAND.                                                   
 * IDABandGetNumJacEvals returns the number of calls made to the  
 *     Jacobian evaluation routine bjac.                          
 * IDABandGetNumResEvals returns the number of calls to the user  
 *     res routine due to finite difference Jacobian evaluation.  
 * IDABandGetLastFlag returns the last error flag set by any of
 *     the IDABAND interface functions.
 *
 * The return value of IDABandGet* is one of:
 *    IDABAND_SUCCESS   if successful
 *    IDABAND_MEM_NULL  if the ida memory was NULL
 *    IDABAND_LMEM_NULL if the idaband memory was NULL
 * -----------------------------------------------------------------
 */

int IDABandGetWorkSpace(void *ida_mem, long int *lenrwB, long int *leniwB);
int IDABandGetNumJacEvals(void *ida_mem, long int *njevalsB);
int IDABandGetNumResEvals(void *ida_mem, long int *nrevalsB);
int IDABandGetLastFlag(void *ida_mem, int *flag);

/* IDABAND return values */

#define IDABAND_SUCCESS    0
#define IDABAND_MEM_NULL  -1 
#define IDABAND_LMEM_NULL -2 
#define IDABAND_ILL_INPUT -3
#define IDABAND_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
