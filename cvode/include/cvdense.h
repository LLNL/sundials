/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2004-05-26 19:54:10 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES dense linear
 * solver, CVDENSE. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvdense_h
#define _cvdense_h

#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * CVDENSE solver constants                                       
 * -----------------------------------------------------------------
 * CVD_MSBJ  : maximum number of steps between dense Jacobian     
 *             evaluations                                        
 *                                                                
 * CVD_DGMAX : maximum change in gamma between dense Jacobian     
 *             evaluations                                        
 * -----------------------------------------------------------------
 */

#define CVD_MSBJ  50   

#define CVD_DGMAX RCONST(0.2)  
 
/*
 * -----------------------------------------------------------------
 * Type : CVDenseJacFn                                            
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function Jac must have the      
 * prototype given below. Its parameters are:                     
 *                                                                
 * N is the problem size.                                         
 *                                                                
 * J is the dense matrix (of type DenseMat) that will be loaded   
 * by a CVDenseJacFn with an approximation to the Jacobian matrix 
 * J = (df_i/dy_j) at the point (t,y).                            
 * J is preset to zero, so only the nonzero elements need to be   
 * loaded. Two efficient ways to load J are:                      
 *                                                                
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < n; j++) {                                    
 *       col_j = DENSE_COL(J,j);                                  
 *       for (i=0; i < n; i++) {                                  
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * 
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < n; j++) {                                    
 *       col_j = (J->data)[j];                                    
 *       for (i=0; i < n; i++) {                                  
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 *                                                                
 * The DENSE_ELEM(A,i,j) macro is appropriate for use in small    
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * t is the current value of the independent variable.            
 *                                                                
 * y is the current value of the dependent variable vector,       
 *      namely the predicted value of y(t).                       
 *                                                                
 * fy is the vector f(t,y).                                       
 *                                                                
 * jac_data is a pointer to user data - the same as the jac_data  
 *          parameter passed to CVDense.                          
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively   
 *     (see cvode.h). The unit roundoff is available as           
 *     UNIT_ROUNDOFF defined in sundialstypes.h                   
 *                                                                
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for      
 * vectors of length N which can be used by a CVDenseJacFn        
 * as temporary storage or work space.                            
 * -----------------------------------------------------------------
 */                                                                
  
typedef void (*CVDenseJacFn)(long int N, DenseMat J, realtype t, 
                             N_Vector y, N_Vector fy, void *jac_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 
 
/*
 * -----------------------------------------------------------------
 * Function : CVDense                                             
 * -----------------------------------------------------------------
 * A call to the CVDense function links the main integrator with  
 * the CVDENSE linear solver.                                     
 *                                                                
 * cvode_mem is the pointer to the integrator memory returned by  
 *              CVodeCreate.                                      
 *                                                                
 * N is the size of the ODE system.                               
 *                                                                
 * The return values of CVDense are:                              
 *    SUCCESS   = 0  if successful                                
 *    LMEM_FAIL = -1 if there was a memory allocation failure     
 * -----------------------------------------------------------------
 */                                                                
  
int CVDense(void *cvode_mem, long int N); 

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVDENSE linear solver                   
 * -----------------------------------------------------------------
 *                                                                
 * CVDenseSetJacFn specifies the dense Jacobian approximation     
 *         routine to be used. A user-supplied djac routine must  
 *         be of type CVDenseJacFn.                               
 *         By default, a difference quotient routine CVDenseDQJac,
 *         supplied with this solver is used.                     
 * CVDenseSetJacData specifies a pointer to user data which is    
 *         passed to the djac routine every time it is called.    
 * -----------------------------------------------------------------
 */                                                                

int CVDenseSetJacFn(void *cvode_mem, CVDenseJacFn djac);
int CVDenseSetJacData(void *cvode_mem, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVDENSE linear solver                
 * -----------------------------------------------------------------
 *                                                                
 * CVDenseGetIntWorkSpace returns the integer workspace used by   
 *     CVDENSE.                                                   
 * CVDenseGetRealWorkSpace returns the real workspace used by     
 *     CVDENSE.                                                   
 * CVDenseGetNumJacEvals returns the number of calls made to the  
 *     Jacobian evaluation routine djac.                          
 * CVDenseGetNumRhsEvals returns the number of calls to the user  
 *     f routine due to finite difference Jacobian evaluation.    
 * -----------------------------------------------------------------
 */                                                                

int CVDenseGetIntWorkSpace(void *cvode_mem, long int *leniwD);
int CVDenseGetRealWorkSpace(void *cvode_mem, long int *lenrwD);
int CVDenseGetNumJacEvals(void *cvode_mem, long int *njevalsD);
int CVDenseGetNumRhsEvals(void *cvode_mem, long int *nfevalsD);

#endif

#ifdef __cplusplus
}
#endif
