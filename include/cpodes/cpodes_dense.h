/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:08 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CPODES dense linear solver, CPDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CPDENSE_H
#define _CPDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_dense.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              C P D E N S E     C O N S T A N T S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * CPDENSE input constants
   * -----------------------------------------------------------------
   * fact_type: Type of factorization used for the constraint Jacobian
   *  CPDENSE_LU: use LU decomposition of constraint Jacobian.
   *  CPDENSE_QR: use QR decomposition of constraint Jacobian.
   *  CPDENSE_SC: use Schur complement
   * -----------------------------------------------------------------
   */

  /* fact_type */

#define CPDENSE_LU   1
#define CPDENSE_QR   2
#define CPDENSE_SC   3

  /* 
   * -----------------------------------------------------------------
   * CPDENSE return values 
   * -----------------------------------------------------------------
   */

#define CPDENSE_SUCCESS           0
#define CPDENSE_MEM_NULL         -1
#define CPDENSE_LMEM_NULL        -2
#define CPDENSE_ILL_INPUT        -3
#define CPDENSE_MEM_FAIL         -4

  /* Additional last_flag values */

#define CPDENSE_JACFUNC_UNRECVR  -5
#define CPDENSE_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Types : CPDenseJacExplFn and CPDenseJacImplFn
   * -----------------------------------------------------------------
   *
   * If the ODE is given in explicit form, a dense Jacobian 
   * approximation function Jac must be of type CPDenseJacExplFn. 
   * Its parameters are:
   *
   * N   is the problem size.
   *
   * Jac is the dense matrix (of type DenseMat) that will be loaded
   *     by a CPDenseJacExplFn with an approximation to the Jacobian matrix
   *     J = (df_i/dy_j) at the point (t,y). J is preset to zero, so only 
   *     the nonzero elements need to be loaded. 
   *
   * t   is the current value of the independent variable.
   *
   * y   is the current value of the dependent variable vector,
   *     namely the predicted value of y(t).
   *
   * fy  is the vector f(t,y).
   *
   * jac_data is a pointer to user data - the same as the jac_data
   *     parameter passed to CVDense.
   *
   * tmp1, tmp2, and tmp3 are pointers to memory allocated for
   * vectors of length N which can be used by a CVDenseJacFn
   * as temporary storage or work space.
   *
   * A CPDenseJacExplFn should return 0 if successful, a positive 
   * value if a recoverable error occurred, and a negative value if 
   * an unrecoverable error occurred.
   *
   * -----------------------------------------------------------------
   *
   * If the ODE is given in implicit form, a dense Jacobian 
   * approximation function djac must be of type CPDenseJacImplFn. 
   * Its parameters are:                     
   *                                                                
   * N   is the problem size, and length of all vector arguments.   
   *                                                                
   * t   is the current value of the independent variable t.        
   *                                                                
   * y   is the current value of the dependent variable vector,     
   *     namely the predicted value of y(t).                     
   *                                                                
   * yp  is the current value of the derivative vector y',          
   *     namely the predicted value of y'(t).                    
   *                                                                
   * f   is the residual vector F(tt,yy,yp).                     
   *                                                                
   * gm  is the scalar in the system Jacobian, proportional to 
   *     the step size h.
   *                                                                
   * jac_data is a pointer to user Jacobian data - the same as the    
   *     jdata parameter passed to CPDenseSetJacFn.                     
   *                                                                
   * Jac is the dense matrix (of type DenseMat) to be loaded by  
   *     an CPDenseJacImplFn routine with an approximation to the   
   *     system Jacobian matrix                                  
   *            J = dF/dy' + gamma*dF/dy                            
   *     at the given point (t,y,y'), where the ODE system is    
   *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
   *     the nonzero elements need to be loaded.  See note below.
   *                                                                
   * tmp1, tmp2, tmp3 are pointers to memory allocated for          
   *     N_Vectors which can be used by an CPDenseJacImplFn routine 
   *     as temporary storage or work space.                     
   *                                                                
   * S CPDenseJacImplFn should return                                
   *     0 if successful,                                           
   *     a positive int if a recoverable error occurred, or         
   *     a negative int if a nonrecoverable error occurred.         
   * In the case of a recoverable error return, the integrator will 
   * attempt to recover by reducing the stepsize (which changes cj).
   *
   * -----------------------------------------------------------------
   *
   * NOTE: The following are two efficient ways to load Jac:         
   * (1) (with macros - no explicit data structure references)      
   *     for (j=0; j < Neq; j++) {                                  
   *       col_j = DENSE_COL(Jac,j);                                 
   *       for (i=0; i < Neq; i++) {                                
   *         generate J_ij = the (i,j)th Jacobian element           
   *         col_j[i] = J_ij;                                       
   *       }                                                        
   *     }                                                          
   * (2) (without macros - explicit data structure references)      
   *     for (j=0; j < Neq; j++) {                                  
   *       col_j = (Jac->data)[j];                                   
   *       for (i=0; i < Neq; i++) {                                
   *         generate J_ij = the (i,j)th Jacobian element           
   *         col_j[i] = J_ij;                                       
   *       }                                                        
   *     }                                                          
   * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
   * efficient in general.  It is only appropriate for use in small 
   * problems in which efficiency of access is NOT a major concern. 
   *                                                                
   * NOTE: If the user's Jacobian routine needs other quantities,   
   *     they are accessible as follows: hcur (the current stepsize)
   *     and ewt (the error weight vector) are accessible through   
   *     CPodeGetCurrentStep and CPodeGetErrWeights, respectively 
   *     (see cpodes.h). The unit roundoff is available as 
   *     UNIT_ROUNDOFF defined in sundials_types.h.
   *                                                                
   * -----------------------------------------------------------------
   */
  
  
  typedef int (*CPDenseJacExplFn)(long int N, realtype t,
                                  N_Vector y, N_Vector fy, 
                                  DenseMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
  typedef int (*CPDenseJacImplFn)(long int N, realtype t, realtype gm,
                                  N_Vector y, N_Vector yp, N_Vector r, 
                                  DenseMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * -----------------------------------------------------------------
   * Type : CPDenseProjJacFn
   * -----------------------------------------------------------------
   *
   *
   * -----------------------------------------------------------------
   */

  typedef int (*CPDenseProjJacFn)(long int Nc, long int Ny, 
                                  realtype t, N_Vector y, N_Vector cy,
                                  DenseMat Jac, void *jac_data,
                                  N_Vector tmp1, N_Vector tmp2); 

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   *      F O R    I M P L I C I T    I N T E G R A T I O N      
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CPDense
   * -----------------------------------------------------------------
   * A call to the CPDense function links the main integrator with
   * the CPDENSE linear solver.
   *
   * cpode_mem is the pointer to the integrator memory returned by
   *           CPodeCreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of CPDense is one of:
   *    CPDENSE_SUCCESS   if successful
   *    CPDENSE_MEM_NULL  if the CPODES memory was NULL
   *    CPDENSE_MEM_FAIL  if there was a memory allocation failure
   *    CPDENSE_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CPDense(void *cpode_mem, long int N);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the CPDENSE linear solver
   * -----------------------------------------------------------------
   *
   * CPDenseSetJacFn specifies the dense Jacobian approximation routine 
   * to be used. A user-supplied djac routine must be of type 
   * CPDenseJacExplFn or CPDenseJacImplFn, depending of the form in which
   * the ODE system is given (explicit or implicit). By default, a
   * difference quotient approximation, supplied with this solver is used.
   * It also specifies a pointer to user data which is passed to the djac 
   * routine every time it is called.
   *
   * The return value of CPDenseSet* is one of:
   *    CPDENSE_SUCCESS   if successful
   *    CPDENSE_MEM_NULL  if the CPODES memory was NULL
   *    CPDENSE_LMEM_NULL if the CPDENSE memory was NULL
   * -----------------------------------------------------------------
   */

  int CPDenseSetJacFn(void *cpode_mem, void *djac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the CPDENSE linear solver
   * -----------------------------------------------------------------
   *
   * CPDenseGetWorkSpace returns the real and integer workspace used
   *                     by CPDENSE.
   * CPDenseGetNumJacEvals returns the number of calls made to the
   *                       Jacobian evaluation routine djac.
   * CPDenseGetNumFctEvals returns the number of calls to the user
   *                       f routine due to finite difference Jacobian
   *                       evaluation.
   * CPDenseGetLastFlag returns the last error flag set by any of
   *                    the CPDENSE interface functions.
   *
   * The return value of CPDenseGet* is one of:
   *    CPDENSE_SUCCESS   if successful
   *    CPDENSE_MEM_NULL  if the CPODES memory was NULL
   *    CPDENSE_LMEM_NULL if the CPDENSE memory was NULL
   * -----------------------------------------------------------------
   */

  int CPDenseGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS);
  int CPDenseGetNumJacEvals(void *cpode_mem, long int *njevals);
  int CPDenseGetNumFctEvals(void *cpode_mem, long int *nfevalsLS);
  int CPDenseGetLastFlag(void *cpode_mem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a CPDENSE return flag
   * -----------------------------------------------------------------
   */

  char *CPDenseGetReturnFlagName(int flag);


  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   *                F O R    P R O J E C T I O N
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CPDenseProj
   * -----------------------------------------------------------------
   * A call to the CPDenseProj function links the main integrator with
   * the CPDENSE linear solver.
   *
   * cpode_mem  the pointer to the integrator memory returned by
   *            CPodeCreate.
   * Nc         the number of constraints
   * Ny         the number of states (size of the ODE system).
   * fact_type  the type of factorization used for the constraint
   *            Jcobian G. Legal values are CPDENSE_LU and CPDENSE_LQ.
   *
   * The return value of CPDense is one of:
   *    CPDENSE_SUCCESS   if successful
   *    CPDENSE_MEM_NULL  if the CPODES memory was NULL
   *    CPDENSE_MEM_FAIL  if there was a memory allocation failure
   *    CPDENSE_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CPDenseProj(void *cpode_mem, long int Nc, long int Ny, int fact_type);

  /*
   * -----------------------------------------------------------------
   * Opyional I/O functions
   * -----------------------------------------------------------------
   */

  int CPDenseProjSetJacFn(void *cpode_mem, void *djacP, void *jacP_data);
  int CPDenseProjGetNumJacEvals(void *cpode_mem, long int *njPevals);
  int CPDenseProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS);

#ifdef __cplusplus
}
#endif

#endif
