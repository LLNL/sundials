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
 * Header file for the CPODES dense linear solver CPLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CPLAPACK_H
#define _CPLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              C P L A P A C K     C O N S T A N T S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * CPLAPACK input constants
   * -----------------------------------------------------------------
   * fact_type: is the type of constraint Jacobian factorization used
   * for the projection. For independent constraints (i.e. Jacobian
   * with full row rank) any of the following three options can be
   * used (although their computational cost varies, depending on 
   * the number of states and constraints)
   *   CPLAPACK_LU:  use LU decomposition of G^T.
   *   CPLAPACK_QR:  use QR decomposition of G^T
   *   CPLAPACK_SC:  use Schur complement
   * If it is known (or suspected) that some constraints are redundant,
   * the following option should be used:
   *   CPLAPACK_QRP: use QR with column pivoting on G^T.  
   * -----------------------------------------------------------------
   */

  /* fact_type */

#define CPLAPACK_LU   1
#define CPLAPACK_QR   2
#define CPLAPACK_SC   3
#define CPLAPACK_QRP  4

  /* 
   * -----------------------------------------------------------------
   * CPLAPACK return values 
   * -----------------------------------------------------------------
   */

#define CPLAPACK_SUCCESS           0
#define CPLAPACK_MEM_NULL         -1
#define CPLAPACK_LMEM_NULL        -2
#define CPLAPACK_ILL_INPUT        -3
#define CPLAPACK_MEM_FAIL         -4

  /* Additional last_flag values */

#define CPLAPACK_JACFUNC_UNRECVR  -5
#define CPLAPACK_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Types : CPLapackDenseJacExplFn and CPLapackDenseJacImplFn
   * -----------------------------------------------------------------
   *
   * If the ODE is given in explicit form, a dense Jacobian 
   * approximation function Jac must be of type CPLapackDenseJacExplFn. 
   * Its parameters are:
   *
   * N   is the problem size.
   *
   * Jac is the dense matrix (of type LapackMat) that will be loaded
   *     by a CPLapackDenseJacExplFn with an approximation to the 
   *     Jacobian matrix J = (df_i/dy_j) at the point (t,y).
   *     Note that Jac is NOT preset to zero!
   *
   * t   is the current value of the independent variable.
   *
   * y   is the current value of the dependent variable vector,
   *     namely the predicted value of y(t).
   *
   * fy  is the vector f(t,y).
   *
   * jac_data is a pointer to user data - the same as the jac_data
   *     parameter passed to CPLapack.
   *
   * tmp1, tmp2, and tmp3 are pointers to memory allocated for
   * vectors of length N which can be used by a CPLapackDenseJacExplFn
   * as temporary storage or work space.
   *
   * A CPLapackDenseJacExplFn should return 0 if successful, a positive 
   * value if a recoverable error occurred, and a negative value if 
   * an unrecoverable error occurred.
   *
   * -----------------------------------------------------------------
   *
   * If the ODE is given in implicit form, a dense Jacobian 
   * approximation function djac must be of type CPLapackDenseJacImplFn.
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
   *     jdata parameter passed to CPLapack.                     
   *                                                                
   * Jac is the dense matrix (of type LapackMat) to be loaded by  
   *     an CPLapackDenseJacImplFn routine with an approximation to the   
   *     system Jacobian matrix                                  
   *            J = dF/dy' + gamma*dF/dy                            
   *     at the given point (t,y,y'), where the ODE system is    
   *     given by F(t,y,y') = 0.
   *     Note that Jac is NOT preset to zero!
   *                                                                
   * tmp1, tmp2, tmp3 are pointers to memory allocated for          
   *     N_Vectors which can be used by an CPLapackDenseJacImplFn routine 
   *     as temporary storage or work space.                     
   *                                                                
   * A CPLapackDenseJacImplFn should return                                
   *     0 if successful,                                           
   *     a positive int if a recoverable error occurred, or         
   *     a negative int if a nonrecoverable error occurred.         
   * In the case of a recoverable error return, the integrator will 
   * attempt to recover by reducing the stepsize (which changes cj).
   *
   * -----------------------------------------------------------------
   *
   * NOTE: The following are two efficient ways to load a dense Jac:         
   * (1) (with macros - no explicit data structure references)      
   *     for (j=0; j < Neq; j++) {                                  
   *       col_j = LAPACK_DENSE_COL(Jac,j);                                 
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
   * A third way, using the LAPACK_DENSE_ELEM(A,i,j) macro, is much less   
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
  
  
  typedef int (*CPLapackDenseJacExplFn)(int N, realtype t,
                                        N_Vector y, N_Vector fy, 
                                        LapackMat Jac, void *jac_data,
                                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
  typedef int (*CPLapackDenseJacImplFn)(int N, realtype t, realtype gm,
                                        N_Vector y, N_Vector yp, N_Vector r, 
                                        LapackMat Jac, void *jac_data,
                                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * -----------------------------------------------------------------
   * Types : CPLapackBandJacExplFn and CPLapackBandJacImplFn
   * -----------------------------------------------------------------
   */

  typedef int (*CPLapackBandJacExplFn)(int N, int mupper, int mlower,
                                       realtype t, N_Vector y, N_Vector fy, 
                                       LapackMat Jac, void *jac_data,
                                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  typedef int (*CPLapackBandJacImplFn)(int N, int mupper, int mlower,
                                       realtype t, realtype gm, 
                                       N_Vector y, N_Vector yp, N_Vector r,
                                       LapackMat Jac, void *jac_data,
                                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

  /*
   * -----------------------------------------------------------------
   * Type : CPLapackDenseProjJacFn
   * -----------------------------------------------------------------
   *
   *
   * -----------------------------------------------------------------
   */

  typedef int (*CPLapackDenseProjJacFn)(int Nc, int Ny, 
                                        realtype t, N_Vector y, N_Vector cy,
                                        LapackMat Jac, void *jac_data,
                                        N_Vector tmp1, N_Vector tmp2); 

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   *      F O R    I M P L I C I T    I N T E G R A T I O N      
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CPLapackDense
   * -----------------------------------------------------------------
   * A call to the CPLapackDense function links the main integrator
   * with the CPLAPACK linear solver using dense Jacobians.
   *
   * cpode_mem is the pointer to the integrator memory returned by
   *           CPodeCreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of CPLapackDense is one of:
   *    CPLAPACK_SUCCESS   if successful
   *    CPLAPACK_MEM_NULL  if the CPODES memory was NULL
   *    CPLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CPLAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CPLapackDense(void *cpode_mem, int N);

  /*
   * -----------------------------------------------------------------
   * Function : CPLapackBand
   * -----------------------------------------------------------------
   * A call to the CPLapackBand function links the main integrator
   * with the CPLAPACK linear solver using banded Jacobians. 
   *
   * cpode_mem is the pointer to the integrator memory returned by
   *           CPodeCreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian approximation.
   *
   * The return value of CPLapackBand is one of:
   *    CPLAPACK_SUCCESS   if successful
   *    CPLAPACK_MEM_NULL  if the CPODES memory was NULL
   *    CPLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CPLAPACK_ILL_INPUT if a required vector operation is missing or
   *                       if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int CPLapackBand(void *cpode_mem, int N, int mupper, int mlower);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the CPLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * CPLapackSetJacFn specifies the Jacobian approximation routine 
   * to be used. When using dense Jacobians, a user-supplied jac 
   * routine must be of type CPLapackDenseJacExplFn or CPLapackDenseJacImplFn, 
   * depending of the form in which the ODE system is given (explicit 
   * or implicit). When using banded Jacobians, a user-supplied jac 
   * routine must be of type CPLapackBandJacExplFn or CPLapackBandJacImplFn, 
   * depending of the form in which the ODE system is given (explicit 
   * or implicit). By default, a difference quotient approximation, 
   * supplied with this solver is used.
   * CPLapackSetJacFn also specifies a pointer to user data which is 
   * passed to the user's jac routine every time it is called.
   *
   * The return value of CPLapackSetJacFn is one of:
   *    CPLAPACK_SUCCESS   if successful
   *    CPLAPACK_MEM_NULL  if the CPODES memory was NULL
   *    CPLAPACK_LMEM_NULL if the CPLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int CPLapackSetJacFn(void *cpode_mem, void *jac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the CPLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * CPLapackGetWorkSpace returns the real and integer workspace used
   *                     by CPLAPACK.
   * CPLapackGetNumJacEvals returns the number of calls made to the
   *                       Jacobian evaluation routine jac.
   * CPLapackGetNumFctEvals returns the number of calls to the user
   *                       f routine due to finite difference Jacobian
   *                       evaluation.
   * CPLapackGetLastFlag returns the last error flag set by any of
   *                    the CPLAPACK interface functions.
   *
   * The return value of CPLapackGet* is one of:
   *    CPLAPACK_SUCCESS   if successful
   *    CPLAPACK_MEM_NULL  if the CPODES memory was NULL
   *    CPLAPACK_LMEM_NULL if the CPLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int CPLapackGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS);
  int CPLapackGetNumJacEvals(void *cpode_mem, long int *njevals);
  int CPLapackGetNumFctEvals(void *cpode_mem, long int *nfevalsLS);
  int CPLapackGetLastFlag(void *cpode_mem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a CPLAPACK return flag
   * -----------------------------------------------------------------
   */

  char *CPLapackGetReturnFlagName(int flag);


  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   *                F O R    P R O J E C T I O N
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CPLapackDenseProj
   * -----------------------------------------------------------------
   * A call to the CPLapackDenseProj function links the main integrator 
   * with the CPLAPACK linear solver using dense Jacobians.
   *
   * cpode_mem  the pointer to the integrator memory returned by
   *            CPodeCreate.
   * Nc         the number of constraints
   * Ny         the number of states (size of the ODE system).
   * fact_type  the type of factorization used for the constraint
   *            Jacobian G. Legal values are CPLAPACK_LU, CPLAPACK_QR,
   *            and CPLAPACK_SC.
   *
   * The return value of CPLapackDenseProj is one of:
   *    CPLAPACK_SUCCESS   if successful
   *    CPLAPACK_MEM_NULL  if the CPODES memory was NULL
   *    CPLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CPLAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CPLapackDenseProj(void *cpode_mem, int Nc, int Ny, int fact_type);

  /*
   * -----------------------------------------------------------------
   * Optional I/O functions
   * -----------------------------------------------------------------
   */

  int CPLapackProjSetJacFn(void *cpode_mem, void *jacP, void *jacP_data);

  int CPLapackProjGetNumJacEvals(void *cpode_mem, long int *njPevals);
  int CPLapackProjGetNumFctEvals(void *cpode_mem, long int *ncevalsLS);


#ifdef __cplusplus
}
#endif

#endif
