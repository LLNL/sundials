/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:13 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the IDA dense linear solver IDALAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _IDALAPACK_H
#define _IDALAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              I D A L A P A C K     C O N S T A N T S
   * =================================================================
   */

  /* 
   * -----------------------------------------------------------------
   * IDALAPACK return values 
   * -----------------------------------------------------------------
   */

#define IDALAPACK_SUCCESS           0
#define IDALAPACK_MEM_NULL         -1
#define IDALAPACK_LMEM_NULL        -2
#define IDALAPACK_ILL_INPUT        -3
#define IDALAPACK_MEM_FAIL         -4

  /* Additional last_flag values */

#define IDALAPACK_JACFUNC_UNRECVR  -5
#define IDALAPACK_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Types : IDALapackDenseJacFn
   * -----------------------------------------------------------------
   *
   * A dense Jacobian approximation function djac must be of type 
   * IDALapackDenseJacFn.
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
   * c_j is the scalar in the system Jacobian, proportional to 
   *     the inverse of the step size h.
   *                                                                
   * jac_data is a pointer to user Jacobian data - the same as the    
   *     jdata parameter passed to IDALapack.                     
   *                                                                
   * Jac is the dense matrix (of type LapackMat) to be loaded by  
   *     an IDALapackDenseJacFn routine with an approximation to the   
   *     system Jacobian matrix                                  
   *            J = dF/dy' + gamma*dF/dy                            
   *     at the given point (t,y,y'), where the ODE system is    
   *     given by F(t,y,y') = 0.
   *     Note that Jac is NOT preset to zero!
   *                                                                
   * tmp1, tmp2, tmp3 are pointers to memory allocated for          
   *     N_Vectors which can be used by an IDALapackDenseJacFn routine 
   *     as temporary storage or work space.                     
   *                                                                
   * A IDALapackDenseJacFn should return                                
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
   *     IDAGetCurrentStep and IDAGetErrWeights, respectively 
   *     (see ida.h). The unit roundoff is available as 
   *     UNIT_ROUNDOFF defined in sundials_types.h.
   *
   * -----------------------------------------------------------------
   */
  
  
  typedef int (*IDALapackDenseJacFn)(int N, realtype t, realtype c_j,
                                     N_Vector y, N_Vector yp, N_Vector r, 
                                     LapackMat Jac, void *jac_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * -----------------------------------------------------------------
   * Types : IDALapackBandJacFn
   * -----------------------------------------------------------------
   */

  typedef int (*IDALapackBandJacFn)(int N, int mupper, int mlower,
                                    realtype t, realtype c_j, 
                                    N_Vector y, N_Vector yp, N_Vector r,
                                    LapackMat Jac, void *jac_data,
                                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : IDALapackDense
   * -----------------------------------------------------------------
   * A call to the IDALapackDense function links the main integrator
   * with the IDALAPACK linear solver using dense Jacobians.
   *
   * ida_mem is the pointer to the integrator memory returned by
   *           IDACreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of IDALapackDense is one of:
   *    IDALAPACK_SUCCESS   if successful
   *    IDALAPACK_MEM_NULL  if the IDA memory was NULL
   *    IDALAPACK_MEM_FAIL  if there was a memory allocation failure
   *    IDALAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int IDALapackDense(void *ida_mem, int N);

  /*
   * -----------------------------------------------------------------
   * Function : IDALapackBand
   * -----------------------------------------------------------------
   * A call to the IDALapackBand function links the main integrator
   * with the IDALAPACK linear solver using banded Jacobians. 
   *
   * ida_mem is the pointer to the integrator memory returned by
   *           IDACreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian approximation.
   *
   * The return value of IDALapackBand is one of:
   *    IDALAPACK_SUCCESS   if successful
   *    IDALAPACK_MEM_NULL  if the IDA memory was NULL
   *    IDALAPACK_MEM_FAIL  if there was a memory allocation failure
   *    IDALAPACK_ILL_INPUT if a required vector operation is missing or
   *                       if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int IDALapackBand(void *ida_mem, int N, int mupper, int mlower);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the IDALAPACK linear solver
   * -----------------------------------------------------------------
   *
   * IDALapackSetJacFn specifies the Jacobian approximation routine 
   * to be used. When using dense Jacobians, a user-supplied jac 
   * routine must be of type IDALapackDenseJacFn. When using banded 
   * Jacobians, a user-supplied jac routine must be of type 
   * IDALapackBandJacFn.
   * By default, a difference quotient approximation, 
   * supplied with this solver is used.
   * IDALapackSetJacFn also specifies a pointer to user data which is 
   * passed to the user's jac routine every time it is called.
   *
   * The return value of IDALapackSetJacFn is one of:
   *    IDALAPACK_SUCCESS   if successful
   *    IDALAPACK_MEM_NULL  if the IDA memory was NULL
   *    IDALAPACK_LMEM_NULL if the IDALAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int IDALapackSetJacFn(void *ida_mem, void *jac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the IDALAPACK linear solver
   * -----------------------------------------------------------------
   *
   * IDALapackGetWorkSpace returns the real and integer workspace used
   *                     by IDALAPACK.
   * IDALapackGetNumJacEvals returns the number of calls made to the
   *                       Jacobian evaluation routine jac.
   * IDALapackGetNumResEvals returns the number of calls to the user
   *                       f routine due to finite difference Jacobian
   *                       evaluation.
   * IDALapackGetLastFlag returns the last error flag set by any of
   *                    the IDALAPACK interface functions.
   *
   * The return value of IDALapackGet* is one of:
   *    IDALAPACK_SUCCESS   if successful
   *    IDALAPACK_MEM_NULL  if the IDA memory was NULL
   *    IDALAPACK_LMEM_NULL if the IDALAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int IDALapackGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS);
  int IDALapackGetNumJacEvals(void *ida_mem, long int *njevals);
  int IDALapackGetNumResEvals(void *ida_mem, long int *nfevalsLS);
  int IDALapackGetLastFlag(void *ida_mem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a IDALAPACK return flag
   * -----------------------------------------------------------------
   */

  char *IDALapackGetReturnFlagName(int flag);

#ifdef __cplusplus
}
#endif

#endif
