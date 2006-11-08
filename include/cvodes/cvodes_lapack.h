/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:12 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the CVODES dense linear solver CVLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CVLAPACK_H
#define _CVLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              C V L A P A C K     C O N S T A N T S
   * =================================================================
   */

  /* 
   * -----------------------------------------------------------------
   * CVLAPACK return values 
   * -----------------------------------------------------------------
   */

#define CVLAPACK_SUCCESS           0
#define CVLAPACK_MEM_NULL         -1
#define CVLAPACK_LMEM_NULL        -2
#define CVLAPACK_ILL_INPUT        -3
#define CVLAPACK_MEM_FAIL         -4

  /* Additional last_flag values */

#define CVLAPACK_JACFUNC_UNRECVR  -5
#define CVLAPACK_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Types : CVLapackDenseJacFn
   * -----------------------------------------------------------------
   *
   * A dense Jacobian approximation function Jac must be of type 
   * CVLapackDenseJacFn. 
   * Its parameters are:
   *
   * N   is the problem size.
   *
   * Jac is the dense matrix (of type LapackMat) that will be loaded
   *     by a CVLapackDenseJacFn with an approximation to the Jacobian 
   *     matrix J = (df_i/dy_j) at the point (t,y). 
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
   *     parameter passed to CVLapack.
   *
   * tmp1, tmp2, and tmp3 are pointers to memory allocated for
   * vectors of length N which can be used by a CVLapackDenseJacFn
   * as temporary storage or work space.
   *
   * A CVLapackDenseJacFn should return 0 if successful, a positive 
   * value if a recoverable error occurred, and a negative value if 
   * an unrecoverable error occurred.
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
   *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
   *     (see cvodes.h). The unit roundoff is available as 
   *     UNIT_ROUNDOFF defined in sundials_types.h.
   *
   * -----------------------------------------------------------------
   */
  
  
  typedef int (*CVLapackDenseJacFn)(int N, realtype t,
                                    N_Vector y, N_Vector fy, 
                                    LapackMat Jac, void *jac_data,
                                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
  /*
   * -----------------------------------------------------------------
   * Types : CVLapackBandJacFn
   * -----------------------------------------------------------------
   */

  typedef int (*CVLapackBandJacFn)(int N, int mupper, int mlower,
                                   realtype t, N_Vector y, N_Vector fy, 
                                   LapackMat Jac, void *jac_data,
                                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : CVLapackDense
   * -----------------------------------------------------------------
   * A call to the CVLapackDense function links the main integrator
   * with the CVLAPACK linear solver using dense Jacobians.
   *
   * cvode_mem is the pointer to the integrator memory returned by
   *           CVodeCreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of CVLapackDense is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODES memory was NULL
   *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CVLAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int CVLapackDense(void *cvode_mem, int N);

  /*
   * -----------------------------------------------------------------
   * Function : CVLapackBand
   * -----------------------------------------------------------------
   * A call to the CVLapackBand function links the main integrator
   * with the CVLAPACK linear solver using banded Jacobians. 
   *
   * cvode_mem is the pointer to the integrator memory returned by
   *           CVodeCreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian approximation.
   *
   * The return value of CVLapackBand is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODES memory was NULL
   *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    CVLAPACK_ILL_INPUT if a required vector operation is missing or
   *                       if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the CVLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * CVLapackSetJacFn specifies the Jacobian approximation routine 
   * to be used. When using dense Jacobians, a user-supplied jac 
   * routine must be of type CVLapackDenseJacFn. When using banded 
   * Jacobians, a user-supplied jac routine must be of type 
   * CVLapackBandJacFn.
   * By default, a difference quotient approximation, supplied with this 
   * solver is used.
   * CVLapackSetJacFn also specifies a pointer to user data which is 
   * passed to the user's jac routine every time it is called.
   *
   * The return value of CVLapackSetJacFn is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODES memory was NULL
   *    CVLAPACK_LMEM_NULL if the CVLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int CVLapackSetJacFn(void *cvode_mem, void *jac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the CVLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * CVLapackGetWorkSpace returns the real and integer workspace used
   *                     by CVLAPACK.
   * CVLapackGetNumJacEvals returns the number of calls made to the
   *                       Jacobian evaluation routine jac.
   * CVLapackGetNumRhsEvals returns the number of calls to the user
   *                       f routine due to finite difference Jacobian
   *                       evaluation.
   * CVLapackGetLastFlag returns the last error flag set by any of
   *                    the CVLAPACK interface functions.
   *
   * The return value of CVLapackGet* is one of:
   *    CVLAPACK_SUCCESS   if successful
   *    CVLAPACK_MEM_NULL  if the CVODES memory was NULL
   *    CVLAPACK_LMEM_NULL if the CVLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int CVLapackGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
  int CVLapackGetNumJacEvals(void *cvode_mem, long int *njevals);
  int CVLapackGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);
  int CVLapackGetLastFlag(void *cvode_mem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a CVLAPACK return flag
   * -----------------------------------------------------------------
   */

  char *CVLapackGetReturnFlagName(int flag);


#ifdef __cplusplus
}
#endif

#endif
