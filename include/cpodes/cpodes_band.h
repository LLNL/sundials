/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:07 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CPODES band linear solver, CPBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CPBAND_H
#define _CPBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_band.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              C P D E N S E     C O N S T A N T S
   * =================================================================
   */

  /* 
   * -----------------------------------------------------------------
   * CPDENSE return values 
   * -----------------------------------------------------------------
   */

#define CPBAND_SUCCESS           0
#define CPBAND_MEM_NULL         -1
#define CPBAND_LMEM_NULL        -2
#define CPBAND_ILL_INPUT        -3
#define CPBAND_MEM_FAIL         -4

  /* Additional last_flag values */

#define CPBAND_JACFUNC_UNRECVR  -5
#define CPBAND_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Type : CPBandJacExplFn and CPBandJacImplFn
   * -----------------------------------------------------------------
   *
   * If the ODE is given in explicit form, a band Jacobian 
   * approximation function Jac must be of type CPBandJacExplFn. 
   * Its parameters are:
   *
   * N is the length of all vector arguments.
   *
   * mupper is the upper half-bandwidth of the approximate banded
   * Jacobian. This parameter is the same as the mupper parameter
   * passed by the user to the CPBand function.
   *
   * mlower is the lower half-bandwidth of the approximate banded
   * Jacobian. This parameter is the same as the mlower parameter
   * passed by the user to the CPBand function.
   *
   * J is the band matrix (of type BandMat) that will be loaded
   * by a CPBandJacFn with an approximation to the Jacobian matrix
   * J = (df_i/dy_j) at the point (t,y).
   * J is preset to zero, so only the nonzero elements need to be
   * loaded. 
   *
   * t is the current value of the independent variable.
   *
   * y is the current value of the dependent variable vector,
   *      namely the predicted value of y(t).
   *
   * fy is the vector f(t,y).
   *
   * jac_data is a pointer to user data - the same as the jac_data
   *          parameter passed to CPBand.
   *
   * tmp1, tmp2, and tmp3 are pointers to memory allocated for
   * vectors of length N which can be used by a CPBandJacFn
   * as temporary storage or work space.
   *
   * A CPBandJacExplFn should return 0 if successful, a positive value
   * if a recoverable error occurred, and a negative value if an 
   * unrecoverable error occurred.
   *
   * -----------------------------------------------------------------
   *
   * If the ODE is given in implicit form, a band Jacobian 
   * approximation function Jac must be of type CPBandJacImplFn. 
   * Its parameters are:
   *  
   * N is the problem size, and length of all vector arguments.   
   *                                                                
   * mupper is the upper bandwidth of the banded Jacobian matrix.   
   *                                                                
   * mlower is the lower bandwidth of the banded Jacobian matrix.   
   *                                                                
   * t is the current value of the independent variable t.        
   *                                                                
   * y is the current value of the dependent variable vector,     
   *    namely the predicted value of y(t).                     
   *                                                                
   * yp is the current value of the derivative vector y',          
   *    namely the predicted value of y'(t).                    
   *                                                                
   * r is the residual vector F(tt,yy,yp).                     
   *                                                                
   * gm  is the scalar in the system Jacobian, proportional to 
   *     the step size h.
   *                                                                
   * jac_data  is a pointer to user Jacobian data - the same as the    
   *    jdata parameter passed to CPBandSetJacFn.                      
   *                                                                
   * J is the band matrix (of type BandMat) to be loaded by    
   *     with an approximation to the system Jacobian matrix
   *            J = dF/dy' + gamma*dF/dy 
   *     at the given point (t,y,y'), where the ODE system is    
   *     given by F(t,y,y') = 0.  Jac is preset to zero, so only 
   *     the nonzero elements need to be loaded.  See note below.
   *                                                                
   * tmp1, tmp2, tmp3 are pointers to memory allocated for          
   *     N_Vectors which can be used by an CPBandJacImplFn routine  
   *     as temporary storage or work space.                     

   * -----------------------------------------------------------------
   *
   * NOTE: the following are three possible ways to load J:
   *
   * (1) (with macros - no explicit data structure references)
   *    for (j=0; j < n; j++) {
   *       col_j = BAND_COL(J,j);
   *       for (i=j-mupper; i <= j+mlower; i++) {
   *         generate J_ij = the (i,j)th Jacobian element
   *         BAND_COL_ELEM(col_j,i,j) = J_ij;
   *       }
   *     }
   *
   * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
   *    for (j=0; j < n; j++) {
   *       col_j = BAND_COL(J,j);
   *       for (k=-mupper; k <= mlower; k++) {
   *         generate J_ij = the (i,j)th Jacobian element, i=j+k
   *         col_j[k] = J_ij;
   *       }
   *     }
   *
   * (3) (without macros - explicit data structure references)
   *     offset = J->smu;
   *     for (j=0; j < n; j++) {
   *       col_j = ((J->data)[j])+offset;
   *       for (k=-mupper; k <= mlower; k++) {
   *         generate J_ij = the (i,j)th Jacobian element, i=j+k
   *         col_j[k] = J_ij;
   *       }
   *     }
   * Caution: J->smu is generally NOT the same as mupper.
   *
   * The BAND_ELEM(A,i,j) macro is appropriate for use in small
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

  typedef int (*CPBandJacExplFn)(long int N, long int mupper, long int mlower,
                                 realtype t, N_Vector y, N_Vector fy, 
                                 BandMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  typedef int (*CPBandJacImplFn)(long int N, long int mupper, long int mlower,
                                 realtype t, realtype gm, 
                                 N_Vector y, N_Vector yp, N_Vector r,
                                 BandMat Jac, void *jac_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /*
   * -----------------------------------------------------------------
   * Function : CPBand
   * -----------------------------------------------------------------
   * A call to the CPBand function links the main CPODES integrator
   * with the CPBAND linear solver.
   *
   * cpode_mem is the pointer to the integrator memory returned by
   *           CPodeCreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian
   *        approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian
   *        approximation.
   *
   * The return value of CPBand is one of:
   *    CPBAND_SUCCESS   if successful
   *    CPBAND_MEM_NULL  if the CPODES memory was NULL
   *    CPBAND_MEM_FAIL  if there was a memory allocation failure
   *    CPBAND_ILL_INPUT if a required vector operation is missing or
   *                     if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int CPBand(void *cpode_mem, long int N,
             long int mupper, long int mlower);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the CPBAND linear solver
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
   * The return value of CPBandSet* is one of:
   *    CPBAND_SUCCESS   if successful
   *    CPBAND_MEM_NULL  if the CPODES memory was NULL
   *    CPBAND_LMEM_NULL if the CPBAND memory was NULL
   * -----------------------------------------------------------------
   */

  int CPBandSetJacFn(void *cpode_mem, void *bjac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the CPBAND linear solver
   * -----------------------------------------------------------------
   *
   * CPBandGetWorkSpace returns the real and integer workspace used
   *                    by CPBAND.
   * CPBandGetNumJacEvals returns the number of calls made to the
   *                      Jacobian evaluation routine bjac.
   * CPBandGetNumRhsEvals returns the number of calls to the user
   *                      f routine due to finite difference Jacobian
   *                      evaluation.
   * CPBandGetLastFlag returns the last error flag set by any of
   *                   the CPBAND interface functions.
   *
   * The return value of CPBandGet* is one of:
   *    CPBAND_SUCCESS   if successful
   *    CPBAND_MEM_NULL  if the CPODES memory was NULL
   *    CPBAND_LMEM_NULL if the CPBAND memory was NULL
   * -----------------------------------------------------------------
   */

  int CPBandGetWorkSpace(void *cpode_mem, long int *lenrwLS, long int *leniwLS);
  int CPBandGetNumJacEvals(void *cpode_mem, long int *njevals);
  int CPBandGetNumRhsEvals(void *cpode_mem, long int *nfevalsLS);
  int CPBandGetLastFlag(void *cpode_mem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a CPBAND return flag
   * -----------------------------------------------------------------
   */

  char *CPBandGetReturnFlagName(int flag);


#ifdef __cplusplus
}
#endif

#endif
