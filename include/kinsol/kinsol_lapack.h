/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:15 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the KINSOL dense linear solver KINLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _KINLAPACK_H
#define _KINLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_nvector.h>

  /*
   * =================================================================
   *              K I N L A P A C K     C O N S T A N T S
   * =================================================================
   */

  /* 
   * -----------------------------------------------------------------
   * KINLAPACK return values 
   * -----------------------------------------------------------------
   */

#define KINLAPACK_SUCCESS           0
#define KINLAPACK_MEM_NULL         -1
#define KINLAPACK_LMEM_NULL        -2
#define KINLAPACK_ILL_INPUT        -3
#define KINLAPACK_MEM_FAIL         -4

  /* Additional last_flag values */

#define KINLAPACK_JACFUNC_UNRECVR  -5
#define KINLAPACK_JACFUNC_RECVR    -6

  /*
   * =================================================================
   *              F U N C T I O N   T Y P E S
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Type : KINLapackDenseJacFn
   * -----------------------------------------------------------------
   * A dense Jacobian approximation function Jac must have the
   * prototype given below. Its parameters are:
   *
   * N is the problem size.
   *
   * Jac is the dense matrix (of type LapackMat) that will be loaded
   * by a KINLapackDenseJacFn with an approximation to the Jacobian matrix
   * Jac = (df_i/dy_j).
   * Two efficient ways to load Jac are:
   *
   * (1) (with macros - no explicit data structure references)
   *     for (j=0; j < n; j++) {
   *       col_j = LAPACK_DENSE_COL(Jac,j);
   *       for (i=0; i < n; i++) {
   *         generate J_ij = the (i,j)th Jacobian element
   *         col_j[i] = J_ij;
   *       }
   *     }
   *
   * (2) (without macros - explicit data structure references)
   *     for (j=0; j < n; j++) {
   *       col_j = (Jac->data)[j];
   *       for (i=0; i < n; i++) {
   *         generate J_ij = the (i,j)th Jacobian element
   *         col_j[i] = J_ij;
   *       }
   *     }
   *
   * The LAPACK_DENSE_ELEM(A,i,j) macro is appropriate for use in small
   * problems in which efficiency of access is NOT a major concern.
   *
   * u   current iterate (unscaled) [input]
   *
   * fu  vector (type N_Vector) containing result of nonlinear
   *     system function evaluated at current iterate:
   *     fval = F(uu) [input]
   *
   * jac_data is a pointer to user data - the same as the jac_data
   *          parameter passed to KINLapackSetJacFn.
   *
   * tmp1, tmp2  available scratch vectors (volatile storage)
   *
   * If successful, the function should return 0 (zero). If an error
   * occurs, then the routine should return a non-zero integer value.
   * -----------------------------------------------------------------
   */

  typedef int (*KINLapackDenseJacFn)(int N, 
                                     N_Vector u, N_Vector fu, 
                                     LapackMat Jac, void *jac_data,
                                     N_Vector tmp1, N_Vector tmp2);

  /*
   * -----------------------------------------------------------------
   * Type : KINLapackBandJacFn
   * -----------------------------------------------------------------
   * A band Jacobian approximation function Jac must have the
   * prototype given below. Its parameters are:
   *
   * N is the problem size.
   *
   * mupper is the upper half-bandwidth of the approximate banded
   * Jacobian. This parameter is the same as the mupper parameter
   * passed by the user to the KINLapackBand function.
   *
   * mlower is the lower half-bandwidth of the approximate banded
   * Jacobian. This parameter is the same as the mlower parameter
   * passed by the user to the KINLapackBand function.
   *
   * Jac is the band matrix (of type LapackMat) that will be loaded
   * by a KINLapackBandJacFn with an approximation to the Jacobian matrix
   * Jac = (df_i/dy_j).
   * Three efficient ways to load J are:
   *
   * (1) (with macros - no explicit data structure references)
   *    for (j=0; j < n; j++) {
   *       col_j = LAPACK_BAND_COL(Jac,j);
   *       for (i=j-mupper; i <= j+mlower; i++) {
   *         generate J_ij = the (i,j)th Jacobian element
   *         LAPACK_BAND_COL_ELEM(col_j,i,j) = J_ij;
   *       }
   *     }
   *
   * (2) (with LAPACK_BAND_COL macro, but without LAPACK_BAND_COL_ELEM macro)
   *    for (j=0; j < n; j++) {
   *       col_j = LAPACK_BAND_COL(Jac,j);
   *       for (k=-mupper; k <= mlower; k++) {
   *         generate J_ij = the (i,j)th Jacobian element, i=j+k
   *         col_j[k] = J_ij;
   *       }
   *     }
   *
   * (3) (without macros - explicit data structure references)
   *     offset = Jac->smu;
   *     for (j=0; j < n; j++) {
   *       col_j = ((Jac->data)[j])+offset;
   *       for (k=-mupper; k <= mlower; k++) {
   *         generate J_ij = the (i,j)th Jacobian element, i=j+k
   *         col_j[k] = J_ij;
   *       }
   *     }
   * Caution: J->smu is generally NOT the same as mupper.
   *
   * The LAPACK_BAND_ELEM(A,i,j) macro is appropriate for use in small
   * problems in which efficiency of access is NOT a major concern.
   *
   * u is the current iterate (unscaled) [input]
   *
   * fu is a vector (type N_Vector) containing result of nonlinear
   *     system function evaluated at current iterate:
   *     fu = F(uu) [input]
   *
   * jac_data is a pointer to user data - the same as the jac_data
   *          parameter passed to KINLapackSetJacFn.
   *
   * tmp1, tmp2  available scratch vectors (volatile storage)
   *
   * If successful, the function should return 0 (zero). If an error
   * occurs, then the routine should return a non-zero integer value.
   * -----------------------------------------------------------------
   */
  
  typedef int (*KINLapackBandJacFn)(int N, int mupper, int mlower,
                                    N_Vector u, N_Vector fu, 
                                    LapackMat Jac, void *jac_data,
                                    N_Vector tmp1, N_Vector tmp2);

  /*
   * =================================================================
   *            E X P O R T E D    F U N C T I O N S 
   * =================================================================
   */

  /*
   * -----------------------------------------------------------------
   * Function : KINLapackDense
   * -----------------------------------------------------------------
   * A call to the KINLapackDense function links the main solver
   * with the KINLAPACK linear solver using dense Jacobians.
   *
   * kinmem is the pointer to the solver memory returned by KINCreate.
   *
   * N is the size of the ODE system.
   *
   * The return value of KINLapackDense is one of:
   *    KINLAPACK_SUCCESS   if successful
   *    KINLAPACK_MEM_NULL  if the KINSOL memory was NULL
   *    KINLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    KINLAPACK_ILL_INPUT if a required vector operation is missing
   * -----------------------------------------------------------------
   */

  int KINLapackDense(void *kinmem, int N);

  /*
   * -----------------------------------------------------------------
   * Function : KINLapackBand
   * -----------------------------------------------------------------
   * A call to the KINLapackBand function links the main solver
   * with the KINLAPACK linear solver using banded Jacobians. 
   *
   * kinmem is the pointer to the solver memory returned by KINCreate.
   *
   * N is the size of the ODE system.
   *
   * mupper is the upper bandwidth of the band Jacobian approximation.
   *
   * mlower is the lower bandwidth of the band Jacobian approximation.
   *
   * The return value of KINLapackBand is one of:
   *    KINLAPACK_SUCCESS   if successful
   *    KINLAPACK_MEM_NULL  if the KINSOL memory was NULL
   *    KINLAPACK_MEM_FAIL  if there was a memory allocation failure
   *    KINLAPACK_ILL_INPUT if a required vector operation is missing or
   *                       if a bandwidth has an illegal value.
   * -----------------------------------------------------------------
   */

  int KINLapackBand(void *kinmem, int N, int mupper, int mlower);

  /*
   * -----------------------------------------------------------------
   * Optional inputs to the KINLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * KINLapackSetJacFn specifies the Jacobian approximation routine 
   * to be used. When using dense Jacobians, a user-supplied jac 
   * routine must be of type KINLapackDenseJacFn. When using banded 
   * Jacobians, a user-supplied jac routine must be of type 
   * KINLapackBandJacFn.
   * By default, a difference quotient approximation, supplied with 
   * this solver is used.
   * KINLapackSetJacFn also specifies a pointer to user data which is 
   * passed to the user's jac routine every time it is called.
   *
   * The return value of KINLapackSetJacFn is one of:
   *    KINLAPACK_SUCCESS   if successful
   *    KINLAPACK_MEM_NULL  if the KINSOL memory was NULL
   *    KINLAPACK_LMEM_NULL if the KINLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int KINLapackSetJacFn(void *kinmem, void *jac, void *jac_data);

  /*
   * -----------------------------------------------------------------
   * Optional outputs from the KINLAPACK linear solver
   * -----------------------------------------------------------------
   *
   * KINLapackGetWorkSpace returns the real and integer workspace used
   *                     by KINLAPACK.
   * KINLapackGetNumJacEvals returns the number of calls made to the
   *                     Jacobian evaluation routine jac.
   * KINLapackGetNumFctEvals returns the number of calls to the user
   *                     f routine due to finite difference Jacobian
   *                     evaluation.
   * KINLapackGetLastFlag returns the last error flag set by any of
   *                     the KINLAPACK interface functions.
   *
   * The return value of KINLapackGet* is one of:
   *    KINLAPACK_SUCCESS   if successful
   *    KINLAPACK_MEM_NULL  if the KINSOL memory was NULL
   *    KINLAPACK_LMEM_NULL if the KINLAPACK memory was NULL
   * -----------------------------------------------------------------
   */

  int KINLapackGetWorkSpace(void *kinmem, long int *lenrwLS, long int *leniwLS);
  int KINLapackGetNumJacEvals(void *kinmem, long int *njevals);
  int KINLapackGetNumFctEvals(void *kinmem, long int *nfevalsLS);
  int KINLapackGetLastFlag(void *kinmem, int *flag);

  /*
   * -----------------------------------------------------------------
   * The following function returns the name of the constant 
   * associated with a KINLAPACK return flag
   * -----------------------------------------------------------------
   */

  char *KINLapackGetReturnFlagName(int flag);


#ifdef __cplusplus
}
#endif

#endif
