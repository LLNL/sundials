/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-03-02 17:59:18 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINDENSE linear solver module header file
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KINDENSE_H
#define _KINDENSE_H

#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Type : KINDenseJacFn
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function Jac must have the
 * prototype given below. Its parameters are:
 *
 * N is the problem size.
 *
 * J is the dense matrix (of type DenseMat) that will be loaded
 * by a KINDenseJacFn with an approximation to the Jacobian matrix
 * J = (df_i/dy_j).
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
 * uu   current iterate (unscaled) [input]
 *
 * fval  vector (type N_Vector) containing result of nonliear
 *     system function evaluated at current iterate:
 *     fy = F(y) [input]
 *
 * jac_data is a pointer to user data - the same as the jac_data
 *          parameter passed to KINDense.
 *
 * vtemp1  available scratch vector (volatile storage)
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINDenseJacFn)(long int N, DenseMat J, 
                             N_Vector uu, N_Vector fval, void *jac_data,
                             N_Vector vtemp1, N_Vector vtemp2);


/*
 * -----------------------------------------------------------------
 * Function : KINDense
 * -----------------------------------------------------------------
 * A call to the KINDense function links the main solver with
 * the KINDENSE linear solver.
 *
 * kinmem pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 * N is the problem size
 *
 * The return value of KINDense is one of:
 *    KINDENSE_SUCCESS   if successful
 *    KINDENSE_MEM_NULL  if the kinsol memory was NULL
 *    KINDENSE_MEM_FAIL  if there was a memory allocation failure
 *    KINDENSE_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int KINDense(void *kinmem, long int N);


/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINDENSE linear solver
 * -----------------------------------------------------------------
 *
 * KINDenseSetJacFn specifies the dense Jacobian approximation
 *                 routine to be used. A user-supplied djac routine
 *                 must be of type KINDenseJacFn. By default, a
 *                 difference quotient routine KINDenseDQJac, supplied
 *                 with this solver is used.                     
 * KINDenseSetJacData specifies a pointer to user data which is
 *                   passed to the djac routine every time it is called.
 *
 * The return value of KINDenseSet* is one of:
 *    KINDENSE_SUCCESS   if successful
 *    KINDENSE_MEM_NULL  if the kinsol memory was NULL
 *    KINDENSE_LMEM_NULL if the kindense memory was NULL
 * -----------------------------------------------------------------
 */

int KINDenseSetJacFn(void *kinmem, KINDenseJacFn djac);
int KINDenseSetJacData(void *kinmem, void *jac_data);


/*
 * -----------------------------------------------------------------
 * Optional outputs from the KINDENSE linear solver
 * -----------------------------------------------------------------
 *
 * KINDenseGetWorkSpace returns the real and integer workspace used
 *                     by KINDENSE.
 * KINDenseGetNumJacEvals returns the number of calls made to the
 *                       Jacobian evaluation routine djac.
 * KINDenseGetNumFuncEvals returns the number of calls to the user
 *                       f routine due to finite difference Jacobian
 *                       evaluation.
 * KINDenseGetLastFlag returns the last error flag set by any of
 *                    the KINDENSE interface functions.
 *
 * The return value of KINDenseGet* is one of:
 *    KINDENSE_SUCCESS   if successful
 *    KINDENSE_MEM_NULL  if the kinsol memory was NULL
 *    KINDENSE_LMEM_NULL if the kindense memory was NULL
 * -----------------------------------------------------------------
 */

int KINDenseGetWorkSpace(void *kinmem, long int *lenrwD, long int *leniwD);
int KINDenseGetNumJacEvals(void *kinmem, long int *njevalsD);
int KINDenseGetNumFuncEvals(void *kinmem, long int *nfevalsD);
int KINDenseGetLastFlag(void *kinmem, int *flag);

/* CVDENSE return values */

#define KINDENSE_SUCCESS    0
#define KINDENSE_MEM_NULL  -1
#define KINDENSE_LMEM_NULL -2
#define KINDENSE_ILL_INPUT -3
#define KINDENSE_MEM_FAIL  -4

#endif

#ifdef __cplusplus
}
#endif
