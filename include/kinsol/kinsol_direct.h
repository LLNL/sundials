/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * Common header file for the direct linear solvers in KINSOL.
 * -----------------------------------------------------------------
 */

#ifndef _KINDLS_H
#define _KINDLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *              K I N D I R E C T     C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * KINDLS return values 
 * -----------------------------------------------------------------
 */

#define KINDLS_SUCCESS           0
#define KINDLS_MEM_NULL         -1
#define KINDLS_LMEM_NULL        -2
#define KINDLS_ILL_INPUT        -3
#define KINDLS_MEM_FAIL         -4

/* Additional last_flag values */

#define KINDLS_JACFUNC_ERR      -5
#define KINDLS_SUNMAT_FAIL      -6

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type: KINDlsJacFn
 * -----------------------------------------------------------------
 *
 * A Jacobian approximation function Jac must be of type KINDlsJacFn.
 * Its parameters are:
 *
 * u        - current iterate (unscaled) [input]
 *
 * fu       - vector (type N_Vector) containing result of nonlinear
 *            system function evaluated at current iterate:
 *            fu = F(u) [input]
 *
 * J        - SUNMatrix that will be loaded by a KINDlsJacFn with
 *            an approximation to the Jacobian matrix 
 *            J = (dF_i/dy_j).
 *
 * user_data  - pointer to user data - the same as the user_data
 *              parameter passed to KINSetFdata.
 *
 * tmp1, tmp2 - available scratch vectors (volatile storage)
 *
 * A KINDlsJacFn should return 0 if successful or a non-zero value
 * otherwise.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: See the relevant SUNMatrix implementation header files and
 *     documentation for mechanisms to inquire about matrix 
 *     dimensions, and for efficient ways to set matrix entries.
 *                                                                
 * -----------------------------------------------------------------
 */
  
typedef int (*KINDlsJacFn)(N_Vector u, N_Vector fu, SUNMatrix J,
                           void *user_data, N_Vector tmp1,
                           N_Vector tmp2);
  

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Required inputs to the KINDLS linear solver interface
 * -----------------------------------------------------------------
 *
 * KINDlsSetLinearSolver specifies the direct SUNLinearSolver object
 * that should be used.  This is required if KINSOL is solving a 
 * problem using the Newton or Picard nonlinear solvers with a 
 * direct linear solver  (i.e. not the fixed point or accelerated 
 * fixed-point solvers).
 *
 * The return value is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINDLS_ILL_INPUT if the arguments are incompatible
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINDlsSetLinearSolver(void *kinmem, 
                                          SUNLinearSolver LS,
                                          SUNMatrix A);
/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINDLS linear solver
 * -----------------------------------------------------------------
 *
 * KINDlsSetJacFn specifies the dense/band/sparse Jacobian 
 * approximation routine to be used for a direct linear solver.
 *
 * By default, a difference quotient approximation is used for
 * dense/band; no default exists for sparse (so this must 
 * be user-supplied).
 *
 * The return value is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINDlsSetJacFn(void *kinmem, KINDlsJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from a KINDLS linear solver
 * -----------------------------------------------------------------
 *
 * KINDlsGetWorkSpace    returns the real and integer workspace used
 *                       by the KINDLS linear solver.
 * KINDlsGetNumJacEvals  returns the number of calls made to the
 *                       Jacobian evaluation routine.
 * KINDlsGetNumFuncEvals returns the number of calls to the user's F
 *                       routine due to finite difference Jacobian
 *                       evaluation.
 * KINDlsGetLastFlag     returns the last error flag set by any of
 *                       the KINDLS interface functions.
 * KINDlsGetReturnFlagName returns the name of the constant
 *                         associated with a KINDLS return flag
 *
 * The return value of KINDlsGet* is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the KINSOL memory was NULL
 *    KINDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINDlsGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int KINDlsGetNumJacEvals(void *kinmem, long int *njevals);
SUNDIALS_EXPORT int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevals);
SUNDIALS_EXPORT int KINDlsGetLastFlag(void *kinmem, long int *flag);
SUNDIALS_EXPORT char *KINDlsGetReturnFlagName(long int flag);

#ifdef __cplusplus
}
#endif

#endif
