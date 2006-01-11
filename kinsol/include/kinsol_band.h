/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:59 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the KINSOL band linear solver, KINBAND.
 * -----------------------------------------------------------------
 */

#ifndef _KINBAND_H
#define _KINBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "sundials_band.h"
#include "sundials_nvector.h"

/*
 * -----------------------------------------------------------------
 * Type : KINBandJacFn
 * -----------------------------------------------------------------
 * A band Jacobian approximation function Jac must have the
 * prototype given below. Its parameters are:
 *
 * N is the problem size.
 *
 * mupper is the upper half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mupper parameter
 * passed by the user to the KINBand function.
 *
 * mlower is the lower half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mlower parameter
 * passed by the user to the KINBand function.
 *
 * J is the band matrix (of type BandMat) that will be loaded
 * by a KINBandJacFn with an approximation to the Jacobian matrix
 * J = (df_i/dy_j).
 * J is preset to zero, so only the nonzero elements need to be
 * loaded. Three efficient ways to load J are:
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
 * u is the current iterate (unscaled) [input]
 *
 * fu is a vector (type N_Vector) containing result of nonlinear
 *     system function evaluated at current iterate:
 *     fu = F(uu) [input]
 *
 * jac_data is a pointer to user data - the same as the jac_data
 *          parameter passed to KINBandSetJacFn.
 *
 * tmp1, tmp2  available scratch vectors (volatile storage)
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINBandJacFn)(long int N, long int mupper, long int mlower,
                            BandMat J, N_Vector u, N_Vector fu, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2);

/*
 * -----------------------------------------------------------------
 * Function : KINBand
 * -----------------------------------------------------------------
 * A call to the KINBand function links the main solver
 * with the KINBAND linear solver.
 *
 * kinmem is the pointer to the integrator memory returned by
 *           KINCreate.
 *
 * N is the problem size
 *
 * mupper is the upper bandwidth of the band Jacobian
 *
 * mlower is the lower bandwidth of the band Jacobian
 *
 * The return value of KINBand is one of:
 *    KINBAND_SUCCESS   if successful
 *    KINBAND_MEM_NULL  if the kinsol memory was NULL
 *    KINBAND_MEM_FAIL  if there was a memory allocation failure
 *    KINBAND_ILL_INPUT if a required vector operation is missing or
 *                      if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

int KINBand(void *kinmem, long int N,
           long int mupper, long int mlower);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINBAND linear solver
 * -----------------------------------------------------------------
 *
 * KINBandSetJacFn specifies the band Jacobian routine to be used. 
 *                A user-supplied bjac routine must be of type 
 *                KINBandJacFn. By default, a difference quotient 
 *                routine KINBandDQJac, supplied with this solver 
 *                is used.
 *                It also specifies a pointer to user data which is
 *                passed to the bjac routine every time it is called.
 *
 * The return value of KINBandSetJacFn is one of:
 *    KINBAND_SUCCESS   if successful
 *    KINBAND_MEM_NULL  if the kinsol memory was NULL
 *    KINBAND_LMEM_NULL if the kinband memory was NULL
 * -----------------------------------------------------------------
 */

int KINBandSetJacFn(void *kinmem, KINBandJacFn bjac, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the KINBAND linear solver
 * -----------------------------------------------------------------
 *
 * KINBandGetWorkSpace returns the real and integer workspace used
 *                    by KINBAND.
 * KINBandGetNumJacEvals returns the number of calls made to the
 *                      Jacobian evaluation routine bjac.
 * KINBandGetNumFuncEvals returns the number of calls to the user
 *                      f routine due to finite difference Jacobian
 *                      evaluation.
 * KINBandGetLastFlag returns the last error flag set by any of
 *                   the KINBAND interface functions.
 *
 * The return value of KINBandGet* is one of:
 *    KINBAND_SUCCESS   if successful
 *    KINBAND_MEM_NULL  if the kinsol memory was NULL
 *    KINBAND_LMEM_NULL if the kinband memory was NULL
 * -----------------------------------------------------------------
 */

int KINBandGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB);
int KINBandGetNumJacEvals(void *kinmem, long int *njevalsB);
int KINBandGetNumFuncEvals(void *kinmem, long int *nfevalsB);
int KINBandGetLastFlag(void *kinmem, int *flag);

/* KINBAND return values */

#define KINBAND_SUCCESS    0
#define KINBAND_MEM_NULL  -1
#define KINBAND_LMEM_NULL -2
#define KINBAND_ILL_INPUT -3
#define KINBAND_MEM_FAIL  -4

#ifdef __cplusplus
}
#endif

#endif
