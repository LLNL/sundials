/*******************************************************************
 *                                                                 *
 * File          : cvband.h                                        *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and          *
 *                 Radu Serban  @ LLNL                             *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * See sundials/cvode/LICENSE or sundials/cvodes/LICENSE           *
 *-----------------------------------------------------------------*
 * This is the header file for the CVODE/CVODES band linear        *
 * solver, CVBAND.                                                 *
 *                                                                 *
 *******************************************************************/
 
#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvband_h
#define _cvband_h

#include <stdio.h>
#include "sundialstypes.h"
#include "band.h"
#include "nvector.h"
 
/******************************************************************
 *                                                                *
 * CVBAND solver constants                                        *
 *----------------------------------------------------------------*
 * CVB_MSBJ  : maximum number of steps between band Jacobian      *
 *             evaluations                                        *
 *                                                                *
 * CVB_DGMAX : maximum change in gamma between band Jacobian      *
 *             evaluations                                        *
 *                                                                *
 ******************************************************************/

#define CVB_MSBJ  50  

#define CVB_DGMAX RCONST(0.2)  
 
/******************************************************************
 *                                                                *           
 * Type : CVBandJacFn                                             *
 *----------------------------------------------------------------*
 * A band Jacobian approximation function Jac must have the       *
 * prototype given below. Its parameters are:                     *
 *                                                                *
 * N is the length of all vector arguments.                       *
 *                                                                *
 * mupper is the upper half-bandwidth of the approximate banded   *
 * Jacobian. This parameter is the same as the mupper parameter   *
 * passed by the user to the CVBand function.                     *
 *                                                                *
 * mlower is the lower half-bandwidth of the approximate banded   *
 * Jacobian. This parameter is the same as the mlower parameter   *
 * passed by the user to the CVBand function.                     *
 *                                                                *
 * J is the band matrix (of type BandMat) that will be loaded     *
 * by a CVBandJacFn with an approximation to the Jacobian matrix  *
 * J = (df_i/dy_j) at the point (t,y).                            *
 * J is preset to zero, so only the nonzero elements need to be   *
 * loaded. Three efficient ways to load J are:                    *
 *                                                                *
 * (1) (with macros - no explicit data structure references)      *
 *    for (j=0; j < n; j++) {                                     *
 *       col_j = BAND_COL(J,j);                                   *
 *       for (i=j-mupper; i <= j+mlower; i++) {                   *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)     *
 *    for (j=0; j < n; j++) {                                     *
 *       col_j = BAND_COL(J,j);                                   *
 *       for (k=-mupper; k <= mlower; k++) {                      *
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k    *
 *         col_j[k] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *  
 * (3) (without macros - explicit data structure references)      *
 *     offset = J->smu;                                           *
 *     for (j=0; j < n; j++) {                                    *
 *       col_j = ((J->data)[j])+offset;                           *
 *       for (k=-mupper; k <= mlower; k++) {                      *
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k    *
 *         col_j[k] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 * Caution: J->smu is generally NOT the same as mupper.           *
 *                                                                *
 * The BAND_ELEM(A,i,j) macro is appropriate for use in small     *
 * problems in which efficiency of access is NOT a major concern. *
 *                                                                *
 * t is the current value of the independent variable.            *
 *                                                                *
 * y is the current value of the dependent variable vector,       *
 *      namely the predicted value of y(t).                       *
 *                                                                *
 * fy is the vector f(t,y).                                       *
 *                                                                *
 * jac_data is a pointer to user data - the same as the jac_data  *
 *          parameter passed to CVBand.                           *
 *                                                                *
 * NOTE: If the user's Jacobian routine needs other quantities,   *
 *     they are accessible as follows: hcur (the current stepsize)*
 *     and ewt (the error weight vector) are accessible through   *
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively   *
 *     (see cvode.h). The unit roundoff is available as           *
 *     UNIT_ROUNDOFF defined in sundialstypes.h                   *
 *                                                                *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for      *
 * vectors of length N which can be used by a CVBandJacFn         *
 * as temporary storage or work space.                            *
 *                                                                *
 ******************************************************************/
  
typedef void (*CVBandJacFn)(long int N, long int mupper, long int mlower,
                            BandMat J, realtype t,
                            N_Vector y, N_Vector fy, void *jac_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 
 
/******************************************************************
 *                                                                *
 * Function : CVBand                                              *
 *----------------------------------------------------------------*
 * A call to the CVBand function links the main CVODE integrator  *
 * with the CVBAND linear solver.                                 *
 *                                                                *
 * cvode_mem is the pointer to the integrator memory returned by  *
 *              CVodeCreate.                                      *
 *                                                                *
 * N is the size of the ODE system.                               *
 *                                                                *
 * mupper is the upper bandwidth of the band Jacobian             *
 *           approximation.                                       *
 *                                                                *
 * mlower is the lower bandwidth of the band Jacobian             *
 *           approximation.                                       *
 *                                                                *
 * The return values of CVBand are:                               *
 *   SUCCESS       = 0  if successful                             *
 *   LMEM_FAIL     = -1 if there was a memory allocation failure. *
 *   LIN_ILL_INPUT = -2 if there was an illegal input.            *
 *                                                                *
 ******************************************************************/

int CVBand(void *cvode_mem, long int N, 
           long int mupper, long int mlower);

/******************************************************************
 * Optional inputs to the CVBAND linear solver                    *
 *----------------------------------------------------------------*
 *                                                                *
 * CVBandSetJacFn specifies the band Jacobian approximation       *
 *         routine to be used. A user-supplied bjac routine must  *
 *         be of type CVBandJacFn.                                *
 *         By default, a difference quotient routine CVBandDQJac, *
 *         supplied with this solver is used.                     *
 * CVBandSetJacData specifies a pointer to user data which is     *
 *         passed to the bjac routine every time it is called.    *
 *                                                                *
 ******************************************************************/

int CVBandSetJacFn(void *cvode_mem, CVBandJacFn bjac);
int CVBandSetJacData(void *cvode_mem, void *jac_data);

/******************************************************************
 * Optional outputs from the CVBAND linear solver                 *
 *----------------------------------------------------------------*
 *                                                                *
 * CVBandGetIntWorkSpace returns the integer workspace used by    *
 *     CVBAND.                                                    *
 * CVBandGetRealWorkSpace returns the real workspace used by      *
 *     CVBAND.                                                    *
 * CVBandGetNumJacEvals returns the number of calls made to the   *
 *     Jacobian evaluation routine bjac.                          *
 * CVBandGetNumRhsEvals returns the number of calls to the user   *
 *     f routine due to finite difference Jacobian evaluation.    *
 *                                                                *
 ******************************************************************/

int CVBandGetIntWorkSpace(void *cvode_mem, long int *leniwB);
int CVBandGetRealWorkSpace(void *cvode_mem, long int *lenrwB);
int CVBandGetNumJacEvals(void *cvode_mem, long int *njevalsB);
int CVBandGetNumRhsEvals(void *cvode_mem, long int *nfevalsB);


/******************************************************************
 *                                                                *           
 * Types : CVBandMemRec, CVBandMem                                *
 *----------------------------------------------------------------*
 * The type CVBandMem is pointer to a CVBandMemRec. This          *
 * structure contains CVBand solver-specific data.                *
 *                                                                *
 ******************************************************************/

typedef struct {

  long int b_n;           /* N = problem dimension                    */

  CVBandJacFn b_jac;      /* jac = Jacobian routine to be called      */

  long int b_ml;          /* b_ml = lower bandwidth of savedJ         */
  
  long int b_mu;          /* b_mu = upper bandwidth of savedJ         */ 
  
  long int b_storage_mu;  /* upper bandwith of M = MIN(N-1,b_mu+b_ml) */
  
  BandMat b_M;            /* M = I - gamma J, gamma = h / l1          */
  
  long int *b_pivots;     /* pivots = pivot array for PM = LU         */
  
  BandMat b_savedJ;       /* savedJ = old Jacobian                    */
  
  long int b_nstlj;       /* nstlj = nst at last Jacobian eval.       */
  
  long int b_nje;         /* nje = no. of calls to jac                */
  
  long int b_nfeB;        /* nfeB = no. of calls to f due to difference
                             quotient band Jacobian approximation     */

  void *b_J_data;         /* J_data is passed to jac                  */
  
} CVBandMemRec, *CVBandMem;


#endif

#ifdef __cplusplus
}
#endif
