/*******************************************************************
 * File          : idasband.h                                      *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/idas/LICENSE                          *
 *-----------------------------------------------------------------*
 * This is the header file for the IDAS band linear solver         *
 * module, IDABAND. It interfaces between the band module and the  *
 * IDAS package when a banded linear solver is appropriate.        *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idasband_h
#define _idasband_h


#include <stdio.h>
#include "idas.h"
#include "sundialstypes.h"
#include "band.h"
#include "nvector.h"

 
/******************************************************************
 *                                                                *           
 * Type : IDABandJacFn                                            *
 *----------------------------------------------------------------*
 * A banded Jacobian approximation function bjac must have the    *
 * prototype given below. Its parameters are:                     *
 *                                                                *
 * Neq is the problem size, and length of all vector arguments.   *
 *                                                                *
 * mupper is the upper bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * mlower is the lower bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * tt  is the current value of the independent variable t.        *
 *                                                                *
 * yy  is the current value of the dependent variable vector,     *
 *        namely the predicted value of y(t).                     *
 *                                                                *
 * yp  is the current value of the derivative vector y',          *
 *        namely the predicted value of y'(t).                    *
 *                                                                *
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * jdata  is a pointer to user Jacobian data - the same as the    *
 *        jdata parameter passed to IDABand.                      *
 *                                                                *
 * resvec is the residual vector F(tt,yy,yp).                     *
 *                                                                *
 * Jac    is the band matrix (of type BandMat) to be loaded by    *
 *        an IDABandJacFn routine with an approximation to the    *
 *        system Jacobian matrix                                  *
 *              J = dF/dy + cj*dF/dy'                             *
 *        at the given point (t,y,y'), where the DAE system is    *
 *        given by F(t,y,y') = 0.  Jac is preset to zero, so only *
 *        the nonzero elements need to be loaded.  See note below.*
 *                                                                *
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          *
 *        N_Vectors which can be used by an IDABandJacFn routine  *
 *        as temporary storage or work space.                     *
 *                                                                *
 *                                                                *
 * NOTE: The following are two efficient ways to load JJ:         *
 *                                                                *
 * (1) (with macros - no explicit data structure references)      *
 *    for (j=0; j < Neq; j++) {                                   *
 *       col_j = BAND_COL(JJ,j);                                  *
 *       for (i=j-mupper; i <= j+mlower; i++) {                   *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)     *
 *    for (j=0; j < Neq; j++) {                                   *
 *       col_j = BAND_COL(JJ,j);                                  *
 *       for (k=-mupper; k <= mlower; k++) {                      *
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k    *
 *         col_j[k] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *
 * NOTE: If the user's Jacobian routine needs other quantities,   *
 *     they are accessible as follows: hcur (the current stepsize)*
 *     and ewt (the error weight vector) are accessible through   *
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively (see  *
 *     ida.h). The unit roundoff is available as                  *
 *     UNIT_ROUNDOFF defined in sundialstypes.h                   *
 *                                                                *
 * A third way, using the BAND_ELEM(A,i,j) macro, is much less    *
 * efficient in general.  It is only appropriate for use in small *
 * problems in which efficiency of access is NOT a major concern. *
 *                                                                *
 * The IDABandJacFn should return                                 *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDAS will attempt   *
 * to recover by reducing the stepsize (which changes cj).        *
 ******************************************************************/
  
typedef int (*IDABandJacFn)(long int Neq, long int mupper, 
                            long int mlower, realtype tt, 
                            N_Vector yy, N_Vector yp, realtype c_j, 
                            void *jdata, N_Vector resvec, BandMat Jac, 
                            N_Vector tmp1, N_Vector tmp2, 
                            N_Vector tmp3);
 
/******************************************************************
 *                                                                *
 * Function : IDABand                                             *
 *----------------------------------------------------------------*
 * A call to the IDABand function links the main IDAS integrator  *
 * with the IDABAND linear solver module.                         *
 *                                                                *
 * ida_mem is the pointer to IDAS memory returned by IDACreate.   *
 *                                                                *
 * mupper is the upper bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * mlower is the lower bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * The return values of IDABand are:                              *
 *    SUCCESS       = 0  if successful                            *
 *    LMEM_FAIL     = -1 if there was a memory allocation failure *
 *    LIN_ILL_INPUT = -2 if the input was illegal or NVECTOR bad. *
 *                                                                *
 * NOTE: The band linear solver assumes a serial implementation   *
 *       of the NVECTOR package. Therefore, IDABand will first    *
 *       test for a compatible N_Vector internal representation   *
 *       by checking (1) the machine environment ID tag and       *
 *       (2) that the functions N_VMake, N_VDispose, N_VGetData,  *
 *       and N_VSetData are implemented.                          *
 *                                                                *
 ******************************************************************/

int IDABand(void *ida_mem, long int Neq, 
            long int mupper, long int mlower);

/******************************************************************
 * Optional inputs to the IDABAND linear solver                   *
 *----------------------------------------------------------------*
 *                                                                *
 * IDABandSetJacFn specifies the dense Jacobian approximation     *
 *         routine to be used. A user-supplied djac routine must  *
 *         be of type IDABandJacFn.                               *
 *         By default, a difference quotient routine IDABandDQJac,*
 *         supplied with this solver is used.                     *
 * IDABandSetJacData specifies a pointer to user data which is    *
 *         passed to the bjac routine every time it is called.    *
 *                                                                *
 ******************************************************************/

int IDABandSetJacFn(void *ida_mem, IDABandJacFn bjac);
int IDABandSetJacData(void *ida_mem, void *jdata);

/******************************************************************
 * Optional outputs from the IDABAND linear solver                *
 *----------------------------------------------------------------*
 *                                                                *
 * IDABandGetIntWorkSpace returns the integer workspace used by   *
 *     IDABAND.                                                   *
 * IDABandGetRealWorkSpace returns the real workspace used by     *
 *     IDABAND.                                                   *
 * IDABandGetNumJacEvals returns the number of calls made to the  *
 *     Jacobian evaluation routine bjac.                          *
 * IDABandGetNumResEvals returns the number of calls to the user  *
 *     res routine due to finite difference Jacobian evaluation.  *
 *                                                                *
 ******************************************************************/

int IDABandGetIntWorkSpace(void *ida_mem, long int *leniwB);
int IDABandGetRealWorkSpace(void *ida_mem, long int *lenrwB);
int IDABandGetNumJacEvals(void *ida_mem, long int *njevalsB);
int IDABandGetNumResEvals(void *ida_mem, long int *nrevalsB);

/******************************************************************
 *                                                                *           
 * Types : IDABandMemRec, IDABandMem                              *
 *----------------------------------------------------------------*
 * The type IDABandMem is pointer to an IDABandMemRec. This       *
 * structure contains IDABand solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {

  long int b_neq;           /* Neq = problem size                           */

  IDABandJacFn b_jac;       /* jac = banded Jacobian routine to be called   */
  
  BandMat b_J;              /* J = dF/dy + cj*dF/dy', banded approximation. */
  
  long int b_mupper;        /* mupper = upper bandwidth of Jacobian matrix. */
  
  long int b_mlower;        /* mlower = lower bandwidth of Jacobian matrix. */
  
  long int b_storage_mu;    /* storage_mu = upper bandwidth with storage for
                               factoring = min(Neq-1, mupper+mlower).       */
  
  long int *b_pivots;       /* pivots = pivot array for PJ = LU             */
  
  long int b_nje;           /* nje = no. of calls to jac                    */
  
  long int b_nreB;          /* nreB = no. of calls to res due to 
                               difference quotient Jacobian evaluation      */

  void *b_jdata;            /* jdata = data structure required by jac.      */
  
} IDABandMemRec, *IDABandMem;

#endif

#ifdef __cplusplus
}
#endif
