/*******************************************************************
 * File          : idasdense.h                                     *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 12 August 2003                                  *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/idas/LICENSE                          *
 *-----------------------------------------------------------------*
 * This is the header file for the IDAS dense linear solver        *
 * module, IDADENSE.                                               *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idasdense_h
#define _idasdense_h


#include <stdio.h>
#include "idas.h"
#include "sundialstypes.h"
#include "dense.h"
#include "nvector.h"

 
/******************************************************************
 *                                                                *           
 * Type : IDADenseJacFn                                           *
 *----------------------------------------------------------------*
 * A dense Jacobian approximation function djac must have the     *
 * prototype given below. Its parameters are:                     *
 *                                                                *
 * Neq is the problem size, and length of all vector arguments.   *
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
 *        jdata parameter passed to IDADense.                     *
 *                                                                *
 * resvec is the residual vector F(tt,yy,yp).                     *
 *                                                                *
 * Jac    is the dense matrix (of type DenseMat) to be loaded by  *
 *        an IDADenseJacFn routine with an approximation to the   *
 *        system Jacobian matrix                                  *
 *              J = dF/dy + c_j*dF/dy'                            *
 *        at the given point (t,y,y'), where the DAE system is    *
 *        given by F(t,y,y') = 0.  Jac is preset to zero, so only *
 *        the nonzero elements need to be loaded.  See note below.*
 *                                                                *
 * tmp1, tmp2, tmp3 are pointers to memory allocated for          *
 *        N_Vectors which can be used by an IDADenseJacFn routine *
 *        as temporary storage or work space.                     *
 *                                                                *
 * NOTE: The following are two efficient ways to load JJ:         *
 * (1) (with macros - no explicit data structure references)      *
 *     for (j=0; j < Neq; j++) {                                  *
 *       col_j = DENSE_COL(JJ,j);                                 *
 *       for (i=0; i < Neq; i++) {                                *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         col_j[i] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 * (2) (without macros - explicit data structure references)      *
 *     for (j=0; j < Neq; j++) {                                  *
 *       col_j = (JJ->data)[j];                                   *
 *       for (i=0; i < Neq; i++) {                                *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         col_j[i] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   *
 * efficient in general.  It is only appropriate for use in small *
 * problems in which efficiency of access is NOT a major concern. *
 *                                                                *
 * NOTE: If the user's Jacobian routine needs other quantities,   *
 *     they are accessible as follows: hcur (the current stepsize)*
 *     and ewt (the error weight vector) are accessible through   *
 *     IDAGetCurrentStep and IDAGetErrWeights, respectively (see  *
 *     ida.h). The unit roundoff is available through a call to   *
 *     UnitRoundoff.                                              *
 *                                                                *
 * The IDADenseJacFn should return                                *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDAS will attempt   *
 * to recover by reducing the stepsize (which changes cj).        *
 ******************************************************************/
  
typedef int (*IDADenseJacFn)(integertype Neq, realtype tt, N_Vector yy, 
                             N_Vector yp, realtype c_j, void *jdata, 
                             N_Vector resvec, DenseMat Jac, 
                             N_Vector tempv1, N_Vector tempv2, 
                             N_Vector tempv3);

/******************************************************************
 *                                                                *
 * Function : IDADense                                            *
 *----------------------------------------------------------------*
 * A call to the IDADense function links the main IDAS integrator *
 * with the IDADENSE linear solver module.                        *
 *                                                                *
 * ida_mem is the pointer to IDAS memory returned by IDACreate.   *
 *                                                                *
 * Neq  is the problem size                                       *
 *                                                                *
 * IDADense returns:                                              *
 *     SUCCESS   = 0   if successful                              *
 *     LMEM_FAIL = -1  if there was a memory allocation failure   *
 *     LIN_ILL_INPUT = -2 if NVECTOR found incompatible           *
 *                                                                *
 * NOTE: The dense linear solver assumes a serial implementation  *
 *       of the NVECTOR package. Therefore, IDADense will first   *
 *       test for a compatible N_Vector internal representation   *
 *       by checking (1) the machine environment ID tag and       *
 *       (2) that the functions N_VMake, N_VDispose, N_VGetData,  *
 *       and N_VSetData are implemented.                          *
 *                                                                *
 ******************************************************************/

int IDADense(void *ida_mem, integertype Neq); 

/******************************************************************
 * Optional inputs to the IDADENSE linear solver                  *
 *----------------------------------------------------------------*
 *                                                                *
 * IDADenseSetJacFn specifies the dense Jacobian approximation    *
 *        routine to be used. A user-supplied djac routine must   *
 *        be of type IDADenseJacFn.                               *
 *        By default, a difference quotient routine IDADenseDQJac,*
 *        supplied with this solver is used.                      *
 * IDADenseSetJacData specifies a pointer to user data which is   *
 *        passed to the djac routine every time it is called.     *
 *                                                                *
 ******************************************************************/

int IDADenseSetJacFn(void *ida_mem, IDADenseJacFn djac);
int IDADenseSetJacData(void *ida_mem, void *jdata);
 
/******************************************************************
 * Optional outputs from the IDADENSE linear solver               *
 *----------------------------------------------------------------*
 *                                                                *
 * IDADenseGetIntWorkSpace returns the integer workspace used by  *
 *     IDADENSE.                                                  *
 * IDADenseGetRealWorkSpace returns the real workspace used by    *
 *     IDADENSE.                                                  *
 * IDADenseGetNumJacEvals returns the number of calls made to the *
 *     Jacobian evaluation routine djac.                          *
 * IDADenseGetNumResEvals returns the number of calls to the user *
 *     res routine due to finite difference Jacobian evaluation.  *
 *                                                                *
 ******************************************************************/

int IDADenseGetIntWorkSpace(void *ida_mem, long int *leniwD);
int IDADenseGetRealWorkSpace(void *ida_mem, long int *lenrwD);
int IDADenseGetNumJacEvals(void *ida_mem, int *njevalsD);
int IDADenseGetNumResEvals(void *ida_mem, int *nrevalsD);

/******************************************************************
 *                                                                *           
 * Types : IDADenseMemRec, IDADenseMem                            *
 *----------------------------------------------------------------*
 * The type IDADenseMem is pointer to an IDADenseMemRec. This     *
 * structure contains IDADense solver-specific data.              *
 *                                                                *
 ******************************************************************/

typedef struct {

  integertype d_neq;     /* Neq = problem dimension              */

  IDADenseJacFn d_jac;   /* jac = Jacobian routine to be called  */
  
  DenseMat d_J;          /* J = dF/dy + cj*dF/dy'                */
  
  integertype *d_pivots; /* pivots = pivot array for PJ = LU     */
  
  int d_nje;             /* nje = no. of calls to jac            */
  
  int d_nreD;            /* nreD = no. of calls to res due to 
                            difference quotient Jacobian evaluation */

  void *d_jdata;         /* jdata is passed to jac               */

} IDADenseMemRec, *IDADenseMem;

#endif

#ifdef __cplusplus
}
#endif
