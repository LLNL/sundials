/*******************************************************************
 *                                                                 *
 * File          : idasdense.h                                     *
 * Programmers   : Alan C. Hindmarsh, Allan G. Taylor, and         *
 *                 Radu Serban @LLNL                               *
 * Version of    : 31 March 2003                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/ida/LICENSE                           *
 *-----------------------------------------------------------------*
 * This is the header file for the IDA dense linear solver         *
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
 * IDADENSE solver optional output indices                        *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * IDADENSE optional output. The symbolic names are used as       *
 * indices into the iopt and ropt arrays passed to IDAMalloc.     *
 * The IDADENSE optional outputs are:                             *
 *                                                                *
 * iopt[DENSE_NJE] : number of Jacobian evaluations, i.e. of      *
 *                   calls made to the dense Jacobian routine     *
 *                   (default or user-supplied).                  *
 *                                                                *
 * iopt[DENSE_LRW] : size (in realtype words) of real workspace   *
 *                   matrices and vectors used by this module.    *
 *                                                                *
 * iopt[DENSE_LIW] : size (in integertype words) of integer       *
 *                   workspace vectors used by this module.       *
 *                                                                *
 ******************************************************************/
 
enum { DENSE_NJE=IDA_IOPT_SIZE, DENSE_LRW, DENSE_LIW };


 
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
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * jdata  is a pointer to user Jacobian data - the same as the    *
 *        jdata parameter passed to IDADense.                     *
 *                                                                *
 * resvec is the residual vector F(tt,yy,yp).                     *
 *                                                                *
 * JJ     is the dense matrix (of type DenseMat) to be loaded by  *
 *        an IDADenseJacFn routine with an approximation to the   *
 *        system Jacobian matrix                                  *
 *              J = dF/dy + cj*dF/dy'                             *
 *        at the given point (t,y,y'), where the DAE system is    *
 *        given by F(t,y,y') = 0.  JJ is preset to zero, so only  *
 *        the nonzero elements need to be loaded.  See note below.*
 *                                                                *
 * tempv1, tempv2, tempv3 are pointers to memory allocated for    *
 *        N_Vectors which can be used by an IDADenseJacFn routine *
 *        as temporary storage or work space.                     *
 *                                                                *
 * Note: The following are two efficient ways to load JJ:         *
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
 * The IDADenseJacFn should return                                *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDA will attempt to *
 * recover by reducing the stepsize (which changes cj).           *
 ******************************************************************/
  
typedef int (*IDADenseJacFn)(integertype Neq, realtype tt, N_Vector yy, 
                             N_Vector yp, realtype cj, void *jdata, 
                             N_Vector resvec, DenseMat JJ, 
                             N_Vector tempv1, N_Vector tempv2, 
                             N_Vector tempv3);

 /******************************************************************
 *                                                                *
 * Function : IDADense                                            *
 *----------------------------------------------------------------*
 * A call to the IDADense function links the main IDA integrator  *
 * with the IDADENSE linear solver module.                        *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by IDAMalloc.    *
 *                                                                *
 * neq  is the problem size                                       *
 *                                                                *
 * djac is the dense Jacobian approximation routine to be used.   *
 *         A user-supplied djac routine must be of type           *
 *         IDADenseJacFn (see above).  Pass NULL for djac if IDA  *
 *         is to use the default difference quotient routine      *
 *         IDADenseDQJac supplied with this module.               *
 *                                                                *
 * jdata is a pointer to user data which is passed to the djac    *
 *         routine every time it is called.                       *
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

int IDADense(void *IDA_mem, integertype neq, 
             IDADenseJacFn djac, void *jdata);

 
/******************************************************************
 *                                                                *
 * Function : IDAReInitDense                                      *
 *----------------------------------------------------------------*
 * A call to the IDAReInitDense function resets the link between  *
 * the main IDA integrator and the IDADENSE linear solver.        *
 * After solving one problem using IDADENSE, call IDAReInit and   *
 * then IDAReInitDense to solve another problem of the same size, *
 * if there is a change in the IDADense parameters djac or jdata. *
 * If there is no change in parameters, it is not necessary to    *
 * call either IDAReInitDense or IDADense for the new problem.    *
 *                                                                *
 * All arguments to IDAReInitDense have the same names and        * 
 * meanings as those of IDADense.  The IDA_mem argument must be   *
 * identical to its value in the previous IDADense call.          *
 *                                                                *
 * The return values of IDAReInitDense are:                       *
 *     SUCCESS   = 0   if successful                              *
 *     LMEM_FAIL = -1  if the IDA_mem argument is NULL            *
 *     LIN_ILL_INPUT = -2 if NVECTOR found incompatible           *
 *                                                                *
 ******************************************************************/

int IDAReInitDense(void *IDA_mem, IDADenseJacFn djac, void *jdata);


#endif

#ifdef __cplusplus
}
#endif
