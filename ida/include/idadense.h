/******************************************************************
 *                                                                *
 * File          : idadense.h                                     *
 * Programmers   : Alan C. Hindmarsh, Allan G. Taylor, and        *
 *                 Radu Serban @LLNL                              *
 * Version of    : 6 March 2002                                   *
 *----------------------------------------------------------------*
 * This is the header file for the IDA dense linear solver        *
 * module, IDADENSE.                                              *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idadense_h
#define _idadense_h


#include <stdio.h>
#include "ida.h"
#include "llnltyps.h"
#include "dense.h"
#include "nvector.h"

/* Return values for IDADense: */

/* SUCCESS = 0 (defined in ida.h) */
enum {IDA_DENSE_FAIL = -1};

 
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
 * iopt[DENSE_LRW] : size (in real words) of real workspace       *
 *                   matrices and vectors used by this module.    *
 *                                                                *
 * iopt[DENSE_LIW] : size (in integer words) of integer           *
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
 * constraints is the vector of inequality constraint options     *
 *             (as passed to IDAMalloc).  Included here to allow  *
 *             for checking of incremented y values in difference *
 *             quotient calculations.                             *
 *                                                                *
 * res    is the residual function for the DAE problem.           *
 *                                                                *
 * rdata  is a pointer to user data to be passed to res, the same *
 *        as the rdata parameter passed to IDAMalloc.             *
 *                                                                *
 * jdata  is a pointer to user Jacobian data - the same as the    *
 *        jdata parameter passed to IDADense.                     *
 *                                                                *
 * resvec is the residual vector F(tt,yy,yp).                     *
 *                                                                *
 * ewt    is the error weight vector.                             *
 *                                                                *
 * hh     is a tentative step size in t.                          *
 *                                                                *
 * uround is the machine unit roundoff.                           *
 *                                                                *
 * JJ     is the dense matrix (of type DenseMat) to be loaded by  *
 *        an IDADenseJacFn routine with an approximation to the   *
 *        system Jacobian matrix                                  *
 *              J = dF/dy + cj*dF/dy'                             *
 *        at the given point (t,y,y'), where the DAE system is    *
 *        given by F(t,y,y') = 0.  JJ is preset to zero, so only  *
 *        the nonzero elements need to be loaded.  See note below.*
 *                                                                *
 * nrePtr is a pointer to the memory location containing the      *
 *        IDA problem data nre = number of calls to res. This     *
 *        Jacobian routine should update this counter by adding   *
 *        on the number of res calls it makes in order to         *
 *        approximate the Jacobian, if any.  For example, if this *
 *        routine calls res a total of Neq times, then it should  *
 *        perform the update *nrePtr += Neq.                      *
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
  
typedef int (*IDADenseJacFn)(integer Neq, real tt, N_Vector yy, N_Vector yp,
                             real cj, N_Vector constraints, ResFn res, void *rdata,
                             void *jdata, N_Vector resvec, N_Vector ewt, real hh,
                             real uround, DenseMat JJ, long int *nrePtr,
                             N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);



 
/******************************************************************
 *                                                                *
 * Function : IDADense                                            *
 *----------------------------------------------------------------*
 * A call to the IDADense function links the main IDA integrator  *
 * with the IDADENSE linear solver module.                        *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by IDAMalloc.    *
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
 * IDADense returns either                                        *
 *     SUCCESS = 0         if successful, or                      *
 *     IDA_DENSE_FAIL = -1 if either IDA_mem was null, or a       *
 *                         malloc failure occurred.               *
 *                                                                *
 * NOTE: The dense linear solver assumes a serial implementation  *
 *       of the NVECTOR package. Therefore, IDADense will first   *
 *       test for a compatible N_Vector internal representation   *
 *       by checking (1) the machine environment ID tag and       *
 *       (2) that the functions N_VMake, N_VDispose, N_VGetData,  *
 *       and N_VSetData are implemented.                          *
 *                                                                *
 ******************************************************************/

int IDADense(void *IDA_mem, IDADenseJacFn djac, void *jdata);


#endif

#ifdef __cplusplus
}
#endif
