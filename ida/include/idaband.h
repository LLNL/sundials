/******************************************************************
 *                                                                *
 * File          : idaband.h                                      *
 * Programmers   : Allan G. Taylor and Alan C. Hindmarsh          *
 * Version of    : 30 August 1999                                 *
 *----------------------------------------------------------------*
 * This is the header file for the IDA band linear solver         *
 * module, IDABAND. It interfaces between the band module and the *
 * IDA package when a banded linear solver is appropriate.        *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idaband_h
#define _idaband_h


#include <stdio.h>
#include "ida.h"
#include "llnltyps.h"
#include "band.h"
#include "nvector.h"

/* Return values for IDABand: */

/* SUCCESS = 0 (defined in ida.h) */
enum {IDA_BAND_FAIL = -1, IDA_BAND_BAD_ARG = -2};

 
/******************************************************************
 *                                                                *
 * IDABAND solver optional output indices                         *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * IDABAND optional output. The symbolic names are used as        *
 * indices into the iopt and ropt arrays passed to IDAMalloc.     *
 * The IDABAND optional outputs are:                              *
 *                                                                *
 * iopt[BAND_NJE] : number of Jacobian evaluations, i.e. of       *
 *                   calls made to the band Jacobian routine      *
 *                   (default or user-supplied).                  *
 *                                                                *
 * iopt[BAND_LRW] : size (in real words) of real workspace        *
 *                   matrices and vectors used by this module.    *
 *                                                                *
 * iopt[BAND_LIW] : size (in integer words) of integer            *
 *                   workspace vectors used by this module.       *
 *                                                                *
 ******************************************************************/
 
enum { BAND_NJE=IDA_IOPT_SIZE, BAND_LRW, BAND_LIW };


 
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
 *        jdata parameter passed to IDABand.                      *
 *                                                                *
 * resvec is the residual vector F(tt,yy,yp).                     *
 *                                                                *
 * ewt    is the error weight vector.                             *
 *                                                                *
 * hh     is a tentative step size in t.                          *
 *                                                                *
 * uround is the machine unit roundoff.                           *
 *                                                                *
 * JJ     is the band matrix (of type BandMat) to be loaded by    *
 *        an IDABandJacFn routine with an approximation to the    *
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
 *        routine calls res a total of M times, then it should    *
 *        perform the update *nrePtr += M.                        *
 *                                                                *
 * tempv1, tempv2, tempv3 are pointers to memory allocated for    *
 *        N_Vectors which can be used by an IDABandJacFn routine  *
 *        as temporary storage or work space.                     *
 *                                                                *
 *                                                                *
 * Note: The following are two efficient ways to load JJ:         *
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
 * A third way, using the BAND_ELEM(A,i,j) macro, is much less    *
 * efficient in general.  It is only appropriate for use in small *
 * problems in which efficiency of access is NOT a major concern. *
 *                                                                *
 * The IDABandJacFn should return                                 *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDA will attempt to *
 * recover by reducing the stepsize (which changes cj).           *
 ******************************************************************/
  
typedef int (*IDABandJacFn)(integer Neq, integer mupper, integer mlower,
                  real tt, N_Vector yy, N_Vector yp, real cj, 
                  N_Vector constraints, ResFn res, void *rdata, void *jdata, 
                  N_Vector resvec, N_Vector ewt, real hh, real uround, 
                  BandMat JJ, long int *nrePtr, N_Vector tempv1,
                  N_Vector tempv2, N_Vector tempv3);

 
/******************************************************************
 *                                                                *
 * Function : IDABand                                             *
 *----------------------------------------------------------------*
 * A call to the IDABand function links the main IDA integrator   *
 * with the IDABAND linear solver module.                         *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by IDAMalloc.    *
 *                                                                *
 * mupper is the upper bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * mlower is the lower bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * bjac is the banded Jacobian approximation routine to be used.  *
 *         A user-supplied bjac routine must be of type           *
 *         IDABandJacFn (see above).  Pass NULL for bjac if IDA   *
 *         is to use the default difference quotient routine      *
 *         IDABandDQJac supplied with this module.                *
 *                                                                *
 * jdata is a pointer to user data which is passed to the bjac    *
 *         routine every time it is called.                       *
 *                                                                *
 * IDABand returns either                                         *
 *    SUCCESS = 0            if successful, or                    *
 *    IDA_BAND_FAIL = -1     if either IDA_mem was NULL or a      *
 *                           malloc failure occurred, or          *
 *    IDA_BAND_BAD_ARG = -2  if mupper or mlower is illegal.      *
 ******************************************************************/

int IDABand(void *IDA_mem, integer mupper, integer mlower,
              IDABandJacFn bjac, void *jdata);



/******************************************************************
 *                                                                *
 * Function : IDABandDQJac                                        *
 *----------------------------------------------------------------*  
 * This routine generates a banded difference quotient            *
 * approximation to the system Jacobian J = dF/dy + cj*dF/dy'     *
 *                                                                *
 ******************************************************************/

int IDABandDQJac(integer Neq, integer mupper, integer mlower, real tt, 
                 N_Vector yy, N_Vector yp, real cj, N_Vector constraints, 
                 ResFn res, void *rdata,  void *jdata, N_Vector resvec, 
                 N_Vector ewt, real hh, real uround, BandMat JJ, 
                 long int *nrePtr,  N_Vector tempv1, N_Vector tempv2,
                 N_Vector tempv3);

#endif

#ifdef __cplusplus
}
#endif
