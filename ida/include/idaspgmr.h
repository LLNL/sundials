/*******************************************************************
 *                                                                 *
 * File          : idaspgmr.h                                      *
 * Programmers   : Alan C. Hindmarsh and Allan G. Taylor           *
 * Version of    : 2 July 2002                                     *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/ida/LICENSE                           *
 *-----------------------------------------------------------------*
 * This is the header file for the IDA Scaled Preconditioned       *
 * GMRES linear solver module, IDASPGMR.                           *
 *                                                                 *
 *******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idaspgmr_h
#define _idaspgmr_h


#include <stdio.h>
#include "ida.h"
#include "sundialstypes.h"
#include "spgmr.h"
#include "nvector.h"


 
/******************************************************************
 *                                                                *
 * IDASPGMR solver optional output indices                        *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * IDASPGMR optional output. The symbolic names are used as       *
 * indices into the iopt and ropt arrays passed to IDAMalloc.     *
 * The IDASPGMR optional outputs are:                             *
 *                                                                *
 * iopt[SPGMR_NPE]  : number of preconditioner evaluations, i.e.  *
 *                    of calls made to user's precond function.   *
 *                                                                *
 * iopt[SPGMR_NLI]  : number of linear iterations.                *
 *                                                                *
 * iopt[SPGMR_NPS]  : number of calls made to user's psolve       *
 *                    function.                                   *
 *                                                                *
 * iopt[SPGMR_NCFL] : number of linear convergence failures.      *
 *                                                                *
 * iopt[SPGMR_LRW]  : size (in realtype words) of real workspace  *
 *                    matrices and vectors used by this module.   *
 *                                                                *
 * iopt[SPGMR_LIW]  : size (in integertype words) of integer      *
 *                    workspace vectors used by this module.      *
 *                                                                *
 ******************************************************************/
 
enum { SPGMR_NPE=IDA_IOPT_SIZE, SPGMR_NLI, SPGMR_NPS, SPGMR_NCFL,
       SPGMR_LRW, SPGMR_LIW };

 
/******************************************************************
 *                                                                *           
 * Type : IDASpgmrPrecondFn                                       *
 *----------------------------------------------------------------*
 * The optional user-supplied functions Precond and PSolve        *
 * together must define the left preconditoner matrix P           *
 * approximating the system Jacobian matrix                       *
 *    J = dF/dy + cj*dF/dy'                                       *
 * (where the DAE system is F(t,y,y') = 0), and solve the linear  *
 * systems P z = r.   Precond is to do any necessary setup        *
 * operations, and PSolve is to compute the solution of P z = r.  *
 *                                                                *
 * The preconditioner setup function Precond is to evaluate and   *
 * preprocess any Jacobian-related data needed by the             *
 * preconditioner solve function PSolve.  This might include      *
 * forming a crude approximate Jacobian, and performing an LU     *
 * factorization on it.  This function will not be called in      *
 * advance of every call to PSolve, but instead will be called    *
 * only as often as necessary to achieve convergence within the   *
 * Newton iteration in IDA.  If the PSolve function needs no      *
 * preparation, the Precond function can be NULL.                 *
 *                                                                *
 * Each call to the Precond function is preceded by a call to     *
 * the system function res with the same (t,y,y') arguments.      *
 * Thus the Precond function can use any auxiliary data that is   *
 * computed and saved by the res function and made accessible     *
 * to Precond.                                                    *
 *                                                                *
 * The error weight vector ewt, step size hh, and unit roundoff   *
 * uround are provided to the Precond function for possible use   *
 * in approximating Jacobian data, e.g. by difference quotients.  *
 *                                                                *
 * A preconditioner setup function Precond must have the          *
 * prototype given below.  Its parameters are as follows:         *
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
 * rr  is the current value of the residual vector F(t,y,y').     *
 *                                                                *
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * res    is the residual function for the DAE problem.           *
 *                                                                *
 * rdata  is a pointer to user data to be passed to res, the same *
 *        as the rdata parameter passed to IDAMalloc.             *
 *                                                                *
 * pdata  is a pointer to user preconditioner data - the same as  *
 *        the pdata parameter passed to IDASpgmr.                 *
 *                                                                *
 * ewt    is the error weight vector.                             *
 *                                                                *
 * constraints  is the constraints vector.                        *
 *                                                                *
 * hh     is a tentative step size in t.                          *
 *                                                                *
 * uround is the machine unit roundoff.                           *
 *                                                                *
 * nrePtr is a pointer to the memory location containing the      *
 * IDA problem data nre = number of calls to res.  This Precond   *
 * routine should update the counter nre by adding on the number  *
 * of res calls it makes in order to compute P, if any.           *
 * Thus if this routine calls res a total of W times, it should   *
 * perform the update *nrePtr += W.                               *
 *                                                                *
 * tempv1, tempv2, tempv3 are pointers to vectors of type         *
 * N_Vector which can be used by an IDASpgmrPrecondFn routine as  *
 * temporary storage or work space.                               *
 *                                                                *
 *                                                                *
 * The IDASpgmrPrecondFn should return                            *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * In the case of a recoverable error return, IDA will attempt to *
 * recover by reducing the stepsize (which changes cj).           *
 ******************************************************************/
  
typedef int (*IDASpgmrPrecondFn)(integertype Neq, realtype tt, N_Vector yy,
         N_Vector yp, N_Vector rr, realtype cj, ResFn res, 
         void *rdata, void *pdata, N_Vector ewt, N_Vector constraints, 
         realtype hh, realtype uround,  long int *nrePtr, N_Vector tempv1, 
         N_Vector tempv2, N_Vector tempv3);


/******************************************************************
 *                                                                *           
 * Type : IDASpgmrPSolveFn                                        *
 *----------------------------------------------------------------*
 * The optional user-supplied function PSolve must compute a      *
 * solution to the linear system P z = r, where P is the left     *
 * preconditioner defined by the user.  If no preconditioning     *
 * is desired, pass NULL for PSolve to IDASpgmr.                  *
 *                                                                *
 * A preconditioner solve function PSolve must have the           *
 * prototype given below.  Its parameters are as follows:         *
 *                                                                *
 * Neq is the problem size, and length of all vector arguments.   *
 *                                                                *
 * tt  is the current value of the independent variable t.        *
 *                                                                *
 * yy  is the current value of the dependent variable vector y.   *
 *                                                                *
 * yp  is the current value of the derivative vector y'.          *
 *                                                                *
 * rr  is the current value of the residual vector F(t,y,y').     *
 *                                                                *
 * cj  is the scalar in the system Jacobian, proportional to 1/hh.*
 *                                                                *
 * res    is the residual function for the DAE problem.           *
 *                                                                *
 * rdata  is a pointer to user data to be passed to res, the same *
 *        as the rdata parameter passed to IDAMalloc.             *
 *                                                                *
 * pdata  is a pointer to user preconditioner data - the same as  *
 *        the pdata parameter passed to IDASpgmr.                 *
 *                                                                *
 * ewt    is the input error weight vector (see delta below).     *
 *                                                                *
 * delta  is an input tolerance for use by PSolve if it uses an   *
 *        iterative method in its solution.   In that case, the   *
 *        the residual vector r - P z of the system should be     *
 *        made less than delta in weighted L2 norm, i.e.,         *
 *            sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .            *
 *                                                                *
 * rvec   is the input right-hand side vector r.                  *
 *                                                                *
 * zvec   is the computed solution vector z.                      *
 *                                                                *
 * nrePtr is a pointer to the memory location containing the      *
 * IDA problem data nre = number of calls to res.  This PSolve    *
 * routine should update the counter nre by adding on the number  *
 * of res calls it makes in order to compute z, if any.           *
 * Thus if this routine calls res a total of W times, it should   *
 * perform the update *nrePtr += W.                               *
 *                                                                *
 * tempv  is an N_Vector which can be used by the PSolve          *
 * routine as temporary storage or work space.                    *
 *                                                                *
 *                                                                *
 * The IDASpgmrPSolveFn should return                             *
 *     0 if successful,                                           *
 *     a positive int if a recoverable error occurred, or         *
 *     a negative int if a nonrecoverable error occurred.         *
 * Following a recoverable error, IDA will attempt to recover by  *
 * updating the preconditioner and/or reducing the stepsize.      *
 ******************************************************************/
  
typedef int (*IDASpgmrPSolveFn)(integertype Neq, realtype tt, N_Vector yy,
             N_Vector yp, N_Vector rr, realtype cj, ResFn res, void *rdata,
             void *pdata, N_Vector ewt, realtype delta, N_Vector rvec,
             N_Vector zvec, long int *nrePtr, N_Vector tempv);

 
/******************************************************************
 *                                                                *
 * Function : IDASpgmr                                            *
 *----------------------------------------------------------------*
 * A call to the IDASpgmr function links the main IDA integrator  *
 * with the IDASPGMR linear solver module.  Its parameters are    *
 * as follows:                                                    *
 *                                                                *
 * IDA_mem   is the pointer to IDA memory returned by IDAMalloc.  *
 *                                                                *
 * precond   is the user's preconditioner setup routine. It is    *
 *           used to evaluate and preprocess any Jacobian-related *
 *           data needed by the psolve routine.  See the          *
 *           description of the type IDASpgmrPrecondFn above.     *
 *           Pass NULL if no such data setup is required.         *
 *                                                                *
 * psolve    is the user's preconditioner solve routine. It is    *
 *           used to solve linear systems P z = r, where P is the *
 *           preconditioner matrix.  See the description of the   *
 *           type IDASpgmrPSolveFn above.  Pass NULL for psolve   * 
 *           if no preconditioning is to be done.  However, a     *
 *           preconditioner of some form is strongly encouraged.  *
 *                                                                *
 * gstype    is the type of Gram-Schmidt orthogonalization to be  *
 *           used.  This must be one of the two enumeration       *
 *           constants MODIFIED_GS or CLASSICAL_GS defined in     *
 *           iterativ.h.  These correspond to using modified or   *
 *           classical Gram-Schmidt algorithms, respectively.     *
 *                                                                *
 * maxl      is the maximum Krylov subspace dimension, an         *
 *           optional input.  Pass 0 to use the default value,    *
 *           MIN(Neq, 5).  Otherwise pass a positive integer.     *
 *                                                                *
 * maxrs     is the maximum number of restarts to be used in the  *
 *           GMRES algorithm, an optional input.  maxrs must be a *
 *           non-negative integer, or -1.  Pass 0 to use the      *
 *           default value, which is 5.  Pass -1 to use the       *
 *           value 0, meaning no restarts.  In any case, maxrs    *
 *           will be restricted to the range 0 to Neq/maxl.       *
 *                                                                *
 * eplifac   is a factor in the linear iteration convergence      *
 *           test constant, an optional input.  Pass 0.0 to use   *
 *           the default, which is 1.0.  Otherwise eplifac must   *
 *           be a positive real number.                           *
 *                                                                *
 * dqincfac  is a factor in the increments to yy used in the      *
 *           difference quotient approximations to matrix-vector  *
 *           products Jv, an optional input.  Pass 0.0 to use     *
 *           the default, which is 1.0.  Otherwise dqincfac must  *
 *           be a positive real number.                           *
 *                                                                *
 * pdata     is a pointer to user preconditioner data.  This      *
 *           pointer is passed to precond and psolve every time   *
 *           these routines are called.                           *
 *                                                                *
 * The return values of IDASpgmr are:                             *
 *    SUCCESS       = 0  if successful                            *
 *    LMEM_FAIL     = -1 if there was a memory allocation failure *
 *    LIN_ILL_INPUT = -2 if there was illegal input.              *
 *                                                                *
 ******************************************************************/

int IDASpgmr(void *IDA_mem, IDASpgmrPrecondFn precond, 
             IDASpgmrPSolveFn psolve, int gstype, int maxl, int maxrs,
             realtype eplifac, realtype dqincfac, void *pdata);


/******************************************************************
 *                                                                *
 * Function : IDAReInitSpgmr                                      *
 *----------------------------------------------------------------*
 * A call to the IDAReInitSpgmr function resets the link between  *
 * the main IDA integrator and the IDASPGMR linear solver.        *
 * After solving one problem using IDASPGMR, call IDAReInit and   *
 * then IDAReInitSpgmr to solve another problem of the same size, *
 * if there is a change in the IDASpgmr parameters precond,       *
 * psolve, gstype, maxrs, eplifac, dqincfac, or pdata, but not in *
 * maxl.  If there is a change in maxl, then IDASpgmr must be     *
 * called again, and the linear solver memory will be reallocated.*
 * If there is no change in parameters, it is not necessary to    *
 * call either IDAReInitSpgmr or IDASpgmr for the new problem.    *
 *                                                                *
 * All arguments to IDAReInitSpgmr have the same names and        * 
 * meanings as those of IDASpgmr.  The IDA_mem argument must be   *
 * identical to its value in the previous IDASpgmr call.          *
 *                                                                *
 * The return values of IDAReInitSpgmr are:                       *
 *    SUCCESS       = 0  if successful                            *
 *    LMEM_FAIL     = -1 if the IDA_mem argument is NULL          *
 *    LIN_ILL_INPUT = -2 if there was illegal input.              *
 *                                                                *
 ******************************************************************/

int IDAReInitSpgmr(void *IDA_mem, IDASpgmrPrecondFn precond, 
             IDASpgmrPSolveFn psolve, int gstype, int maxl, int maxrs,
             realtype eplifac, realtype dqincfac, void *pdata);


#endif

#ifdef __cplusplus
}
#endif
