/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-01-24 23:54:35 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * This is the public header file for the scaled preconditioned
 * Bi-CGSTAB linear solver module, IDASPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPBCG_H
#define _IDASPBCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "spbcg.h"
#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * Type : IDASpbcgPrecSetupFn
 * -----------------------------------------------------------------
 * The optional user-supplied functions PrecSetup and PrecSolve
 * together must define the left preconditoner matrix P
 * approximating the system Jacobian matrix
 *    J = dF/dy + c_j*dF/dy'
 * (where the DAE system is F(t,y,y') = 0), and solve the linear
 * systems P z = r. PrecSetup is to do any necessary setup
 * operations, and PrecSolve is to compute the solution of
 * P z = r.
 *
 * The preconditioner setup function PrecSetup is to evaluate and
 * preprocess any Jacobian-related data needed by the
 * preconditioner solve function PrecSolve. This might include
 * forming a crude approximate Jacobian, and performing an LU
 * factorization on it. This function will not be called in
 * advance of every call to PrecSolve, but instead will be called
 * only as often as necessary to achieve convergence within the
 * Newton iteration. If the PrecSolve function needs no
 * preparation, the PrecSetup function can be NULL.
 *
 * Each call to the PrecSetup function is preceded by a call to
 * the system function res with the same (t,y,y') arguments.
 * Thus the PrecSetup function can use any auxiliary data that is
 * computed and saved by the res function and made accessible
 * to PrecSetup.
 *
 * A preconditioner setup function PrecSetup must have the
 * prototype given below. Its parameters are as follows:
 *
 * tt  is the current value of the independent variable t.
 *
 * yy  is the current value of the dependent variable vector,
 *     namely the predicted value of y(t).
 *
 * yp  is the current value of the derivative vector y',
 *     namely the predicted value of y'(t).
 *
 * rr  is the current value of the residual vector F(t,y,y').
 *
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.
 *
 * prec_data is a pointer to user preconditioner data - the same as
 *           the pdata parameter passed to IDASpbcg.
 *
 * tmp1, tmp2, tmp3 are pointers to vectors of type N_Vector
 *                  which can be used by an IDASpbcgPrecSetupFn
 *                  routine as temporary storage or work space.
 *
 * NOTE: If the user's preconditioner needs other quantities,
 *       they are accessible as follows: hcur (the current step size)
 *       and ewt (the error weight vector) are accessible through
 *       IDAGetCurrentStep and IDAGetErrWeights, respectively (see
 *       ida.h). The unit roundoff is available as UNIT_ROUNDOFF
 *       defined in sundialstypes.h.
 *
 * The IDASpbcgPrecSetupFn should return:
 *     0 if successful,
 *     a positive int if a recoverable error occurred, or
 *     a negative int if a nonrecoverable error occurred.
 * In the case of a recoverable error return, the integrator will
 * attempt to recover by reducing the step size (which changes cj).
 * -----------------------------------------------------------------
 */

typedef int (*IDASpbcgPrecSetupFn)(realtype tt,
                                   N_Vector yy, N_Vector yp, N_Vector rr,
                                   realtype c_j, void *prec_data,
                                   N_Vector tmp1, N_Vector tmp2,
                                   N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type : IDASpbcgPrecSolveFn
 * -----------------------------------------------------------------
 * The optional user-supplied function PrecSolve must compute a
 * solution to the linear system P z = r, where P is the left
 * preconditioner defined by the user. If no preconditioning
 * is desired, pass NULL for PrecSolve to IDASpbcg.
 *
 * A preconditioner solve function PrecSolve must have the
 * prototype given below. Its parameters are as follows:
 *
 * tt is the current value of the independent variable t.
 *
 * yy is the current value of the dependent variable vector y.
 *
 * yp is the current value of the derivative vector y'.
 *
 * rr is the current value of the residual vector F(t,y,y').
 *
 * rvec is the input right-hand side vector r.
 *
 * zvec is the computed solution vector z.
 *
 * c_j is the scalar in the system Jacobian, proportional to 1/hh.
 *
 * delta is an input tolerance for use by PrecSolve if it uses an
 *       iterative method in its solution. In that case, the
 *       residual vector r - P z of the system should be
 *       made less than delta in weighted L2 norm, i.e.,
 *            sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .
 *       Note: The error weight vector ewt can be obtained
 *       through a call to the routine IDAGetErrWeights.
 *
 * prec_data is a pointer to user preconditioner data - the same as
 *           the pdata parameter passed to IDASpbcg.
 *
 * tmp is an N_Vector which can be used by the PrecSolve
 *     routine as temporary storage or work space.
 *
 * The IDASpbcgPrecSolveFn should return:
 *     0 if successful,
 *     a positive int if a recoverable error occurred, or
 *     a negative int if a nonrecoverable error occurred.
 * Following a recoverable error, the integrator will attempt to
 * recover by updating the preconditioner and/or reducing the
 * step size.
 * -----------------------------------------------------------------
 */

typedef int (*IDASpbcgPrecSolveFn)(realtype tt,
                                   N_Vector yy, N_Vector yp, N_Vector rr,
                                   N_Vector rvec, N_Vector zvec,
                                   realtype c_j, realtype delta, void *prec_data,
                                   N_Vector tmp);

/*
 * -----------------------------------------------------------------
 * Type : IDASpbcgJacTimesVecFn
 * -----------------------------------------------------------------
 * The user-supplied function jtimes is to generate the product
 * J*v for given v, where J is the Jacobian matrix
 *    J = dF/dy + c_j*dF/dy'
 * or an approximation to it, and v is a given vector.
 * It should return 0 if successful and a nonzero int otherwise.
 *
 * A function jtimes must have the prototype given below. Its
 * parameters are as follows:
 *
 *   tt   is the current value of the independent variable.
 *
 *   yy   is the current value of the dependent variable vector,
 *        namely the predicted value of y(t).
 *
 *   yp   is the current value of the derivative vector y',
 *        namely the predicted value of y'(t).
 *
 *   rr   is the current value of the residual vector F(t,y,y').
 *
 *   v    is the N_Vector to be multiplied by J.
 *
 *   Jv   is the output N_Vector containing J*v.
 *
 *   c_j  is the scalar in the system Jacobian, proportional
 *        to 1/hh.
 *
 *   jac_data is a pointer to user Jacobian data, the same as the
 *            pointer passed to CVSpbcg.
 *
 *   tmp1, tmp2 are two N_Vectors which can be used by Jtimes for
 *              work space.
 * -----------------------------------------------------------------
 */

typedef int (*IDASpbcgJacTimesVecFn)(realtype tt,
                                     N_Vector yy, N_Vector yp, N_Vector rr,
                                     N_Vector v, N_Vector Jv,
                                     realtype c_j, void *jac_data,
                                     N_Vector tmp1, N_Vector tmp2);

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcg
 * -----------------------------------------------------------------
 * A call to the IDASpbcg function links the main integrator with
 * the IDASPBCG linear solver module. Its parameters are as
 * follows:
 *
 * IDA_mem   is the pointer to memory block returned by IDACreate.
 *
 * maxl      is the maximum Krylov subspace dimension, an
 *           optional input. Pass 0 to use the default value.
 *           Otherwise pass a positive integer.
 *
 * The return values of IDASpbcg are:
 *    IDASPBCG_SUCCESS    if successful
 *    IDASPBCG_MEM_NULL   if the ida memory was NULL
 *    IDASPBCG_MEM_FAIL   if there was a memory allocation failure
 *    IDASPBCG_ILL_INPUT  if there was illegal input.
 * -----------------------------------------------------------------
 */

int IDASpbcg(void *ida_mem, int maxl);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the IDASPBCG linear solver
 * -----------------------------------------------------------------
 * IDASpbcgSetPrecSolveFn specifies the PrecSolve function.
 *                        Default is NULL.
 * IDASpbcgSetPrecSetupFn specifies the PrecSetup function.
 *                        Default is NULL.
 * IDASpbcgSetPrecData specifies a pointer to user preconditioner
 *                     data. This pointer is passed to PrecSetup and
 *                     PrecSolve every time these routines are called.
 *                     Default is NULL.
 * IDASpbcgSetJacTimesVecFn specifies the jtimes function.
 *                          Default is to use an internal finite
 *                          difference approximation routine.
 * IDASpbcgSetJacData specifies a pointer to user Jacobian data.
 *                    This pointer is passed to jtimes every time this
 *                    routine is called.
 *                    Default is NULL.
 * IDASpbcgSetEpsLin specifies the factor in the linear iteration
 *                   convergence test constant.
 *                   Default is 0.05.
 * IDASpbcgSetIncrementFactor specifies a factor in the increments
 *                            to yy used in the difference quotient
 *                            approximations to matrix-vector products Jv.
 *                            Default is 1.0.
 *
 * The return value of IDASpbcgSet* is one of:
 *    IDASPBCG_SUCCESS   if successful
 *    IDASPBCG_MEM_NULL  if the ida memory was NULL
 *    IDASPBCG_LMEM_NULL if the idaspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int IDASpbcgSetPrecSolveFn(void *ida_mem, IDASpbcgPrecSolveFn psolve);
int IDASpbcgSetPrecSetupFn(void *ida_mem, IDASpbcgPrecSetupFn pset);
int IDASpbcgSetPrecData(void *ida_mem, void *prec_data);
int IDASpbcgSetJacTimesVecFn(void *ida_mem, IDASpbcgJacTimesVecFn jtimes);
int IDASpbcgSetJacData(void *ida_mem, void *jac_data);
int IDASpbcgSetEpsLin(void *ida_mem, realtype eplifac);
int IDASpbcgSetIncrementFactor(void *ida_mem, realtype dqincfac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the IDASPBCG linear solver
 * -----------------------------------------------------------------
 * IDASpbcgGetWorkSpace returns the real and integer workspace used
 *                      by IDASPBCG.
 * IDASpbcgGetNumPrecEvals returns the number of preconditioner
 *                         evaluations, i.e. the number of calls made
 *                         to PrecSetup with jok==FALSE.
 * IDASpbcgGetNumPrecSolves returns the number of calls made to
 *                          PrecSolve.
 * IDASpbcgGetNumLinIters returns the number of linear iterations.
 * IDASpbcgGetNumConvFails returns the number of linear
 *                         convergence failures.
 * IDASpbcgGetNumJtimesEvals returns the number of calls to jtimes
 * IDASpbcgGetNumResEvals returns the number of calls to the user
 *                        res routine due to finite difference Jacobian
 *                        times vector evaluation.
 * IDASpbcgGetLastFlag returns the last error flag set by any of
 *                     the IDASPBCG interface functions.
 *
 * The return value of IDASpbcgGet* is one of:
 *    IDASPBCG_SUCCESS   if successful
 *    IDASPBCG_MEM_NULL  if the ida memory was NULL
 *    IDASPBCG_LMEM_NULL if the idaspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int IDASpbcgGetWorkSpace(void *ida_mem, long int *lenrwSG, long int *leniwSG);
int IDASpbcgGetNumPrecEvals(void *ida_mem, long int *npevals);
int IDASpbcgGetNumPrecSolves(void *ida_mem, long int *npsolves);
int IDASpbcgGetNumLinIters(void *ida_mem, long int *nliters);
int IDASpbcgGetNumConvFails(void *ida_mem, long int *nlcfails);
int IDASpbcgGetNumJtimesEvals(void *ida_mem, long int *njvevals);
int IDASpbcgGetNumResEvals(void *ida_mem, long int *nrevalsSG); 
int IDASpbcgGetLastFlag(void *ida_mem, int *flag);

/* IDASPBCG return values */

#define IDASPBCG_SUCCESS     0
#define IDASPBCG_MEM_NULL   -1 
#define IDASPBCG_LMEM_NULL  -2 
#define IDASPBCG_ILL_INPUT  -3
#define IDASPBCG_MEM_FAIL   -4

#ifdef __cplusplus
}
#endif

#endif
