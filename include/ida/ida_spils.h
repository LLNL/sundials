/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *         Alan Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 * Header file for the Scaled and Preconditioned Iterative Linear 
 * Solver interface in IDA.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPILS_H
#define _IDASPILS_H

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  IDASPILS Constants
  ===============================================================*/

#define IDASPILS_SUCCESS     0
#define IDASPILS_MEM_NULL   -1 
#define IDASPILS_LMEM_NULL  -2 
#define IDASPILS_ILL_INPUT  -3
#define IDASPILS_MEM_FAIL   -4
#define IDASPILS_PMEM_NULL  -5
#define IDASPILS_SUNLS_FAIL -6

  
/*===============================================================
  IDASPILS user-supplied function prototypes
  ===============================================================*/
  
/*---------------------------------------------------------------
  Type: IDASpilsPrecSetupFn

  The optional user-supplied functions PrecSetup and PrecSolve
  together must define the left preconditoner matrix P
  approximating the system Jacobian matrix
     J = dF/dy + c_j*dF/dy'
  (where the DAE system is F(t,y,y') = 0), and solve the linear
  systems P z = r.   PrecSetup is to do any necessary setup
  operations, and PrecSolve is to compute the solution of
  P z = r.
 
  The preconditioner setup function PrecSetup is to evaluate and
  preprocess any Jacobian-related data needed by the
  preconditioner solve function PrecSolve.  This might include
  forming a crude approximate Jacobian, and performing an LU
  factorization on it.  This function will not be called in
  advance of every call to PrecSolve, but instead will be called
  only as often as necessary to achieve convergence within the
  Newton iteration.  If the PrecSolve function needs no
  preparation, the PrecSetup function can be NULL when passed to
  IDASpilsSetPreconditioner().
 
  Each call to the PrecSetup function is preceded by a call to
  the system function res with the same (t,y,y') arguments.
  Thus the PrecSetup function can use any auxiliary data that is
  computed and saved by the res function and made accessible
  to PrecSetup.
 
  A preconditioner setup function PrecSetup must have the
  prototype given below.  Its parameters are as follows:
 
  tt  is the current value of the independent variable t.
 
  yy  is the current value of the dependent variable vector,
      namely the predicted value of y(t).
 
  yp  is the current value of the derivative vector y',
      namely the predicted value of y'(t).
 
  rr  is the current value of the residual vector F(t,y,y').
 
  c_j is the scalar in the system Jacobian, proportional to 1/hh.
 
  user_data is a pointer to user data, the same as the user_data
      parameter passed to IDASetUserData.
 
  NOTE: If the user's preconditioner needs other quantities,
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through
      IDAGetCurrentStep and IDAGetErrWeights, respectively (see
      ida.h). The unit roundoff is available as
      UNIT_ROUNDOFF defined in sundials_types.h
 
  The IDASpilsPrecSetupFn should return
      0 if successful,
      a positive int if a recoverable error occurred, or
      a negative int if a nonrecoverable error occurred.
  In the case of a recoverable error return, the integrator will
  attempt to recover by reducing the stepsize (which changes cj).
  ---------------------------------------------------------------*/
typedef int (*IDASpilsPrecSetupFn)(realtype tt, N_Vector yy,
                                   N_Vector yp, N_Vector rr,
				   realtype c_j, void *user_data);

  
/*---------------------------------------------------------------
  Type: IDASpilsPrecSolveFn

  The optional user-supplied function PrecSolve must compute a
  solution to the linear system P z = r, where P is the left
  preconditioner defined by the user.  If no preconditioning
  is desired, pass NULL for PrecSolve to 
  IDASpilsSetPreconditioner() (or do not call that routine at 
  all).
 
  A preconditioner solve function PrecSolve must have the
  prototype given below.  Its parameters are as follows:
 
  tt is the current value of the independent variable t.
 
  yy is the current value of the dependent variable vector y.
 
  yp is the current value of the derivative vector y'.
 
  rr is the current value of the residual vector F(t,y,y').
 
  rvec is the input right-hand side vector r.
 
  zvec is the computed solution vector z.
 
  c_j is the scalar in the system Jacobian, proportional to 1/hh.
 
  delta is an input tolerance for use by PrecSolve if it uses an
      iterative method in its solution.   In that case, the
      the residual vector r - P z of the system should be
      made less than delta in weighted L2 norm, i.e.,
             sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .
      Note: the error weight vector ewt can be obtained
      through a call to the routine IDAGetErrWeights.
 
  user_data is a pointer to user data, the same as the user_data
      parameter passed to IDASetUserData.
 
  The IDASpilsPrecSolveFn should return
      0 if successful,
      a positive int if a recoverable error occurred, or
      a negative int if a nonrecoverable error occurred.
  Following a recoverable error, the integrator will attempt to
  recover by updating the preconditioner and/or reducing the
  stepsize.
  ---------------------------------------------------------------*/
typedef int (*IDASpilsPrecSolveFn)(realtype tt, N_Vector yy,
                                   N_Vector yp, N_Vector rr,
				   N_Vector rvec, N_Vector zvec,
				   realtype c_j, realtype delta,
                                   void *user_data);

/*---------------------------------------------------------------
 Type: IDASpilsJacTimesSetupFn

 The user-supplied Jacobian-times-vector product setup function 
 JacTimesSetup and the user-supplied Jacobian-times-vector 
 product function JTimes together must generate the product
 J*v for v, where J is the Jacobian matrix
     J = dF/dy + c_j*dF/dy'
 or an approximation to it, and v is a given vector. 

 Each call to the JacTimesSetup function is preceded by a call 
 to the residual res with the same (t,y) arguments.  Thus the 
 JacTimesSetup function can use any auxiliary data that is 
 computed and saved by the res function and made accessible to 
 JacTimesSetup.

 A function JacTimesSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
          namely the predicted value of y(t).

 fy      is the vector f(t,y).

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the IDASetUserData function.

 Returned value:
 The value to be returned by the JacTimesSetup function is a flag
 indicating whether it was successful.  This value should be
   0   if successful,
   > 0 for a recoverable error (step will be retried),
   < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*IDASpilsJacTimesSetupFn)(realtype tt, N_Vector yy,
                                       N_Vector yp, N_Vector rr,
                                       realtype c_j, void *user_data);


/*---------------------------------------------------------------
  Type: IDASpilsJacTimesVecFn

  The user-supplied function jtimes is to generate the product
  J*v for given v, where J is the Jacobian matrix
     J = dF/dy + c_j*dF/dy'
   or an approximation to it, and v is a given vector.
  It should return 0 if successful and a nonzero int otherwise.
 
  A function jtimes must have the prototype given below. Its
  parameters are as follows:
 
    tt   is the current value of the independent variable.
 
    yy   is the current value of the dependent variable vector,
         namely the predicted value of y(t).
 
    yp   is the current value of the derivative vector y',
         namely the predicted value of y'(t).
 
    rr   is the current value of the residual vector F(t,y,y').
 
    v    is the N_Vector to be multiplied by J.
 
    Jv   is the output N_Vector containing J*v.
 
    c_j  is the scalar in the system Jacobian, proportional
         to 1/hh.
 
    user_data is a pointer to user data, the same as the
         pointer passed to IDASetUserData.
 
    tmp1, tmp2 are two N_Vectors which can be used by Jtimes for
          work space.
  ---------------------------------------------------------------*/
typedef int (*IDASpilsJacTimesVecFn)(realtype tt, N_Vector yy,
                                     N_Vector yp, N_Vector rr,
				     N_Vector v, N_Vector Jv,
				     realtype c_j, void *user_data,
				     N_Vector tmp1, N_Vector tmp2);


/*===============================================================
  IDASPILS Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  Required inputs for the IDASPILS linear solver interface:

  IDASpilsSetLinearSolver specifies the iterative SUNLinearSolver 
  object that IDA should use.  

  The return value is one of:
     IDASPILS_SUCCESS   if successful
     IDASPILS_MEM_NULL  if the IDA memory was NULL
     IDASPILS_ILL_INPUT if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASpilsSetLinearSolver(void *ida_mem, 
                                            SUNLinearSolver LS);

  
/*---------------------------------------------------------------
  Optional inputs to the IDASPILS linear solver -- ALL of these 
  must be called AFTER the corresponding iterative linear solver 
  object has been attached to IDA).

  IDASpilsSetPreconditioner specifies the PrecSetup and PrecSolve
            functions.
            Default is NULL for both arguments.
  IDASpilsSetJacTimes specifies the jtsetup and jtimes functions.
            Default is to use an internal finite difference
            approximation routine for the jtimes (with no setup).
  IDASpilsSetEpsLin specifies the factor in the linear iteration
            convergence test constant.
            Default is 0.05
  IDASpilsSetIncrementFactor specifies a factor in the increments
            to yy used in the difference quotient approximations
            to matrix-vector products Jv.
            Default is 1.0
                                                                 
  The return value of IDASpilsSet* is one of:
     IDASPILS_SUCCESS   if successful
     IDASPILS_MEM_NULL  if the ida memory was NULL
     IDASPILS_LMEM_NULL if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASpilsSetPreconditioner(void *ida_mem,
                                              IDASpilsPrecSetupFn pset, 
					      IDASpilsPrecSolveFn psolve);
SUNDIALS_EXPORT int IDASpilsSetJacTimes(void *ida_mem,
                                        IDASpilsJacTimesSetupFn jtsetup,
                                        IDASpilsJacTimesVecFn jtimes);

SUNDIALS_EXPORT int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac);
SUNDIALS_EXPORT int IDASpilsSetIncrementFactor(void *ida_mem,
                                               realtype dqincfac);


/*---------------------------------------------------------------
  Optional outputs from the IDASPILS linear solver interface:
                                                                 
  IDASpilsGetWorkSpace returns the real and integer workspace used 
      by IDASPILS.                                                  

  IDASpilsGetNumPrecEvals returns the number of preconditioner   
      evaluations, i.e. the number of calls made to PrecSetup    
      with jok==SUNFALSE.                                           

  IDASpilsGetNumPrecSolves returns the number of calls made to   
      PrecSolve.                                                 

  IDASpilsGetNumLinIters returns the number of linear iterations.

  IDASpilsGetNumConvFails returns the number of linear           
      convergence failures.                                      

  IDASpilsGetNumJtimesEvals returns the number of calls to jtimes

  IDASpilsGetNumResEvals returns the number of calls to the user 
      res routine due to finite difference Jacobian times vector 
      evaluation.                                                

  IDASpilsGetLastFlag returns the last error flag set by any of
      the IDASPILS interface functions.
                                                                 
  The return value of IDASpilsGet* is one of:
     IDASPILS_SUCCESS   if successful
     IDASPILS_MEM_NULL  if the ida memory was NULL
     IDASPILS_LMEM_NULL if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASpilsGetWorkSpace(void *ida_mem,
                                         long int *lenrwLS,
                                         long int *leniwLS);
SUNDIALS_EXPORT int IDASpilsGetNumPrecEvals(void *ida_mem,
                                            long int *npevals);
SUNDIALS_EXPORT int IDASpilsGetNumPrecSolves(void *ida_mem,
                                             long int *npsolves);
SUNDIALS_EXPORT int IDASpilsGetNumLinIters(void *ida_mem,
                                           long int *nliters);
SUNDIALS_EXPORT int IDASpilsGetNumConvFails(void *ida_mem,
                                            long int *nlcfails);
SUNDIALS_EXPORT int IDASpilsGetNumJtimesEvals(void *ida_mem,
                                              long int *njvevals);
SUNDIALS_EXPORT int IDASpilsGetNumResEvals(void *ida_mem,
                                           long int *nrevalsLS); 
SUNDIALS_EXPORT int IDASpilsGetLastFlag(void *ida_mem,
                                        long int *flag);

/*---------------------------------------------------------------
  The following function returns the name of the constant 
  associated with an IDASPILS return flag
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT char *IDASpilsGetReturnFlagName(long int flag);

#ifdef __cplusplus
}
#endif

#endif
