/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-12-07 19:44:38 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE/CVODES scaled
 * preconditioned Bi-CGSTAB linear solver, CVSPBCG.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _CVSPBCG_H
#define _CVSPBCG_H

#include "spbcg.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * CVSPBCG solver constants
 * -----------------------------------------------------------------
 * CVSPBCG_MAXL   : default value for the maximum Krylov
 *                  dimension
 *
 * CVSPBCG_MSBPRE : maximum number of steps between
 *                  preconditioner evaluations
 *
 * CVSPBCG_DGMAX  : maximum change in gamma between
 *                  preconditioner evaluations
 *
 * CVSPBCG_DELT   : default value for factor by which the
 *                  tolerance on the nonlinear iteration is
 *                  multiplied to get a tolerance on the linear
 *                  iteration
 * -----------------------------------------------------------------
 */

#define CVSPBCG_MAXL   5
#define CVSPBCG_MSBPRE 50
#define CVSPBCG_DGMAX  RCONST(0.2)
#define CVSPBCG_DELT   RCONST(0.05)

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgPrecSetupFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup function PrecSetup and
 * the user-supplied preconditioner solve function PrecSolve
 * together must define left and right preconditoner matrices
 * P1 and P2 (either of which may be trivial), such that the
 * product P1*P2 is an approximation to the Newton matrix
 * M = I - gamma*J. Here J is the system Jacobian J = df/dy,
 * and gamma is a scalar proportional to the integration step
 * size h. The solution of systems P z = r, with P = P1 or P2,
 * is to be carried out by the PrecSolve function, and PrecSetup
 * is to do any necessary setup operations.
 *
 * The user-supplied preconditioner setup function PrecSetup
 * is to evaluate and preprocess any Jacobian-related data
 * needed by the preconditioner solve function PrecSolve.
 * This might include forming a crude approximate Jacobian,
 * and performing an LU factorization on the resulting
 * approximation to M. This function will not be called in
 * advance of every call to PrecSolve, but instead will be called
 * only as often as necessary to achieve convergence within the
 * Newton iteration. If the PrecSolve function needs no
 * preparation, the PrecSetup function can be NULL.
 *
 * For greater efficiency, the PrecSetup function may save
 * Jacobian-related data and reuse it, rather than generating it
 * from scratch. In this case, it should use the input flag jok
 * to decide whether to recompute the data, and set the output
 * flag *jcurPtr accordingly.
 *
 * Each call to the PrecSetup function is preceded by a call to
 * the CVRhsFn f with the same (t,y) arguments. Thus the PrecSetup
 * function can use any auxiliary data that is computed and
 * saved by the f function and made accessible to PrecSetup.
 *
 * A function PrecSetup must have the prototype given below.
 * Its parameters are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *          namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE  means that Jacobian data, if saved from
 *                  the previous PrecSetup call, can be reused
 *                  (with the current value of gamma).
 *         A Precset call with jok == TRUE can only occur after
 *         a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         to be set by PrecSetup as follows:
 *         Set *jcurPtr = TRUE if Jacobian data was recomputed.
 *         Set *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                        but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * P_data  is a pointer to user data - the same as the P_data
 *         parameter passed to CVSpbcg.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *                      for N_Vectors which can be used by
 *                      CVSpbcgPrecSetupFn as temporary storage or
 *                      work space.
 *
 * NOTE: If the user's preconditioner needs other quantities,
 *       they are accessible as follows: hcur (the current stepsize)
 *       and ewt (the error weight vector) are accessible through
 *       CVodeGetCurrentStep and CVodeGetErrWeights, respectively).
 *       The unit roundoff is available as UNIT_ROUNDOFF defined in
 *       sundialstypes.h.
 *
 * Returned value:
 * The value to be returned by the PrecSetup function is a flag
 * indicating whether it was successful. This value should be
 *     0 if successful,
 *   > 0 for a recoverable error (step will be retried),
 *   < 0 for an unrecoverable error (integration is halted).
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  realtype gamma, void *P_data,
                                  N_Vector tmp1, N_Vector tmp2,
                                  N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve function PrecSolve
 * is to solve a linear system P z = r in which the matrix P is
 * one of the preconditioner matrices P1 or P2, depending on the
 * type of preconditioning chosen.
 *
 * A function PrecSolve must have the prototype given below.
 * Its parameters are as follows:
 *
 * t      is the current value of the independent variable.
 *
 * y      is the current value of the dependent variable vector.
 *
 * fy     is the vector f(t,y).
 *
 * r      is the right-hand side vector of the linear system.
 *
 * z      is the output vector computed by PrecSolve.
 *
 * gamma  is the scalar appearing in the Newton matrix.
 *
 * delta  is an input tolerance for use by PSolve if it uses
 *        an iterative method in its solution. In that case,
 *        the residual vector Res = r - P z of the system
 *        should be made less than delta in weighted L2 norm,
 *        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
 *        Note: the error weight vector ewt can be obtained
 *        through a call to the routine CVodeGetErrWeights.
 *
 * lr     is an input flag indicating whether PrecSolve is to use
 *        the left preconditioner P1 or right preconditioner
 *        P2: lr = 1 means use P1, and lr = 2 means use P2.
 *
 * P_data is a pointer to user data - the same as the P_data
 *        parameter passed to CVSpbcg.
 *
 * tmp    is a pointer to memory allocated for an N_Vector
 *        which can be used by PSolve for work space.
 *
 * Returned value:
 * The value to be returned by the PrecSolve function is a flag
 * indicating whether it was successful. This value should be
 *   0 if successful,
 *   positive for a recoverable error (step will be retried),
 *   negative for an unrecoverable error (integration is halted).
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
                                  N_Vector r, N_Vector z,
                                  realtype gamma, realtype delta,
                                  int lr, void *P_data, N_Vector tmp);

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgJacTimesVecFn
 * -----------------------------------------------------------------
 * The user-supplied function jtimes is to generate the product
 * J*v for given v, where J is the Jacobian df/dy, or an
 * approximation to it, and v is a given vector. It should return
 * 0 if successful and a nonzero int otherwise.
 *
 * A function jtimes must have the prototype given below. Its
 * parameters are as follows:
 *
 *   v        is the N_Vector to be multiplied by J.
  *
 *   Jv       is the output N_Vector containing J*v.
 *
 *   t        is the current value of the independent variable.
 *
 *   y        is the current value of the dependent variable
 *            vector.
 *
 *   fy       is the vector f(t,y).
 *
 *   jac_data is a pointer to user Jacobian data, the same as the
 *            pointer passed to CVSpbcg.
 *
 *   tmp      is a pointer to memory allocated for an N_Vector
 *            which can be used by Jtimes for work space.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                    N_Vector y, N_Vector fy,
                                    void *jac_data, N_Vector tmp);

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcg
 * -----------------------------------------------------------------
 * A call to the CVSpbcg function links the main CVODE integrator
 * with the CVSPBCG linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * pretype   is the type of user preconditioning to be done.
 *           This must be one of the four enumeration constants
 *           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
 *           in iterative.h. These correspond to no preconditioning,
 *           left preconditioning only, right preconditioning
 *           only, and both left and right preconditioning,
 *           respectively.
 *
 * maxl      is the maximum Krylov dimension. This is an
 *           optional input to the CVSPBCG solver. Pass 0 to
 *           use the default value CVSPBCG_MAXL=5.
 *
 * The return value of CVSpbcg is one of:
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_MEM_FAIL  if there was a memory allocation failure
 *    CVSPBCG_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

int CVSpbcg(void *cvode_mem, int pretype, int maxl);

/*
 * -----------------------------------------------------------------
 * Function: CVSpbcgSetPrecType
 * -----------------------------------------------------------------
 * CVSpbcgSetPrecType resets the type of preconditioner, pretype,
 *                    from the value set in a prior call to CVSpbcg.
 *                    This must be one of PREC_NONE, PREC_LEFT,
 *                    PREC_RIGHT, or PREC_BOTH.
 * -----------------------------------------------------------------
 */

int CVSpbcgSetPrecType(void *cvode_mem, int pretype);

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPBCG linear solver
 * -----------------------------------------------------------------
 * CVSpbcgSetDelt specifies the factor by which the tolerance on
 *                the nonlinear iteration is multiplied to get a
 *                tolerance on the linear iteration. This is an
 *                optional input to the CVSPBCG solver.
 *                Default value is 0.05.
 * CVSpbcgSetPrecSetupFn specifies the PrecSetup function.
 *                       Default is NULL.
 * CVSpbcgSetPrecSolveFn specifies the PrecSolve function.
 *                       Default is NULL.
 * CVSpbcgSetPrecData specifies a pointer to user preconditioner
 *                    data. This pointer is passed to PrecSetup and
 *                    PrecSolve every time these routines are called.
 *                    Default is NULL.
 * CVSpbcgSetJacTimesVecFn specifies the jtimes function.
 *                         Default is to use an internal finite
 *                         difference approximation routine.
 * CVSpbcgSetJacData specifies a pointer to user Jacobian data.
 *                   This pointer is passed to jtimes every time this
 *                   routine is called.
 *                   Default is NULL.
 *
 * The return value of CVSpbcgSet* is one of:
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_LMEM_NULL if the cvspbcg memory was NULL
 *    CVSPBCG_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

int CVSpbcgSetDelt(void *cvode_mem, realtype delt);
int CVSpbcgSetPrecSetupFn(void *cvode_mem, CVSpbcgPrecSetupFn pset);
int CVSpbcgSetPrecSolveFn(void *cvode_mem, CVSpbcgPrecSolveFn psolve);
int CVSpbcgSetPrecData(void *cvode_mem, void *P_data);
int CVSpbcgSetJacTimesVecFn(void *cvode_mem, CVSpbcgJacTimesVecFn jtimes);
int CVSpbcgSetJacData(void *cvode_mem, void *jac_data);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSPBCG linear solver
 * -----------------------------------------------------------------
 * CVSpbcgGetWorkSpace returns the real and integer work space used
 *                     by CVSPBCG.
 * CVSpbcgGetNumPrecEvals returns the number of preconditioner
 *                        evaluations, i.e. the number of calls made
 *                        to PrecSetup with jok==FALSE.
 * CVSpbcgGetNumPrecSolves returns the number of calls made to
 *                         PrecSolve.
 * CVSpbcgGetNumLinIters returns the number of linear iterations.
 * CVSpbcgGetNumConvFails returns the number of linear
 *                        convergence failures.
 * CVSpbcgGetNumJtimesEvals returns the number of calls to jtimes.
 * CVSpbcgGetNumRhsEvals returns the number of calls to the user
 *                       f routine due to finite difference Jacobian
 *                       times vector evaluation.
 * CVSpbcgGetLastFlag returns the last error flag set by any of
 *                    the CVSPBCG interface functions.
 *
 * The return value of CVSpbcgGet* is one of:
 *    CVSPBCG_SUCCESS   if successful
 *    CVSPBCG_MEM_NULL  if the cvode memory was NULL
 *    CVSPBCG_LMEM_NULL if the cvspbcg memory was NULL
 * -----------------------------------------------------------------
 */

int CVSpbcgGetWorkSpace(void *cvode_mem, long int *lenrwSG, long int *leniwSG);
int CVSpbcgGetNumPrecEvals(void *cvode_mem, long int *npevals);
int CVSpbcgGetNumPrecSolves(void *cvode_mem, long int *npsolves);
int CVSpbcgGetNumLinIters(void *cvode_mem, long int *nliters);
int CVSpbcgGetNumConvFails(void *cvode_mem, long int *nlcfails);
int CVSpbcgGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
int CVSpbcgGetNumRhsEvals(void *cvode_mem, long int *nfevalsSG); 
int CVSpbcgGetLastFlag(void *cvode_mem, int *flag);

/* CVSPBCG return values */

#define CVSPBCG_SUCCESS    0
#define CVSPBCG_MEM_NULL  -1
#define CVSPBCG_LMEM_NULL -2
#define CVSPBCG_ILL_INPUT -3
#define CVSPBCG_MEM_FAIL  -4

#endif

#ifdef __cplusplus
}
#endif
