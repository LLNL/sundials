/******************************************************************
 *                                                                *
 * File          : kinspgmr.h                                     *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 27 June 2002                                   *
 *----------------------------------------------------------------*
 * This is the header file for the KINSol scaled, preconditioned  *
 * GMRES linear solver, KINSpgmr.                                 *
 *                                                                *
 * Note: The type integertype must be large enough to store the   *
 * value of the linear system size Neq.                           *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _kinspgmr_h
#define _kinspgmr_h


#include <stdio.h>
#include "kinsol.h"   /*  for KINSOL_IOPT_SIZE, etc.  */
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"

 
/******************************************************************
 *                                                                *
 * KINSpgmr solver statistics indices                             *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * KINSpgmr-specific statistic. The symbolic names are used as    *
 * indices into the iopt and ropt arrays and values of both arrays*
 * are set in this module                                         *
 *                                                                *
 * The KINSpgmr statistics are:                                   *
 *                                                                *
 * iopt[SPGMR_NLI]  (output) number of linear iterations.         *
 *                                                                *
 * iopt[SPGMR_NPE]  (output) number of preconditioner evaluations *
 *                                                                *
 * iopt[SPGMR_NPS]  (output) number of calls made to user's psolve*
 *                    function.                                   *
 *                                                                *
 * iopt[SPGMR_NCFL] (output) number of linear convergence failures*
 *                                                                *
 ******************************************************************/
 
enum { SPGMR_NLI=KINSOL_IOPT_SIZE, SPGMR_NPE, SPGMR_NPS, SPGMR_NCFL };


/******************************************************************
 *                                                                *
 * KINSpgmr solver constants                                      *
 *----------------------------------------------------------------*
 * KINSPGMR_MAXL    : default value for the maximum Krylov        *
 *                   dimension is MIN(N, KINSPGMR_MAXL)           *
 *                                                                * 
 * KINSPGMR_MSBPRE  : maximum number of steps between             *
 *                   preconditioner evaluations                   *
 *                                                                *
 ******************************************************************/

#define KINSPGMR_MAXL    10

#define KINSPGMR_MSBPRE  10 

   /* Constants for error returns from KINSpgmr. */

#define KIN_MEM_NULL      -1
#define KINSPGMR_MEM_FAIL -2
#define SPGMR_MEM_FAIL    -3
 
 
/******************************************************************
 *                                                                *           
 * Type : KINSpgmrPrecondFn                                       *
 *----------------------------------------------------------------*
 * The user-supplied preconditioner setup function precondset and *
 * the user-supplied preconditioner solve function precondsolve   *
 * together must define the right preconditoner matrix P chosen   *
 * so as to provide an easier system for the Krylov solver        *
 * to solve. precondset is called to provide any matrix data      *
 * required by the subsequent call(s) to precondsolve. The data is*
 * stored in the memory allocated to P_data and the structuring of*
 * that memory is up to the user.     More specifically,          *
 * the user-supplied preconditioner setup function precondset     *
 * is to evaluate and preprocess any Jacobian-related data        *
 * needed by the preconditioner solve function precondsolve.      *
 * This might include forming a crude approximate Jacobian,       *
 * and performing an LU factorization on the resulting            *
 * approximation to J.  This function will not be called in       *
 * advance of every call to precondsolve, but instead will be     *
 * called only as often as necessary to achieve convergence       *
 * within the Newton iteration in KINSol.  If the precondsolve    *
 * function needs no preparation, the precondset function can be  *
 * NULL.                                                          *
 *                                                                *
 * precondset should not modify the contents of the arrays        *
 * uu or fval as those arrays are used elsewhere in the           *
 * iteration process.                                             *
 *                                                                *
 * Each call to the precondset function is preceded by a call to  *
 * the system function func. Thus the precondset function can use *
 * any auxiliary data that is computed by the func function and   *
 * saved in a way accessible to precondset.                       *
 *                                                                *
 * The two scaling arrays, fscale and uscale, and unit roundoff   *
 * uround are provided to the precondset function for possible use*
 * in approximating Jacobian data, e.g. by difference quotients.  *
 * These arrays should also not be altered                        *
 *                                                                *
 * A function precondset must have the prototype given below.     *
 * Its parameters are as follows:                                 *
 *                                                                *
 * Neq     is the length of all vector arguments.                 *
 *                                                                *
 * uu      an N_Vector giving the current iterate for the system. *
 *                                                                *
 * uscale  an N_Vector giving the diagonal entries of the uu-     *
 *         scaling matrix.                                        *
 *                                                                *
 * fval    an N_Vector giving the current function value          *
 *                                                                *
 * fscale  an N_Vector giving the diagonal entries of the func-   *
 *         scaling matrix.                                        *
 *                                                                *
 * vtemp1, vtemp2  two N_Vector temporaries for use by preconset  *
 *                                                                *
 * func    the function func defines the system being solved:     *
 *         func(uu) = 0., and its name is passed initially to     *
 *         KINSol in the call to KINMalloc                        *
 *                                                                *
 * uround  is the machine unit roundoff.                          *
 *                                                                *
 * nfePtr  is a pointer to the memory location containing the     *
 *           KINSol problem data nfe = number of calls to func.   *
 *           The precondset routine should update this counter by *
 *           adding on the number of func calls made in order to  *
 *           approximate the Jacobian, if any.  For example, if   *
 *           the routine calls func a total of W times, then the  *
 *           update is *nfePtr += W.                              *
 *                                                                *
 * P_data is a pointer to user data - the same as the P_data      *
 *           parameter passed to KINSpgmr.                        *
 *                                                                *
 *                                                                *
 * Returned value:                                                *
 * The value to be returned by the precondset function is a flag  *
 * indicating whether it was successful.  This value should be    *
 *   0  if successful,                                            *
 *   1  if failure, in which case KINSol stops                    *
 *                                                                *
 *                                                                *
 ******************************************************************/

typedef int (*KINSpgmrPrecondFn)(integertype Neq, 
                                 N_Vector uu, N_Vector uscale ,
                                 N_Vector fval, N_Vector fscale,
                                 N_Vector vtemp1, N_Vector vtemp2,
                                 SysFn func, realtype uround,
                                 long int *nfePtr, void *P_data);


/******************************************************************
 *                                                                *           
 * Type : KINSpgmrPrecondSolveFn                                  *
 *----------------------------------------------------------------*
 * The user-supplied preconditioner solve function precondsolve   *
 * is to solve a linear system P x = r in which the matrix P is   *
 * the (right) preconditioner matrix P.                           *
 *                                                                *
 * precondset should not modify the contents of the iterate       *
 * array uu  or the current function value array  fval as those   *
 * are used elsewhere in the iteration process.                   *
 *                                                                *
 * A function precondsolve must have the prototype given below.   *
 * Its parameters are as follows:                                 *
 *                                                                *
 * Neq     system size, and length of all vector arguments.       *
 *                                                                *
 * uu      an N_Vector giving the current iterate for the system. *
 *                                                                *
 * uscale  an N_Vector giving the diagonal entries of the uu-     *
 *         scaling matrix.                                        *
 *                                                                *
 * fval    an N_Vector giving the current function value          *
 *                                                                *
 * fscale  an N_Vector giving the diagonal entries of the func-   *
 *         scaling matrix.                                        *
 *                                                                *
 * vtem    an N_Vector to hold the RHS vector r on input and      *
 *         the result vector x on return                          *
 *                                                                *
 * ftem    an N_Vector work array, usually set on input as vtemp  *
 *                                                                *
 * func    the function func defines the system being solved:     *
 *         func(uu) = 0.                                          *
 *                                                                *
 * uround  is the machine unit roundoff.                          *
 *                                                                *
 * nfePtr is a pointer to the memory location containing the      *
 *          KINSol problem data nfe = number of calls to func. The*
 *          precondsolve routine should update this counter by    *
 *          adding  on the number of func calls made in order to  *
 *          carry out the solution, if any.  For example, if the  *
 *          routine calls func a total of W times, then the update*
 *          is  *nfePtr += W.                                     *
 *                                                                *
 * P_data is a pointer to user data - the same as the P_data      *
 *          parameter passed to KINSpgmr.                         *
 *                                                                *
 * Returned value:                                                *
 * The value to be returned by the precondsolve function is a flag*
 * indicating whether it was successful.  This value should be    *
 *   0 if successful,                                             *
 *   1 if failure, in which case KINSol stops                     *
 *                                                                *
 ******************************************************************/
  
typedef int (*KINSpgmrPrecondSolveFn)(integertype Neq,
                                      N_Vector uu, N_Vector uscale, 
                                      N_Vector fval, N_Vector fscale, 
                                      N_Vector vv, N_Vector ftem,
                                      SysFn func, realtype u_round,
                                      long int *nfePtr, void *P_data);


/********************************************************************
 *                                                                  *
 * type : KINSpgmruserAtimesFn                                      *
 *                                                                  *
 *    The user-supplied A times v routine (optional) where A is the *
 *    Jacobian matrix dF/du, or an approximation to it, and v is a  *
 *    given  vector.  This routine computes the product z = J v.    *
 *    It should return 0 if successful and a nonzero int otherwise. *
 *                                                                  *
 *   f_data is a pointer to the structure used to handle data for   *
 *     the user-supplied system function and also to contain data   *
 *     for evaluation of the jacobian of func                       *
 *                                                                  *
 *   v is the N_Vector to be multiplied by J                        *
 *     (preconditioned and unscaled as received)                    *
 *                                                                  *
 *   z is the N_Vector resulting from the application of J to v     *
 *                                                                  * 
 *   new_uu   is an input flag indicating whether or not the uu     *
 *     vector has been changed since the last call to this function.*
 *     If this function computes and saves Jacobian data, then      *
 *     this computation can be skipped if new_uu = FALSE.           *
 *                                                                  *
 *   uu  is the N_Vector equal to the current iterate u             *
 *                                                                  *
 ********************************************************************/

typedef int (*KINSpgmruserAtimesFn)(void *f_data, N_Vector v, 
                                    N_Vector z, booleantype *new_uu, 
                                    N_Vector uu);
  
 
/******************************************************************
 *                                                                *
 * Function : KINSpgmr                                            *
 *----------------------------------------------------------------*
 * A call to the KINSpgmr function links the main KINSol solver   *
 * with the KINSpgmr linear solver. Among other things, it sets   *
 * the generic names linit, lsetup, lsolve, and lfree to the      *
 * specific names for this package:                               *
 *                   KINSpgmrInit                                 *
 *                          KINSpgmrSetup                         *
 *                                  KINSpgmrSolve                 *
 *                                              KINSpgmrFree      *
 *                                                                *
 * kin_mem is the pointer to KINSol memory returned by            *
 *             KINSolMalloc.                                      *
 *                                                                *
 * maxl      is the maximum Krylov dimension. This is an          *
 *             optional input to the KINSpgmr solver. Pass 0 to   *
 *             use the default value MIN(Neq, KINSPGMR_MAXL=10).  *
 *                                                                *
 * maxlrst   is the maximum number of linear solver restarts      *
 *            allowed. Values outside the range 0 to 2*Neq/maxl   *
 *            will be restricted to that range. 0, meaning no     *
 *            restarts, is a safe starting value.                 *
 *                                                                *
 * msbpre    is the maximum number of steps calling the solver    *
 *           precondsolve without calling the preconditioner      *
 *           precondset. (The default is KINSPGMR_MSBPRE = 10)    *
 *                                                                *
 * precondset  is the user's preconditioner routine. It is used to*
 *             evaluate and preprocess any Jacobian-related data  *
 *             needed by the precondsolve routine.  See the       *
 *             documentation for the type KINSpgmrPrecondFn for   *
 *             full details.  Pass NULL if no such setup of       *
 *             Jacobian data is required.  A precond routine is   *
 *             NOT required, but rather provided when needed by   *
 *             user's precondsolve routine                        *
 *                                                                *
 * precondsolve  is the user's preconditioner solve routine. It   *
 *             is used to solve Px=b, where P is a preconditioner *
 *             matrix.  See the documentation for the type        *
 *             KINSpgmrPrecondSolveFn for full details.  The only *
 *              case in which psolve is allowed to be NULL is when*
 *             no preconditioning is to be done. The NULL is taken*
 *             as a flag that preconditioning is not desired.     *
 *                                                                *
 * userAtimes  is an optional routine supplied by the user to     *
 *             perform the matrix-vector multiply J v, where J is *
 *             an approximate Jacobian matrix for that iteration. *
 *             Enter NULL if no such routine is required. If one  *
 *             is supplied, conforming to the definitions given   *
 *             in this file, enter its function name.             *
 *                                                                *
 * P_data      is a pointer to user preconditioner data. This     *
 *             pointer is passed to precondset and precondsolve   *
 *             every time these routines are called.              *
 *                                                                *
 *       KINSpgmr returns an integer that is one of the following *
 *       defined constants:                                       *
 *            SUCCESS, KIN_MEM_NULL, KINSPGMR_MEM_FAIL, and       *
 *            SPGMR_MEM_FAIL.                                     *
 *                                                                *
 ******************************************************************/
  
int KINSpgmr(void *kin_mem, int maxl, int maxlrst, int msbpre,
             KINSpgmrPrecondFn precondset, 
             KINSpgmrPrecondSolveFn precondsolve,
             KINSpgmruserAtimesFn userAtimes,
             void *P_data);
  
#endif
#ifdef __cplusplus
}
#endif
