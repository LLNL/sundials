/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2004-07-27 23:52:40 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/kinsol/LICENSE
 * -----------------------------------------------------------------
 * KINSOL solver module header file
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _kinsol_h
#define _kinsol_h

#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * Enumeration for inputs to KINSetEtaForm (eta choice)
 * -----------------------------------------------------------------
 * ETACONSTANT : use constant value for eta (default value is 0.1
 *               but a different value can be specified via a call
 *               to KINSetEtaConstValue)
 *
 * ETACHOICE1 : use choice #1 as given in Eisenstat and Walker's
 *              paper of SIAM J.Sci.Comput.,17 (1996), pp 16-32,
 *              wherein eta is defined to be:
 *
 *              eta(k+1) = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
 *                       ---------------------------------------------
 *                                       ||F(u_k)||_L2
 *
 *                                                      1+sqrt(5)
 *              eta_safe = eta(k)^ealpha where ealpha = ---------
 *                                                          2
 *
 * ETACHOICE2 : use choice #2 as given in Eisenstat and Walker's
 *              paper wherein eta is defined to be:
 *
 *                                  [ ||F(u_k+1)||_L2 ]^ealpha
 *              eta(k+1) = egamma * [ --------------- ]
 *                                  [  ||F(u_k)||_L2  ]
 *
 *              where egamma = [0,1] and ealpha = (1,2]
 *
 *              eta_safe = egamma*(eta(k)^ealpha)
 *
 *              Note: The default values of the scalar coefficients
 *              egamma and ealpha (both required) are egamma = 0.9
 *              and ealpha = 2.0, but the routine KINSetEtaParams
 *              can be used to specify different values.
 *
 * When using either ETACHOICE1 or ETACHOICE2, if eta_safe > 0.1
 * then the following safeguard is applied:
 *
 *  eta(k+1) = MAX {eta(k+1), eta_safe}
 *
 * The following safeguards are always applied when using either
 * ETACHOICE1 or ETACHOICE2 so that eta_min <= eta <= eta_max:
 *
 *  eta(k+1) = MAX {eta(k+1), eta_min}
 *  eta(k+1) = MIN {eta(k+1), eta_max}
 *
 * where eta_min = 1.0e-4 and eta_max = 0.9 (see KINForcingTerm).
 * -----------------------------------------------------------------
 */

enum { ETACHOICE1 = 1, ETACHOICE2 = 2, ETACONSTANT = 3 }; 

/*
 * -----------------------------------------------------------------
 * Enumeration for global strategy
 * -----------------------------------------------------------------
 * Choices are INEXACT_NEWTON and LINESEARCH.
 * -----------------------------------------------------------------
 */

enum { INEXACT_NEWTON = 1 , LINESEARCH = 2 };
 
/*
 * -----------------------------------------------------------------
 * Type : SysFn
 * -----------------------------------------------------------------
 * The user-supplied subroutine implementing the nonlinear system
 * function (vector-valued function) F must take as input the
 * dependent variable vector uu (type N_Vector), and set fval (type
 * N_Vector) equal to F(uu) before returning. Additional workspace
 * is allocated by the user and referenced by the f_data memory
 * pointer.
 * 
 * Note: The user-defined routine (internally referenced by a
 * a pointer (type SysFn) named func) should have a void return
 * value type.
 * -----------------------------------------------------------------
 */

typedef void (*SysFn)(N_Vector uu, N_Vector fval, void *f_data );

/*
 * -----------------------------------------------------------------
 * User-Callable Routines
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCreate
 * -----------------------------------------------------------------
 * KINCreate allocates and initializes an internal memory block for
 * the KINSOL solver module.
 *
 * If successful, KINCreate returns a pointer to the initialized
 * memory block which should be passed to KINMalloc. If an
 * error occurs, then KINCreate returns a NULL pointer.
 * -----------------------------------------------------------------
 */

void *KINCreate(void);

/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions (KINSOL)
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *     Function Name    |      Optional Input  [Default Value]
 *                      |
 * -----------------------------------------------------------------
 *                      |
 * KINSetFdata          | pointer to user-allocated memory that is
 *                      | passed to the user-supplied subroutine
 *                      | implementing the nonlinear system function
 *                      | F(u)
 *                      | [NULL]
 *                      |
 * KINSetErrFile        | pointer (type FILE) indicating where all
 *                      | warning/error messages should be sent
 *                      | [stderr]
 *                      |
 * KINSetInfoFile       | pointer (type FILE) specifying where
 *                      | informative (non-error) messages should
 *                      | be sent (see KINSetPrintLevel)
 *                      | [stdout]
 *                      |
 * KINSetPrintLevel     | level of verbosity of output:
 *                      |
 *                      |  0  no statistical information is
 *                      |     displayed (default level)
 *                      |
 *                      |  1  for each nonlinear iteration display
 *                      |     the following information: the scaled
 *                      |     norm (L2) of the system function
 *                      |     evaluated at the current iterate, the
 *                      |     scaled norm of the Newton step (only if
 *                      |     using INEXACT_NEWTON), and the number
 *                      |     of function evaluations performed
 *                      |     thus far
 *                      |
 *                      |  2  display level 1 output and the
 *                      |     following values for each iteration:
 *                      |
 *                      |       fnorm (L2) = ||fscale*func(u)||_L2
 *                      |       (only for INEXACT_NEWTON)
 *                      |
 *                      |       scaled fnorm (for stopping) =
 *                      |       ||fscale*ABS(func(u))||_L-infinity
 *                      |       (for INEXACT_NEWTON and LINESEARCH)
 *                      |
 *                      |  3  display level 2 output plus additional
 *                      |     values used by the global strategy
 *                      |     (only if using LINESEARCH), and
 *                      |     statistical information for the linear
 *                      |     solver
 *                      | [0]
 *                      |
 * KINSetNumMaxIters    | maximum number of nonlinear iterations
 *                      | [MXITER_DEFAULT] (defined in kinsol.c)
 *                      |
 * KINSetNoPrecInit     | flag controlling whether or not the
 *                      | KINSol routine makes an initial call
 *                      | to the preconditioner setup routine (pset)
 *                      | (possible values are TRUE and FALSE)
 *                      | [FALSE]
 *                      |
 * KINSetMaxPrecCalls   | maximum number of nonlinear iterations
 *                      | that may be performed between calls to
 *                      | the preconditioner setup routine (pset)
 *                      | [MSBPRE] (defined in kinsol.c)
 *                      |
 * KINSetEtaForm        | flag indicating which method to use to
 *                      | compute the value of the eta coefficient
 *                      | used in the calculation of the linear
 *                      | solver convergence tolerance:
 *                      |
 *                      |  eps = (eta+uround)*||fscale*func(u)||_L2
 *                      |
 *                      | the linear solver tests for convergence by
 *                      | checking if the following inequality has
 *                      | been satisfied:
 *                      |
 *                      |  ||fscale*(func(u)+J(u)*p)||_L2 <= eps
 *                      |
 *                      | where J(u) is the system Jacobian
 *                      | evaluated at the current iterate, and p
 *                      | denotes the Newton step
 *                      |
 *                      | choices for computing eta are as follows:
 *                      |
 *                      |  ETACHOICE1  (refer to KINForcingTerm)
 *                      |
 *                      |  eta = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
 *                      |        ---------------------------------------------
 *                      |                        ||F(u_k)||_L2
 *                      | 
 *                      |  ETACHOICE2  (refer to KINForcingTerm)
 *                      |
 *                      |                [ ||F(u_k+1)||_L2 ]^alpha
 *                      |  eta = gamma * [ --------------- ]
 *                      |                [  ||F(u_k)||_L2  ]
 *                      |
 *                      |  where gamma = [0,1] and alpha = (1,2]
 *                      |
 *                      |  ETACONSTANT  use a constant value for eta
 *                      | [ETACHOICE1]
 *                      |
 * KINSetEtaConstValue  | constant value of eta - use with
 *                      | ETACONSTANT option
 *                      | [0.1]
 *                      |
 * KINSetEtaParams      | values of eta_gamma (egamma) and eta_alpha
 *                      | (ealpha) coefficients - use with ETACHOICE2
 *                      | option
 *                      | [0.9 and 2.0]
 *                      |
 * KINSetNoMinEps       | flag controlling whether or not the value
 *                      | of eps is bounded below by 0.01*fnormtol
 *                      | (see KINSetFuncNormTol)
 *                      |
 *                      |  FALSE  constrain value of eps by setting
 *                      |         to the following:
 *                      |
 *                      |          eps = MAX{0.01*fnormtol, eps}
 *                      |
 *                      |  TRUE  do not constrain value of eps
 *                      | [FALSE]
 *                      |
 * KINSetMaxNewtonStep  | maximum scaled length of Newton step
 *                      | (reset to value of one if user-supplied
 *                      | value is less than one)
 *                      | [1000*||uscale*u_0||_L2]
 *                      |
 * KINSetRelErrFunc     | real scalar equal to realative error in
 *                      | computing F(u) (used in difference-
 *                      | quotient approximation of matrix-vector
 *                      | product J(u)*v)
 *                      | [(uround)^1/2]
 *                      |
 * KINSetFuncNormTol    | real scalar used as stopping tolerance on
 *                      | ||fscale*ABS(func(u))||_L-infinity (see
 *                      | KINStop and KINInitialStop)
 *                      | [(uround)^1/3]
 *                      |
 * KINSetScaledStepTol  | real scalar used as stopping tolerance on
 *                      | the maximum scaled step length:
 *                      |
 *                      |  ||    u_k+1 - u_k    ||
 *                      |  || ----------------- ||_L-infinity
 *                      |  || ABS(u_k+1)+uscale ||
 *                      |
 *                      | (see KINStop)
 *                      | [(uround)^2/3]
 *                      |
 * KINSetConstraints    | pointer to an array (type N_Vector) of
 *                      | constraints on the solution vector u
 *                      | 
 *                      | if constraints[i] =
 *                      |
 *                      |   0  u[i] not constrained
 *                      |
 *                      |  +1  u[i] constrained to be >= 0
 *                      |  -1  u[i] constrained to be <= 0
 *                      |
 *                      |  +2  u[i] constrained to be > 0
 *                      |  -2  u[i] constrained to be < 0
 *                      |
 *                      | if a NULL pointer is given, then no
 *                      | constraints are applied to vector u
 *                      | [NULL]
 *                      |
 * -----------------------------------------------------------------
 * Note: If successful, these functions return SUCCESS. If an
 * argument has an illegal value, then an error message is printed
 * to the file specified by errfp and an error code is returned.
 * -----------------------------------------------------------------
 */

int KINSetFdata(void *kinmem, void *f_data);
int KINSetErrFile(void *kinmem, FILE *errfp);
int KINSetInfoFile(void *kinmem, FILE *infofp);
int KINSetPrintLevel(void *kinmemm, int printfl);
int KINSetNumMaxIters(void *kinmem, long int mxiter);
int KINSetNoPrecInit(void *kinmem, booleantype noPrecInit);
int KINSetMaxPrecCalls(void *kinmem, long int msbpre);
int KINSetEtaForm(void *kinmem, int etachoice);
int KINSetEtaConstValue(void *kinmem, realtype eta);
int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha);
int KINSetNoMinEps(void *kinmem, booleantype noMinEps);
int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep);
int KINSetRelErrFunc(void *kinmem, realtype relfunc);
int KINSetFuncNormTol(void *kinmem, realtype fnormtol);
int KINSetScaledStepTol(void *kinmem, realtype scsteptol);
int KINSetConstraints(void *kinmem, N_Vector constraints);

/*
 * -----------------------------------------------------------------
 * KINSet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSet* subroutines are the
 * following (all but SUCCESS are prefixed by "KINS_"):
 *
 * SUCCESS : means the associated variable was successfully set [0]
 *
 * NO_MEM : means a NULL KINSOL memory block pointer was given
 *          (must call the KINCreate and KINMalloc memory
 *          allocation subroutines prior to calling KINSol) [-1]
 *
 * ILL_INPUT : means the supplied parameter was invalid (check error
 *             message) [-2]
 * -----------------------------------------------------------------
 */

/* error return values for KINSet* functions */
/* Note: SUCCESS = 0*/

enum { KINS_NO_MEM = -1, KINS_ILL_INPUT = -2 };

/*
 * -----------------------------------------------------------------
 * Function : KINMalloc
 * -----------------------------------------------------------------
 * KINMalloc allocates additional memory for vector storage and
 * sets a couple problem-specific KINSOL variables.
 *
 * Note: Additional vectors must be initialized by the user and
 * passed to the KINSol routine.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  func  name of user-supplied subroutine implementing the
 *        nonlinear function F(u)
 *
 *  tmpl  implementation-specific template vector (type N_Vector)
 *        (created using either N_VNew_Serial or N_VNew_Parallel)
 *
 * If successful, KINMalloc returns SUCCESS. If an error occurs,
 * then KINMalloc prints an error message and returns an error
 * code.
 * -----------------------------------------------------------------
 */

int KINMalloc(void *kinmem, SysFn func, N_Vector tmpl);

/*
 * -----------------------------------------------------------------
 * KINMalloc Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINMalloc subroutine are the
 * following (all but SUCCESS are prefixed by "KINM_"):
 *
 * SUCCESS : means the necessary system memory was successfully
 *           allocated [0]
 *
 * NO_MEM : means a NULL KINSOL memory block pointer was given
 *          (must call the KINCreate routine before calling
 *          KINMalloc) [-1]
 *
 * ILL_INPUT : means the name of a user-supplied subroutine
 *             implementing the nonlinear system function F(u)
 *             was not given [-2]
 *
 * MEM_FAIL : means an error occurred during memory allocation
 *            (either insufficient system resources are available
 *            or the vector kernel has not yet been initialized via
 *            a call to NV_SpecInit_Serial or
 *            NV_SpecInit_Parallel) [-3]
 * -----------------------------------------------------------------
 */

/* error return values for KINMalloc */
/* Note: SUCCESS = 0 */

enum { KINM_NO_MEM = -1, KINM_ILL_INPUT = -2, KINM_MEM_FAIL = -3 };

/*
 * -----------------------------------------------------------------
 * Function : KINResetSysFunc
 * -----------------------------------------------------------------
 * KINResetSysFunc changes the user-provided routine which
 * defines the nonlinear problem to be solved.
 *
 * KINResetSysFunc is provided for a user who wishes to solve
 * several problems of the same size but with different
 * functions.
 * -----------------------------------------------------------------
 */

int KINResetSysFunc(void *kinmem, SysFn func);

/*
 * -----------------------------------------------------------------
 * KINResetSysFunc Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINResetSysFunc subroutine are
 * the following (all but SUCCESS are prefixed by "KINS_"):
 *
 * SUCCESS : means the associated variable was successfully set [0]
 *
 * NO_MEM : means a NULL KINSOL memory block pointer was given
 *          (must call the KINCreate and KINMalloc memory
 *          allocation subroutines prior to calling KINSol) [-1]
 *
 * ILL_INPUT : means the supplied parameter was invalid (check error
 *             message) [-2]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINSol
 * -----------------------------------------------------------------
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system. If the initial guess (initial value assigned to vector u)
 * doesn't violate any user-defined constraints, then the subroutine
 * attempts to solve the system F(u) = 0 using a nonlinear Krylov
 * subspace projection method. The Newton-Krylov iterations are
 * stopped if either of the following conditions is satisfied:
 *
 *  ||F(u)||_L-infinity <= 0.01*fnormtol
 *
 *  ||u[i+1] - u[i]||_L-infinity <= scsteptol
 *
 * However, if the current iterate satisfies the second stopping
 * criterion, it doesn't necessarily mean an approximate solution
 * has been found since the algorithm may have stalled, or the
 * user-specified step tolerance (scsteptol) may be too large.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  uu  vector set to initial guess by user before calling KINSol,
 *      but which upon return contains an approximate solution of
 *      the nonlinear system F(u) = 0
 *
 *  strategy  global strategy applied to Newton step if unsatisfactory
 *            (INEXACT_NEWTON or LINESEARCH)
 *
 *  u_scale  vector containing diagonal elements of scaling matrix
 *           for vector u chosen so that the components of
 *           u_scale*u (as a matrix multiplication) all have
 *           about the same magnitude when u is close to a root
 *           of F(u)
 *
 *  f_scale  vector containing diagonal elements of scaling matrix
 *           for F(u) chosen so that the components of
 *           f_scale*F(u) (as a matrix multiplication) all have
 *           roughly the same magnitude when u is not too near a
 *           root of F(u)
 *
 * Note: The components of vectors u_scale and f_scale should be
 * positive.
 *
 * If successful, KINSol returns a positive value (SUCCESS,
 * KINSOL_INITIAL_GUESS_OK or KINSOL_STEP_LT_STPTOL) and vector uu
 * contains an approximate solution of the given nonlinear system.
 * If an error occurs, then an error message is printed and an
 * error code is returned.
 * -----------------------------------------------------------------
 */

int KINSol(void *kinmem, N_Vector uu, int strategy,
	   N_Vector u_scale, N_Vector f_scale);

/*
 * -----------------------------------------------------------------
 * KINSol Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSol subroutine are the
 * following (all but SUCCESS are prefixed by "KINSOL_"):
 *
 * SUCCESS : means ||fscale*ABS(func(u))||_L-infinity <= 0.01*fnormtol
 *           and the current iterate uu is probably an approximate
 *           solution of the nonlinear system F(u) = 0 [0]
 *
 * INITIAL_GUESS_OK : means the initial user-supplied guess already
 *                    satisfies the stopping criterion given above
 *                    [1]
 *
 * STEP_LT_STPTOL : means the following inequality has been
 *                  satisfied (stopping tolerance on scaled
 *                  step length):
 *
 *                    ||    u_k+1 - u_k    ||
 *                    || ----------------- ||_L-infinity <= scsteptol
 *                    || ABS(u_k+1)+uscale ||
 *
 *                  so the current iterate (denoted above by u_k+1)
 *                  may be an approximate solution of the given
 *                  nonlinear system, but it is also quite possible
 *                  that the algorithm is "stalled" (making
 *                  insufficient progress) near an invalid solution,
 *                  or the real scalar scsteptol is too large [2]
 *
 * LNSRCH_NONCONV : means the line search algorithm was unable to
 *                  find an iterate sufficiently distinct from the
 *                  current iterate
 *
 *                  failure to satisfy the sufficient decrease
 *                  condition could mean the current iterate is
 *                  "close" to an approximate solution of the given
 *                  nonlinear system, the finite-difference
 *                  approximation of the matrix-vector product
 *                  J(u)*v is inaccurate, or the real scalar
 *                  scsteptol is too large [-5]
 *
 * MAXITER_REACHED : means the maximum number of nonlinear iterations
 *                   has been reached [-6]
 *
 * MXNEWT_5X_EXCEEDED : means five consecutive steps have been taken
 *                      that satisfy the following inequality:
 *
 *                       ||uscale*p||_L2 > 0.99*mxnewtstep
 *
 *                      where p denotes the current step and
 *                      mxnewtstep is a real scalar upper bound
 *                      on the scaled step length
 *
 *                      such a failure may mean ||fscale*func(u)||_L2
 *                      asymptotes from above to a finite value, or
 *                      the real scalar mxnewtstep is too small [-7]
 *
 * LINESEARCH_BCFAIL : means the line search algorithm (implemented
 *                     in KINLineSearch) was unable to satisfy the
 *                     beta-condition for MXNBCF + 1 nonlinear
 *                     iterations (not necessarily consecutive),
 *                     which may indicate the algorithm is making
 *                     poor progress [-8]
 *
 * KRYLOV_FAILURE : means the Krylov subspace projection method
 *                  failed to converge [-9]
 *
 * PRECONDSET_FAILURE : means the user-supplied routine pset
 *                      (used to compute the preconditioner)
 *                      encountered an unrecoverable error [-10]
 *
 * PRECONDSOLVE_FAILURE : means the user-supplied routine psolve
 *                        (used to solve the preconditioned linear
 *                        system) encountered an unrecoverable
 *                        error [-11]
 *
 * NO_MEM : means a NULL KINSOL memory block pointer was given
 *          (must call the KINCreate and KINMalloc memory
 *          allocation subroutines prior to calling KINSol) [-1]
 *
 * NO_MALLOC : means additional system memory has not yet been
 *             allocated for vector storage (forgot to call the
 *             KINMalloc routine) [-2]
 *
 * INPUT_ERROR : means at least one input parameter was invalid
 *               (check error output) [-3]
 *
 * LSOLV_NO_MEM : means system memory has not yet been allocated
 *                for the linear solver (check linear solver error
 *                code) [-4]
 * -----------------------------------------------------------------
 */

/* KINSol return values */
  
enum { SUCCESS = 0, KINSOL_NO_MEM = -1, KINSOL_NO_MALLOC = -2,
       KINSOL_INPUT_ERROR = -3, KINSOL_LSOLV_NO_MEM = -4,
       KINSOL_INITIAL_GUESS_OK = 1,KINSOL_STEP_LT_STPTOL = 2,
       KINSOL_LNSRCH_NONCONV = -5, KINSOL_MAXITER_REACHED = -6,
       KINSOL_MXNEWT_5X_EXCEEDED = -7, KINSOL_LINESEARCH_BCFAIL = -8,
       KINSOL_KRYLOV_FAILURE = -9, KINSOL_PRECONDSET_FAILURE = -10,
       KINSOL_PRECONDSOLVE_FAILURE = -11 };
 
/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions (KINSOL)
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSOL solver:
 *
 *       Function Name       |      Returned Value
 *                           |
 * -----------------------------------------------------------------
 *                           |
 * KINGetIntWorkSpace        | integer workspace size (total number
 *                           | of long int-sized blocks of memory
 *                           | allocated by KINSOL for vector
 *                           | storage)
 *                           |
 * KINGetRealWorkSpace       | real workspace size (total number of
 *                           | realtype-sized blocks of memory
 *                           | allocated by KINSOL for vector
 *                           | storage)
 *                           |
 * KINGetNumFuncEvals        | total number evaluations of the
 *                           | nonlinear system function F(u)
 *                           | (number of direct calls made to the
 *                           | user-supplied subroutine by KINSOL
 *                           | module member functions)
 *                           |
 * KINGetNumNonlinSolvIters  | total number of nonlinear iterations
 *                           | performed
 *                           |
 * KINGetNumBetaCondFails    | total number of beta-condition
 *                           | failures (see KINLineSearch)
 *                           |
 *                           | KINSOL halts if the number of such
 *                           | failures exceeds the value of the
 *                           | constant MXNBCF (defined in kinsol.c)
 *                           |
 * KINGetNumBacktrackOps     | total number of backtrack operations
 *                           | (step length adjustments) performed
 *                           | by the line search algorithm (see
 *                           | KINLineSearch)
 *                           |
 * KINGetFuncNorm            | scaled norm of the nonlinear system
 *                           | function F(u) evaluated at the
 *                           | current iterate:
 *                           |
 *                           |  ||fscale*func(u)||_L2
 *                           |
 * KINGetStepLength          | scaled norm (or length) of the step
 *                           | used during the previous iteration:
 *                           |
 *                           |  ||uscale*p||_L2
 * -----------------------------------------------------------------
 */

int KINGetIntWorkSpace(void *kinmem, long int *leniw);
int KINGetRealWorkSpace(void *kinmem, long int *lenrw);
int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters);
int KINGetNumFuncEvals(void *kinmem, long int *nfevals);
int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails); 
int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr);
int KINGetFuncNorm(void *kinmem, realtype *fnorm);
int KINGetStepLength(void *kinmem, realtype *steplength);

/*
 * -----------------------------------------------------------------
 * KINGet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSet* subroutines are the
 * following (NO_MEM is prefixed by "KING_"):
 *
 * OKAY : means the information was successfully retrieved [0]
 * 
 * NO_MEM : means a NULL KINSOL memory block pointer was given
 *          (must call the KINCreate and KINMalloc memory
 *          allocation subroutines prior to calling KINSol) [-1]
 * -----------------------------------------------------------------
 */

/* KINGet* return values */

enum { OKAY = 0, KING_NO_MEM = -1 };

/*
 * -----------------------------------------------------------------
 * Function : KINFree
 * -----------------------------------------------------------------
 * KINFree frees system memory resources reserved for the KINSOL
 * solver module.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 * -----------------------------------------------------------------
 */

void KINFree(void *kinmem);

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_linit)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_linit initializes solver-specific data structures (including
 * variables used as counters or for storing statistical information),
 * but system memory allocation should be done by the subroutine
 * that actually initializes the environment for liner solver
 * package. If the linear system is to be preconditioned, then the
 * variable setupNonNull (type booleantype) should be set to TRUE
 * (predefined constant) and the kin_lsetup routine should be
 * appropriately defined.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 * If the necessary variables have been successfully initialized,
 * then the kin_linit function should return LINIT_OK. Otherwise,
 * the subroutine should indicate a failure has occurred by
 * returning LINIT_ERR and printing an informative error message
 * to the kin_mem->kin_errfp pointer.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * kin_linit Return Values
 * -----------------------------------------------------------------
 * The possible return values for the kin_linit subroutine are the
 * following:
 *
 * LINIT_OK : means the kin_linit routine successfully initialized
 *            the linear solver memory block [0]
 *
 * LINIT_ERR : means an error occurred [-1]
 * -----------------------------------------------------------------
 */

/* kin_linit return values */

enum { LINIT_OK = 0, LINIT_ERR = -1 };

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lsetup)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_lsetup interfaces with the user-supplied pset subroutine (the
 * preconditioner setup routine), and updates relevant variable
 * values (see KINSpgmrSetup). Simply stated, the kin_lsetup routine
 * prepares the linear solver for a subsequent call to the
 * user-supplied kin_lsolve function.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 * If successful, the kin_lsetup routine should return 0 (zero). If
 * a recoverable error occurs, then the subroutine should return a
 * positive integer value. However, if an unrecoverable error
 * occurs, then the function should return a negative integer value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lsolve)(KINMem kin_mem, N_Vector xx,
 *                              N_Vector bb, realtype *res_norm)
 * -----------------------------------------------------------------
 * kin_lsolve interfaces with the subroutine implementing the
 * numerical method to be used to solve the linear system J*xx = bb,
 * and must increment the relevant counter variable values in
 * addition to computing certain values used by the global strategy
 * and forcing term routines (see KINInexactNewton, KINLineSearch,
 * KINForcingTerm, and KINSpgmrSolve).
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 *  xx  vector (type N_Vector) set to initial guess by kin_lsolve
 *      routine prior to calling the linear solver, but which upon
 *      return contains an approximate solution of the linear
 *      system J*xx = bb, where J denotes the system Jacobian
 *
 *  bb  vector (type N_Vector) set to -func(u) (negative of the
 *      value of the system function evaluated at the current
 *      iterate) by KINLinSolDrv before kin_lsolve is called
 *
 *  res_norm  holds the value of the L2-norm (Euclidean norm) of
 *            the residual vector upon return
 *
 * If successful, the kin_lsolve routine should return 0 (zero). If
 * a recoverable error occurs, then the subroutine should return a
 * positive integer value. However, if an unrecoverable error
 * occurs, then the function should return a negative integer value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : int (*kin_lfree)(KINMem kin_mem)
 * -----------------------------------------------------------------
 * kin_lfree is called by KINFree and should free (deallocate) all
 * system memory resources allocated for the linear solver module
 * (see KINSpgmrFree).
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINMalloc
 *
 * If successful, the kin_lfree routine should return 0 (zero). If
 * an error occurs, then the subroutine should return a negative
 * integer value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Linear Solver Return Values
 * -----------------------------------------------------------------
 * The subset of possible linear solver error codes that are
 * reported by KINSOL subroutines interfacing with the linear
 * solver module are the following:
 *
 * SUCCESS : means the routine exited normally [0]
 *
 * LIN_NO_MEM : means a NULL KINSOL memory block pointer was
 *              given (must call KINCreate and KINMalloc prior
 *              to initializing the linear solver module) [-2]
 *
 * LIN_ILL_INPUT : means at least one input parameter was invalid
 *                 (check error message(s)) [-1]
 *
 * LIN_NO_LMEM : means system memory has not yet been allocated for
 *               the linear solver module [-3]
 * -----------------------------------------------------------------
 */

/* Note: SUCCESS = 0  */

enum { LIN_ILL_INPUT = -1, LIN_NO_MEM = -2, LIN_NO_LMEM = -3 };

#endif

#ifdef __cplusplus
}
#endif
