/*******************************************************************
 * File          : kinsol.h                                        *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/kinsol/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the interface file for the main KINSol solver           *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _kinsol_h
#define _kinsol_h

#include <stdio.h>
#include "sundialstypes.h"
#include "nvector.h"

/******************************************************************
 * Enumeration for inputs to KINSetEtaForm                        *
 *----------------------------------------------------------------*
 * Eta choice                                                     *
 *                                                                *
 * ETACONSTANT: constant eta, default of 0.1 or user supplied     *
 *               choice, for which see KINSetEtaConstValue        *
 *                                                                *
 * ETACHOICE1:  which uses choice 1 of Eisenstat and Walker's     *
 *              paper of SIAM J.Sci.Comput.,17 (1996), pp 16-32,  *
 *              wherein eta is:                                   *
 *                                                                *
 *  eta(k) =                                                      *
 *    ABS( norm(func(uu(k))) - norm(func(uu(k-1))+J(uu(k-1))*p) ) *
 *                       / norm(func(uu(k-1)))                    *
 *                                                                *
 * ETACHOICE2: which uses choice 2 of Eisenstat and Walker,       *
 *              wherein eta is                                    *
 *                                                                *
 *  eta(k) =                                                      *
 *    egamma * ( norm(func(uu(k))) / norm(func(u(k-1))) )^ealpha  *
 *                                                                *
 *             In choice 2, egamma and ealpha, both required,     *
 *             are either defaults (egamma = 0.9, ealpha = 2)     *
 *             or from  user input, see KINSetEtaParams           *
 *                                                                *
 * For eta(k) determined by either Choice 1 or Choice 2, a value  *
 * eta_safe is determined, and the safeguard                      *
 *              eta(k) <- max(eta_safe,eta(k))                    *
 * is applied to prevent eta from becoming too small too quickly. *
 * For Choice 1,                                                  *
 *              eta_safe = eta(k-1)^((1.+sqrt(5.))/2.)            *
 * and for Choice 2,                                              *
 *              eta_safe = egamma*eta(k-1)^ealpha.                *
 * These safeguards are turned off if they drop below 0.1.        *
 * Also, eta is never allowed to be less than eta_min = 1.e-4     *
 *----------------------------------------------------------------*
 * Global strategy                                                *
 *                                                                *
 * Choices are INEXACT_NEWTON or LINESEARCH                       *
 *----------------------------------------------------------------*/

enum { ETACHOICE1, ETACHOICE2, ETACONSTANT }; 

enum { INEXACT_NEWTON , LINESEARCH };
 
/******************************************************************
 * Type : SysFn                                                   *
 *----------------------------------------------------------------*        
 * The func function which defines the system to be solved :      *
 *              func(uu) = 0      must have type SysFn.           *
 * func takes as input the problem size Neq and the dependent     *
 * variable vector uu. The function stores the result of func(uu) *
 * in fval. The necessary work space, besides uu and fval,  is    *
 * provided by the pointer f_data.                                *
 * The uu argument is of type N_Vector.                           *
 * A SysFn function does not have a return value.                 *
 ******************************************************************/

typedef void (*SysFn)(N_Vector uu, N_Vector fval, void *f_data );

/*================================================================*
 *                                                                *
 *          U S E R - C A L L A B L E   R O U T I N E S           *
 *                                                                *
 *================================================================*/

/*----------------------------------------------------------------*
 * Function : KINCreate                                           *
 *----------------------------------------------------------------*
 * KINCreate creates an internal memory block for a problem to    *
 * be solved by KINSOL.                                           *
 *                                                                *
 * If successful, KINCreate returns a pointer to initialized      *
 * problem memory. This pointer should be passed to KINMalloc.    *
 * If an initialization error occurs, KINCreate prints an error   *
 * message to standard err and returns NULL.                      *
 *----------------------------------------------------------------*/

void *KINCreate(void);

/*----------------------------------------------------------------*
 * Solver optional input specification functions                  *
 *----------------------------------------------------------------*
 * The following functions can be called to set optional inputs   *
 * to values other than the defaults given below:                 *
 *                                                                *
 * Function             |  Optional input / [ default value ]     *
 *                      |                                         * 
 * -------------------------------------------------------------- *
 *                      |                                         * 
 * KINSetFdata          | a pointer to user data that will be     *
 *                      | passed to the user's func function      *
 *                      | every time func is called.              *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * KINSetErrFile        | the file pointer for an error file      *
 *                      | where all KINSOL warning and error      *
 *                      | messages will be written. This parameter*
 *                      | can be stdout (standard output), stderr *
 *                      | (standard error), a file pointer        *
 *                      | (corresponding to a user error file     *
 *                      | opened for writing) returned by fopen.  *
 *                      | If not called, then all messages will   *
 *                      | be written to standard output.          *
 *                      | [NULL]                                  *
 *                      |                                         * 
 * KINSetPrintLevel     | allows user to select from 4 levels  of *
 *                      | output:                                 *
 *                      |                                         * 
 *                      |  0 no statistics printed                *
 *                      |                                         * 
 *                      |  1 output the nonlinear iteration       *
 *                      |    count, the  scaled norm of func(uu), *
 *                      |    and number of func calls.            *
 *                      |                                         * 
 *                      |  2 same as 1 with addition of global    *
 *                      |    strategy statistics:                 *
 *                      |    f1 = 0.5*norm(fscale*func(uu))**2    *
 *                      |  f1new = 0.5*norm(fscale*func(unew))**2 *
 *                      |                                         * 
 *                      |  3 same as 2 with addition of further   *
 *                      |    Krylov iteration statistics.         *
 *                      |                                         *
 *                      | [0]                                     * 
 *                      |                                         * 
 * KINSetNumMaxIters    | maximum allowable number of nonlinear   *
 *                      | iterations.                             *
 *                      | [MXITER_DEFAULT]                        *
 *                      |                                         * 
 * KINSetNoPrecInit     | flag to control the intial call to the  *
 *                      | preconditioner setup routine.           *
 *                      |  FALSE - force the initial call         *
 *                      |  TRUE  - prevent the initial call       *
 *                      | Use the choice TRUE only after a series *
 *                      | of calls with a FALSE value.            *
 *                      | [FALSE]                                 *
 *                      |                                         * 
 * KINSetMaxPrecCalls   | maximum number of steps calling the     *
 *                      | preconditioner solve without calling    *
 *                      | the preconditioner setup                *
 *                      | [10]                                    *
 *                      |                                         * 
 * KINSetEtaForm        | flag indicating which of three methods  *
 *                      | to use for computing eta, the coeff. in *
 *                      | the linear solver convergence tolerance *
 *                      | eps, given by                           *
 *                      |    eps = (eta+u_round)*norm(func(uu))   *
 *                      | Here, all norms are the scaled L2 norm. *
 *                      | The linear solver attempts to produce a *
 *                      | step p such that                        *
 *                      |    norm(func(u) + J(uu)*p) <= eps.      *
 *                      | Two of the methods for computing eta    *
 *                      | calculate it based on the convergence   *
 *                      | process in the routine KINForcingTerm.  *
 *                      | The third method does not require       *
 *                      | calculation; a constant eta is selected *
 *                      | The allowed values are (see above)      *
 *                      | ETACONSTANT, ETACHOICE1, or ETACHOICE2  *
 *                      | [ETACHOICE1]                            *
 *                      |                                         * 
 * KINSetEtaConstValue  | constant value of eta for ETACONSTANT   *
 *                      | [0.1]                                   *
 *                      |                                         * 
 * KINSetEtaParams      | parameter values for eta in the case    *
 *                      | ETACHOICE2: egamma and ealpha           *
 *                      | [0.9  2.0]                              *
 *                      |                                         * 
 * KINSetNoMinEps       | flag to control lower bound on eps,     *
 *                      | the linear solver convergence tolerance *
 *                      |  FALSE - use standard eps min testing   *
 *                      |  TRUE  - remove protection against eps  *
 *                      |          becoming too small             *
 *                      | [FALSE]                                 *
 *                      |                                         * 
 * KINSetMaxNewtonStep  | maximum allowable length of a Newton    *
 *                      | step. Default value is calculated from  *
 *                      | 1000*max(norm(uscale*uu(0),norm(uscale))*
 *                      |                                         * 
 * KINSetRelErrFunc     | relative error in computing func(uu)    *
 *                      | [roundoff unit]                         *
 *                      |                                         * 
 * KINSetFuncNormTol    | a real (scalar) value containing the    *
 *                      | stopping tolerance on                   *
 *                      |    maxnorm( fscale * func(uu) )         *
 *                      | The default value is                    *
 *                      |    (uround) to the 1/3 power            *
 *                      | where uround is the unit roundoff.      *
 *                      |                                         * 
 * KINSetScaledStepTol  | a real (scalar) value containing the    *
 *                      | stopping tolerance on the maximum       *
 *                      | scaled step  uu(k) - uu(k-1).           *
 *                      | The default value is                    *
 *                      |     (uround) to the 2/3 power           *
 *                      |                                         * 
 * KINSetConstraints    | pointer to an array (type N_Vector) of  *
 *                      | constraints on uu.  If the pointer      *
 *                      | passed is NULL, then NO constraints     *
 *                      | are applied to uu.  If constraints(i)   *
 *                      | = +2 or -2, then uu(i) will be          *
 *                      | constrained to be > 0.0 or < 0.0,       *
 *                      | respectively.  If constraints(i) = +1   *
 *                      | or -1, then uu(i) will be constrained   *
 *                      | to be >= 0.0 or <= 0.0, respectively.   *
 *                      | A zero value in constraints(i) implies  *
 *                      | there is no constraint on uu(i).        *
 *                      |                                         * 
 * -------------------------------------------------------------- *
 * If successful, these functions return SUCCESS. If an argument  *
 * has an illegal value, they print an error message to the       *
 * file specified by errfp and return one of the error flags      *  
 * defined below.                                                 *
 *----------------------------------------------------------------*/

int KINSetFdata(void *kinmem, void *f_data);
int KINSetErrFile(void *kinmem, FILE *errfp);
int KINSetPrintLevel(void *kinmemm, int printfl);
int KINSetNumMaxIters(void *kinmem, long int mxiter);
int KINSetNoPrecInit(void *kinmem, booleantype noPrecInit);
int KINSetMaxPrecCalls(void *kinmem, long int msbpre);
int KINSetEtaForm(void *kinmem, int etachoice);
int KINSetEtaConstValue(void *kinmem, realtype eta);
int KINSetEtaParams(void *kinmem, realtype egamma, realtype eaplpha);
int KINSetNoMinEps(void *kinmem, booleantype noMinEps);
int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep);
int KINSetRelErrFunc(void *kinmem, realtype relfunc);
int KINSetFuncNormTol(void *kinmem, realtype fnormtol);
int KINSetScaledStepTol(void *kinmem, realtype scsteptol);
int KINSetConstraints(void *kinmem, N_Vector constraints);

/* Error return values for KINSet* functions */
/* SUCCESS = 0*/
enum {KINS_NO_MEM = -1, KINS_ILL_INPUT = -2};

/******************************************************************
 * Function   : KINMalloc                                         *
 * -------------------------------------------------------------- *
 *                                                                *
 *            This function allocates main memory for the KINSol  *
 *            package. It also allocates several N_Vectors used   *
 *            by the package. Other N_Vectors are also to be      *
 *            allocated by the user and supplied to KINSol.       *
 *                                                                *
 * nvspec  is a pointer to a vector specification structure       *
 *                                                                *
 * If successful, KINMalloc returns SUCCESS.                      *
 * If an initialization error occurs, KINMalloc prints an error   *
 * message to the file specified by errfp and returns an error    *
 * flag                                                           *
 *****************************************************************/

int KINMalloc(void *kinmem, SysFn func, NV_Spec nvspec);

/* Error return values for KINMalloc */
/* SUCCESS = 0 */
enum {KINM_NO_MEM = -1, KINM_ILL_INPUT = -2, KINM_MEM_FAIL = -3};

/******************************************************************
 * Function: KINResetSysFunc                                      *
 *----------------------------------------------------------------*
 * KINResetSysFunc changes the user-provided routine which        *
 *   defines the nonlinear problem to be solved.                  *
 *   KINResetSysFunc is provided for user who wish to solve       *
 *   several problems of the same size but with different         *
 *   functions.                                                   *
 ******************************************************************/

int KINResetSysFunc(void *kinmem, SysFn func);

/******************************************************************
 * Function : KINSol                                              *
 *----------------------------------------------------------------*
 * KINSol initializes memory for a problem previously allocated   *
 * by a call to KINMalloc. It also checks the initial value of uu *
 * (the initial guess) against the constraints and checks if the  *
 * initial guess is a solution of the system. It then attempts to *
 * solve the system func(uu) = 0, where the function func is      *
 * supplied by the user. The Newton-Krylov iterations are stopped *
 * if either func(uu) is smaller in norm than fnormtol(see below),*
 * or the scaled difference between successive iterates is        *
 * smaller in norm than scsteptol (see below).  However, the      *
 * second termination may mean the iterations have stalled at a   *
 * point that is not near a root.                                 *
 * The input arguments for KINSol and their meanings are as       *
 * follows:                                                       *
 *                                                                *
 * kinmem  pointer to KINSol memory block returned by the         *
 *         preceding KINCreate call                               *
 *                                                                *
 * uu      is the solution vector for the system func(uu) = 0.    *
 *         uu is to be set to an initial value if a nonzero       *
 *         vector starting value is desired.                      *
 *                                                                *
 * func    is the system function for the system:   func(uu) = 0. *          
 *                                                                *
 * globalstrategy  is a variable which indicates which global     *
 *         strategy to apply the computed increment delta in the  *
 *         solution uu.  Choices are INEXACT_NEWTON, LINESEARCH.  *
 *                                                                *
 * uscale  is an array (type N_Vector) of diagonal elements of the*
 *         scaling matrix for uu. The elements of uscale must be  *
 *         positive values. The scaling matrix uscale should be   *
 *         chosen so that uscale * uu (as a matrix multiplication)*
 *         should have all its components with roughly the same   *
 *         magnitude when uu is close to a root of func.          *
 *                                                                *
 * fscale  is an array (type N_Vector) of diagonal elements of the*
 *         scaling matrix for func. The elements of fscale must be*
 *         positive values.  The scaling matrix fscale should be  *
 *         chosen so that fscale * func(uu) (as a matrix          *
 *         multiplication) should have all its components with    *
 *         roughly the same magnitude when uu is NOT too near a   *
 *         root of func.                                          *
 *****************************************************************/

int KINSol(void *kinmem, N_Vector u,
           int strategy, N_Vector u_scale, N_Vector f_scale);

/******************************************************************
 *         KINSOL termination codes                               *
 ******************************************************************     
 *  KINSol returns an integer-valued termination code             *
 *                                                                *
 *  The termination values SUCCESS and KINSOL_***** are:          *
 *                                                                *
 *  SUCCESS :    means maxnorm(fscale*func(uu) <= fnormtol, where *
 *               maxnorm() is the maximum norm function N_VMaxNorm*
 *               uu is probably an approximate root of func.      *
 *                                                                *
 *  INITIAL_GUESS_OK: means the initial guess uu has been found   *
 *               to already satisfy the system to the desired     *
 *               accuracy. No calculation was performed other     *
 *               than testing uu.                                 *
 *                                                                *
 *  STEP_LT_STPTOL:  means the scaled distance between the last   *
 *               two steps is less than scsteptol.  uu may be an  *
 *               approximate root of func, but it is also possible*
 *               that the algorithm is making very slow progress  *
 *               and is not near a root or that scsteptol is too  *
 *               large                                            *
 *                                                                *
 *  LNSRCH_NONCONV: means the LineSearch module failed to reduce  *
 *               norm(func) sufficiently on the last global step  *
 *               Either uu is close to a root of f and no more    *
 *               accuracy is possible, or the finite-difference   *
 *               approximation to j*v is inaccurate, or scsteptol *
 *               is too large. Check the outputs ncfl and nni: if *
 *               ncfl is close to nni, it may be the case that the*
 *               Krylov iteration is converging very slowly. In   *
 *               this case, the user may want to use precondition-*
 *               ing and/or increase the maxl value in the        *
 *               KINSpgmr input list (that is,  increase the max  *
 *               dimension of the Krylov subspace by setting maxl *
 *               to nonzero (thus not using the default value of  *
 *               KINSPGMR_MAXL, or if maxl is being set, increase *
 *               its value                                        *
 *                                                                *
 *  MAXITER_REACHED: means that the maximum allowable number of   *
 *               nonlinear iterations has been reached. This is by*
 *               default 200, but may be changed through optional *
 *               input iopt[MXITER].                              *
 *                                                                *
 *  MXNEWT_5X_EXCEEDED: means 5 consecutive steps of length mxnewt*
 *               (maximum Newton stepsize limit) have been taken. *
 *               Either norm(f) asymptotes from above to a finite *
 *               value in some direction, or mxnewt is too small. *
 *               Mxnewt is computed internally (by default) as    *
 *               mxnewt = 1000*max(norm(uscale*uu0),1), where     *
 *               uu0 is the initial guess for uu, and norm() is   *
 *               the Euclidean norm. Mxnewt can be  set by the    *
 *               user through optional input ropt[MXNEWTSTEP].    *
 *                                                                *
 *  LINESEARCH_BCFAIL: means that more than the allowed maximum   *
 *               number of failures (MXNBCF) occurred when trying *
 *               to satisfy the beta condition in the linesearch  *
 *               algorithm. It is likely that the iteration is    *
 *               making poor progress.                            *
 *                                                                *
 *  KRYLOV_FAILURE: means there was a failure of the Krylov       *
 *               iteration process to converge                    *
 *                                                                *
 *  PRECONDSET_FAILURE: means there was a nonrecoverable          *
 *               error in PrecondSet causing the iteration to halt*
 *                                                                *
 *  PRECONDSOLVE_FAILURE: means there was a nonrecoverable        *
 *            error in PrecondSolve causing the iteration to halt.*
 *                                                                *
 *  NO_MEM:    the KINSol memory pointer received was NULL        *
 *                                                                *
 *  NO_MALLOC: the KINSol memory was not allocated                *
 *                                                                *
 *  INPUT_ERROR: one or more input parameters or arrays was in    *
 *               eror. See the listing in errfp for further info  *
 *                                                                *
 *  LSOLV_NO_MEM: The linear solver memory pointer (lmem) was     *
 *             received as NULL. The return value from the linear *
 *             solver needs to be checked and the cause found.    *
 ******************************************************************/

/* KINSol return values */
  
enum { SUCCESS = 0,
       KINSOL_NO_MEM = -1, KINSOL_NO_MALLOC = -2,
       KINSOL_INPUT_ERROR = -3, KINSOL_LSOLV_NO_MEM = -4,
       KINSOL_INITIAL_GUESS_OK = 1,KINSOL_STEP_LT_STPTOL = 2,
       KINSOL_LNSRCH_NONCONV = 3, KINSOL_MAXITER_REACHED = 4,
       KINSOL_MXNEWT_5X_EXCEEDED = 5, KINSOL_LINESEARCH_BCFAIL = 6,
       KINSOL_KRYLOV_FAILURE = 7, KINSOL_PRECONDSET_FAILURE = 8,
       KINSOL_PRECONDSOLVE_FAILURE = 9 };
 
/*----------------------------------------------------------------*
 * Solver optional output extraction functions                    *
 *----------------------------------------------------------------*
 * The following functions can be called to get optional outputs  *
 * and statistics related to the KINSOL solver.                   *
 * -------------------------------------------------------------- *
 * KINGetIntWorkSpace returns the KINSOL integer workspace size   *
 * KINGetRealWorkSpace returns the KINSOL real workspace size     *
 * KINGetNumFuncEvals returns the number of calls to the user's   *
 *       func function                                            *
 * KINGetNumNonlinSolvIters returns the number of nonlinear       *
 *       iterations performed                                     *
 * KINGetNumBetaCondFails returns the total numebr of times the   *
 *       beta condition could not be met in the linesearch        *
 *       algorithm. The nonlinear iteration is halted if this     *
 *       value ever exceeds MXNBCF (10).                          *
 * KINGetNumBacktrackOps returns the number of backtrack          *
 *       operations done in the linesearch algorithm              *
 * KINGetFuncNorm returns the scaled norm at a given iteration:   *
 *       norm(fscale(func(uu))                                    *
 * KINGetStepLength returns the last step length in the global    *
 *       strategy routine (KINLineSearch or KINInexactNewton).    *
 *****************************************************************/

int KINGetIntWorkSpace(void *kinmem, long int *leniw);
int KINGetRealWorkSpace(void *kinmem, long int *lenrw);
int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters);
int KINGetNumFuncEvals(void *kinmem, long int *nfevals);
int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails); 
int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr);
int KINGetFuncNorm(void *kinmem, realtype *fnorm);
int KINGetStepLength(void *kinmem, realtype *steplength);

/* KINGet* return values */
enum { OKAY=0, KING_NO_MEM=-1 };

/******************************************************************
 * Function : KINFree                                             *
 *----------------------------------------------------------------*
 * KINFree frees the problem memory kinsol_mem allocated by       *
 * KINMalloc.  Its only argument is the pointer kinsol_mem        *
 * returned by KINMalloc   .                                      *
 *                                                                *
 ******************************************************************/

void KINFree(void *kinmem);
 
/******************************************************************
 * Types : struct KINMemRec, KINMem                               *
 *----------------------------------------------------------------*
 * The type KINMem is type pointer to struct KINMemRec. This      *
 * structure contains fields to keep track of problem status.     *
 ******************************************************************/

typedef struct KINMemRec {

  realtype kin_uround;    /* machine unit roundoff */

  /* Problem Specification Data */

  SysFn kin_func;              /* = func(uu)= 0.                              */
  void *kin_f_data;            /* ptr to func work space                      */
  realtype kin_fnormtol;       /* ptr to func norm tolerance                  */
  realtype kin_scsteptol;      /* ptr to scaled step tolerance                */
  int kin_globalstrategy;      /* choices are INEXACT_NEWTON & LINESEARCH     */
  int kin_printfl;             /* print (output) option selected              */
  long int kin_mxiter;         /* max number of nonlinear iterations          */
  long int kin_msbpre;         /* max number of iterations without calling the
                                  preconditioner.  this variable set by the
                                  setup routine for the linear solver         */
  int kin_etaflag;             /* eta computation choice                      */
  booleantype kin_noMinEps;
  booleantype kin_precondflag; /* precondition is in use                      */
  booleantype kin_setupNonNull;/*  preconditioning setup routine is non-null
                                   and preconditioning is in use              */
  booleantype kin_constraintsSet; /* if set, user has given a valid constraints
                                     array and constraints will be used       */
  booleantype kin_precondcurrent; /* if set, the preconditioner is current ,
                                     else not                                 */
  booleantype kin_callForcingTerm;/* if set, call the routine KINForcingTerm  */
  realtype kin_mxnewtstep;     /* max allowable step length of a Newton step  */
  realtype kin_sqrt_relfunc;   /* relative error bound for func(uu) 
                                     (sqrt of error used in the code)         */
  realtype kin_stepl;          /* step length of current step (w/scaling)     */
  realtype kin_stepmul;        /* scalar quantity the step was scaled by      */
  realtype kin_eps;            /* current eps value for the iteration         */
  realtype kin_eta;            /* current eta value for the iteration         */
  realtype kin_eta_gamma;      /* gamma value for use in eta calculation      */
  realtype kin_eta_alpha;      /* alpha value for use in eta calculation      */
  booleantype kin_noPrecInit;
  realtype kin_pthrsh;         /* threshold value for calling preconditioner  */

  /* Counters */

  long int  kin_nni;           /* number of nonlinear iterations              */
  long int  kin_nfe;           /* number of func references/calls             */
  long int  kin_nnilpre;       /* nni value at last precond call              */
  long int  kin_nbcf;          /* number of times the beta condition could not 
                                  be met in LineSearch                        */
  long int  kin_nbktrk;        /* number of backtracks                        */
  long int  kin_ncscmx;        /* number of consecutive steps of size
                                  mxnewtstep taken during the last nonlinear
                                  iteration                                   */

  /* Vectors of length Neq */

  N_Vector kin_uu;         /* pointer to user supplied solution vector and
                              current iterate during most of the process  */
  N_Vector kin_unew;       /* pointer to the newly calculated iterate;
                              also bb vector in linear solver             */
  N_Vector kin_fval;       /* pointer to returned func vector             */
  N_Vector kin_uscale;     /* pointer to user-supplied scaling vector 
                              for uu                                      */
  N_Vector kin_fscale;     /* pointer to user-supplied scaling vector 
                              for func                                    */
  N_Vector kin_pp;         /* pointer to the incremental change vector 
                              for uu in this iteration; also vector xx
                              in the linear solver                        */
  N_Vector kin_constraints;/* pointer to user supplied constraints vector */ 
  N_Vector kin_vtemp1;     /* scratch vector  # 1                         */
  N_Vector kin_vtemp2;     /* scratch vector  # 2                         */

  /*  Space requirements for KINSOL */ 

  long int kin_lrw1;        /* no. of realtype words in 1 N_Vector          */ 
  long int kin_liw1;        /* no. of integer words in 1 N_Vector           */ 
  long int kin_lrw;         /* no. of realtype words in KINSOL work vectors */
  long int kin_liw;         /* no. of integer words in KINSOL work vectors  */

  /* Linear Solver Data */
 
  /* Linear Solver functions to be called */

  int (*kin_linit)(struct KINMemRec *kin_mem);

  int (*kin_lsetup)(struct KINMemRec *kin_mem);

  int (*kin_lsolve)(struct KINMemRec *kin_mem, N_Vector xx, N_Vector bb, 
                    realtype *res_norm );
  
  int (*kin_lfree)(struct KINMemRec *kin_mem);
  
  
  /* Linear Solver specific memory */
  
  void *kin_lmem;           

  /* Saved Values */

  realtype kin_fnorm;    /* current value for the norm of func(uu)            */
  realtype kin_f1norm;   /* current value for the expression:  fnorm*fnorm/2  
                            NOTE: the value f1normp is computed in
                            KINLineSearch or KINInexactNewton and supplied
                            to the calling routine to set this value.
                            f1normp remains a variable of scope to KINSol.    */
  realtype kin_res_norm; /* current value for the norm of the residual        */
  realtype kin_sfdotJp;  /* scaled f value dotted with J p, where p is the
                            computed increment in the computed solution.
                            This term is used in the global strategy routine
                            and in KINForcingTerm.                            */
  realtype kin_sJpnorm;  /* norm of scaled J p, as above, also used in
                            the KINForcingTerm routine                        */
  
  /* In the last two descriptions above, J is the Jacobian matrix evaluated at
     the last iterate u, and p is the current increment from the last iterate u.
     In both cases, p is scaled by a factor lambda (rl in KINLineSearch) to 
     represent the actual step chosen by the global strategy routine.         */


  booleantype kin_MallocDone;

  /* Message File */
  
  FILE *kin_errfp;  /* where KINSol error/warning/info messages are sent */

  /* Pointer to vector specification strucutre */
  
  NV_Spec kin_nvspec;

} *KINMem;


/******************************************************************
 * Communication between user and a KINSOL Linear Solver           *
 *----------------------------------------------------------------*
 * Return values of the linear solver specification routine.      *
 * The values of these are given in the enum statement below.     *
 *                                                                *
 *    SUCCESS      : The routine was successful.                  *
 *                                                                *
 *    LIN_NO_MEM   : KINSOL memory = NULL.                        *
 *                                                                *
 *    LMEM_FAIL    : A memory allocation failed.                  *
 *                                                                *
 *    LIN_ILL_INPUT: Some input was illegal (see message).        *
 *                                                                *
 *    LIN_NO_LMEM  : The linear solver's memory = NULL.           *
 *----------------------------------------------------------------*/

/* SUCCESS = 0  */
enum {LMEM_FAIL=-1, LIN_ILL_INPUT=-2, LIN_NO_MEM=-3, LIN_NO_LMEM=-4};

/******************************************************************
 * Communication between kinsol.c and a KINSol Linear Solver      *
 *----------------------------------------------------------------*
 * (1) kin_linit return values                                    *
 *                                                                *
 * LINIT_OK    : The kin_linit routine succeeded.                 *
 *                                                                *
 * LINIT_ERR   : The kin_linit routine failed. Each linear solver *
 *               init routine should print an appropriate error   *
 *               message to (kin_mem->kin_errfp).                 *
 ******************************************************************/

/* kin_linit return values */

#define LINIT_OK        0
#define LINIT_ERR      -1

/*******************************************************************
 * int (*kin_linit)(KINMem kin_mem);                               *
 *-----------------------------------------------------------------*
 * The purpose of kin_linit is to perform any needed               *
 * initializations of solver-specific memory, such as              *
 * counters/statistics. Actual memory allocation should be done    *
 * by the routine that initializes the linear solver package.      *
 * The kin_linit routine should set  setupNonNull to be TRUE if    * 
 * the setup operation for the linear solver is non-empty and      *
 * FALSE if the setup operation does nothing. An LInitFn should    *
 * return LINIT_OK (= 0) if it has successfully initialized the    *
 * KINSol linear solver and LINIT_ERR (= -1) otherwise. These      *
 * constants are defined above. If an error does occur, an         *
 * appropriate message should be sent to (kin_mem->kin_errfp).     *
 *******************************************************************/

/*******************************************************************
 * int (*kin_lsetup)(KINMem kin_mem);                              *
 *-----------------------------------------------------------------*
 * The job of kin_lsetup is to prepare the linear solver for       *
 * subsequent calls to kin_lsolve.                                 *
 *                                                                 *
 * kin_mem - problem memory pointer of type KINMem. See the big    *
 *          typedef earlier in this file.                          *
 *                                                                 *
 * The kin_lsetup routine should return 0 if successful,           *
 * a positive value for a recoverable error, and a negative value  *
 * for an unrecoverable error.                                     *
 *******************************************************************/

/*******************************************************************
 * int (*kin_lsolve)(KINMem kin_mem, N_Vector bb, N_Vector xx,     *
 *                   realtype *res_norm);                          *
 *-----------------------------------------------------------------*
 * kin_lsolve must solve the linear equation J x = b, where        *
 * J is an approximate Jacobian matrix, x is the approximate system*
 * solution,  and the RHS vector b is input. The solution is to be *
 * returned in the vector b. kin_lsolve returns a positive value   *
 * for a recoverable error and a negative value for an             *
 * unrecoverable error. Success is indicated by a 0 return value.  *
 *******************************************************************/

/*******************************************************************
 * void (*kin_lfree)(KINMem kin_mem);                              *
 *-----------------------------------------------------------------*
 * kin_lfree should free up any memory allocated by the linear     *
 * solver. This routine is called once a problem has been          *
 * completed and the linear solver is no longer needed.            *
 *******************************************************************/

#endif
#ifdef __cplusplus
}
#endif
