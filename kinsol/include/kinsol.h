/*******************************************************************
 *                                                                 *
 * File          : kinsol.h                                        *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 25 July 2002                                    *
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
 *                                                                *
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
 *                                                                *
 ******************************************************************/

typedef void (*SysFn)(integertype Neq, N_Vector uu, N_Vector fval,
                      void *f_data );

/******************************************************************
 *                                                                *
 * Function   : KINMalloc                                         *
 *                                                                *
 *            This function allocates main memory for the KINSol  *
 *            package. It also allocates several vectors of size  *
 *            Neq used by the package. Other N_Vectors are also   *
 *            to be allocated by the user and supplied to KINSol. *
 *                                                                *
 * Neq     size of nonlinear system, and size of vectors being    *
 *         handled by the current allocation call to KINMalloc.   *
 *                                                                *
 * msgfp   pointer to a FILE used to receive error messages from  *
 *         KINMalloc.  Pass NULL to direct messages to stdout.    *
 *                                                                *
 * machEnv is a pointer to machine environment-specific           *
 *         information. Pass NULL for the sequential case         *
 *         (see nvector.h).                                       *
 *                                                                *
 * If successful, KINMalloc returns a pointer to initialized      *
 * problem memory. This pointer should be passed to KINSol. If    *
 * an initialization error occurs, KINMalloc prints an error      *
 * message to the file specified by msgfp and returns NULL.       *
 *                                                                *
 *****************************************************************/

void *KINMalloc(integertype Neq, FILE *msgfp, M_Env machEnv);


#define KINSOL_IOPT_SIZE 10
#define KINSOL_ROPT_SIZE 8
#define OPT_SIZE        40
 
/******************************************************************
 *                                                                *
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
 * Neq     is the number of equations in the algebraic system     *
 *                                                                *
 * kinmem  pointer to KINSol memory block returned by the         *
 *         preceding KINMalloc call                               *
 *                                                                *
 * uu      is the solution vector for the system func(uu) = 0.    *
 *         uu is to be set to an initial value if a nonzero       *
 *         vector starting value is desired.                      *
 *                                                                *
 * func    is the system function for the system:   func(uu) = 0. *          
 *                                                                *
 * globalstrategy  is a variable which indicates which global     *
 *         strategy to apply the computed increment delta in the  *
 *         solution uu.  Choices are INEXACT_NEWTON, LINESEARCH,  *
 *         as given in the enum statement below.                  *
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
 *                                                                *
 * fnormtol  is a real (scalar) value containing the stopping     *
 *         tolerance on maxnorm( fscale * func(uu) ) .            *
 *         If fnormtol is input as 0, then a default value of     *
 *         (uround) to the 1/3 power will be used.                *
 *         uround is the unit roundoff for the machine            *
 *         in use for the calculation. (See UnitRoundoff in       *
 *         sundialsmath module.)                                  *
 *                                                                *
 * scsteptol  is a real (scalar) value containing the stopping    *
 *         tolerance on the maximum scaled step  uu(k) - uu(k-1). *
 *         If scsteptol is input as 0., then a default value of   *
 *         (uround) to the 2/3 power will be used.                *
 *         uround is the unit roundoff for the machine            *
 *         in use for the calculation. (see UnitRoundoff in       *
 *         sundialsmath module                                        *
 *                                                                *
 * constraints  is a pointer to an array (type N_Vector) of       *
 *         constraints on uu .  If the pointer passed in is NULL, *
 *         then NO constraints are applied to uu . A NULL pointer *
 *         also stops application of the constraint on the        *
 *         maximum relative change in uu, controlled by the input *
 *         variable relu which is input via ropt[RELU].           *
 *         A positive value in constraints[i] implies that the    *
 *         i-th* component of uu is to be constrained > 0.        *
 *         A negative value in constraints[i] implies that the    *
 *         i-th component of uu is to be constrained < 0.         *
 *         A zero value in constraints[i] implies there is no     *
 *         constraint on uu[i].                                   *
 *                                                                *
 * optIn   is a flag indicating whether optional inputs from the  *
 *         user in the arrays iopt and ropt are to be used        *
 *         pass FALSE to ignore all optional inputs and TRUE      *
 *         to use all optional inputs that are present.           *
 *         Either choice does NOT affect outputs in other         *
 *         positions of iopt or ropt.                             *
 *                                                                *
 * iopt    is the user-allocated array (of size OPT_SIZE) that    *
 *         will hold optional integer inputs and outputs.         *
 *         The user can pass NULL if he/she does not              *
 *         wish to use optional integer inputs or outputs.        *
 *         If optIn is TRUE, the user should preset to 0 those    *
 *         locations for which default values are to be used.     *
 *         See Optional Inputs and Outputs, below.                *
 *                                                                *
 * ropt    is the user-allocated array (of size OPT_SIZE) that    *
 *         will hold optional real inputs and outputs.            *
 *         The user can pass NULL if he/she does not              *
 *         wish to use optional real inputs or outputs.           *
 *         If optIn is TRUE, the user should preset to 0.0 the    *
 *         optional input locations for which default values are  *
 *         to be used.                                            *
 *         See Optional Inputs and Outputs, below.                *
 *                                                                *
 * f_data  is a pointer to work space for use by the user-supplied*
 *         function func. The space allocated to f_data is        *
 *         allocated by the user's program before the call to     *
 *         KINMalloc.                                             *
 *                                                                *
 *****************************************************************/

enum { INEXACT_NEWTON , LINESEARCH };  /* globalstrategy */
 
/******************************************************************
 *                                                                *
 * Optional Inputs and Outputs                                    *
 *----------------------------------------------------------------*
 * The user should declare two arrays for optional input and      *
 * output, an iopt array for optional integer input and output    *
 * and an ropt array for optional real input and output. These    *
 * arrays should both be of size OPT_SIZE.                        *
 * So the user's declaration should look like:                    *
 *                                                                *
 * long int iopt[OPT_SIZE];                                       *
 * realtype ropt[OPT_SIZE];                                       *
 *                                                                *
 * The following definitions are indices into the iopt and ropt   *
 * arrays. A brief description of the contents of these positions *
 * follows.                                                       *
 *                                                                *
 * iopt[PRINTFL]    (input)  allows user to select from 4 levels  *
 *                  of output to FILE msgfp.                      *
 *                  =0 no statistics printed   (DEFAULT)          *
 *                  =1 output the nonlinear iteration count, the  *
 *                     scaled norm of func(uu), and number of     *
 *                     func calls.                                *
 *                  =2 same as 1 with the addition of global      *
 *                     strategy statistics:                       *
 *                     f1 = 0.5*norm(fscale*func(uu))**2   and    *
 *                     f1new = 0.5*norm(fscale*func(unew))**2 .   *
 *                  =3 same as 2 with the addition of further     *
 *                     Krylov iteration statistics.               *
 *                                                                *
 * iopt[MXITER]     (input) maximum allowable number of nonlinear *
 *                   iterations. The default is MXITER_DEFAULT.   *
 *                                                                *
 * iopt[PRECOND_NO_INIT] (input) Set to 1 to prevent the initial  *
 *                    call to the routine precondset upon a given *
 *                    call to KINSol. Set to 0 or leave unset to  *
 *                    force the initial call to precondset.       *
 *                    Use the choice of 1 only after beginning the*
 *                    first of a series of calls with a 0 value.  *
 *                    If a value other than 0 or 1 is encountered,*
 *                    this element is set to the default, 0.      *
 *                                                                *
 * iopt[ETACHOICE]   (input) a flag indicating which of three     *
 *                    methods to use for computing eta, the       *
 *                    coefficient in the linear solver            *
 *                    convergence tolerance eps, given by         *
 *                      eps = (eta+u_round)*norm(func(uu)) .      *
 *                    Here, all norms are the scaled L2 norm.     *
 *                    The linear solver attempts to produce a step*
 *                    p such that norm(func(u) + J(uu)*p) <= eps. *
 *                    Two of the methods for computing eta        *
 *                    calculate a value based on the convergence  *
 *                    process in the routine KINForcingTerm.      *
 *                    The third method does not require           *
 *                    calculation; a constant eta is selected.    *
 *                                                                *
 *                    The default if iopt[ETACHOICE] is not       *
 *                    specified is ETACHOICE1 (see below).        *
 *                                                                *
 *                    The allowed values and their meanings are:  *
 *               ETACONSTANT: constant eta, default of 0.1 or user*
 *                  supplied choice, for which see ropt[ETACONST].*
 *                                                                *
 *               ETACHOICE1: [default] which uses choice 1 of     *
 *                  Eisenstat and Walker's paper of SIAM J. Sci.  *
 *                  Comput.,17 (1996), pp 16-32, wherein eta is:  *
 *                          eta(k) =                              *
 *   ABS( norm(func(uu(k))) - norm(func(uu(k-1))+J(uu(k-1))*p) )  *
 *                       / norm(func(uu(k-1)))                    *
 *                                                                *
 *               ETACHOICE2:  which uses choice 2 of              *
 *                  Eisenstat and Walker, wherein eta is:         *
 *                  eta(k) = egamma *                             *
 *              ( norm(func(uu(k))) / norm(func(u(k-1))) )^ealpha *
 *                                                                *
 *                  In choice 2, egamma and ealpha, both required,*
 *                  are either defaults (egamma = 0.9, ealpha = 2)*
 *                  or from  user input, see ropt[ETAALPHA] and   *
 *                  ropt[ETAGAMMA] below.                         *
 *                                                                *
 *                  For eta(k) determined by either Choice 1 or   *
 *                  Choice 2, a value eta_safe is determined, and *
 *                  the safeguard   eta(k) <- max(eta_safe,eta(k))*
 *                  is applied to prevent eta(k) from becoming too*
 *                  small too quickly.                            *
 *                   For Choice 1,                                *
 *                     eta_safe = eta(k-1)^((1.+sqrt(5.))/2.)     *
 *                  and for Choice 2,                             *
 *                     eta_safe = egamma*eta(k-1)^ealpha.         *
 *                  (These safeguards are turned off if they drop *
 *                  below 0.1.  Also, eta is never allowed to be  *
 *                  less than eta_min = 1.e-4.)                   *
 *                                                                *
 * iopt[NO_MIN_EPS] (input) Set to 1 or greater to remove         *
 *                  protection agains eps becoming too small.     *
 *                  This option is useful for debugging linear    *
 *                  and nonlinear solver interactions. Set to 0   *
 *                  for standard eps minimum value testing.       *
 *                                                                *
 * iopt[NNI]        (output) total number of nonlinear iterations.*
 *                                                                *
 * iopt[NFE]        (output) total number of calls to the user-   *
 *                   supplied system function func.               *
 *                                                                *
 * iopt[NBCF]       (output) total number of times the beta       *
 *                   condition could not be met in the linesearch *
 *                   algorithm. The nonlinear iteration is halted *
 *                   if this value ever exceeds MXNBCF (10).      *
 *                                                                *
 * iopt[NBKTRK]     (output) total number of backtracks in the    *
 *                   linesearch algorithm.                        *
 *                                                                *
 * ropt[MXNEWTSTEP] (input) maximum allowable length of a Newton  *
 *                   step. The default value is calculated from   *
 *                   1000*max(norm(uscale*uu(0),norm(uscale)).    *
 *                                                                *
 * ropt[RELFUNC]    (input) relative error in computing func(uu)  *
 *                   if known. Default is the machine epsilon.    *
 *                                                                *
 * ropt[RELU]       (input) a scalar constraint which restricts   *
 *                   the update of uu to  del(uu)/uu < ropt[RELU] *
 *                   The default is no constraint on the relative *
 *                   step in uu.                                  *
 *                                                                *
 * ropt[ETAGAMMA]   (input) the coefficient egamma in the eta     *
 *                   computation. (See iopt[ETACHOICE] above.)    *
 *                                                                *
 * ropt[ETAALPHA]   (input) the coefficient ealpha in the eta     *
 *                   computation. (See iopt[ETACHOICE] above.)    *
 *                                                                *
 * ropt[ETACONST]   (input)  a user specified constant value for  *
 *                   eta, used in lieu of that computed by routine*
 *                   KINForcingTerm. (See iopt[ETACHOICE] above.) *
 *                                                                *
 * ropt[FNORM]      (output) the scaled norm at a given iteration:*
 *                   norm(fscale(func(uu)).                       *
 *                                                                *
 * ropt[STEPL]      (output) last step length in the global       *
 *                   strategy routine (KINLineSearch or           *
 *                   KINInexactNewton).                           *
 *                                                                *
 *****************************************************************/

/******************************************************************
 *                                                                *
 *         KINSOL termination codes                               *
 ******************************************************************     
 *  KINSol returns an integer-valued termination code             *
 *                                                                *
 *  The termination values KINSOL_***** are now given:            *
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
 *                                                                *
 *  INPUT_ERROR: one or more input parameters or arrays was in    *
 *               eror. See the listing in msgfp for further info  *
 *                                                                *
 *  LSOLV_NO_MEM: The linear solver memory pointer (lmem) was     *
 *             received as NULL. The return value from the linear *
 *             solver needs to be checked and the cause found.    *
 *                                                                *
 ******************************************************************/


int KINSol(void *kinmem, integertype Neq, N_Vector uu, SysFn func, 
           int globalstrategy, N_Vector uscale, N_Vector fscale,
           realtype fnormtol, realtype scsteptol, N_Vector constraints, 
           booleantype optIn, long int iopt[], realtype ropt[], void *f_data);


/* KINSol return values */
  
enum { KINSOL_NO_MEM=-1, KINSOL_INPUT_ERROR=-2, KINSOL_LSOLV_NO_MEM=-3, 
       KINSOL_SUCCESS=1, KINSOL_INITIAL_GUESS_OK=2,KINSOL_STEP_LT_STPTOL=3, 
       KINSOL_LNSRCH_NONCONV=4, KINSOL_MAXITER_REACHED=5, 
       KINSOL_MXNEWT_5X_EXCEEDED=6, KINSOL_LINESEARCH_BCFAIL=7,
       KINSOL_KRYLOV_FAILURE = 8, KINSOL_PRECONDSET_FAILURE=9, 
       KINSOL_PRECONDSOLVE_FAILURE=10};
 
 
 
/******************************************************************
 *                                                                *
 * Function : KINFree                                             *
 *----------------------------------------------------------------*
 * KINFree frees the problem memory kinsol_mem allocated by       *
 * KINMalloc.  Its only argument is the pointer kinsol_mem        *
 * returned by KINMalloc   .                                      *
 *                                                                *
 ******************************************************************/

void KINFree(void *kin_mem);
 
 
/* iopt indices */

enum { PRINTFL=0 , MXITER, PRECOND_NO_INIT, NNI ,NFE ,NBCF, NBKTRK,
       ETACHOICE, NO_MIN_EPS};

/* ropt indices */

enum { MXNEWTSTEP=0 , RELFUNC , RELU , FNORM , STEPL,
       ETACONST, ETAGAMMA, ETAALPHA };

enum {ETACHOICE1 = 0, ETACHOICE2, ETACONSTANT}; 
  /* 3 methods to determine eta
     check iopt[ETACHOICE] against these three constants
     Note that setting ETACHOICE1 to 0 implies that it is
     the default , see iopt, ropt conventions */

/******************************************************************
 *                                                                *
 * Types : struct KINMemRec, KINMem                               *
 *----------------------------------------------------------------*
 * The type KINMem is type pointer to struct KINMemRec. This      *
 * structure contains fields to keep track of problem status.     *
 *                                                                *
 ******************************************************************/


typedef struct KINMemRec {

  realtype kin_uround;    /* machine unit roundoff */

  /* Problem Specification Data */

  integertype  kin_Neq;        /*  system size                                */
  SysFn kin_func;              /* = func(uu)= 0.                              */
  void *kin_f_data;            /* ptr to func work space                      */
  realtype kin_fnormtol;       /* ptr to func norm tolerance                  */
  realtype kin_scsteptol;      /* ptr to scaled step tolerance                */
  int kin_globalstrategy;      /* choices are INEXACT_NEWTON & LINESEARCH     */
  long int *kin_iopt;          /* pointer to integer optional INPUTs and
                                  OUTPUTs                                     */
  realtype *kin_ropt;          /* pointer to real optional INPUTs and OUTPUTs */
  int kin_printfl;             /* print (output) option selected              */
  int kin_mxiter;              /* max number of nonlinear iterations          */
  int kin_msbpre;              /* max number of iterations without calling the
                                  preconditioner.  this variable set by the
                                  setup routine for the linear solver         */
  int kin_etaflag;             /* eta computation choice                      */
  booleantype kin_ioptExists;  /* logical set when the array iopt is non-null */
  booleantype kin_roptExists;  /* logical set when the array ropt is non-null */
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
  realtype kin_relu;           /* scalar bound on (del(uu)/uu)                */
  realtype kin_stepl;          /* previous steplength in LineSearch or
                                  InexactNewton algorithm                     */
  realtype kin_eps;            /* current eps value for the iteration         */
  realtype kin_eta;            /* current eta value for the iteration         */
  realtype kin_eta_gamma;      /* gamma value for use in eta calculation      */
  realtype kin_eta_alpha;      /* alpha value for use in eta calculation      */
  realtype kin_pthrsh;         /* threshold value for calling preconditioner  */

  /* Counters */

  long int  kin_nni;       /* number of nonlinear iterations              */
  long int  kin_nfe;       /* number of func references/calls             */
  long int  kin_nnilpre;   /* nni value at last precond call              */
  long int  kin_nbcf;      /* number of times the beta condition could not 
                              be met in LineSearch                        */
  long int  kin_nbktrk;    /* number of backtracks                        */
  long int  kin_ncscmx;    /* number of consecutive max stepl occurences
                              in the global strategy                      */

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


  /* Linear Solver Data */
 
  /* Linear Solver functions to be called */

  int (*kin_linit)(struct KINMemRec *kin_mem, booleantype *setupNonNull);

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


  /* Arrays for Optional Input and Optional Output */


  /* Message File */
  
  FILE *kin_msgfp;  /* where KINSol error/warning/info messages are sent */

  /* Pointer to Machine Environment-Specific Information */
  
  M_Env kin_machenv;

} *KINMem;



/******************************************************************
 *                                                                *
 * Communication between kinsol.c and a KINSol Linear Solver      *
 *----------------------------------------------------------------*
 * (1) kin_linit return values                                    *
 *                                                                *
 * LINIT_OK    : The kin_linit routine succeeded.                 *
 *                                                                *
 * LINIT_ERR   : The kin_linit routine failed. Each linear solver *
 *               init routine should print an appropriate error   *
 *               message to (kin_mem->msgfp).                     *
 *                                                                *
 *                                                                *
 ******************************************************************/

/* kin_linit return values */

#define LINIT_OK        0
#define LINIT_ERR      -1

/*******************************************************************
 *                                                                 *
 * int (*kin_linit)(KINMem kin_mem, booleantype *setupNonNull);    *
 *-----------------------------------------------------------------*
 * The purpose of kin_linit is to perform any needed               *
 * initializations of solver-specific memory, such as              *
 * counters/statistics. Actual memory allocation should be done    *
 * by the routine that initializes the linear solver package.      *
 * The kin_linit routine should set  *setupNonNull to be TRUE if   * 
 * the setup operation for the linear solver is non-empty and      *
 * FALSE if the setup operation does nothing. An LInitFn should    *
 * return LINIT_OK (= 0) if it has successfully initialized the    *
 * KINSol linear solver and LINIT_ERR (= -1) otherwise. These      *
 * constants are defined above. If an error does occur, an         *
 * appropriate message should be sent to (kin_mem->msgfp).         *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
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
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * int (*kin_lsolve)(KINMem kin_mem, N_Vector bb, N_Vector xx,     *
 *                   realtype *res_norm);                          *
 *-----------------------------------------------------------------*
 * kin_lsolve must solve the linear equation J x = b, where        *
 * J is an approximate Jacobian matrix, x is the approximate system*
 * solution,  and the RHS vector b is input. The solution is to be *
 * returned in the vector b. kin_lsolve returns a positive value   *
 * for a recoverable error and a negative value for an             *
 * unrecoverable error. Success is indicated by a 0 return value.  *
 *                                                                 *
 *******************************************************************/

/*******************************************************************
 *                                                                 *
 * void (*kin_lfree)(KINMem kin_mem);                              *
 *-----------------------------------------------------------------*
 * kin_lfree should free up any memory allocated by the linear     *
 * solver. This routine is called once a problem has been          *
 * completed and the linear solver is no longer needed.            *
 *                                                                 *
 *******************************************************************/


#endif
#ifdef __cplusplus
}
#endif
