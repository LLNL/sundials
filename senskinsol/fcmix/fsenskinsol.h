/********************************************************************
 *                                                                  *
 * File          : fsenskinsol.h                                    *
 * Programmers   : Allan G. Taylor and Alan C. Hindmarsh, and       *
 *                 Keith Grant @ LLNL                               *
 *                                                                  *
 * Version of    : 14 Nov 2000 - Changed SENSKFUN to KFUN           *
 *               : 13 Sep 2000                                      *
 *               : 15 Jan 1999 (standard KINSOL version)            *
 *------------------------------------------------------------------*
 *  This is the header file for the FSENSKINSOL Interface Package   *
 *  See below for usage details                                     *
 *                                                                  *
 *******************************************************************/

/***************************************************************************
 *
 *               FSENSKINSOL Interface Package
 *
 * The FSENSKINSOL Interface Package is a package of C functions which
 * supports the use of the SensKINSOL solver for the solution of the
 * parameter sensitivities of nonlinear systems F(u)=0, in a mixed Fortran/C
 * setting. While SensKINSOL is written in C, it is assumed here that the
 * user's calling program and user-supplied problem-defining routines are
 * written in Fortran. This package provides the necessary interface, to
 * both the serial and the parallel (MPI) versions of SensKINSOL.
 *
 * The user-callable functions, with the corresponding SensKINSOL functions,
 * are as follows:
 *
 *    FPSENSKINMALLOC and FSSENSKINMALLOC are interfaces to SensKINMalloc
 *    for parallel and serial environments, respectively.
 *
 *    FSENSKINSPGMR00, FSENSKINSPGMR10, FSENSKINSPGMR20, FSENSKINSPGMR01,
 *    FSENSKINSPGMR11, and FSENSKINSPGMR21 interface to SensKINSpgmr for
 *    the various options. NOTE: only the first of these six routines are
 *    found in fsenskinsolp.c (parallel version) or fsenskinsols.c (serial
 *    version). To allow resolving of references to be done more readily, 
 *    the last five routines are in separate files, fsenskinspgmr10.c, etc.
 *
 *    FSENSKINSOL      interfaces to SensKINSol  or
 *    FSENSKININIT     interfaces to SensKINInit
 *
 *    FSENSKINFREE     interfaces to SensKINFree
 *
 *    FSENSKINLININIT  interfaces to SensKINLinInit
 *
 *    FSENSKINDIFF     interfaces to SensKINDiff
 *
 *    FSENSKINLINSOLVE interfaces to SensKINLinSolve
 *
 * The user-supplied functions, each with the corresponding interface
 * function which calls it (and its type within SensKINSOL), are as follows:
 *
 *    KFUN is called by the interface function SensKINfunc of type SysFn
 *    Note that versions of KFUN for sensitivity analysis require five
 *    arguments compared to the three for the standard FKINSOL interface
 *
 *    FATIMES is called by the interface function KINUAtimes of type
 *    KINSpgmruserAtimesFn
 *
 *    KPSOL is called by the interface function KINPSol of type
 *    KINSpgmrSolveFn
 *
 *    KPRECO is called by the interface function KINPreco of type
 *    KINSpgmrPrecondFn
 *
 *
 * In contrast to the case of direct use of SensKINSOL, the names of all
 * user-supplied routines here are fixed, in order to maximize portability
 * for the resulting mixed-language program.
 *
 *
 * Important note on portability:
 *
 * This package uses generic names in the actual code. For simplicity, the
 * uppercase, no underscore translation is used universally in this
 * documentation [e.g. KFUN is used to denote a given Fortran-callable
 * subroutine that is actually in the code as K_FUN, the generic name].
 * The actual names are determined by parameters set in fcmixpar.h.
 *
 ***************************************************************************
 *
 *               Usage of the FSENSKINSOL Interface Package
 *
 * The usage of FSENSKINSOL requires calls to six to nine interface
 * functions, depending on the method options selected, and one to three
 * user-supplied routines which define the problem to be solved. These
 * function calls and user routines are summarized separately below.
 * Some details are omitted, and the user is referred to the SensKINSOL
 * manual for more complete documentation. Information on the arguments of
 * any given user-callable interface routine, or of a given user-supplied
 * function called by an interface function, can be found in the
 * documentation on the corresponding function in the KINSOL package.
 *
 *
 * (1) User-supplied system routine: KFUN
 *
 * The user must in all cases supply the following Fortran routine
 *
 *    SUBROUTINE KFUN (NEQ, UU, FVAL, NP, P)
 *    DIMENSION UU(*), FVAL(*), P(*)
 *
 * It must set the FVAL array to f(u), the system function, as a function
 * of the array UU = u.  Here UU and FVAL are arrays representing vectors,
 * which are distributed vectors in the parallel case, and NEQ is the
 * length of these arrays (meaning, in the parallel case, the local length
 * NLOCAL on the current processor).
 *
 *
 * (2) Optional user-supplied Jacobian-vector product routine: FATIMES
 *
 * As an option, the user may supply a routine that computes the product
 * of the system Jacobian and a given vector. This has the form
 *
 *    SUBROUTINE FATIMES(V, Z, NEWU, UU, IER)
 *    DIMENSION V(*), Z(*), UU(*)
 *
 * This must set the array Z to the product J*V, where J is the Jacobian
 * matrix J = df/du, and V is a given array.  Here UU is an array containing
 * the current value of the unknown vector u.  NEWU is an input integer
 * indicating whether UU has changed since FATIMES was last called
 * (1 = yes, 0 = no).  If FATIMES computes and saves Jacobian data, then
 * no such computation is necessary when NEWU = 0.  Here V, Z, and UU are
 * arrays of length NEQ, the problem size, or the local length of all
 * distributed vectors in the parallel case.  FATIMES should return IER = 0
 * if successful, or a nonzero IER otherwise.
 *
 *
 * (3) Initialization:  FKINITMPI and FPSENSKINMALLOC
 *
 *
 * (3.1) To initialize the use of the MPI (Message Passing Interface)
 * library by SensKINSOL, the user must make the following call
 * (parallel/MPI version only):
 *
 *    CALL FKINITMPI (NLOCAL, NGLOBAL, IER)
 *
 * The arguments are:
 *
 *    NLOCAL  = local size of vectors on this processor
 *
 *    NGLOBAL = the system size, and the global size of vectors (the sum
 *              of all values of NLOCAL)
 *
 *    IER     = return completion flag. Values are
 *                0 = success,
 *               -1 = failure.
 *
 * The routine FKINITMPI, part of the standard KINSOL FORTRAN interface,
 * is found in the file fkinsolp.c rather than in fsenskinsolp.c
 *
 * (3.2) To set various problem and solution parameters and allocate
 * internal memory, make the following call:
 *
 *    CALL FPSENSKINMALLOC(NEQ, NP, IDIFFRHS, DIFFSCL, P, PBAR, IER)
 *
 * for the parallel case or
 *
 *    CALL FSSENSKINMALLOC(NEQ, NP, IDIFFRHS, DIFFSCL, P, PBAR, IER)
 *
 * for the serial case.
 *
 *
 * The arguments are:
 *
 *    NEQ      = The (global) problem size
 *
 *    NP       = The number of sensitivity parameters (np >= 0)
 *
 *    IDIFFRHS = Integer flag determining the form of the RHS 
 *               difference equation
 *                  = 0 : Forward difference
 *                  = 1 : Backward difference
 *                  = 2 : Centered difference
 *
 *    DIFFSCL  = Positive scale factor for the RHS difference increment.
 *               The nominal value for diffscl is 1.0
 *
 *    P        = Array of sensitivity parameters of length np
 *
 *    PBAR     = Array of sensitivity parameter nominal magnitudes > 0
 *
 *    IER      = return completion flag. Values are
 *                  =  0 : success,
 *                  = -1 : failure.
 *
 * See printed message for details in case of failure.
 *
 *
 * (4) Specification of linear system solution method.
 *
 * The solution method in SensKINSOL involves the solution of linear systems
 * related to the Jacobian J = df/du of the nonlinear system.  SensKINSOL
 * presently includes only one choice for the treatment of these systems,
 * and the user must call a routine with a specific name to make this
 * choice.
 *
 *
 * (4.1) SPGMR treatment of the linear systems.
 *
 * For the Scaled Preconditioned GMRES solution of the linear systems,
 * the user must make one of the following six setup calls:
 *
 *    CALL FSENSKINSPGMR00(MAXL, MAXLRST, MSBPRE)
 *
 * if no preconditioning or user-supplied ATimes routine is to be used;
 *
 *
 * There are five additional routines that can be called instead of
 * FSENSKINSPGMR00 for the cases where preconditioning or user-supplied
 * atimes routines are to be specified. Their arguments and call are
 * documented here, but they are to be found in files like fsenskinspgmr01.c
 *
 *    CALL FSENSKINSPGMR01(MAXL, MAXLRST, MSBPRE)  if no preconditioning
 *    is used but a user-supplied atimes routine is called.
 *
 *    CALL FSENSKINSPGMR10(MAXL, MAXLRST, MSBPRE)   if preconditioning is
 *    used but no preconditioning setup routine. also, no atimes routine.
 *
 *    CALL FSENSKINSPGMR11(MAXL, MAXLRST, MSBPRE)  if preconditioning is used
 *    but no preconditioning setup routine. A user-supplied atimes routine IS
 *    supplied (FATIMES).
 *
 *    CALL FSENSKINSPGMR20(MAXL, MAXLRST, MSBPRE)  if preconditioning is used
 *    and a preconditioning setup routine IS used. No user-supplied atimes
 *    routine.
 *
 *    CALL FSENSKINSPGMR21(MAXL, MAXLRST, MSBPRE)  if preconditioning is used
 *    and a precondiioning setup routine IS used. A user-supplied atimes
 *    routine is supplied as well.
 *
 *    The arguments are:
 *
 *       MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 *
 *       MAXLRST  = maximum number of linear system restarts; 0 indicates
 *                  default.
 *
 *       MSBPRE   = maximum number of preconditioning solve calls without
 *                  calling the preconditioning setup routine; 0 indicates
 *                  default; applies only for FSENSKINSPGMR20 and
 *                  FSENSKINSPGMR21.
 *
 *
 * In the four cases FSENSKINSPGMR10, FSENSKINSPGMR11, FSENSKINSPGMR20,
 * and FSENSKINSPGMR21, the user program must include the following routine
 * for solution of the preconditioner linear system:
 *
 *    SUBROUTINE KPSOL (NEQ, UU, USCALE, FVAL, FSCALE, VTEM, FTEM,
 *                      UROUND, NFE, IER)
 *    DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)
 *
 *
 * Typically this routine will use only NEQ, UU, FVAL, VTEM and FTEM.
 * It must solve the preconditioner linear system Pz = r, where r = VTEM is
 * input, and store the solution z in VTEM as well.  Here P is the right
 * preconditioner. If scaling is being used, the routine supplied must also
 * account for scaling on either coordinate or function value. NEQ is the
 * (global) problem size.
 *
 *
 * In the two cases FSENSKINSPGMR20 and FSENSKINSPGMR21, the user program
 * must also include the following routine for the evaluation and
 * preprocessing of the preconditioner:
 *
 *    SUBROUTINE KPRECO (NEQ, UU, USCALE, FVAL, FSCALE, VTEMP1, VTEMP2,
 *                       UROUND, NFE, IER)
 *    DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEMP1(*), VTEMP2(*)
 *
 *
 * It must perform any evaluation of Jacobian-related data and preprocessing
 * needed for the solution of the preconditioner linear systems by KPSOL.
 * The variables UU through FSCALE are for use in the preconditioning setup
 * process. Typically, the system function KFUN is called, so that FVAL will
 * have been updated. UU is the current solution iterate. VTEMP1 and VTEMP2
 * are available for work space. If scaling is being used, USCALE and FSCALE
 * are available for those operatins requiring scaling.  NEQ is the (global)
 * problem size.
 *
 * On return, set IER = 0 if KPRECO was successful, set IER positive if an
 * error occurred.
 *
 *
 * (5) The solver: FSENSKINSOL (or FSENSKININIT)
 *
 * Carrying out the solving of the nonlinear system is accomplished by
 * making calls as follows:
 *
 *    CALL FSENSKINSOL (NEQ, UU, GLOBALSTRAT, USCALE, FSCALE, FNORMTOL,
 *       SCSTEPTOL, CONSTRAINTS, OPTIN, IOPT, ROPT, IER)
 *
 *
 * The arguments are:
 *
 *    NEQ         = (INTEGER) number of equations (unknowns) in the
 *                  nonlinear system
 *
 *    UU          = Array containing the initial guess when called,
 *                  returns the solution
 *
 *    GLOBALSTRAT = (INTEGER)a number defining the global strategy choice:
 *                      0 = INEXACTNEWTON,
 *                      1 = LINESEARCH .
 *
 *    USCALE      = Array of scaling factors for the UU vector
 *
 *    FSCALE      = Array of scaling factors for the FVAL (function) vector
 *
 *    FNORMTOL    = Tolerance on the norm of f(u) to accept convergence.
 *
 *    SCSTEPTOL   = Tolerance on minimum scaled step size
 *
 *    CONSTRAINTS = Array of constraint values, by element of the
 *                  solution UU
 *
 *    OPTIN       = Integer used as a flag to indicate whether possible input
 *                  values in IOPT are to be used for input 0 = NO, 1 = YES.
 *
 *    IOPT        = INTEGER array of input and output values
 *
 *    ROPT        = Array of real input and output values
 *
 *    IER         = Integer error flag as returned by SensKINSOL.
 *                  1 or 2 =  success; other = error condition. See  
 *                  SensKINSOL documentation for full details.
 *
 * If solution vector UU of the nonlinear system F(UU) = 0 is already 
 * known, the user has the option of calling FSENSKININIT instead of 
 * FSENSKINSOL. FSENSKININIT has the same calling sequence as FSENSKININIT, 
 * but simply checks that the "initial guess" is an accurate solution and 
 * does the KINSOL framework initialization required for the subsequent 
 * linear sensitivity solutions.
 *
 *
 * (6) Sensitivity Solutions: FSENSKINLININIT, FSENSKINDIFF,
 *     and FSENSKINLINSOLVE
 *
 *
 * (6.1) Following obtaining the solution UU of the nonlinear system
 * F(UU) = 0 by calling  FSENSKINSOL or verifying that UU is a solution 
 * by calling FSENSKININIT, SensKINSol can calculate the sensitivity of
 * UU to the parameters "P" passed to FPSENSKINMALLOC or FSSENSKINMALLOC.
 * A sensitivity solver initialization is first done by
 *
 *    CALL FSENSKINLININIT (IER)
 *
 *    IER     = Integer error flag as returned by SensKINLinInit.
 *              0 = success; nonzero = error condition. See 
 *              SensKINSOL documentation or sens_kinsol.h for 
 *              full details.
 *
 *
 * (6.2) Before each sensitivity solution is done, the right hand side 
 * difference method and/or the RHS difference increment scaling can 
 * optionally be reset by 
 *
 *    CALL FSENSKINDIFF (IDIFFRHS, DIFFSCL, IER)
 *
 * with
 *
 *    IDIFFRHS = Integer flag determining the form of the RHS 
 *               difference equation
 *                 = 0 : Forward difference
 *                 = 1 : Backward difference
 *                 = 2 : Centered difference
 *
 *    DIFFSCL = Positive scale factor for the RHS difference increment.
 *              The nominal value for diffscl is 1.0
 *
 *    IER     = Integer error flag as returned by SensKINDIFF. 
 *              0 = success; nonzero = error condition. See
 *              SensKINSOL documentation or sens_kinsol.h for
 *              full details.
 *
 *
 * (6.3) For each parameter "P", the sensitivity is calculated by
 *
 *    CALL FSENSKINLINSOLVE (IP, WW, IER)
 *
 * with
 *
 *    IP      = Index of the current parameter P(IP) (count from 1).
 *
 *    WW      = Array of sensitivites of the solution UU of
 *              the system F(UU,P) = 0 to variation of P(IP).
 *              WW is returned by FSENSKINLINSOLVE.
 *
 *    IER     = Integer error flag as returned by SensKINLinSolve.
 *              0 = success; nonzero = error condition. See
 *              SensKINSOL documentation or sens_kinsol.h for
 *              full details.
 *
 *
 * (7) Memory freeing: FSENSKINFREE and FKFREEMPI
 *
 * To the free the internal memory created by the calls to FKINITMPI
 * (parallel/MPI version only) and FSENSKINMALLOC, make the following calls,
 * in this order:
 *
 *    CALL FSENSKINFREE
 *    CALL FKFREEMPI  (parallel/MPI version only)
 *
 * The routine FKFREEMPI, part of the standard KINSOL FORTRAN interface, is
 * found in the file fkinsolp.c rather than in fsenskinsolp.c
 *
 * * * - - - * * * - - - * * * - - - * * * - - - * * * - - - * * * - - - * * *
 *
 * Informational output:
 *
 * Some of the optional outputs are NFE, NNI, NLI, NPS, NCFL, NPE, stored
 * in the array iopt at NFE+1, NNI+1, SPGMR_NLI+1, SPGMR_NPS+1, SPGMR_NPE+1,
 * SPGMR_NPS+1, SPGMR_NCFL+1, respectively. (See the KINSOL and KINSPGMR
 * header files for descriptions and information on other outputs. The +1
 * on each index is due to the differences in labeling arrays between C
 * and Fortran.)
 *
 ***************************************************************************/

#ifndef _fsenskinsol_h
#define _fsenskinsol_h

#include "fcmixpar.h"


/* generic names are translated through the define statements below for a
  specific platform/compiler */

#if (CRAY)

#define  F_PSENSKINMALLOC  FPSENSKINMALLOC
#define  F_SSENSKINMALLOC  FSSENSKINMALLOC
#define  F_SENSKINSPGMR00  FSENSKINSPGMR00
#define  F_SENSKINSPGMR01  FSENSKINSPGMR01
#define  F_SENSKINSPGMR10  FSENSKINSPGMR10
#define  F_SENSKINSPGMR11  FSENSKINSPGMR11
#define  F_SENSKINSPGMR20  FSENSKINSPGMR20
#define  F_SENSKINSPGMR21  FSENSKINSPGMR21
#define  F_SENSKINSOL      FSENSKINSOL
#define  F_SENSKININIT     FSENSKININIT
#define  F_SENSKINLININIT  FSENSKINLININIT
#define  F_SENSKINDIFF     FSENSKINDIFF
#define  F_SENSKINLINSOLVE FSENSKINLINSOLVE
#define  F_SENSKINFREE     FSENSKINFREE
#define  K_FUN             KFUN

#elif  (UNDERSCORE)

#define  F_PSENSKINMALLOC  fpsenskinmalloc_
#define  F_SSENSKINMALLOC  fssenskinmalloc_
#define  F_SENSKINSPGMR00  fsenskinspgmr00_
#define  F_SENSKINSPGMR01  fsenskinspgmr01_
#define  F_SENSKINSPGMR10  fsenskinspgmr10_
#define  F_SENSKINSPGMR11  fsenskinspgmr11_
#define  F_SENSKINSPGMR20  fsenskinspgmr20_
#define  F_SENSKINSPGMR21  fsenskinspgmr21_
#define  F_SENSKINSOL      fsenskinsol_
#define  F_SENSKININIT     fsenskininit_
#define  F_SENSKINLININIT  fsenskinlininit_
#define  F_SENSKINDIFF     fsenskindiff_
#define  F_SENSKINLINSOLVE fsenskinlinsolve_
#define  F_SENSKINFREE     fsenskinfree_
#define  K_FUN             kfun_

#else

#define  F_PSENSKINMALLOC  fpsenskinmalloc
#define  F_SSENSKINMALLOC  fssenskinmalloc
#define  F_SENSKINSPGMR00  fsenskinspgmr00
#define  F_SENSKINSPGMR01  fsenskinspgmr01
#define  F_SENSKINSPGMR10  fsenskinspmgr10
#define  F_SENSKINSPGMR11  fsenskinspgmr11
#define  F_SENSKINSPGMR20  fsenskinspgmr20
#define  F_SENSKINSPGMR21  fsenskinspgmr21
#define  F_SENSKINSOL      fsenskinsol
#define  F_SENSKININIT     fsenskininit
#define  F_SENSKINLININIT  fsenskinlininit
#define  F_SENSKINDIFF     fsenskindiff
#define  F_SENSKINLINSOLVE fsenskinlinsolve
#define  F_SENSKINFREE     fsenskinfree
#define  K_FUN             kfun

#endif

/* KINSOL header files  */

#include "llnltyps.h"  /* definitions of types real and integer     */
#include "kinsol.h"    /* definition of type SysFn                  */
#include "nvector.h"   /* definition of type N_Vector, machEnvType  */


/* Prototypes: Functions called by the solver */

void SensKINfunc(integer Neq, N_Vector uu, N_Vector fval, void *f_data);
void K_FUN (integer *nlocal, real *udata, real *fdata, integer *np, real *p);


/* Declarations for global variables, shared among various routines */

void *SENSKIN_smem;           /* Sensitivity memory structure */

#endif
