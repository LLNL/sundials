 /******************************************************************
 *                                                                 *
 * File          : fkinsol.h                                       *
 * Programmers   : Allan G. Taylor, Alan C. Hindmarsh, and         * 
 *                 Radu Serban @ LLNL                              *
 * Version of    : 27 January 2004                                 *
 *-----------------------------------------------------------------*
 *  This is the header file for the FKINSOL Interface Package      *
 *  See below for usage details                                    *
 *                                                                 *
 *******************************************************************/

/***************************************************************************

                  FKINSOL Interface Package

 The FKINSOL Interface Package is a package of C functions which support the
 use of the KINSOL solver for the solution of nonlinear systems f(u) = 0, 
 in a mixed Fortran/C setting.  While KINSOL is written in C, it is assumed 
 here that the user's calling program and user-supplied problem-defining 
 routines are written in Fortran. This package provides the necessary interface
 to KINSOL for both the serial and the parallel NVECTOR implementations.

 The user-callable functions, with the corresponding KINSOL functions,
 are as follows:
  FNVINITS and FNVINITP interface to NV_SpecInit_Serial and
       NV_SpecInitParallel (defined by nvector_serial
       and nvector_parallel, respectively)
  FKINMALLOC  interfaces to KINMalloc 
  FKINSPGMR   interfaces to KINSpgmr
  FKINSOL   interfaces to KINSol
  FKINFREE  interfaces to KINFree
  FNVFREES and FNVFREEP interface to NV_SpecFree_Serial and
       NV_SpecFree_parallel(defined by nvector_serial
       and nvector_parallel, respectively)

 The user-supplied functions, each with the corresponding interface function
 which calls it (and its type within KINSOL), are as follows:
  FKFUN   is called by the interface function FKINfunc of type SysFn
  FJTIMES is called by the interface function FKINJtimes of type
          KINSpgmrJacTimesVecFn
  FKPSOL  is called by the interface function FKINPSol of type KINSpgmrPrecSolveFn
  FKPSET  is called by the interface function FKINPSet of type KINSpgmrPrecSetupFn
 In contrast to the case of direct use of KINSOL, the names of all 
 user-supplied routines here are fixed, in order to maximize portability for
 the resulting mixed-language program.

 Important note on portability:
 This package uses generic names in the actual code. For simplicity, the
 uppercase, no underscore translation is used universally in this documentation
 [e.g. FKFUN is used to denote a given Fortran-callable subroutine that is 
 actually in the code as FK_FUN, the generic name]. The actual names are
 determined by parameters set at the configuration stage

 =========================================================================

                  Usage of the FKINSOL Interface Package

 The usage of FKINSOL requires calls to six interface functions, and one to 
 three user-supplied routines which define the problem to be solved.  These
 function calls and user routines are summarized separately below.

 Some details are omitted, and the user is referred to the KINSOL manual
 for more complete documentation.  Information on the arguments of any
 given user-callable interface routine, or of a given user-supplied
 function called by an interface function, can be found in the
 documentation on the corresponding function in the KINSOL package.

 The number labels on the instructions below end with s for instructions
 that apply to the serial version of CVODE only, and end with p for
 those that apply to the parallel version only.


 (1) User-supplied system routine: FKFUN
 The user must in all cases supply the following Fortran routine
       SUBROUTINE FKFUN (UU, FVAL)
       DIMENSION UU(*), FVAL(*)
 It must set the FVAL array to f(u), the system function, as a function 
 of the array UU = u.  Here UU and FVAL are arrays representing vectors,
 which are distributed vectors in the parallel case, and NEQ is the
 problem size.

 (2) Optional user-supplied Jacobian-vector product routine: FKJTIMES
 As an option, the user may supply a routine that computes the product
 of the system Jacobian and a given vector.  This has the form
      SUBROUTINE FKJTIMES(V, Z, NEWU, UU, IER)
      DIMENSION V(*), Z(*), UU(*)
 This must set the array Z to the product J*V, where J is the Jacobian
 matrix J = df/du, and V is a given array.  Here UU is an array containing
 the current value of the unknown vector u.  NEWU is an input integer 
 indicating whether UU has changed since FKJTIMES was last called 
 (1 = yes, 0 = no).  If FKJTIMES computes and saves Jacobian data, then 
 no such computation is necessary when NEWU = 0.  Here V, Z, and UU are 
 arrays of length NEQ, the problem size, or the local length of all 
 distributed vectors in the parallel case.  FKJTIMES should return IER = 0 
 if successful, or a nonzero IER otherwise.

 (3) Initialization:  FNVINITS / FNVINITP and FKINMALLOC

 (3.1s) To initialize the serial machine environment, the user must make
 the following call:
       CALL FNVINITS (NEQ, IER)
 The arguments are:
 NEQ     = size of vectors
 IER     = return completion flag. Values are 0 = success, -1 = failure.

 (3.1p) To initialize the parallel machine environment, the user must make 
 the following call:
       CALL FNVINITP (NLOCAL, NGLOBAL, IER)
 The arguments are:
 NLOCAL  = local size of vectors on this processor
 NGLOBAL = the system size, and the global size of vectors (the sum 
           of all values of NLOCAL)
 IER     = return completion flag. Values are 0 = success, -1 = failure.

 (3.2) To set various problem and solution parameters and allocate
 internal memory, make the following call:
       CALL FKINMALLOC(MSBPRE, FNORMTOL, SCSTEPTOL, CONSTRAINTS,
                       OPTIN, IOPT, ROPT, IER)
 The arguments are:
 MSBPRE      = maximum number of preconditioning solve calls without calling
               the preconditioning setup routine; 0 indicates default.
 FNORMTOL    = tolerance on the norm of f(u) to accept convergence
 SCSTEPTOL   = tolerance on minimum scaled step size
 CONSTRAINTS = array of constraint values, on components of the solution UU
 INOPT       = integer used as a flag to indicate whether possible input values
               in IOPT are to be used for input: 0 = NO, 1 = YES.
 IOPT        = array for integer optional inputs and outputs
               (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 ROPT        = array of real optional inputs and outputs
 IER         = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
               See printed message for details in case of failure.

 (4) Specification of linear system solution method.
 The solution method in KINSOL involves the solution of linear systems 
 related to the Jacobian J = df/du of the nonlinear system.  KINSOL presently 
 includes only one choice for the treatment of these systems, and the user
 must call a routine with a specific name to make this choice.

 (4.1) SPGMR treatment of the linear systems.
 For the Scaled Preconditioned GMRES solution of the linear systems,
 the user must make the call:
      CALL FKINSPGMR(MAXL, MAXLRST, IER)

 In the above routine, the arguments are as follows:
 MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 MAXLRST  = maximum number of linear system restarts; 0 indicates default.
 IER      = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
            See printed message for details in case of failure.

 If the user program includes the KJTIMES routine for the evaluation of the 
 Jacobian vector product, the following call must be made
      CALL FKINSPGMRSETJAC(FLAG, IER)
 The argument FLAG=0 specifies using the internal finite differences
 approximation to the Jacobian vector product, while FLAG=1 specifies that 
 FKJTIMES is provided.

 Usage of the user-supplied routine KPSOL for solution of the preconditioner 
 linear system is specified by calling
      CALL FKINSPGMRSETPSOL(FLAG, IER)
 where FLAG=0 indicates no KPSOL (default) and FLAG=1 specifies using FKPSOL.
 The user-supplied routine KPSOL must be of the form:

      SUBROUTINE FKPSOL (UU, USCALE, FVAL, FSCALE, VTEM, FTEM, IER)
      DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)

 Typically this routine will use only UU, FVAL, VTEM and FTEM.
 It must solve the preconditioner linear system Pz = r, where r = VTEM is 
 input, and store the solution z in VTEM as well.  Here P is the right 
 preconditioner. If scaling is being used, the routine supplied must also 
 account for scaling on either coordinate or function value.

 Usage of the user-supplied routine FKSET for construction of the preconditioner 
 is specified by calling
      CALL FKINSPGMRSETPPSET(FLAG, IER)
 where FLAG=0 indicates no FKPSET (default) and FLAG=1 specifies using FKPSET.
 The user-supplied routine FKPSET must be of the form:

      SUBROUTINE FKPSET (UU, USCALE, FVAL, FSCALE, VTEMP1, VTEMP2, IER)
      DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEMP1(*), VTEMP2(*)

 It must perform any evaluation of Jacobian-related data and preprocessing 
 needed for the solution of the preconditioner linear systems by FKPSOL.  The 
 variables UU through FSCALE are for use in the preconditioning setup process.
 Typically, the system function KFUN is called, so that FVAL will have been
 updated. UU is the current solution iterate. VTEMP1 and VTEMP2 are available
 for work space. If scaling is being used, USCALE and FSCALE are available for
 those operatins requiring scaling.  NEQ is the (global) problem size.

 On return, set IER = 0 if FKPSET was successful, set IER positive if an error 
 occurred.

 (5) The solver: FKINSOL
 Solving the nonlinear system is accomplished by making the following call:
      CALL FKINSOL (UU, GLOBALSTRAT, USCALE, FSCALE, IER)
 The arguments are:
 UU          = array containing the initial guess on input, and the solution
               on return
 GLOBALSTRAT = (INTEGER) a number defining the global strategy choice:
               0 = InexactNewton, 1 = LineSearch
 USCALE      = array of scaling factors for the UU vector
 FSCALE      = array of scaling factors for the FVAL (function) vector
 IER         = INTEGER error flag as returned by KINSOL: 1 means success,
               2 means initial guess satisfies f(u) = 0 (approx.),
               3 means apparent stalling (small step),
               a value < 0 or > 3 means other error or failure.
               See KINSOL documentation for detailed information.

 (6) Memory freeing: FKINFREE and FNVFREES / FNVFREEP
 To the free the internal memory created by the calls to FKINMALLOC and
 either FNVINITS or FNVINITP, make the following calls, in this order:
      CALL FKINFREE
      CALL FNVFREES (serial case) or CALL FNVFREEP (parallel case)

 (7) Optional inputs and outputs: IOPT/ROPT

 The optional inputs available by way of IOPT and ROPT have the following
 names, locations, and descriptions.  For further details see the KINSOL
 documentation.  A zero value results in the default.

 PRINTFL =         IOPT(1) = optional output print flag
 MXITER  =         IOPT(2) = maximum Newton iterations
 PRECOND_NO_INIT = IOPT(3) = flag to suppress initial preconditioner setup call
 ETACHOICE =       IOPT(8) = choice of forcing term (1 = Choice 1, 2 = Choice 2,
                             3 = constant)
 NO_MIN_EPS =      IOPT(9) = flag to suppress minimum tolerance (eps)
 MXNEWTSTEP = ROPT(1) = max size of Newton step
 RELFUNC =    ROPT(2) = relative error in computing f(u)
 RELU =       ROTP(3) = control on relative change in components of u per step
 ETACONST = ROPT(6), ETAGAMMA = ROPT(7), ETAALPHA = ROPT(8): constants in
            optional choices of forcing terms.

 The optional outputs available by way of IOPT and ROPT have the following
 names, locations, and descriptions.  For further details see the KINSOL
 documentation.
 
 NNI  =   IOPT(4) = number of Newton iterations
 NFE  =   IOPT(5) = number of f evaluations
 NBCF =   IOPT(6) = number of Linesearch beta condition failures
 NBKTRK = IOPT(7) = number of Linesearch backtracks
 FNORM =  ROPT(4) = final scaled norm of f(u)
 STEPL =  ROPT(5) = scaled last step length

The following optional outputs are specific to the SPGMR module:
 NLI  = IOPT(11) = number of linear (Krylov) iterations
 NPE  = IOPT(12) = number of preconditioner evaluations
 NPS  = IOPT(13) = number of preconditioner solves
 NCFL = IOPT(14) = number of linear convergence failures

*******************************************************************************/

#ifndef _fkinsol_h
#define _fkinsol_h

#include "nvector.h"

/* Generic names are translated through the define statements below for a
  specific platform/compiler */

#if defined(SUNDIALS_UNDERSCORE_NONE)

#define FKIN_MALLOC        fkinmalloc
#define FKIN_SPGMR         fkinspgmr
#define FKIN_SPGMRSETJAC   fkinspgmrsetjac
#define FKIN_SPGMRSETPSOL  fkinspgmrsetpsol
#define FKIN_SPGMRSETPSET  fkinspgmrsetpset
#define FKIN_SOL           fkinsol
#define FKIN_FREE          fkinfree
#define FK_FUN             fkfun
#define FK_PSET            fkpreco
#define FK_PSOL            fkpsol
#define FK_JTIMES          fkjtimes

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#define FKIN_MALLOC        fkinmalloc__
#define FKIN_SPGMR         fkinspgmr__
#define FKIN_SPGMRSETJAC   fkinspgmrsetjac__
#define FKIN_SPGMRSETPSOL  fkinspgmrsetpsol__
#define FKIN_SPGMRSETPSET  fkinspgmrsetpset__
#define FKIN_SOL           fkinsol__
#define FKIN_FREE          fkinfree__
#define FK_FUN             fkfun__
#define FK_PSET            fkpset__
#define FK_PSOL            fkpsol__
#define FK_JTIMES          fkjtimes__

#else

#define FKIN_MALLOC        fkinmalloc_
#define FKIN_SPGMR         fkinspgmr_
#define FKIN_SPGMRSETJAC   fkinspgmrsetjac_
#define FKIN_SPGMRSETPSOL  fkinspgmrsetpsol_
#define FKIN_SPGMRSETPSET  fkinspgmrsetpset_
#define FKIN_SOL           fkinsol_
#define FKIN_FREE          fkinfree_
#define FK_FUN             fkfun_
#define FK_PSET            fkpset_
#define FK_PSOL            fkpsol_
#define FK_JTIMES          fkjtimes_

#endif


/* KINSOL header files  */

#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "kinsol.h"        /* definition of type SysFn                      */
#include "nvector.h"       /* definition of type N_Vector, machEnvType      */


/* Prototypes: Functions called by the solver */

void FKINfunc(N_Vector uu, N_Vector fval, void *f_data);

int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *P_data,
             N_Vector vtemp1, N_Vector vtemp2);

int FKINPSol(N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale, 
             N_Vector vv, void *P_data,
             N_Vector vtemp);

int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *J_data);


/* Declarations for global variables, shared among various routines */

extern NV_Spec F2C_nvspec;

void *KIN_mem;
int *KIN_iopt;
realtype *KIN_ropt;

#endif
