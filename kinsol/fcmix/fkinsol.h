 /******************************************************************
 *                                                                  *
 * File          : fkinsol.h                                        *
 * Programmers   : Allan G. Taylor and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 21 December 2001                                 *
 *------------------------------------------------------------------*
 *  This is the header file for the FKINSOL Interface Package       *
 *  See below for usage details                                     *
 *                                                                  *
 *******************************************************************/

/***************************************************************************

                  FKINSOL Interface Package

The FKINSOL Interface Package is a package of C functions which support
the use of the KINSOL solver for the solution of nonlinear systems f(u)=0, 
in a mixed Fortran/C setting.  While KINSOL is written in C, it is assumed 
here that the user's calling program and user-supplied problem-defining 
routines are written in Fortran. This package provides the necessary interface,
to both the serial and the parallel (MPI) versions of KINSOL.

The user-callable functions, with the corresponding KINSOL functions,
are as follows:
  FKINITMPI interfaces to PVecInitMPI  (FKINITMPI is located in fkinsolp.c)
  FPKINMALLOC and FSKINMALLOC  interface to KINMalloc 
  FKINSPGMR00, FKINSPGMR10, FKINSPGMR20, FKINSPGMR01, FKINSPGMR11, and
       FKINSPGMR21 interface to KINSpgmr for the various options
       NOTE: To allow resolving of references to be done more readily, 
       these are in separate files.
  FKINSOL   interfaces to KINSol
  FKINFREE  interfaces to KINFree
  FKFREEMPI interfaces to PVecFreeMPI (FKFREEMPI is located in fkinsolp.c)

FKINITMPI and FKFREEMPI are used only in the Fortran interface with the 
parallel (MPI) version of KINSOL.

The user-supplied functions, each with the corresponding interface function
which calls it (and its type within KINSOL), are as follows:
  KFUN   is called by the interface function KINfunc of type SysFn
  FATIMES is called by the interface function KINUAtimes of type
          KINSpgmruserAtimesFn
  KPSOL  is called by the interface function KINPSol of type KINSpgmrSolveFn
  KPRECO is called by the interface function KINPreco of type KINSpgmrPrecondFn
In contrast to the case of direct use of KINSOL, the names of all 
user-supplied routines here are fixed, in order to maximize portability for
the resulting mixed-language program.

Important note on portability:
This package uses generic names in the actual code. For simplicity, the
uppercase, no underscore translation is used universally in this documentation
[e.g. KFUN is used to denote a given Fortran-callable subroutine that is 
actually in the code as K_FUN, the generic name]. The actual names are
determined by parameters set in fcmixpar.h.

****************************************************************************

                 Usage of the FKINSOL Interface Package

The usage of FKINSOL requires calls to four to six interface functions, 
depending on the method options selected, and one to three
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.

Some details are omitted, and the user is referred to the KINSOL manual
for more complete documentation.  Information on the arguments of any
given user-callable interface routine, or of a given user-supplied
function called by an interface function, can be found in the
documentation on the corresponding function in the KINSOL package.


(1) User-supplied system routine: KFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE KFUN (NEQ, UU, FVAL)
      DIMENSION UU(*), FVAL(*)
It must set the FVAL array to f(u), the system function, as a function 
of the array UU = u.  Here UU and FVAL are arrays representing vectors,
which are distributed vectors in the parallel case, and NEQ is the
length of these arrays (meaning, in the parallel case, the local length 
NLOCAL on the current processor).

(2) Optional user-supplied Jacobian-vector product routine: FATIMES
As an option, the user may supply a routine that computes the product
of the system Jacobian and a given vector.  This has the form
      SUBROUTINE FATIMES(V, Z, NEWU, UU, IER)
      DIMENSION V(*), Z(*), UU(*)
This must set the array Z to the product J*V, where J is the Jacobian
matrix J = df/du, and V is a given array.  Here UU is an array containing
the current value of the unknown vector u.  NEWU is an input integer 
indicating whether UU has changed since FATIMES was last called 
(1 = yes, 0 = no).  If FATIMES computes and saves Jacobian data, then 
no such computation is necessary when NEWU = 0.  Here V, Z, and UU are 
arrays of length NEQ, the problem size, or the local length of all 
distributed vectors in the parallel case.  FATIMES should return IER = 0 
if successful, or a nonzero IER otherwise.

(3) Initialization:  FKINITMPI and FPKINMALLOC

(3.1) To initialize the use of the MPI (Message Passing Interface) library
by KINSOL, the user must make the following call (parallel/MPI version only):
      CALL FKINITMPI (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors on this processor
NGLOBAL = the system size, and the global size of vectors (the sum 
          of all values of NLOCAL)
IER     = return completion flag. Values are 0 = success, -1 = failure.

(3.2) To set various problem and solution parameters and allocate
internal memory, make the following call:
      CALL FPKINMALLOC(NEQ, IER)  or  CALL FSKINMALLOC(NEQ, IER)
in the parallel or serial case, respectively.  The arguments are:
NEQ    = the (global) problem size
IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
         See printed message for details in case of failure.

(4) Specification of linear system solution method.
The solution method in KINSOL involves the solution of linear systems 
related to the Jacobian J = df/du of the nonlinear system.  KINSOL presently 
includes only one choice for the treatment of these systems, and the user
must call a routine with a specific name to make this choice.

(4.1) SPGMR treatment of the linear systems.
For the Scaled Preconditioned GMRES solution of the linear systems,
the user must make one of the following six setup calls:
      CALL FKINSPGMR00(MAXL, MAXLRST, MSBPRE)
       if no preconditioning or user-supplied ATimes routine is to be used;

There are five additional routines that can be called instead of FKINSPGMR00 
for the cases where preconditioning or user-supplied atimes routines are to be
specified. Their arguments and call are documented here, but they are to be
found in files like fkinspgmr01.c

      CALL FKINSPGMR01(MAXL, MAXLRST, MSBPRE)  if no preconditioning
          but a user-supplied atimes routine
                                                       
      CALL FKINSPGMR10(MAXL, MAXLRST, MSBPRE)   if preconditioning is
         used but no preconditioning setup routine. also, no atimes routine.

      CALL FKINSPGMR11(MAXL, MAXLRST, MSBPRE)  if preconditioning is used but
         no preconditioning setup routine. A user-supplied atimes routine IS
	 supplied (FATIMES).

      CALL FKINSPGMR20(MAXL, MAXLRST, MSBPRE)  if preconditioning is used and
         a preconditioning setup routine IS used. No user-supplied atimes
	 routine.

      CALL FKINSPGMR21(MAXL, MAXLRST, MSBPRE)  if preconditioning is use dand
         a precondiioning setup routine IS used. A user-supplied atimes routine
	 is supplied as well.
                                                      
The arguments are:
MAXL     = maximum Krylov subspace dimension; 0 indicates default.
MAXLRST  = maximum number of linear system restarts; 0 indicates default.
MSBPRE   = maximum number of preconditioning solve calls without calling the
            preconditioning setup routine; 0 indicates default; applies only
	    for FKINSPGMR20 and FKINSPGMR21.

In the four cases FKINSPGMR10, FKINSPGMR11, FKINSPGMR20, and FKINSPGMR21,
the user program must include the following routine for solution of the 
preconditioner linear system:

      SUBROUTINE KPSOL (NEQ, UU, USCALE, FVAL, FSCALE, VTEM, FTEM, 
                        UROUND, NFE, IER)
      DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)

Typically this routine will use only NEQ, UU, FVAL, VTEM and FTEM.
It must solve the preconditioner linear system Pz = r, where r = VTEM is 
input, and store the solution z in VTEM as well.  Here P is the right 
preconditioner. If scaling is being used, the routine supplied must also 
account for scaling on either coordinate or function value.
NEQ is the (global) problem size.

In the two cases FKINSPGMR20 and FKINSPGMR21, the user program must also 
include the following routine for the evaluation and preprocessing of the 
preconditioner:

      SUBROUTINE KPRECO (NEQ, UU, USCALE, FVAL, FSCALE, VTEMP1, VTEMP2,
                         UROUND, NFE, IER)
      DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEMP1(*), VTEMP2(*)

It must perform any evaluation of Jacobian-related data and preprocessing 
needed for the solution of the preconditioner linear systems by KPSOL.  The 
variables UU through FSCALE are for use in the preconditioning setup process.
Typically, the system function KFUN is called, so that FVAL will have been
updated. UU is the current solution iterate. VTEMP1 and VTEMP2 are available
for work space. If scaling is being used, USCALE and FSCALE are available for
those operatins requiring scaling.  NEQ is the (global) problem size.

On return, set IER = 0 if KPRECO was successful, set IER positive if an error 
occurred.


(5) The solver: FKINSOL
Carrying out the solving of the nonlinear system is accomplished by making 
calls as follows:
      CALL FKINSOL (NEQ, UU, GLOBALSTRAT, USCALE, FSCALE, FNORMTOL,
      SCSTEPTOL, CONSTRAINTS, OPTIN, IOPT,ROPT, IER)
The arguments are:
NEQ   = (INTEGER) number of equations (unknowns) in the nonlinear system
UU    = array containing the initial guess when called, returns the solution
GLOBALSTRAT = (INTEGER)a number defining the global strategy choice:
         0 = InexactNewton, 1 = LineSearch .
USCALE = array of scaling factors for the UU vector
FSCALE = array of scaling factors for the FVAL (function) vector
FNORMTOL = tolerance on the norm of f(u) to accept convergence.
SCSTEPTOL = tolerance on minimum scaled step size
CONSTRAINTS = array of constraint values, by element of the solution UU
OPTIN    = integer used as a flag to indicate whether possible input values
           in IOPT are to be used for input: 0 = NO, 1 = YES.
IOPT     = array for integer optional inputs and outputs
           (declare as INTEGER*4 or INTEGER*8 according to C type long int)
ROPT     = array of real optional inputs and outputs
IER      = INTEGER error flag as returned by KINSOL . See KINSOL documentation
           for further information.

(6) Memory freeing: FKINFREE and FKFREEMPI
To the free the internal memory created by the calls to FKINITMPI (parallel/
MPI version only) and FKINMALLOC, make the following calls, in this order:
      CALL FKINFREE
      CALL FKFREEMPI  (parallel/MPI version only)
 * * * - - - * * * - - - * * * - - - * * * - - - * * * - - - * * * - - - * * *
Informational output:
     Some of the optional outputs are NFE, NNI, NLI, NPS, NCFL, NPE, stored 
in the array iopt at NFE+1, NNI+1, SPGMR_NLI+1, SPGMR_NPS+1, SPGMR_NPE+1, 
SPGMR_NPS+1, SPGMR_NCFL+1, respectively. (See the KINSOL and KINSPGMR header
files for descriptions and information on other outputs. The +1 on each
index is due to the differences in labeling arrays between C and Fortran.)


****************************************************************************/
#ifndef _fkinsol_h
#define _fkinsol_h

#include "fcmixpar.h"


/* generic names are translated through the define statements below for a
  specific platform/compiler */

#if (CRAY)

#define K_PRECO        KPRECO
#define K_PSOL         KPSOL
#define F_KINITMPI     FKINITMPI
#define F_KFREEMPI     FKFREEMPI
#define F_PKINMALLOC   FPKINMALLOC
#define F_SKINMALLOC   FSKINMALLOC
#define F_KINSPGMR00   FKINSPGMR00
#define F_KINSPGMR01   FKINSPGMR01
#define F_KINSPGMR10   FKINSPGMR10
#define F_KINSPGMR11   FKINSPGMR11
#define F_KINSPGMR20   FKINSPGMR20
#define F_KINSPGMR21   FKINSPGMR21
#define F_KINSOL       FKINSOL
#define F_KINFREE      FKINFREE
#define K_FUN          KFUN
#define F_ATIMES       FATIMES
#define F_DQJAC        FDQJAC

#elif  (UNDERSCORE)

#define  K_PRECO       kpreco_
#define  K_PSOL        kpsol_
#define  F_KINITMPI    fkinitmpi_
#define  F_KFREEMPI    fkfreempi_
#define  F_PKINMALLOC  fpkinmalloc_
#define  F_SKINMALLOC  fskinmalloc_
#define  F_KINSPGMR00  fkinspgmr00_
#define  F_KINSPGMR01  fkinspgmr01_
#define  F_KINSPGMR10  fkinspgmr10_
#define  F_KINSPGMR11  fkinspgmr11_
#define  F_KINSPGMR20  fkinspgmr20_
#define  F_KINSPGMR21  fkinspgmr21_
#define  F_KINSOL      fkinsol_
#define  F_KINFREE     fkinfree_
#define  K_FUN         kfun_
#define  F_ATIMES      fatimes_
#define  F_DQJAC       fdqjac_

#else

#define  K_PRECO       kpreco
#define  K_PSOL        kpsol
#define  F_KINITMPI    fkinitmpi
#define  F_KFREEMPI    fkfreempi
#define  F_PKINMALLOC  fpkinmalloc
#define  F_SKINMALLOC  fskinmalloc
#define  F_KINSPGMR00  fkinspgmr00
#define  F_KINSPGMR01  fkinspgmr01
#define  F_KINSPGMR10  fkinspmgr10
#define  F_KINSPGMR11  fkinspgmr11
#define  F_KINSPGMR20  fkinspgmr20
#define  F_KINSPGMR21  fkinspgmr21
#define  F_KINSOL      fkinsol
#define  F_KINFREE     fkinfree
#define  K_FUN         kfun
#define  F_ATIMES      fatimes
#define  F_DQJAC       fdqjac

#endif

/* KINSOL header files  */

#include "llnltyps.h"  /* definitions of types real and integer             */
#include "kinsol.h"    /* definition of type SysFn                          */
#include "nvector.h"   /* definition of type N_Vector, machEnvType          */


/* Prototypes: Functions called by the solver */


void KINfunc(integer Neq, N_Vector uu, N_Vector fval, void *f_data);


int KINPreco(integer Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
	     N_Vector vtemp1, N_Vector vtemp2,
	     SysFn func, real u_round,
	     long int *nfePtr, void *P_data);


int KINPSol(integer Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
	     N_Vector vtem, N_Vector ftem,
	     SysFn func, real u_round,
	     long int *nfePtr, void *P_data);

int KINUAtimes(void *f_data, N_Vector v, N_Vector z, boole *new_uu,
               N_Vector uu);

/* Declarations for global variables, shared among various routines */

void *KIN_kmem;
machEnvType KIN_machEnv; /* typedef machEnvType is defined (differently) in 
                            the serial and parallel versions of nvector.h */
                             


#endif
