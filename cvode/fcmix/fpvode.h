/* File fpvode.h: Header file for the FPVODE Interface Package
   Version of 1 March 2002 */

#ifndef _fpvode_h
#define _fpvode_h

/***************************************************************************

                  FPVODE Interface Package

The FPVODE Interface Package is a package of C functions which support
the use of the PVODE solver (MPI version), for the solution of ODE
systems dy/dt = f(t,y), in a mixed Fortran/C setting.  While PVODE is
written in C, it is assumed here that the user's calling program and
user-supplied problem-defining routines are written in Fortran.

The user-callable functions, with the corresponding CVODE/PVODE functions,
are as follows:
  FPVINITMPI interfaces to PVecInitMPI
  FPVMALLOC  interfaces to CVodeMalloc
  FPVREINIT  interfaces to CVReInit
  FCVDIAG    interfaces to CVDiag
  FCVSPGMR00, FCVSPGMR01, FCVSPGMR10, FCVSPGMR11, FCVSPGMR20, FCVSPGMR21,
    FCVREINSPGMR00, FCVREINSPGMR01, FCVREINSPGMR10, FCVREINSPGMR11,
    FCVREINSPGMR20, FCVREINSPGMR21,
             interface to CVSpgmr for the various options
  FCVODE     interfaces to CVode
  FCVDKY     interfaces to CVodeDky
  FCVFREE    interfaces to CVodeFree
  FPVFREEMPI interfaces to PVecFreeMPI

The user-supplied functions, each listed with the corresponding interface
function which calls it (and its type within CVODE), are as follows:
  PVFUN    is called by the interface function CVf of type RhsFn
  PVPSOL   is called by the interface function CVPSol of type CVSpgmrSolveFn
  PVPRECO  is called by the interface function CVPreco of type CVSpgmrPrecondFn
  PVJTIMES is called by the interface function CVJtimes of type CVSpgmrJtimesFn
In contrast to the case of direct use of PVODE, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program.

Important note on portability.
In this package, the names of the interface functions, and the names of
the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header file fpvode.h.  Those mapping definitions depend in turn on a
pair of parameters, CRAY and UNDERSCORE, defined in the header file
fcmixpar.h, which is machine-dependent.  The names into which the dummy
names are mapped are either upper or lower case, and may or may not have
an underscore appended, depending on these parameters.  Check, and if 
necessary modify, the file fcmixpar.h for a given machine environment.

****************************************************************************

                 Usage of the FPVODE Interface Package

The usage of FPVODE requires calls to five to seven interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.

Some details are omitted, and the user is referred to the user documents
on PVODE and CVODE for more complete documentation.  Information on the
arguments of any given user-callable interface routine, or of a given
user-supplied function called by an interface function, can be found in
the documentation on the corresponding function in the CVODE package.


(1) User-supplied right-hand side routine: PVFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE PVFUN (NLOC, T, Y, YDOT)
      DIMENSION Y(*), YDOT(*)
It must set the YDOT array to f(t,y), the right-hand side of the ODE
system, as function of T = t and the array Y = y.  Here Y and YDOT
are distributed vectors, and NLOC is the length of the segment local 
to the current processor.

(2) Optional user-supplied Jacobian-vector product routine: PVJTIMES
As an option, the user may supply a routine that computes the product
of the system Jacobian J = df/dy and a given vector v.  If supplied, it
must have the following form:
      SUBROUTINE PVJTIMES (NLOC, V, FJV, T, Y, FY, VNRM, EWT, H, UROUND,
     1                     NFE, WORK, IER)
      DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), WORK(*)
Typically this routine will use only NLOC, T, Y, V, and FJV.  It must
compute the product vector Jv, where the vector v is stored in V, and store
the product in FJV.  On return, set IER = 0 if PVJTIMES was successful,
and nonzero otherwise.

(3) Initialization:  FPVINITMPI, FPVMALLOC, FPVREINIT

(3.1) To initialize the use of the MPI (Message Passing Interface) library
by PVODE, the user must make the following call:
      CALL FPVINITMPI (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors on this processor
NGLOBAL = the ODE problem size, and the global size of vectors (the sum 
          of all values of NLOCAL)
IER     = return completion flag. Values are 0 = success, -1 = failure.
Note: If MPI was initialized by the user, the communicator must be
set to MPI_COMM_WORLD.  If not, this routine initializes MPI and sets
the communicator equal to MPI_COMM_WORLD.

(3.2) To set various problem and solution parameters and allocate
internal memory, make the following call:
      CALL FPVMALLOC(NEQ, T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1               IOPT, ROPT, IER)
The arguments are:
NEQ    = the global problem size
T0     = initial value of t
Y0     = array of initial conditions
METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
ITMETH = nonlinear iteration method: 1 = functional iteration, 2 = Newton iter.
IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array
RTOL   = relative tolerance (scalar)
ATOL   = absolute tolerance (scalar or array)
INOPT  = optional input flag: 0 = none, 1 = inputs used
IOPT   = array of length 40 for integer optional inputs and outputs
         (declare as INTEGER*4 or INTEGER*8 according to C type long int)
ROPT   = array of length 40 for real optional inputs and outputs
         The optional inputs are MAXORD, MXSTEP, MXHNIL, SLDET, H0, HMAX,
         HMIN, stored in IOPT(1), IOPT(2), IOPT(3), IOPT(14), ROPT(1),
         ROPT(2), ROPT(3), respectively.  If any of these optional inputs
         are used, set the others to zero to indicate default values.
         The optional outputs are NST, NFE, NSETUPS, NNI, NCFN, NETF, QU, QCUR,
         LENRW, LENIW, NOR, HU, HCUR, TCUR, TOLSF, stored in IOPT(4) .. IOPT(13),
         IOPT(15), ROPT(4) .. ROPT(7), resp.  See the CVODE manual for details. 
IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
         See printed message for details in case of failure.

(3.3) To re-initialize the PVODE solver for the solution of a new problem
of the same size as one already solved, make the following call:
      CALL FPVREINIT(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1               IOPT, ROPT, IER)
The arguments have the same names and meanings as those of FPVMALLOC,
except that NEQ has been omitted from the argument list (being unchanged
for the new problem).  FPVREINIT performs the same initializations as
FPVMALLOC, but does no memory allocation, using instead the existing
internal memory created by the previous FPVMALLOC call.  The call to
specify the linear system solution method may or may not be needed;
see paragraph (4.2) below.

(4) Specification of linear system solution method.
In the case of a stiff system, the implicit BDF method involves the solution
of linear systems related to the Jacobian J = df/dy of the ODE system.
PVODE presently includes two choices for the treatment of these systems, 
and the user must call a routine with a specific name to make the
desired choice.

(4.1) Diagonal approximate Jacobian.
This choice is appropriate when the Jacobian can be well approximated by
a diagonal matrix.  The user must make the call:
      CALL FCVDIAG (IER)
IER is an error return flag: 0 = success, -1 = memory failure.
There is no additional user-supplied routine.  Optional outputs specific
to the approximate diagonal Jacobian case are LRW and LIW, stored in
IOPT(16) and IOPT(17), respectively.  (See the CVODE manual for descriptions.)

(4.2) SPGMR treatment of the linear systems.
For the Scaled Preconditioned GMRES solution of the linear systems,
the user must make one of the following six calls:
      CALL FCVSPGMR00 (IGSTYPE, MAXL, DELT, IER)              if no
                   preconditioning is to be done and PVJTIMES is not supplied;
      CALL FCVSPGMR10 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)    if the
           preconditioner involves no data setup and PVJTIMES is not supplied;
      CALL FCVSPGMR20 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)    if the
           preconditioner involves data setup and PVJTIMES is not supplied;
      CALL FCVSPGMR01 (IGSTYPE, MAXL, DELT, IER)              if no
           preconditioning is to be done but PVJTIMES is supplied;
      CALL FCVSPGMR11 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)    if the
           preconditioner involves no data setup but PVJTIMES is supplied;
      CALL FCVSPGMR21 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)    if the
           preconditioner involves data setup and PVJTIMES is supplied.
(In the two-digit suffix on the name above, the first digit is the number of
preconditioner routines supplied, and the second digit is 1 or 0 according as
PVJTIMES is supplied or not.)

In all cases, the arguments are:
IPRETYPE = preconditioner type: 1 = left only, 2 = right only, 3 = both sides.
IGSTYPE  = Gram-schmidt process type: 0 = modified G-S, 1 = classical G-S.
MAXL     = maximum Krylov subspace dimension; 0 indicates default.
DELT     = linear convergence tolerance factor; 0.0 indicates default.
IER      = error return flag: 0 = success, -1 = memory allocation failure,
           -2 = illegal input.

     In the cases FCVSPGMR10, FCVSPGMR11, FCVSPGMR20, and FCVSPGMR21, the user
program must include the following routine for solution of the preconditioner
linear system:
      SUBROUTINE PVPSOL (NLOC, T, Y, FY, VT, GAMMA, EWT, DELTA, NFE, R, LR, Z, IER)
      DIMENSION Y(*), FY(*), VT(*), EWT(*), R(*), Z(*), 
Typically this routine will use only NLOC, T, Y, GAMMA, R, LR, and Z.  It
must solve the preconditioner linear system Pz = r, where r = R is input, 
and store the solution z in Z.  Here P is the left preconditioner if LR = 1
and the right preconditioner if LR = 2.  The preconditioner (or the product
of the left and right preconditioners if both are nontrivial) should be an 
approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).

     In the cases FCVSPGMR20 and FCVSPGMR21, the user program must also include
the following routine for the evaluation and preprocessing of the preconditioner:
      SUBROUTINE PVPRECO (NLOC, T, Y, FY, JOK, JCUR, GAMMA, EWT, H, UROUND, 
     1                   NFE, V1, V2, V3, IER)
      DIMENSION Y(*), FY(*), EWT(*), V1(*), V2(*), V3(*) 
Typically this routine will use only NLOC, T, Y, JOK, and GAMMA. It must
perform any evaluation of Jacobian-related data and preprocessing needed
for the solution of the preconditioner linear systems by PVPSOL.
The JOK argument allows for Jacobian data to be saved and reused:  If 
JOK = 0, this data should be recomputed from scratch.  If JOK = 1, a saved
copy of it may be reused, and the preconditioner constructed from it.
On return, set JCUR = 1 if Jacobian data was computed, and 0 otherwise.
Also on return, set IER = 0 if PVPRECO was successful, set IER positive if a 
recoverable error occurred, and set IER negative if a non-recoverable error
occurred.

     Optional outputs specific to the SPGMR case are NPE, NLI, NPS, NCFL,
LRW, and LIW, stored in IOPT(16) ... IOPT(21), respectively.  (See the CVODE
manual for descriptions.)

     If a sequence of problems is being solved using the SPGMR linear solver,
then following the call to FPVREINIT, a call to the FCVSPGMR* routine may or
may not be needed.  First, if the choice among the six SPGMR options is the
same and the input arguments are the same, no FCVSPGMR* call is needed.
If a different choice of options is desired, or there is a change in input
arguments other than MAXL, then the user program should call one of the
routines FCVREINSPGMR00, FCVREINSPGMR01, FCVREINSPGMR10, FCVREINSPGMR11,
FCVREINSPGMR20, or FCVREINSPGMR21.  In this case, the FCVREINSPGMR routine
reinitializes the SPGMR linear solver, but without reallocating its memory.
The arguments of each FCVREINSPGMR routine have the same names and meanings
as the corresponding FCVSPGMR routine.  Finally, if the value of MAXL is
being changed, then a call to one of the six FCVSPGMR* routines must be made,
where again a different choice of that routine is allowed.  

(5) The integrator: FCVODE
Carrying out the integration is accomplished by making calls as follows:
      CALL FCVODE (TOUT, T, Y, ITASK, IER)
The arguments are:
TOUT  = next value of t at which a solution is desired (input)
T     = value of t reached by the solver on output
Y     = array containing the computed solution on output
ITASK = task indicator: 0 = normal mode (overshoot TOUT and interpolate)
        1 = one-step mode (return after each internal step taken)
IER   = completion flag: 0 = success, values -1 ... -8 are various
        failure modes (see CVODE manual).
The current values of the optional outputs are available in IOPT and ROPT.

(6) Computing solution derivatives: FCVDKY
To obtain a derivative of the solution, of order up to the current method
order, make the following call:
      CALL FCVDKY (T, K, DKY)
The arguments are:
T   = value of t at which solution derivative is desired
K   = derivative order (0 .le. K .le. QU)
DKY = array containing computed K-th derivative of y on return

(7) Memory freeing: FCVFREE and FPVFREEMPI
To the free the internal memory created by the calls to FPVINITMPI and
FPVMALLOC, make the following calls, in this order:
      CALL FCVFREE
      CALL FPVFREEMPI

****************************************************************************/

#include "fcmixpar.h"   /* parameters for function name definitions */

/* Definitions of interface function names */

#if (CRAY)

#define FPV_INITMPI FPVINITMPI
#define FPV_MALLOC  FPVMALLOC
#define FPV_REINIT  FPVREINIT
#define FCV_DIAG    FCVDIAG
#define FCV_SPGMR00 FCVSPGMR00
#define FCV_SPGMR10 FCVSPGMR10
#define FCV_SPGMR20 FCVSPGMR20
#define FCV_SPGMR01 FCVSPGMR01
#define FCV_SPGMR11 FCVSPGMR11
#define FCV_SPGMR21 FCVSPGMR21
#define FCV_REINSPGMR00 FCVREINSPGMR00
#define FCV_REINSPGMR10 FCVREINSPGMR10
#define FCV_REINSPGMR20 FCVREINSPGMR20
#define FCV_REINSPGMR01 FCVREINSPGMR01
#define FCV_REINSPGMR11 FCVREINSPGMR11
#define FCV_REINSPGMR21 FCVREINSPGMR21
#define FCV_CVODE   FCVODE
#define FCV_DKY     FCVDKY
#define FCV_FREE    FCVFREE
#define FPV_FREEMPI FPVFREEMPI
#define FCV_FUN     PVFUN
#define FCV_PSOL    PVPSOL
#define FCV_PRECO   PVPRECO
#define FCV_JTIMES  PVJTIMES

#elif (UNDERSCORE)

#define FPV_INITMPI fpvinitmpi_
#define FPV_MALLOC  fpvmalloc_
#define FPV_REINIT  fpvreinit_
#define FCV_DIAG    fcvdiag_
#define FCV_SPGMR00 fcvspgmr00_
#define FCV_SPGMR10 fcvspgmr10_
#define FCV_SPGMR20 fcvspgmr20_
#define FCV_SPGMR01 fcvspgmr01_
#define FCV_SPGMR11 fcvspgmr11_
#define FCV_SPGMR21 fcvspgmr21_
#define FCV_REINSPGMR00 fcvreinspgmr00_
#define FCV_REINSPGMR10 fcvreinspgmr10_
#define FCV_REINSPGMR20 fcvreinspgmr20_
#define FCV_REINSPGMR01 fcvreinspgmr01_
#define FCV_REINSPGMR11 fcvreinspgmr11_
#define FCV_REINSPGMR21 fcvreinspgmr21_
#define FCV_CVODE   fcvode_
#define FCV_DKY     fcvdky_
#define FCV_FREE    fcvfree_
#define FPV_FREEMPI fpvfreempi_
#define FCV_FUN     pvfun_
#define FCV_PSOL    pvpsol_
#define FCV_PRECO   pvpreco_
#define FCV_JTIMES  pvjtimes_

#else

#define FPV_INITMPI fpvinitmpi
#define FPV_MALLOC  fpvmalloc
#define FPV_REINIT  fpvreinit
#define FCV_DIAG    fcvdiag
#define FCV_SPGMR00 fcvspgmr00
#define FCV_SPGMR10 fcvspgmr10
#define FCV_SPGMR20 fcvspgmr20
#define FCV_SPGMR01 fcvspgmr01
#define FCV_SPGMR11 fcvspgmr11
#define FCV_SPGMR21 fcvspgmr21
#define FCV_REINSPGMR00 fcvreinspgmr00
#define FCV_REINSPGMR10 fcvreinspgmr10 
#define FCV_REINSPGMR20 fcvreinspgmr20 
#define FCV_REINSPGMR01 fcvreinspgmr01 
#define FCV_REINSPGMR11 fcvreinspgmr11 
#define FCV_REINSPGMR21 fcvreinspgmr21 
#define FCV_CVODE   fcvode
#define FCV_DKY     fcvdky
#define FCV_FREE    fcvfree
#define FPV_FREEMPI fpvfreempi
#define FCV_FUN     pvfun
#define FCV_PSOL    pvpsol
#define FCV_PRECO   pvpreco
#define FCV_JTIMES  pvjtimes

#endif


/* CVODE header files  */

#include "llnltyps.h" /* definitions of types real and integer             */
#include "cvode.h"    /* definition of type RHSFn                          */
#include "nvector.h"  /* definition of type N_Vector, machEnvType          */


/* Prototypes: Functions Called by the CVODE Solver */

void CVf(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);

int CVPreco(integer N, real tn, N_Vector y, N_Vector fy, boole jok,
                   boole *jcurPtr, real gamma, N_Vector ewt, real h,
                   real uround, long int *nfePtr, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int CVPSol(integer N, real tn, N_Vector y, N_Vector fy, N_Vector vtemp,
                  real gamma, N_Vector ewt, real delta, long int *nfePtr,
                  N_Vector r, int lr, void *P_data, N_Vector z);

int CVJtimes(integer N, N_Vector v, N_Vector Jv, RhsFn f, 
             void *f_data, real t, N_Vector y, N_Vector fy,
             real vnrm, N_Vector ewt, real h, real uround, 
             void *jac_data, long int *nfePtr, N_Vector work);


/* Declarations for global variables, shared among various routines */

void *CV_cvodemem;
N_Vector CV_yvec;
machEnvType PV_machEnv;


#endif
