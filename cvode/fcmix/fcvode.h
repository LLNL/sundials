/* File fcvode.h: Header file for the FCVODE Interface Package
   Version of 26 June 2002 */

#ifndef _fcvode_h
#define _fcvode_h

/***************************************************************************

                  FCVODE Interface Package

The FCVODE Interface Package is a package of C functions which support
the use of the CVODE solver, for the solution of ODE
systems dy/dt = f(t,y), in a mixed Fortran/C setting.  While CVODE is
written in C, it is assumed here that the user's calling program and
user-supplied problem-defining routines are written in Fortran.
This package provides the necessary interface to CVODE for both the 
serial and the parallel NVECTOR implementations.

The user-callable functions, with the corresponding CVODE functions,
are as follows:
  FMENVINITS and FMENVINITP interface to M_EnvInit_Serial and
       M_EnvInitParallel (defined by nvector_serial
       and nvector_parallel, respectively)

  FCVMALLOC  interfaces to CVodeMalloc

  FCVREINIT  interfaces to CVReInit

  FCVDIAG    interfaces to CVDiag

  FCVDENSE0, FCVDENSE1, FCVREINDENSE0, FCVREINDENSE1
             interface to CVDense for the various option

  FCVBAND0, FCVBAND1, FCVREINBAND0, FCVREINBAND1
             interface to CVBand for the various options

  FCVSPGMR00, FCVSPGMR01, FCVSPGMR10, FCVSPGMR11, FCVSPGMR20, FCVSPGMR21,
  FCVREINSPGMR00, FCVREINSPGMR01, FCVREINSPGMR10, FCVREINSPGMR11,
  FCVREINSPGMR20, FCVREINSPGMR21,
             interface to CVSpgmr for the various options

  FCVODE     interfaces to CVode

  FCVDKY     interfaces to CVodeDky

  FCVFREE    interfaces to CVodeFree

  FMENVFREES and FMENVFREEP interface to M_EnvFree_Serial and
       M_EnvFree_parallel(defined by nvector_serial
       and nvector_parallel, respectively)

The user-supplied functions, each listed with the corresponding interface
function which calls it (and its type within CVODE), are as follows:
  CVFUN    is called by the interface function CVf of type RhsFn
  CVDJAC   is called by the interface function CVDenseJac of type CVDenseJacFn
  CVBJAC   is called by the interface function CVBandJac of type CVBandJacFn
  CVPSOL   is called by the interface function CVPSol of type CVSpgmrSolveFn
  CVPRECO  is called by the interface function CVPreco of type CVSpgmrPrecondFn
  CVJTIMES is called by the interface function CVJtimes of type CVSpgmrJtimesFn
In contrast to the case of direct use of CVODE, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program.

Important note on portability.
In this package, the names of the interface functions, and the names of
the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header file fcvode.h.  Those mapping definitions depend in turn on a
pair of parameters, CRAY and UNDERSCORE, defined in the header file
fcmixpar.h, which is machine-dependent.  The names into which the dummy
names are mapped are either upper or lower case, and may or may not have
an underscore appended, depending on these parameters.  Check, and if 
necessary modify, the file fcmixpar.h for a given machine environment.

=========================================================================

                 Usage of the FPVODE Interface Package

The usage of FCVODE requires calls to five to seven interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.

Some details are omitted, and the user is referred to the user documents
on CVODE for more complete documentation.  Information on the
arguments of any given user-callable interface routine, or of a given
user-supplied function called by an interface function, can be found in
the documentation on the corresponding function in the CVODE package.


(1) User-supplied right-hand side routine: CVFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE CVFUN (NEQ, T, Y, YDOT)
      DIMENSION Y(*), YDOT(*)
It must set the YDOT array to f(t,y), the right-hand side of the ODE
system, as function of T = t and the array Y = y.  Here Y and YDOT
are distributed vectors, and NEQ is the problem size.

(2) Optional user-supplied dense Jacobian approximation routine: CVDJAC
As an option when using the DENSE linear solver, the user may supply a
routine that computes a dense approximation of the system Jacobian 
J = df/dy. If supplied, it must have the following form:
      SUBROUTINE CVDJAC (NEQ, JAC, T, Y, FY, EWT, H, UROUND,
     1                   NFE, WORK1, WORK2, WORK3)
      DIMENSION JAC(*), Y(*), FY(*), EWT(*), WORK1(*), WORK2(*), WORK3(*)
Typically this routine will use only NEQ, T, Y, and JAC. It must compute
the Jacobian and store it, column-wise in JAC.

(3) Optional user-supplied band Jacobian approximation routine: CVBJAC
As an option when using the BAND linear solver, the user may supply a
routine that computes a band approximation of the system Jacobian 
J = df/dy. If supplied, it must have the following form:
      SUBROUTINE CVBJAC (NEQ, MU, ML, JAC, T, Y, FY, EWT, H, UROUND,
     1                   NFE, WORK1, WORK2, WORK3)
      DIMENSION JAC(*), Y(*), FY(*), EWT(*), WORK1(*), WORK2(*), WORK3(*)
Typically this routine will use only NEQ, MU, ML, T, Y, and JAC. 
It must compute the Jacobian and store it, column-wise in JAC.

(4) Optional user-supplied Jacobian-vector product routine: CVJTIMES
As an option when using the SPGMR linear solver, the user may supply a 
routine that computes the product of the system Jacobian J = df/dy and 
a given vector v.  If supplied, it must have the following form:
      SUBROUTINE CVJTIMES (NEQ, V, FJV, T, Y, FY, VNRM, EWT, H, UROUND,
     1                     NFE, WORK, IER)
      DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), WORK(*)
Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
compute the product vector Jv, where the vector v is stored in V, and store
the product in FJV.  On return, set IER = 0 if CVJTIMES was successful,
and nonzero otherwise.

(5) Initialization:  FMENVINITS / FMENVINITP , FCVMALLOC, FCVREINIT

(5.1s) To initialize the serial machine environment, the user must make
the following call:
       CALL FMENVINITS (NEQ, IER)
The arguments are:
NEQ     = size of vectors
IER     = return completion flag. Values are 0 = success, -1 = failure.

(5.1p) To initialize the parallel machine environment, the user must make 
the following call:
       CALL FMENVINITP (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors on this processor
NGLOBAL = the system size, and the global size of vectors (the sum 
          of all values of NLOCAL)
IER     = return completion flag. Values are 0 = success, -1 = failure.

(5.2) To set various problem and solution parameters and allocate
internal memory, make the following call:
      CALL FCVMALLOC(NEQ, T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
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

(5.3) To re-initialize the CVODE solver for the solution of a new problem
of the same size as one already solved, make the following call:
      CALL FCVREINIT(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1               IOPT, ROPT, IER)
The arguments have the same names and meanings as those of FCVMALLOC,
except that NEQ has been omitted from the argument list (being unchanged
for the new problem).  FCVREINIT performs the same initializations as
FCVMALLOC, but does no memory allocation, using instead the existing
internal memory created by the previous FCVMALLOC call.  The call to
specify the linear system solution method may or may not be needed;
see paragraph (4.2) below.

(6) Specification of linear system solution method.
In the case of a stiff system, the implicit BDF method involves the solution
of linear systems related to the Jacobian J = df/dy of the ODE system.
CVODE presently includes four choices for the treatment of these systems, 
and the user must call a routine with a specific name to make the
desired choice.

(6.1) Diagonal approximate Jacobian.
This choice is appropriate when the Jacobian can be well approximated by
a diagonal matrix.  The user must make the call:
      CALL FCVDIAG (IER)
IER is an error return flag: 0 = success, -1 = memory failure.
There is no additional user-supplied routine.  Optional outputs specific
to the approximate diagonal Jacobian case are LRW and LIW, stored in
IOPT(16) and IOPT(17), respectively.  (See the CVODE manual for descriptions.)

(6.2) DENSE treatment of the linear system.
The user must make one of the following two calls:
      CALL FCVDENSE0(IER)
          if CVDJAC is not supplied 

      CALL FCVDENSE1(IER)
          if CVDJAC is supplied 

In all cases, the argument is:
IER = error return flag: 0 = success , -1 = memory allocation failure,
           -2 = illegal input. 

  In the case FCVDENSE1, the user program must include the CVDJAC routine 
for the evaluation of the dense approximation to the Jacobian.

(6.3) BAND treatment of the linear system
The user must make one of the following two calls:
      CALL FCVBAND0(MU, ML, IER)
          if CVBJAC is not supplied

      CALL FCVBAND1(MU, ML, IER)
          if CVBJAC is supplied

In all cases, the arguments are:
MU  = upper bandwidth
ML  = lower bandwidth
IER = error return flag: 0 = success , -1 = memory allocation failure,
           -2 = illegal input.     

  In the case FCVBAND1, the user program must include the CVBJAC routine 
for the evaluation of the banded approximation to the Jacobian.

(6.4) SPGMR treatment of the linear systems.
For the Scaled Preconditioned GMRES solution of the linear systems,
the user must make one of the following six calls:
      CALL FCVSPGMR00 (IGSTYPE, MAXL, DELT, IER)              
          if no preconditioning is to be done and CVJTIMES is not supplied;

      CALL FCVSPGMR10 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
          if the preconditioner involves no data setup and CVJTIMES is not supplied;

      CALL FCVSPGMR20 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
           if the preconditioner involves data setup and CVJTIMES is not supplied;

      CALL FCVSPGMR01 (IGSTYPE, MAXL, DELT, IER)
           if the preconditioning is to be done but CVJTIMES is supplied;

      CALL FCVSPGMR11 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
           if the preconditioner involves no data setup but CVJTIMES is supplied;

      CALL FCVSPGMR21 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
           if the preconditioner involves data setup and CVJTIMES is supplied.

(In the two-digit suffix on the name above, the first digit is the number of
preconditioner routines supplied, and the second digit is 1 or 0 according as
CVJTIMES is supplied or not.)

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
      SUBROUTINE CVPSOL (NEQ, T, Y, FY, VT, GAMMA, EWT, DELTA, NFE, R, LR, Z, IER)
      DIMENSION Y(*), FY(*), VT(*), EWT(*), R(*), Z(*), 
Typically this routine will use only NEQ, T, Y, GAMMA, R, LR, and Z.  It
must solve the preconditioner linear system Pz = r, where r = R is input, 
and store the solution z in Z.  Here P is the left preconditioner if LR = 1
and the right preconditioner if LR = 2.  The preconditioner (or the product
of the left and right preconditioners if both are nontrivial) should be an 
approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).

     In the cases FCVSPGMR20 and FCVSPGMR21, the user program must also include
the following routine for the evaluation and preprocessing of the preconditioner:
      SUBROUTINE CVPRECO (NEQ, T, Y, FY, JOK, JCUR, GAMMA, EWT, H, UROUND, 
     1                   NFE, V1, V2, V3, IER)
      DIMENSION Y(*), FY(*), EWT(*), V1(*), V2(*), V3(*) 
Typically this routine will use only NEQ, T, Y, JOK, and GAMMA. It must
perform any evaluation of Jacobian-related data and preprocessing needed
for the solution of the preconditioner linear systems by CVPSOL.
The JOK argument allows for Jacobian data to be saved and reused:  If 
JOK = 0, this data should be recomputed from scratch.  If JOK = 1, a saved
copy of it may be reused, and the preconditioner constructed from it.
On return, set JCUR = 1 if Jacobian data was computed, and 0 otherwise.
Also on return, set IER = 0 if CVPRECO was successful, set IER positive if a 
recoverable error occurred, and set IER negative if a non-recoverable error
occurred.

     Optional outputs specific to the SPGMR case are NPE, NLI, NPS, NCFL,
LRW, and LIW, stored in IOPT(16) ... IOPT(21), respectively.  (See the CVODE
manual for descriptions.)

     If a sequence of problems is being solved using the SPGMR linear solver,
then following the call to FCVREINIT, a call to the FCVSPGMR* routine may or
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

(7) The integrator: FCVODE
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

(8) Computing solution derivatives: FCVDKY
To obtain a derivative of the solution, of order up to the current method
order, make the following call:
      CALL FCVDKY (T, K, DKY)
The arguments are:
T   = value of t at which solution derivative is desired
K   = derivative order (0 .le. K .le. QU)
DKY = array containing computed K-th derivative of y on return

(9) Memory freeing: FCVFREE and FMENVFREES / FMENVFREEP
To the free the internal memory created by the calls to FMENVINITS
 or FMENVINITP and FCVMALLOC, make the following calls, in this order:
      CALL FCVFREE
      CALL FMENVFREES or CALL FMENVFREEP  


****************************************************************************/

#include "fcmixpar.h"   /* parameters for function name definitions */

/* Definitions of interface function names */

#if (CRAY)

#define FCV_MALLOC  FCVMALLOC
#define FCV_REINIT  FCVREINIT
#define FCV_DIAG    FCVDIAG
#define FCV_DENSE0  FCVDENSE0
#define FCV_DENSE1  FCDENSE1
#define FCV_REINDENSE0 FCVREINDENSE0
#define FCV_REINDENSE1 FCVREINDENSE1
#define FCV_BAND0   FCVBAND0
#define FCV_BAND1   FCVBAND1
#define FCV_REINBAND0 FCVREINBAND0
#define FCV_REINBAND1 FCVREINBAND1
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
#define FCV_FUN     CVFUN
#define FCV_DJAC    CVDJAC
#define FCV_BJAC    CVBJAC
#define FCV_PSOL    CVPSOL
#define FCV_PRECO   CVPRECO
#define FCV_JTIMES  CVJTIMES

#elif (UNDERSCORE)

#define FCV_MALLOC  fcvmalloc_
#define FCV_REINIT  fcvreinit_
#define FCV_DIAG    fcvdiag_
#define FCV_DENSE0  fcvdense0_
#define FCV_DENSE1  fcvdense1_
#define FCV_REINDENSE0 fcvreindense0_
#define FCV_REINDENSE1 fcvreindense1_
#define FCV_BAND0   fcvband0_
#define FCV_BAND1   fcvband1_
#define FCV_REINBAND0 fcvreinband0_
#define FCV_REINBAND1 fcvreinband1_
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
#define FCV_FUN     cvfun_
#define FCV_DJAC    cvdjac_
#define FCV_BJAC    cvbjac_
#define FCV_PSOL    cvpsol_
#define FCV_PRECO   cvpreco_
#define FCV_JTIMES  cvjtimes_

#else

#define FCV_MALLOC  fcvmalloc
#define FCV_REINIT  fcvreinit
#define FCV_DIAG    fcvdiag
#define FCV_DENSE0  fcvdense0
#define FCV_DENSE1  fcvdense1
#define FCV_REINDENSE0 fcvreindense0
#define FCV_REINDENSE1 fcvreindense1
#define FCV_BAND0   fcvband0
#define FCV_BAND1   fcvband1
#define FCV_REINBAND0 fcvreinband0
#define FCV_REINBAND1 fcvreinband1
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
#define FCV_FUN     cvfun
#define FCV_DJAC    cvdjac
#define FCV_BJAC    cvbjac
#define FCV_PSOL    cvpsol
#define FCV_PRECO   cvpreco
#define FCV_JTIMES  cvjtimes

#endif


/* CVODE header files  */

#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "cvode.h"         /* definition of type RHSFn                      */
#include "nvector.h"       /* definition of type N_Vector, machEnvType      */
#include "dense.h"         /* definition of DenseMat                        */
#include "band.h"          /* definition of BandMat                         */

/* Prototypes: Functions Called by the CVODE Solver */

void CVf(integertype N, realtype t, N_Vector y, N_Vector ydot, void *f_data);

void CVDenseJac(integertype N, DenseMat J, RhsFn f, void *f_data,
                realtype t, N_Vector y, N_Vector fy, N_Vector ewt,
                realtype h, realtype uround, void *jac_data,
                long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);

void CVBandJac(integertype N, integertype mupper, integertype mlower,
               BandMat J, RhsFn f, void *f_data, realtype t,
               N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
               realtype uround, void *jac_data, long int *nfePtr,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int CVPreco(integertype N, realtype tn, N_Vector y, N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma, N_Vector ewt, realtype h,
            realtype uround, long int *nfePtr, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int CVPSol(integertype N, realtype tn, N_Vector y, N_Vector fy, N_Vector vtemp,
           realtype gamma, N_Vector ewt, realtype delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z);

int CVJtimes(integertype N, N_Vector v, N_Vector Jv, RhsFn f, 
             void *f_data, realtype t, N_Vector y, N_Vector fy,
             realtype vnrm, N_Vector ewt, realtype h, realtype uround, 
             void *jac_data, long int *nfePtr, N_Vector work);


/* Declarations for global variables, shared among various routines */

void *CV_cvodemem;
N_Vector CV_yvec;

#endif
