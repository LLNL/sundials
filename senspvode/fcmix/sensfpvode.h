/******************************************************************
 *                                                                *
 * File          : sensfpvode.h                                   *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * Header file for the SensFPVODE Interface Package               *
 *                                                                *
 ******************************************************************/


#ifndef _sensfpvode_h
#define _sensfpvode_h

/***************************************************************************

                  SENSFPVODE Interface Package

The SENSFPVODE Interface Package is a package of C functions which,
together with the FPVODE interface package, support the use of the
SensPVODE solver (MPI version), for the solution and sensitivity
analysis of ODE systems dy/dt = f(t,y,p), in a mixed Fortran/C
setting.  In this context, y is a vector that contains the state
variables, and p is an array of real-valued parameters.  While
SensPVODE is written in C, it is assumed here that the user's calling
program and user-supplied problem-defining routines are written in
Fortran.

The user-callable functions in this package are listed below, together
with the corresponding SensPVODE functions:
  SFPVMALLOC  interfaces to SensCVodeMalloc
  SFPVREINIT  interfaces to SensCVReInit
  SFCVSPGMR0, SFCVSPGMR1, and SFCVSPGMR2 interface to SensCVSpgmr
  SFCVFREE    interfaces to SensCVodeFree

The user-supplied functions, each with the corresponding interface function
which calls it (and its type within CVODE), are as follows:
  PVFUN   is called by the interface function SensCVf of type RhsFn
  PVPSOL  is called by the interface function SCVPSol of type CVSpgmrSolveFn
  PVPRECO is called by the interface function SCVPreco of type CVSpgmrPrecondFn
In contrast to the case of direct use of SensPVODE, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program.

Important note on the sensitivity version of user-supplied functions.
To perform sensitivity analysis in the mixed Fortan/C setting, it is
necessary to extend the call sequence to each of the user-supplied
functions (PVFUN, PVPSOL, PVPRECO) by appending the array of real-valued
parameters as an argument. For example, 

      SUBROUTINE PVFUN(NLOC, T, Y, YDOT, P) 

has P appended to the call sequence so that PVFUN can access the real
parameters through P and thereby evaluate y' = f(t,y,p).

Important note on portability.
In this package, the names of the interface functions, and the names of
the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header file sensfpvode.h.  Those mapping definitions depend in turn on a
pair of parameters, CRAY and UNDERSCORE, defined in the header file
fcmixpar.h, which is machine-dependent.  The names into which the dummy
names are mapped are either upper or lower case, and may or may not have
an underscore appended, depending on these parameters.  Check, and if 
necessary modify, the file fcmixpar.h for a given machine environment.

****************************************************************************

                 Usage of the SENSFPVODE Interface Package

The usage of SENSFPVODE requires calls to five to seven interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.

Some details are omitted, and the user is referred to the SensPVODE manual
for more complete documentation.  Information on the arguments of any
given user-callable interface routine, or of a given user-supplied
function called by an interface function, can be found in the
documentation on the corresponding function in the SensPVODE package.


(1) User-supplied right-hand side routine: PVFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE PVFUN (NLOC, T, Y, YDOT, P)
      DIMENSION Y(*), YDOT(*), P(*)
It must set the YDOT array to f(t,y,p), the right-hand side of the ODE
system, as function of T = t, the array Y = y, and the array P = p.  
Here Y and YDOT are distributed vectors, and NLOC is the length of the
segment local to the current processor.

(2) Initialization:  FPVINITMPI, SFPVMALLOC, SFPVREINIT

(2.1) To initialize the use of the MPI (Message Passing Interface) library
by SensPVODE, the user must make the following call:
      CALL FPVINITMPI (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors y, etc. on this processor
NGLOBAL = the number of ODEs contained in y' = f(t,y,p), and the sum of
          all values of NLOCAL.
IER     = return completion flag. Values are 0 = success, -1 = failure.
Note: If MPI was initialized by the user, the communicator must be
set to MPI_COMM_WORLD.  If not, this routine initializes MPI and sets
the communicator equal to MPI_COMM_WORLD.

(2.2) To set various problem and solution parameters and allocate
internal memory, make the following call:
      CALL SFPVMALLOC(NY, NS, NTOTAL, T0, Y0, METH, ITMETH, IATOL, RTOL, 
     1                ATOL, INOPT, IOPT, ROPT, IER, P, PBAR, RHOMAX)
The arguments are:
NY     = the number of ODEs contained in y' = f(t,y,p)
NS     = the number of sensitivity vectors to be computed
NTOTAL = the total number of ODEs to be solved. Normally NTOTAL is equal
         to NY*(NS+1), but this is not required.
T0     = initial value of t
Y0     = is a vector of length NTOTAL that contains the initial values for
         the NY state variables and NY*NS sensitivity variables at time t = T0.
METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
ITMETH = nonlinear iteration method: 1 = functional iteration, 2 = Newton iter.
IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array
RTOL   = relative tolerance (scalar)
ATOL   = absolute tolerance (scalar or array)
INOPT  = optional input flag: 0 = none, 1 = inputs used
IOPT   = array of length 40 for integer optional inputs and outputs
         (declare as INTEGER*4 or INTEGER*8 according to C type long int)
ROPT   = array of length 40 for real optional inputs and outputs
         The optional inputs are MAXORD, MXSTEP, MXHNIL, H0, HMAX, HMIN,
         stored in IOPT(1) .. IOPT(3) and ROPT(1) .. ROPT(3), respectively.
         If any of these optional inputs are used, set the others to zero
         to indicate default values.
         The optional outputs are NST, NFE, NSETUPS, NNI, NCFN, NETF, QU, QCUR,
         LENRW, LENIW, HU, HCUR, TCUR, TOLSF, stored in IOPT(4) .. IOPT(13) 
         and ROPT(4) .. ROPT(7), resp.  See the CVODE manual for details. 
IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
         See printed message for details in case of failure.
P      = an array of real-valued parameters. The length of P is (at least) NS.
PBAR   = an array of nonzero, real values that are used to scale the NS
         sensitivity vectors. Normally PBAR(I) = P(I), if P(I) is nonzero.
RHOMAX = a real value used for selecting a finite difference formula to 
         estimate the scaled sensitivity derivatives. RHOMAX = 0.0 is a
         reasonable default value; it selects a centered difference formula
         that uses 2 evaluations of f(t,y,p) to estimate each scaled 
         sensitivity derivative.

(2.3) To reinitialize the SensPVODE solver for the solution of a new
problem of the same size as one already solved, make the following call:
      CALL SFPVREINIT(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1                IOPT, ROPT, IER, P, PBAR, RHOMAX)
The arguments have the same names and meanings as those of SFPVMALLOC,
except that NY, NS, and NTOTAL have been omitted from the argument
list (being unchanged for the new problem). SFPVREINIT performs the
same initializations as SFPVMALLOC, but does no memory allocation,
using instead the existing internal memory created by the previous 
SFPVMALLOC call.
     
(3) Specification of linear system solution method (SPGMR)
In the case of a stiff system, the implicit BDF method involves the solution
of linear systems related to the Jacobian J = df/dy of the ODE system.
SensPVODE uses the Scaled Preconditioned GMRES (SPGMR) method for solving
these linear systems. To do so, the user must make one of the following
three calls:
      CALL SFCVSPGMR0 (IGSTYPE, MAXL, DELT)             if no preconditioning
                                                        is to be done;
      CALL SFCVSPGMR1 (IPRETYPE, IGSTYPE, MAXL, DELT)   if the preconditioner
                                                        involves no data setup;
      CALL SFCVSPGMR2 (IPRETYPE, IGSTYPE, MAXL, DELT)   if the preconditioner
                                                        involves data setup.
The arguments are:
IPRETYPE = preconditioner type: 1 = left only, 2 = right only, 3 = both sides.
IGSTYPE  = Gram-Schmidt process type: 0 = modified G-S, 1 = classical G-S.
MAXL     = maximum Krylov subspace dimension; 0 indicates default.
DELT     = linear convergence tolerance factor; 0.0 indicates default.
     In the second and third cases, the user program must include the 
following routine for solution of the preconditioner linear system:
      SUBROUTINE PVPSOL (NLOC, T, Y, FY, VT, GAMMA, EWT, DELTA, NFE, R, LR, Z, IER, PAR)
      DIMENSION Y(*), FY(*), VT(*), EWT(*), R(*), Z(*), PAR(*)
Note that PAR contains the array of real parameters p in y' = f(t,y,p).
Typically PVPSOL will use only NLOC, T, Y, GAMMA, R, LR, Z, and PAR. 
It must solve the preconditioner linear system Pz = r, where r = R is input, 
and store the solution z in Z.  Here P is the left preconditioner if LR = 1
and the right preconditioner if LR = 2.  The preconditioner (or the product
of the left and right preconditioners if both are nontrivial) should be an 
approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).
     In the third case (SFCVSPGMR2), the user program must also include the
following routine for the evaluation and preprocessing of the preconditioner:
      SUBROUTINE PVPRECO (NLOC, T, Y, FY, JOK, JCUR, GAMMA, EWT, H, UROUND, 
     1                   NFE, V1, V2, V3, IER, PAR)
      DIMENSION Y(*), FY(*), EWT(*), V1(*), V2(*), V3(*), PAR(*)
Again, note that PAR contains the array of real parameters p in y' = f(t,y,p).
Typically PVPRECO will use only NLOC, T, Y, JOK, and GAMMA. It must
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
LRW, and LIW, stored in IOPT(14) ... IOPT(19), respectively.  (See the CVODE
manual for descriptions.)

(4) The integrator: FCVODE
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

(5) Computing solution derivatives: FCVDKY
To obtain a derivative of the solution, of order up to the current method
order, make the following call:
      CALL FCVDKY (T, K, DKY)
The arguments are:
T   = value of t at which solution derivative is desired
K   = derivative order (0 .le. K .le. QU)
DKY = array containing computed K-th derivative of y on return

(6) Memory freeing: SFCVFREE and FPVFREEMPI
To the free the internal memory created by the calls to FPVINITMPI and
SFPVMALLOC, make the following calls, in this order:
      CALL SFCVFREE
      CALL FPVFREEMPI

****************************************************************************/

#include "fcmixpar.h"   /* parameters for function name definitions */
#include "fpvode.h"

/* Definitions of interface function names */

#if (CRAY)

#define SFPV_MALLOC  SFPVMALLOC
#define SFPV_REINIT  SFPVREINIT
#define SFCV_FREE    SFCVFREE
#define SFCV_FUN     PVFUN
#define SFCV_SPGMR0  SFCVSPGMR0
#define SFCV_SPGMR1  SFCVSPGMR1
#define SFCV_SPGMR2  SFCVSPGMR2
#define SFCV_PSOL    PVPSOL
#define SFCV_PRECO   PVPRECO

#elif (UNDERSCORE)

#define SFPV_MALLOC  sfpvmalloc_
#define SFPV_REINIT  sfpvreinit_
#define SFCV_FREE    sfcvfree_
#define SFCV_FUN     pvfun_
#define SFCV_SPGMR0  sfcvspgmr0_
#define SFCV_SPGMR1  sfcvspgmr1_
#define SFCV_SPGMR2  sfcvspgmr2_
#define SFCV_PSOL    pvpsol_
#define SFCV_PRECO   pvpreco_

#else

#define SFPV_MALLOC  sfpvmalloc
#define SFPV_REINIT  sfpvreinit
#define SFCV_FREE    sfcvfree
#define SFCV_FUN     pvfun
#define SFCV_SPGMR0  sfcvspgmr0
#define SFCV_SPGMR1  sfcvspgmr1
#define SFCV_SPGMR2  sfcvspgmr2
#define SFCV_PSOL    pvpsol
#define SFCV_PRECO   pvpreco

#endif


/* Prototypes: Functions Called by the CVODE Solver */

void SensCVf(integer N, real t, N_Vector y, N_Vector ydot, void *p_save);

int SCVPreco(integer N, real tn, N_Vector y, N_Vector fy, boole jok,
	     boole *jcurPtr, real gamma, N_Vector ewt, real h,
	     real uround, long int *nfePtr, void *P_data,
	     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int SCVPSol(integer N, real t, N_Vector y, N_Vector fy, N_Vector vtemp,
           real gamma, N_Vector ewt, real delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z);


/************************************************************************
 *                                                                      *
 * Type : pData                                                         *
 *----------------------------------------------------------------------*
 * The type pData is a type for the block of data that contains a       *
 * pointer to the array of real parameters p used in evaluating         *
 * y' = f(t,y,p).                                                       *
 * This block of data is created by a user call to SFPVMALLOC and is    *
 * not normally seen by the user.                                       *
 *                                                                      *
 * p is a pointer to the array of real parameters used in y' = f(t,y,p) *
 *                                                                      *
 ************************************************************************/

typedef struct {
  real *p;
} *pData;

#endif
