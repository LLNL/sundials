/******************************************************************
 *                                                                *
 * File          : sensfpvbbd.h                                   *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 25 August 2000                                 *
 *----------------------------------------------------------------*
 * Header file for the SENSFPVBBD Interface Package               *
 *                                                                *
 ******************************************************************/


#ifndef _sensfpvbbd_h
#define _sensfpvbbd_h

/***************************************************************************

                  SENSFPVBBD Interface Package

The SENSFPVBBD Interface Package is a package of C functions which,
together with the SENSFPVODE Interface Package, support the use of the
SensPVODE solver (MPI version) with the PVBBDPRE preconditioner
module for the solution and sensitivity analysis of ODE systems in a
mixed Fortran/C setting.  The ODE systems have the form dy/dt = f(t,y,p),
where y is a vector that contains the state variables and p is an
array of real-valued parameters.  For the linear systems that arise 
within SensPVODE, the routines in PVBBDPRE provide a preconditioner 
that is block-diagonal with banded blocks.  
While SensPVODE and PVBBDPRE are written in C, it is assumed here that
the user's calling program and user-supplied problem-defining routines
are written in Fortran.

The user-callable functions in this package, with the corresponding
SensPVODE and PVBBDPRE functions, are as follows: 
  SFPVBBDIN  interfaces to PVBBDAlloc and SensCVSpgmr 
  FPVBBDOPT  accesses optional outputs
  FPVBBDF    interfaces to PVBBDFree

In addition to the Fortran right-hand side function PVFUN, the
following are the user-supplied functions required by this package,
each with the corresponding interface function which calls it (and its
type within PVBBDPRE):
  PVLOCFN  is called by the interface function SPVgloc of type PVLocalFn
  PVCOMMF  is called by the interface function PVcfn of type PVCommFn
(The names of all user-supplied routines here are fixed, in order to
maximize portability for the resulting mixed-language program.)

Important note on the sensitivity version of user-supplied functions.
To perform sensitivity analysis in the mixed Fortan/C setting, it is
necessary to extend the call sequence to the user-supplied function
PVLOCFN by appending the array of real-valued parameters as an
argument. For example, 

      SUBROUTINE PVLOCFN (NLOC, T, YLOC, GLOC, P)

has P appended to the call sequence so that PVLOCFN can access the real
parameters through P if this information is needed to compute a
function g(t,y,p) which approximates the right-hand side function f(t,y,p).

Important note on portability.
In this package, the names of the interface functions, and the names of
the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header file sensfpvbbd.h.  Those mapping definitions depend in turn on a
pair of parameters, CRAY and UNDERSCORE, defined in the header file
fcmixpar.h, which is machine-dependent.  The names into which the dummy
names are mapped are either upper or lower case, and may or may not have
an underscore appended, depending on these parameters.  Check, and if 
necessary modify, the file fcmixpar.h for a given machine environment.

****************************************************************************

              Usage of the SENSFPVODE/SENSFPVBBD Interface Packages

The usage of combined interface packages SENSFPVODE and SENSFPVBBD requires
calls to seven to nine interface functions, and three user-supplied
routines which define the problem to be solved and indirectly define
the preconditioner.  These function calls and user routines are
summarized separately below.

Some details are omitted, and the user is referred to the SensPVODE 
User Guide for more complete information.

(1) User-supplied right-hand side routine: PVFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE PVFUN (N, T, Y, YDOT, P)
      DIMENSION Y(*), YDOT(*), P(*)
It must set the YDOT array to f(t,y,p), the right-hand side of the ODE 
system, as function of T = t, the array Y = y, and the array P = p.
Here Y and YDOT are distributed vectors, and NLOC is the length of the
segment local to the current processor.

(2) User-supplied routines to define preconditioner: PVLOCFN and PVCOMMF

The routines in the PVBBDPRE module provide a preconditioner matrix
for SensPVODE that is block-diagonal with banded blocks.  The blocking
corresponds to the distribution of the dependent variable vector y
among the processors.  Each preconditioner block is generated from the
Jacobian of the local part (on the current processor) of a given
function g(t,y,p) approximating f(t,y,p).  The blocks are generated by a
difference quotient scheme on each processor independently, utilizing
the assumed banded structure with given half-bandwidths.

(2.1) Local approximate function PVLOCFN.
The user must supply a subroutine of the form
      SUBROUTINE PVLOCFN (NLOC, T, YLOC, GLOC, P)
      DIMENSION YLOC(*), GLOC(*), P(*)
to compute the function g(t,y,p) which approximates the right-hand
side function f(t,y,p).  This function is to be computed locally,
i.e. without inter-processor communication.  (The case where g is
mathematically identical to f is allowed.)  It takes as input the
local vector length NLOC, the independent variable value T = t, and
the local real dependent variable array YLOC, and P = parameter array
p. It is to compute the local part of g(t,y,p) and store this in the
real array GLOC.

(2.2) Communication function PVCOMMF.
The user must also supply a subroutine of the form
      SUBROUTINE PVCOMMF (NLOC, T, YLOC)
      DIMENSION YLOC(*)
which is to perform all inter-processor communication necessary to
evaluate the approximate right-hand side function g described above.
This function takes as input the local vector length NLOC, the
independent variable value T = t, and the local real dependent
variable array YLOC.  It is expected to save communicated data in 
work space defined by the user, and made available to PVLOCFN.
Each call to the PVCOMMF is preceded by a call to PVFUN with the same
(t,y,p) arguments.  Thus PVCOMMF can omit any communications done by PVFUN
if relevant to the evaluation of g.

(3) Initialization:  FPVINITMPI, SFPVMALLOC, and SFPVBBDIN

(3.1) To initialize the use of the MPI (Message Passing Interface) library
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

(3.2) To set various problem and solution parameters and allocate
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

(3.3) To specify the SPGMR linear system solver, and to allocate memory 
and initialize data associated with the SPGMR method and the PVBBDPRE
preconditioner, make the following call:
      CALL SFPVBBDIN (MUDQ, MLDQ, MU, ML, DQRELY, IPRETYPE, IGSTYPE, MAXL,
     1               DELT, IER)
The arguments are:
MUDQ,MLDQ = upper and lower half-bandwidths to be used in the computation
            of the local Jacobian blocks by difference quotients.
            These may be smaller than the true half-bandwidths of the
            Jacobian of the local block of g, when smaller values may
            provide greater efficiency.
MU, ML    = upper and lower half-bandwidths of the band matrix that 
            is retained as an approximation of the local Jacobian block.
            These may be smaller than MUDQ and MLDQ.
DQRELY    = relative increment factor in y for difference quotients
            (optional). 0.0 indicates the default, sqrt(unit roundoff).
IPRETYPE  = preconditioner type: 1 = left only, 2 = right only, 3 = both sides.
IGSTYPE   = Gram-schmidt process type: 0 = modified G-S, 1 = classical G-S.
MAXL      = maximum Krylov subspace dimension; 0 indicates default.
DELT      = linear convergence tolerance factor; 0.0 indicates default.
IER       = return completion flag.  Values are 0 = success, and -1 = failure.

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
        failure modes (see CVODE User Guide).
The current values of the optional outputs are available in IOPT and ROPT.

(5) Optional outputs: FPVBBDOPT
Optional outputs specific to the SPGMR solver are NPE, NLI, NPS, NCFL,
LRW, and LIW, stored in IOPT(14) ... IOPT(19), respectively.
To obtain the optional outputs associated with the PVBBDPRE module, make
the following call:
      CALL FPVBBDOPT (LENRPW, LENIPW, NGE)
The arguments returned are:
LENRPW = length of real preconditioner work space, in real words.
         This size is local to the current processor.
LENIPW = length of integer preconditioner work space, in integer words.
         This size is local to the current processor.
NGE    = number of g(t,y,p) evaluations (calls to PVLOCFN) so far.

(6) Computing solution derivatives: FCVDKY
To obtain a derivative of the solution (optionally), of order up to
the current method order, make the following call:
      CALL FCVDKY (T, K, DKY)
The arguments are:
T   = value of t at which solution derivative is desired
K   = derivative order (0 .le. K .le. QU)
DKY = array containing computed K-th derivative of y on return

(7) Memory freeing: FPVBBDF, SFCVFREE, and FPVFREEMPI
To the free the internal memory created by the calls to  SFPVBBDIN, 
FPVINITMPI, and SFPVMALLOC, make the following calls, in this order:
      CALL FPVBBDF
      CALL SFCVFREE
      CALL FPVFREEMPI

****************************************************************************/

#include "fcmixpar.h"   /* parameters for function name definitions  */
#include "fpvbbd.h"     /* sensitivity version of PVBBDPRE interface */

/* Definitions of interface function names */

#if (CRAY)

#define SFPV_BBDIN   SFPVBBDIN
#define SFPV_GLOCFN  PVLOCFN

#elif (UNDERSCORE)

#define SFPV_BBDIN   sfpvbbdin_
#define SFPV_GLOCFN  pvlocfn_

#else

#define SFPV_BBDIN   sfpvbbdin
#define SFPV_GLOCFN  pvlocfn

#endif


/* CVODE header files  */

#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definition of type N_Vector                       */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "senscvspgmr.h" /* SensCVSpgmr prototype      			   */
#include "sensitivity.h" /* sensitivity data types and prototype           */
			     

/* Prototypes: Functions Called by the PVBBDPRE Module */

void SPVgloc(integer Nloc, real t, real *yloc, real *gloc, void *s_data);

#endif
