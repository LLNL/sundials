/*************************************************************************
 *                                                                       *
 * File          : fkinbbd.h                                             *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL           *
 * Version of    : 17 January 2001                                       *
 *-----------------------------------------------------------------------*
 *                                                                       *
 * This is the Fortran interface include file for the BBD preconditioner *
 * (KINBBDPRE)                                                           *
 *                                                                       *
 ************************************************************************/
#ifndef _fkinbbd_h
#define _fkinbbd_h

/***************************************************************************

                  FKINBBD Interface Package

The FKINBBD Interface Package is a package of C functions which, together
with the FKINSOL Interface Package, support the use of the KINSOL solver 
(parallel MPI version) with the KINBBDPRE preconditioner module,
for the solution of nonlinear systems in a mixed Fortran/C setting.  The
combination of KINSOL and KINBBDPRE solves systems f(u) = 0. with
the SPGMR (scaled preconditioned GMRES) method for the linear systems
that arise, and with a preconditioner that is block-diagonal with
banded blocks.  While KINSOL and KINBBDPRE are written in C, it is
assumed here that the user's calling program and user-supplied
problem-defining routines are written in Fortran.

The user-callable functions in this package, with the corresponding
KINSOL and KINBBDPRE functions, are as follows: 
  FKINBBDINIT0  interfaces to KBBDAlloc and KINSpgmr, internal ATimes routine
  FKINBBDINIT1  interfaces to KBBDAlloc and KINSpgmr, external ATimes routine
                (stored in the file fkinbbdinit1.c)
  FKINBBDOPT    accesses optional outputs
  FKINBBDFREE   interfaces to KBBDFree

In addition to the Fortran system function KFUN, and optional Jacobian
vector product routine FATIMES, the following are the user-supplied 
functions required by this package, each with the corresponding
interface function which calls it (and its type within KINBBDPRE):
  KLOCFN is called by the interface function KINgloc of type KINLocalFn
  KCOMMFN is called by the interface function KINgcomm of type KINCommFn
(The names of all user-supplied routines here are fixed, in order to
maximize portability for the resulting mixed-language program.)

Important note on portability:
The names used in this interface package make use of the preprocessor to
expand them appropriately for different machines/platforms.
Later in this file, each name is expanded appropriately for those platforms.
An example:  F_KINSOL is replaced with either  KINSOL, kinsol_ or kinsol,
depending on the platform.  See fcmixpar.h . This documentation, however, will
simply assume the first form, for convenience. That is, KINSOL, etc.
****************************************************************************

              Usage of the FKINSOL/FKINBBD Interface Packages

The usage of combined interface packages FKINSOL and FKINBBD requires calls
to seven or eight interface functions, and three or four user-supplied
routines which define the problem to be solved and indirectly define
the preconditioner.  These function calls and user routines are
summarized separately below.

Some details are omitted, and the user is referred to the KINSOL User Guide
for more complete information.

(1) User-supplied system function routine: KFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE KFUN (NLOC, UU, FVAL)
      DIMENSION UU(*), FVAL(*)
It must set the FVAL array to f(u), the system function, as a function of
the array UU = u.  Here UU and FVAL are distributed vectors, and NLOC is
the length of the local part of these vectors on the current processor.

(2) Optional user-supplied Jacobian-vector product routine: FATIMES
As an option, the user may supply a routine that computes the product
of the system Jacobian and a given vector.  This has the form
      SUBROUTINE FATIMES(V, Z, NEWU, UU, IER)
      DIMENSION V(*), Z(*), UU(*)
This must set the array Z to the product J*V, where J is the Jacobian
matrix J = dF/du, and V is a given array.  Here UU is an array containing
the current value of the unknown vector u.  NEWU is an input integer 
indicating whether UU has changed since FATIMES was last called 
(1 = yes, 0 = no).  If FATIMES computes and saves Jacobian data, then 
no such computation is necessary when NEWU = 0.  Here V, Z, and UU are 
arrays of length NLOC, the local length of all distributed vectors.  
FATIMES should return IER = 0 if successful, or a nonzero IER otherwise.

(3) User-supplied routines to define preconditoner: KLOCFN and KCOMMFN

The routines in the KINBBDPRE (kinbbdpre.c) module provide a preconditioner 
matrix for KINSOL that is block-diagonal with banded blocks.  The blocking
corresponds to the distribution of the dependent variable vector u
among the processors.  Each preconditioner block is generated from the
Jacobian of the local part (on the current processor) of a given
function g(u) approximating f(u).  The blocks are generated by a
difference quotient scheme on each processor independently, utilizing
the assumed banded structure with given half-bandwidths.

(3.1) Local approximate function KLOCFN.
The user must supply a subroutine of the form
      SUBROUTINE KLOCFN (NLOC, ULOC, GLOC)
      DIMENSION ULOC(*), GLOC(*)
to compute the function g(u) which approximates the system function f(u). 
This function is to be computed locally, i.e. without  inter-processor
communication.  (The case where g is mathematically identical to f is allowed).
It takes as input the local vector length NLOC and the local real solution
array ULOC.  It is to compute the local part of g(u) and store this in the real array GLOC.

(3.2) Communication function KCOMMFN.
The user must also supply a subroutine of the form
      SUBROUTINE KCOMMFN (NLOC, ULOC)
      DIMENSION ULOC(*)
which is to perform all inter-processor communication necessary to
evaluate the approximate system function g described above.
This function takes as input the local vector length NLOC, and the local real
dependent variable array ULOC.  It is expected to save communicated data in 
work space defined by the user, and made available to KLOCFN.
Each call to the KCOMMFN function is preceded by a call to KFUN with the same
arguments.  Thus KCOMMFN can omit any communications done by KFUN
if relevant to the evaluation of g.

(4) Initialization:  FKINITMPI, FPKINMALLOC, and FKINBBDINIT0/FKINBBDINIT1.

(4.1) To initialize the use of the MPI (Message Passing Interface) library
by KINSOL, the user must make the following call:
      CALL FKINITMPI (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors on this processor
NGLOBAL = the ODE problem size, and the global size of vectors (the sum 
          of all values of NLOCAL)
IER     = return completion flag. Values are 0 = success, -1 = failure.

(4.2) To allocate internal memory for KINSOL, make the following call:
      CALL FPKINMALLOC(NEQ, IER)  or  CALL FSKINMALLOC(NEQ, IER)
in the parallel or serial case, respectively.  The arguments are:
NEQ    = problem size (global).
IER    = return completion flag.  Values are 0 = success, and -1 = failure.
         See printed message for details in case of failure.

(4.3) To specify the SPGMR linear system solver, and to allocate memory 
and initialize data associated with the SPGMR method and the BBD
preconditioner, make the following call:
      CALL FKINBBDINIT0 (MAXL, MAXLRST, MSBPRE, MU, ML, IER)
if not supplying a routine FATIMES, or
      CALL FKINBBDINIT1 (MAXL, MAXLRST, MSBPRE, MU, ML, IER) 
if supplying a routine FATIMES for Jacobian-vector products.

The arguments are:
MAXL     = maximum Krylov subspace dimension; 0 indicates default.
MAXLRST  = maximum number of linear solver restarts.
MSBPRE   = maximum number of preconditioner solve calls without a precondtioner
           setup call.
MU, ML   = upper and lower half-bandwidths to be used in the computation
           of the local Jacobian blocks.  These may be smaller than the
           true half-bandwidths of the Jacobian of the local block of g,
           when smaller values may provide greater efficiency.
IER      = return completion flag.  Values are 0 = success, and -1 = failure.

(5) Solver: FKINSOL
Carrying out the solving of the nonlinear system is accomplished by making 
calls as follows:
      CALL FKINSOL (NEQ, UU, GLOBALSTRAT, USCALE, FSCALE, FNORMTOL,
      SCSTEPTOL, CONSTRAINTS, OPTIN, IOPT,ROPT, IER)
The arguments are:
NEQ   = (integer) global number of unknowns in the nonlinear system
UU    = array containing the initial guess when called, returns the solution
GLOBALSTRAT = (integer) a number defining the global strategy choice:
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
IER      = integer error flag as returned by KINSOL . See KINSOL documentation
           for further information.

(6) Optional outputs: FKINBBDOPT
Optional outputs specific to the SPGMR solver are NPE, NLI, NPS, NCFL,
LRW, and LIW, stored in IOPT(14) ... IOPT(19), respectively.
To obtain the optional outputs associated with the KINBBDPRE module, make
the following call:
      CALL FKINBBDOPT (LENRPW, LENIPW, NGE)
The arguments returned are:
LENRPW = length of real preconditioner work space, in real words.
         This size is local to the current processor.
LENIPW = length of integer preconditioner work space, in integer words.
         This size is local to the current processor.
NGE    = number of g(u) evaluations (calls to KLOCFN) so far.


(7) Memory freeing: FKBBDFREE, FKINFREE, and FKFREEMPI
To the free the internal memory created by the calls to  
FKINBBDINIT0/FKINBBDINIT1, FKINITMPI, and FPKINMALLOC, make the 
following calls, in this order:
      CALL FKINBBDFREE
      CALL FKINFREE
      CALL FKFREEMPI

****************************************************************************/


#include "fcmixpar.h"

/* definitions of interface function names */

#if (CRAY)

#define F_KINBBDINIT0    FKINBBDINIT0
#define F_KINBBDINIT1    FKINBBDINIT1
#define K_COMMFN         KCOMMFN
#define K_LOCFN          KLOCFN
#define F_KINBBDOPT      FKINBBDOPT
#define F_KINBBDFREE     FKINBBDFREE


#elif  (UNDERSCORE)

#define F_KINBBDINIT0    fkinbbdinit0_
#define F_KINBBDINIT1    fkinbbdinit1_
#define K_COMMFN         kcommfn_
#define K_LOCFN          klocfn_
#define F_KINBBDOPT      fkinbbdopt_
#define F_KINBBDFREE     fkinbbdfree_

#else

#define F_KINBBDINIT0   fkinbbdinit0
#define F_KINBBDINIT1   fkinbbdinit1
#define K_COMMFN        kcommfn
#define K_LOCFN         klocfn
#define F_KINBBDOPT     fkinbbdopt
#define F_KINBBDFREE    fkinbbdfree

#endif



/* KINSOL header files  */

#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"   /* definition of type N_Vector                      */


/* Prototypes: Functions Called by the KINBBDPRE Module */

void KINgloc(integer Nloc, N_Vector uu, N_Vector gval, void *f_data);

void KINgcomm(integer Nloc, real *uloc, void *f_data);


/* Declarations for global variables, shared among various routines */

void *KBBD_Data;


#endif
