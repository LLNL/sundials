/***************************************************************************
 *                                                                         *
 * File        : fkinbbd.h                                                 *
 * Programmers : Allan G Taylor, Alan C. Hindmarsh, and Radu Serban @ LLNL *
 * Version of  : 27 January 2004                                           *
 *-------------------------------------------------------------------------*
 *                                                                         *
 * This is the Fortran interface include file for the BBD preconditioner   *
 * (KINBBDPRE)                                                             *
 *                                                                         *
 ***************************************************************************/
#ifndef _fkinbbd_h
#define _fkinbbd_h

/*******************************************************************************

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
   FKINBBDINIT  interfaces to KBBDPrecAlloc and KINSpgmr
   FKINBBDOPT   accesses optional outputs
   FKINBBDFREE  interfaces to KBBDPrecFree

 In addition to the Fortran system function KFUN, and optional Jacobian
 vector product routine KJTIMES, the following are the user-supplied 
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
 An example:  F_KINSOL is replaced with either  FKINSOL, fkinsol_ or fkinsol,
 depending on the platform.  See fcmixpar.h. This documentation, however,
 will simply assume the first form, for convenience, that is, FKINSOL, etc.

 ==============================================================================

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
      SUBROUTINE KFUN (NEQ, UU, FVAL)
      DIMENSION UU(*), FVAL(*)
 It must set the FVAL array to f(u), the system function, as a function of
 the array UU = u.  Here UU and FVAL are distributed vectors, and NEQ is
 the problem dimension.

 (2) Optional user-supplied Jacobian-vector product routine: KJTIMES
 As an option, the user may supply a routine that computes the product
 of the system Jacobian and a given vector.  This has the form
      SUBROUTINE KJTIMES(V, Z, NEWU, UU, IER)
      DIMENSION V(*), Z(*), UU(*)
 This must set the array Z to the product J*V, where J is the Jacobian
 matrix J = dF/du, and V is a given array.  Here UU is an array containing
 the current value of the unknown vector u.  NEWU is an input integer 
 indicating whether UU has changed since KJTIMES was last called 
 (1 = yes, 0 = no).  If KJTIMES computes and saves Jacobian data, then 
 no such computation is necessary when NEWU = 0.  Here V, Z, and UU are 
 arrays of length NLOC, the local length of all distributed vectors.  
 KJTIMES should return IER = 0 if successful, or a nonzero IER otherwise.

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
 array ULOC.  It is to compute the local part of g(u) and store this in the 
 realtype array GLOC.

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

 (4) Initialization:  FNVINITP, FKINMALLOC, and FKINBBDINIT.

 (4.1) To initialize the parallel machine environment, the user must make 
 the following call:
       CALL FNVINITP (NLOCAL, NGLOBAL, IER)
 The arguments are:
 NLOCAL  = local size of vectors on this processor
 NGLOBAL = the system size, and the global size of vectors (the sum 
           of all values of NLOCAL)
 IER     = return completion flag. Values are 0 = success, -1 = failure.

 (4.2) To allocate internal memory for KINSOL, make the following call:
       CALL FKINMALLOC (MSBPRE, FNORMTOL, SCSTEPTOL, CONSTRAINTS,
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


 (4.3) To specify the SPGMR linear system solver, and to allocate memory 
 and initialize data associated with the SPGMR method and the BBD
 preconditioner, make the following call:
      CALL FKINBBDINIT (NLOCAL, MAXL, MAXLRST, MU, ML, IER)

 The arguments are:
 NLOCAL   = local size of vectors on this processor
 MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 MAXLRST  = maximum number of linear solver restarts.
 MU, ML   = upper and lower half-bandwidths to be used in the computation
            of the local Jacobian blocks.  These may be smaller than the
            true half-bandwidths of the Jacobian of the local block of g,
            when smaller values may provide greater efficiency.
 IER      = return completion flag.  Values are 0 = success, and -1 = failure.

 (5) Solver: FKINSOL
 Solving of nonlinear system is accomplished by making the following call:
      CALL FKINSOL (UU, GLOBALSTRAT, USCALE, FSCALE, IER)
 The arguments are:
 UU          = array containing the initial guess when called, returns the solution
 GLOBALSTRAT = (INTEGER) a number defining the global strategy choice:
               0 = InexactNewton, 1 = LineSearch .
 USCALE      = array of scaling factors for the UU vector
 FSCALE      = array of scaling factors for the FVAL (function) vector
 IER         = integer error flag as returned by KINSOL . See KINSOL documentation
               for further information.

 (6) Optional outputs: FKINBBDOPT
 In addition to the optional inputs and outputs available with the FKINSOL
 interface package, there are optional outputs specific to the KINBBDPRE
 module.  These are accessed by making the following call:
      CALL FKINBBDOPT (LENRPW, LENIPW, NGE)
 The arguments returned are:
 LENRPW = length of real preconditioner work space, in realtype words.
          This size is local to the current processor.
 LENIPW = length of integer preconditioner work space, in integertype words.
          This size is local to the current processor.
 NGE    = number of g(u) evaluations (calls to KLOCFN) so far.

 (7) Memory freeing: FKINBBDFREE, FKINFREE, and FNVFREEP
 To the free the internal memory created by the calls to  
 FKINBBDINIT, FNVINITP, and FKINMALLOC, make the 
 following calls, in this order:
      CALL FKINBBDFREE
      CALL FKINFREE
      CALL FNVFREEP

*******************************************************************************/

/* definitions of interface function names */

#if defined(SUNDIALS_UNDERSCORE_NONE)

#define FKIN_BBDINIT    fkinbbdinit
#define FKIN_BBDOPT     fkinbbdopt
#define FKIN_BBDFREE    fkinbbdfree
#define FK_COMMFN       fkcommfn
#define FK_LOCFN        fklocfn

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#define FKIN_BBDINIT    fkinbbdinit__
#define FKIN_BBDOPT     fkinbbdopt__
#define FKIN_BBDFREE    fkinbbdfree__
#define FK_COMMFN       fkcommfn__
#define FK_LOCFN        fklocfn__

#else

#define FKIN_BBDINIT    fkinbbdinit_
#define FKIN_BBDOPT     fkinbbdopt_
#define FKIN_BBDFREE    fkinbbdfree_
#define FK_COMMFN       fkcommfn_
#define FK_LOCFN        fklocfn_

#endif


/* KINSOL header files  */

#include "sundialstypes.h"  /* definitions of types realtype and integertype */
#include "nvector.h"        /* definition of type N_Vector                   */


/* Prototypes: Functions Called by the KINBBDPRE Module */

void FKINgloc(integertype Nloc, N_Vector uu, N_Vector gval, void *f_data);

void FKINgcomm(integertype Nloc, N_Vector uu, void *f_data);


/* Declaration for global variable shared among various routines */

void *KBBD_Data;


#endif
