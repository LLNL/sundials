/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the Fortran interface include file for the BBD
 * preconditioner (ARKBBDPRE)
 *--------------------------------------------------------------*/

/*===============================================================
                   FARKBBD Interface Package

 The FARKBBD Interface Package is a package of C functions which,
 together with the FARKODE Interface Package, support the use of
 the ARKODE solver and MPI-parallel N_Vector module, along with
 the ARKBBDPRE preconditioner module, for the solution of ODE
 systems in a mixed Fortran/C setting.  The combination of ARKODE
 and ARKBBDPRE solves the linear systems arising from the solution
 of the implicit portions of the ODE system
       dy/dt = fe(t,y) + fi(t,y)
 using a Krylov iterative linear solver via the ARKSPILS interface,
 and with a preconditioner that is block-diagonal with banded
 blocks.  While ARKODE and ARKBBDPRE are written in C, it is
 assumed here that the user's calling program and user-supplied
 problem-defining routines are written in Fortran.

 The user-callable functions in this package, with the
 corresponding ARKODE and ARKBBDPRE functions, are as follows:

   Fortran               ARKODE
   --------------        ---------------------------
   FARKBBDINIT           ARKBBDPrecInit
   FARKBBDREINIT         ARKBBDPrecReInit
   FARKBBDOPT            (accesses optional outputs)
   --------------        ---------------------------

 In addition to the Fortran implicit right-hand side function
 FARKIFUN, the user-supplied functions used by this package are
 listed below, each with the corresponding interface function
 which calls it (and its type within ARKBBDPRE or ARKODE):

   Fortran           ARKODE            Type
   --------------    -----------       -----------------
   FARKGLOCFN        FARKgloc          ARKLocalFn
   FARKCOMMFN        FARKcfn           ARKCommFn
   FARKJTSETUP(*)    FARKJTSetup       ARKSpilsJTSetupFn
   FARKJTIMES(*)     FARKJtimes        ARKSpilsJtimesFn
   --------------    -----------       -----------------
   (*) = optional

 Important notes on portability:

 The names of all user-supplied routines here are fixed, in
 order to maximize portability for the resulting mixed-language
 program.

 Additionally, the names of the interface functions, and the
 names of the Fortran user routines called by them, appear as
 dummy names which are mapped to actual values by a series of
 definitions in the header file farkbbd.h.

 ================================================================

            Usage of the FARKODE/FARKBBD Interface Packages

 The usage of the combined interface packages FARKODE and
 FARKBBD requires calls to a variety of interface functions, and
 three or more user-supplied routines which define the problem to
 be solved and indirectly define the preconditioner.  These
 function calls and user routines are summarized separately below.

 Some details are omitted, and the user is referred to the ARKODE
 user document for more complete information.

 (1) User-supplied implicit right-hand side routine: FARKIFUN

     If any portion of the ODE system should be treated
     implicitly (and hence would require a linear solver at all),
     the user must supply the following Fortran routine:

       SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)

     Sets the YDOT array to fi(T,Y), the implicit portion of the
     right-hand side of the ODE system, as function of time T and
     the state variable array Y.

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       YDOT -- array containing state derivatives [realtype,
               output]
       IPAR -- array containing integer user data that was passed
               to FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.

 (2) User-supplied routines to define preconditoner: FARKGLOCFN
     and FARKCOMMFN

     The routines in the ARKBBDPRE module provide a
     preconditioner matrix for ARKODE that is block-diagonal
     with banded blocks.  The blocking corresponds to the
     distribution of the dependent variable vector y among the
     processors.  Each preconditioner block is generated from the
     Jacobian of the local part (on the current processor) of a
     given function g(t,y) approximating fi(t,y).  The blocks are
     generated by a difference quotient scheme on each processor
     independently, utilizing an assumed banded structure with
     given half-bandwidths.  A separate pair of half-bandwidths
     defines the band matrix retained.

 (2.1) Local approximate function FARKGLOCFN.

     The user must supply a subroutine of the form

       SUBROUTINE FARKGLOCFN(NLOC, T, YLOC, GLOC, IPAR, RPAR, IER)

     Sets the GLOC array to the function g(T,YLOC) which
     approximates the right-hand side function fi(T,YLOC).  This
     function is to be computed locally, i.e. without
     inter-processor communication (the case where g is
     mathematically identical to fi is allowed).

     The arguments are:
       NLOC -- local problem size [long int, input]
       T    -- current time [realtype, input]
       YLOC -- array containing local state variables
               [realtype, input]
       GLOC -- array containing local state derivatives
               [realtype, output]
       IPAR -- array containing integer user data that was passed
               to FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.

 (2.2) Communication function FARKCOMMFN.

     The user must also supply a subroutine of the form:

       SUBROUTINE FARKCOMMFN(NLOC, T, YLOC, IPAR, RPAR, IER)

     This performs all inter-processor communication necessary to
     evaluate g described above.  It is expected to save
     communicated data in work space defined by the user, and
     made available to FARKLOCFN.  Each call to the FARKCOMMFN is
     preceded by a call to FARKIFUN with the same (T,YLOC)
     arguments.  Thus FARKCOMMFN can omit any communications done
     by FARKIFUN if relevant to the evaluation of g.

     The arguments are:
       NLOC -- local problem size [long int, input]
       T    -- current time [realtype, input]
       YLOC -- array containing local state variables
               [realtype, input]
       IPAR -- array containing integer user data that was passed
               to FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.
 
 (3) Optional user-supplied Jacobian-vector setup and product 
     functions: FARKJTSETUP and FARKJTIMES

     As an option, the user may supply a routine that computes
     the product of the system Jacobian  J = dfi(t,y)/dy and a
     given vector v.  If supplied, a 'setup' routine to prepare
     any user data structures must exist, and have the form:
 
       SUBROUTINE FARKJTSETUP(T, Y, FY, H, IPAR, RPAR, IER)

     Typically this routine will use only T and Y.  It must perform any 
     relevant preparations for subsequent calls to the user-provided
     FARKJTIMES routine (see below).

     The arguments are:
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.

     The accompanying Jacobian matrix-vector product routine must
     have the following form:

       SUBROUTINE FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)

     Typically this routine will use only NEQ, T, Y, V, and FJV.
     It must compute the product vector J*v, and store the product in FJV.

     The arguments are:
       V    -- vector to multiply [realtype, input]
       FJV  -- product vector [realtype, output]
       T    -- current time [realtype, input]
       Y    -- state variables [realtype, input]
       FY   -- state derivatives [realtype, input]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed
               to FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WORK -- array containing temporary workspace of same size
               as Y [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.

 (4) Initialization:  FNVINITP, generic iterative linear solver
     initialization, FARKMALLOC, FARKSPILSINIT, and FARKBBDINIT.

 (4.1) To initialize the parallel vector specification, the user
     must make the following call:

        CALL FNVINITP(COMM, 4, NLOCAL, NGLOBAL, IER)

     where the second argument is an int containing the ARKODE 
     solver ID (4). The other arguments are:
        COMM = the MPI communicator [int, input]
        NLOCAL = local vector size on this processor 
           [long int, input]
        NGLOBAL = system size, and the global size of vectors 
           (the sum of all values of NLOCAL) [long int, input]
        IER = return completion flag [int, ouptut]. 
                  0 = success, 
                 -1 = failure.

 (4.2) To initialize a generic iterative linear solver structure for
     solving linear systems arising from implicit or IMEX treatment
     of the IVP, the user must make one of the following calls:

          CALL FSUNPCGINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPBCGSINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPFGMRINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPGMRINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPTFQMRINIT(4, PRETYPE, MAXL, IER)

     In each of these, one argument is an int containing the ARKODE solver
     ID (4).

     The other arguments are:

        PRETYPE = type of preconditioning to perform (0=none, 1=left,
           2=right, 3=both) [int, input]
        MAXL = maximum Krylov subspace dimension [int, input]
        IER = return completion flag [int, output]:
	          0 = success,
		 -1 = failure.

 (4.3) To set various problem and solution parameters and
     allocate internal memory, make the following call:

       CALL FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT,
      &                ROUT, IPAR, RPAR, IER)

     The arguments are:
        T0 = initial value of t [realtype, input]
	Y0 = array of initial conditions [realtype, input]
	IMEX = flag denoting integration method [int, input]:
                  0 = implicit,
                  1 = explicit,
                  2 = imex
        IATOL = type for absolute tolerance ATOL [int, input]:
                  1 = scalar,
                  2 = array,
                  3 = user-supplied function; the user must
                      supply a routine FARKEWT to compute the
		      error weight vector.
        RTOL = scalar relative tolerance [realtype, input]
	ATOL = scalar/array absolute tolerance [realtype, input]
	IOUT = array of length at least 29 for integer optional 
               outputs [long int, output]
	ROUT = array of length at least 6 for real optional 
               outputs [realtype, output]
	IPAR = array of user integer data [long int, in/out]
	RPAR = array with user real data [realtype, in/out]
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for details).

     The user data arrays IPAR and RPAR are passed unmodified to
     all subsequent calls to user-provided routines. Changes to
     either array inside a user-provided routine will be
     propagated. Using these two arrays, the user can dispense
     with COMMON blocks to pass data betwen user-provided
     routines.

 (4.4) Create the ARKSPILS interface to attach the generic 
     iterative linear solver to ARKode, by making the following call:
    
       CALL FARKSPILSINIT(IER)

     The arguments are:
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured

 (4.5) To allocate memory and initialize data associated with the
     ARKBBDPRE preconditioner, make the following call:

       CALL FARKBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)

      The arguments are:
        NLOCAL = local vector size on this process 
             [long int, input]
        MUDQ = upper half-bandwidth to be used in the computation
             of the local Jacobian blocks by difference 
             quotients.  These may be smaller than the true 
             half-bandwidths of the Jacobian of the local block 
             of g, when smaller values may provide greater 
             efficiency [long int, input]
        MLDQ = lower half-bandwidth to be used in the computation
             of the local Jacobian blocks by difference 
             quotients [long int, input]
        MU = upper half-bandwidth of the band matrix that is
             retained as an approximation of the local Jacobian
             block (may be smaller than MUDQ) [long int, input]
        ML = lower half-bandwidth of the band matrix that is
             retained as an approximation of the local Jacobian
             block (may be smaller than MLDQ) [long int, input]
        DQRELY = relative increment factor in y for difference 
             quotients [realtype, input]
                    0.0 = default (sqrt(unit roundoff))
        IER = return completion flag [int, output]:
                    0 = success
                   <0 = an error occurred

 (4.6) To specify whether the Krylov linear solver should use the
     supplied FARKJTSETUP and FARKJTIMES routines, or the internal 
     finite difference approximation, make the call

        CALL FARKSPILSSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKJTSETUP and FARKJTIMES 
     are provided (FLAG=0 specifies to use and internal finite 
     difference approximation to this product).  The int return 
     flag IER=0 if successful, and nonzero otherwise.
 
 (5) Re-initialization: FARKREINIT, FARKBBDREINIT

     If a sequence of problems of the same size is being solved
     using the Krylov linear solver in combination with the
     ARKBBDPRE preconditioner, then the ARKODE package can be
     reinitialized for the second and subsequent problems so as
     to avoid further memory allocation.  First, in place of the
     call to FARKMALLOC, make the following call:

       CALL FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)

     The arguments have the same names and meanings as those of
     FARKMALLOC.  FARKREINIT performs the same initializations
     as FARKMALLOC, but does no memory allocation, using instead
     the existing internal memory created by the previous
     FARKMALLOC call.  The subsequent call to specify the linear
     system solution method may or may not be needed.

     If there is no change in any of the linear solver or 
     preconditioner arguments, then no additional calls are 
     necessary.  

     Following the call to FARKREINIT, if there is no change in 
     any of the linear solver arguments, but the user wishes to 
     modify the values of MUDQ, MLDQ or DQRELY from the previous 
     call to FARKBBDINIT, then a user may call:

       CALL FARKBBDREINIT(MUDQ, MLDQ, DQRELY, IER)

     The arguments have the same names and meanings as those of 
     FARKBBDINIT.

     However, if there is a change in any of the linear solver 
     arguments or other preconditioner arguments, then a call to
     FSUNPCGINIT, FSUNSPBCGSINIT, FSUNSPFGMRINIT, FSUNSPGMRINIT, 
     or FSUNSPTFQMRINIT is required; in this case the linear 
     solver memory is reallocated.  Following this call, the 
     ARKSPILS interface must also be reconstructed using another
     call to FARKSPILSINIT (interface memory is freed and 
     reallocated), as well as a subsequent call to FARKBBDINIT.
 
 (6) The integrator: FARKODE

     Carrying out the integration is accomplished by making calls
     as follows:

       CALL FARKODE(TOUT, T, Y, ITASK, IER)

     The arguments are:
       TOUT = next value of t at which a solution is desired
           [realtype, input]
       T = value of t reached by the solver [realtype, output]
       Y = state variable array [realtype, output]
       ITASK = task indicator [int, input]:
                 1 = normal mode (overshoot TOUT and interpolate)
                 2 = one-step mode (return after each internal
                     step taken)
                 3 = normal tstop mode (like 1, but integration
                     never proceeds past TSTOP, which must be
                     specified through a call to FARKSETRIN using
                     the key 'STOP_TIME')
                 4 = one step tstop (like 2, but integration
                     never goes past TSTOP)
       IER = completion flag [int, output]:
                  0 = success,
                  1 = tstop return,
                  2 = root return,
                  values -1 ... -10 are failure modes (see
                    ARKODE manual).
     The current values of the optional outputs are immediately
     available in the IOUT and ROUT arrays.

 (7) Optional outputs: FARKBBDOPT

     Optional outputs specific to the ARKSpils linear solver
     interface are:
        LENRWLS = IOUT(14) from ARKSpilsGetWorkSpace
        LENIWLS = IOUT(15) from ARKSpilsGetWorkSpace
        LSTF    = IOUT(16) from ARKSpilsGetLastFlag
        NFELS   = IOUT(17) from ARKSpilsGetNumRhsEvals
        NJTV    = IOUT(18) from ARKSpilsGetNumJtimesEvals
        NPE     = IOUT(19) from ARKSpilsGetNumPrecEvals
        NPS     = IOUT(20) from ARKSpilsGetNumPrecSolves
        NLI     = IOUT(21) from ARKSpilsGetNumLinIters
        NCFL    = IOUT(22) from ARKSpilsGetNumConvFails
     See the ARKODE manual for descriptions.

     To obtain the optional outputs associated with the ARKBBDPRE
     module, make the following call:

       CALL FARKBBDOPT(LENRWBBD, LENIWBBD, NGEBBD)

     The arguments returned are:
       LENRWBP = length of real preconditioner work space, in
           realtype words (this size is local to the current
	   processor if run in parallel) [long int, output]
       LENIWBP = length of integer preconditioner work space, in
           integer words (processor-local) [long int, output]
       NGEBBD = number of g(t,y) evaluations (calls to ARKLOCFN)
           so far [long int, output]

 (8) Computing solution derivatives: FARKDKY

     To obtain a derivative of the solution, of order up to the
     method order, make the following call:

       CALL FARKDKY(T, K, DKY, IER)

     The arguments are:
       T = time at which solution derivative is desired, within
           the interval [TCUR-HU,TCUR], [realtype, input].
       K = derivative order (0 .le. K .le. QU) [int, input]
       DKY = array containing computed K-th derivative of y
           [realtype, output]
       IER = return flag [int, output]:
                    0 = success
                   <0 = illegal argument.
 
 (9) Memory freeing: FARKFREE

     To free the internal memory created by the calls to
     FARKMALLOC, FNVINITP, FARKSPILSINIT and FARKBBDINIT,
     make the call:

       CALL FARKFREE()

===============================================================*/

#ifndef _FARKBBD_H
#define _FARKBBD_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files  */
/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_BBDINIT    SUNDIALS_F77_FUNC(farkbbdinit,   FARKBBDINIT)
#define FARK_BBDREINIT  SUNDIALS_F77_FUNC(farkbbdreinit, FARKBBDREINIT)
#define FARK_BBDOPT     SUNDIALS_F77_FUNC(farkbbdopt,    FARKBBDOPT)
#define FARK_GLOCFN     SUNDIALS_F77_FUNC(farkglocfn,    FARKGLOCFN)
#define FARK_COMMFN     SUNDIALS_F77_FUNC(farkcommfn,    FARKCOMMFN)

#else

#define FARK_BBDINIT    farkbbdinit_
#define FARK_BBDREINIT  farkbbdreinit_
#define FARK_BBDOPT     farkbbdopt_
#define FARK_GLOCFN     farkglocfn_
#define FARK_COMMFN     farkcommfn_

#endif

/* Prototypes of exported functions */
void FARK_BBDINIT(long int *Nloc, long int *mudq,
		  long int *mldq, long int *mu,
		  long int *ml, realtype* dqrely, int *ier);
void FARK_BBDREINIT(long int *mudq, long int *mldq,
		    realtype* dqrely, int *ier);
void FARK_BBDOPT(long int *lenrwbbd, long int *leniwbbd,
		 long int *ngebbd);

/* Prototypes: Functions Called by the ARKBBDPRE Module */
int FARKgloc(long int Nloc, realtype t, N_Vector yloc,
	     N_Vector gloc, void *user_data);
int FARKcfn(long int Nloc, realtype t, N_Vector y, 
	    void *user_data);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
