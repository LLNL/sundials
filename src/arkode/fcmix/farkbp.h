/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * This is the Fortran interface include file for the BAND
 * preconditioner (ARKBANDPRE).
 *--------------------------------------------------------------*/

/*===============================================================
                  FARKBP Interface Package

 The FARKBP Interface Package is a package of C functions which,
 together with the FARKODE Interface Package, support the use of 
 the ARKODE solver and N_VSerial vector module with the ARKBANDPRE
 preconditioner module, for the solution of ODE systems in a mixed
 Fortran/C setting.  The combination of ARKODE and ARKBANDPRE 
 solves the linear systems arising from the solution of the 
 implicit portions of the ODE system
       dy/dt = fe(t,y) + fi(t,y)  
 using a Krylov iterative linear solver (SPGMR, SPTFQMR, SPBCG, 
 SPFGMR or PCG), and with a banded preconditioner.  This 
 preconditioner can be constructed using either a user-supplied 
 routine, or automatically via finite differencing.
 
 The user-callable functions in this package, with the 
 corresponding ARKODE and ARKBBDPRE functions, are as follows: 

   Fortran              ARKODE
   -------------        ---------------------------
   FARKBPINIT           ARKBandPrecInit
   FARKBPOPT            (accesses optional outputs)
   -------------        ---------------------------
 
 In addition to the Fortran implicit right-hand side function 
 FARKIFUN, the user may (optionally) supply a routine FARKJTIMES 
 that is called by the C interface function FARKJtimes of type 
 ARKSpilsJtimesFn.  

 
 Important notes on portability:

 The names of all user-supplied routines here are fixed, in 
 order to maximize portability for the resulting mixed-language 
 program.

 Additionally, the names of the interface functions, and the 
 names of the Fortran user routines called by them, appear as 
 dummy names which are mapped to actual values by a series of 
 definitions in the header file farkbp.h.
 
 ================================================================
 
           Usage of the FARKODE/FARKBP Interface Packages
 
 The usage of the combined interface packages FARKODE and FARKBP 
 requires calls to seven to ten interface functions, and one or 
 two user-supplied routines which define the problem to be solved
 and indirectly define the preconditioner.  These function calls 
 and user routines are summarized separately below.
 
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
 
 (2) Optional user-supplied Jacobian-vector product: FARKJTIMES

     As an option, the user may supply a routine that computes 
     the product of the system Jacobian  J = dfi(t,y)/dy and a 
     given vector v.  If supplied, it must have the following 
     form:

       SUBROUTINE FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, 
      &                      WORK, IER)

     Typically this routine will use only NEQ, T, Y, V, and FJV. 
     It must compute the product vector J*v where the vector V, 
     and store the product in FJV.  

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

 (3) Initialization:  FNVINITS, FARKMALLOC, linear solver 
     initialization, and FARKBPINIT.
 
 (3.1) To initialize the serial vector specification, the user 
      must make the following call:

        CALL FNVINITS(4, NEQ, IER)

     where the first argument is an int containing the ARKODE 
     solver ID (4). The other arguments are:
        NEQ = size of vectors [long int, input]
	IER = return completion flag [int, output]:
	          0 = success, 
		 -1 = failure.
 
 (3.2) To set various problem and solution parameters and 
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
	IOUT = array of length 22 for integer optional outputs
	   [long int, output]
	ROUT = array of length 6 for real optional outputs 
	   [realtype, output]
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

 (3.3) Attach one of the 3 SPILS linear solvers.
    
 (3.3A) To specify the SPGMR (Scaled Preconditioned GMRES) 
     linear solver make the following call:

       CALL FARKSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)

     The arguments are:
        IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
	IGSTYPE = Gram-schmidt process type [int, input]: 
              1 = modified G-S
              2 = classical G-S.
	MAXL = maximum Krylov subspace dimension [int; input]; 
	      0 = default
	DELT = convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     If a sequence of problems of the same size is being solved 
     using the SPGMR linear solver, then following the call to 
     FARKREINIT, a call to the FARKSPGMRREINIT routine is needed
     if any of IPRETYPE, IGSTYPE, DELT is being changed.  In that
     case, call FARKSPGMRREINIT as follows:

       CALL FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKSPGMR.  If 
     MAXL is being changed, then the user should call FARKSPGMR 
     instead.
 
 (3.3B) To specify the SPBCG (Scaled Preconditioned Bi-CGSTAB) 
     linear solver make the following call:

       CALL FARKSPBCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 
              0 = default.
       DELT = convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     If a sequence of problems of the same size is being solved 
     using the SPBCG linear solver, then following the call to 
     FARKREINIT, a call to the FARKSPBCGREINIT routine is needed
     if any of its arguments is being changed.  The call is:

       CALL FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKSPBCG.

 (3.3C) To specify the SPTFQMR (Scaled Preconditioned TFQMR) 
     linear solver make the following call:

       CALL FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 
              0 = default.
       DELT = convergence tolerance factor [realtype, input]
	      0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     If a sequence of problems of the same size is being solved 
     using the SPTFQMR linear solver, then following the call to
     FARKREINIT, a call to the FARKSPTFQMRREINIT routine is 
     needed if any of its arguments is being changed.  The call 
     is:

       CALL FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)

     The arguments have the same meanings as for FARKSPTFQMR.

 (3.3D) To specify the SPFGMR (Scaled Preconditioned Flexible GMRES) 
     linear solver make the following call:

       CALL FARKSPFGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)

     The arguments are:
        IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
	IGSTYPE = Gram-schmidt process type [int, input]: 
              1 = modified G-S
              2 = classical G-S.
	MAXL = maximum Krylov subspace dimension [int; input]; 
	      0 = default
	DELT = convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     If a sequence of problems of the same size is being solved 
     using the SPFGMR linear solver, then following the call to 
     FARKREINIT, a call to the FARKSPFGMRREINIT routine is needed
     if any of IPRETYPE, IGSTYPE, DELT is being changed.  In that
     case, call FARKSPFGMRREINIT as follows:

       CALL FARKSPFGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKSPFGMR.  If 
     MAXL is being changed, then the user should call FARKSPFGMR 
     instead.
 
 (3.3E) To specify the PCG (Preconditioned Conjugate Gradient) 
     linear solver make the following call:

       CALL FARKPCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = use preconditioning
       MAXL = maximum Krylov subspace dimension [int, input]; 
              0 = default.
       DELT = convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     If a sequence of problems of the same size is being solved 
     using the PCG linear solver, then following the call to 
     FARKREINIT, a call to the FARKPCGREINIT routine is needed
     if any of its arguments is being changed.  The call is:

       CALL FARKPCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKPCG.

 (3.4) To allocate memory and initialize data associated with the
      ARKBANDPRE preconditioner, make the following call:

        CALL FARKBPINIT(NEQ, MU, ML, IER)

      The arguments are:
        NEQ = problem size [long int, input]
        MU = upper half-bandwidth of the band matrix that is 
	     retained as an approximation of the Jacobian 
	     [long int, input]
        ML = lower half-bandwidth of the band matrix approximant
	     to the Jacobian [long int, input]
        IER = return completion flag [int, output]:
                    0 = success
                   <0 = an error occurred

 (3.5) To specify whether the Krylov linear solver should use the
     supplied FARKJTIMES or the internal finite difference 
     approximation, make the call

        CALL FARKSPILSSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKJTIMES is provided 
     (FLAG=0 specifies to use and internal finite difference 
     approximation to this product).  The int return flag IER=0
     if successful, and nonzero otherwise.
 
 (4) The integrator: FARKODE

     Carrying out the integration is accomplished by making calls
     as follows:

       CALL FARKODE(TOUT, T, Y, ITASK, IER)

     The arguments are:
       TOUT = next value of t at which a solution is desired 
           [realtype, input]
       T = value of t reached by the solver [realtype, output]
       Y = state variable array on output [realtype, output]
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
 
 (5) Optional outputs: FARKBPOPT

     Optional outputs specific to the SP* linear solvers are:
        LENRWLS = IOUT(14) from ARKSp*GetWorkSpace (real space)
        LENIWLS = IOUT(15) from ARKSp*GetWorkSpace (int space)
        LSTF    = IOUT(16) from ARKSp*GetLastFlag
        NFELS   = IOUT(17) from ARKSp*GetRhsEvals
        NJTV    = IOUT(18) from ARKSp*GetJtimesEvals
        NPE     = IOUT(19) from ARKSp*GetPrecEvals
        NPS     = IOUT(20) from ARKSp*GetPrecSolves
        NLI     = IOUT(21) from ARKSp*GetLinIters
        NCFL    = IOUT(22) from ARKSp*GetConvFails
     See the ARKODE manual for descriptions.

     To obtain the optional outputs associated with the 
     ARKBANDPRE module, make the following call:

       CALL FARKBPOPT(LENRWBP, LENIWBP, NFEBP)

     The arguments returned are:
       LENRWBP = length of real preconditioner work space, in 
           realtype words (this size is local to the current 
	   processor if run in parallel) [long int, output]
       LENIWBP = length of integer preconditioner work space, in
           integer words (processor-local) [long int, output]
       NFEBP = number of fi(t,y) evaluations [long int, output]
 
 (6) Computing solution derivatives: FARKDKY

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
 
 (7) Memory freeing: FARKFREE

     To free the internal memory created by the calls to 
     FARKMALLOC, FNVINITS or FNVINITP, and FARKBPINIT, make the
     call:

       CALL FARKFREE()
 
===============================================================*/

#ifndef _FARKBP_H
#define _FARKBP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files  */
/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_BPINIT    SUNDIALS_F77_FUNC(farkbpinit, FARKBPINIT)
#define FARK_BPOPT     SUNDIALS_F77_FUNC(farkbpopt, FARKBPOPT)

#else

#define FARK_BPINIT    farkbpinit_
#define FARK_BPOPT     farkbpopt_

#endif

/* Prototypes of exported function */
void FARK_BPINIT(long int *N, 
		 long int *mu, 
		 long int *ml, 
		 int *ier);
void FARK_BPOPT(long int *lenrwbp, 
		long int *leniwbp, 
		long int *nfebp);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
