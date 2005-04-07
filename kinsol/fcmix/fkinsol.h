/*
 * -----------------------------------------------------------------
 * $Revision: 1.25.2.3 $
 * $Date: 2005-04-07 00:18:39 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file for the FKINSOL Interface Package.
 * See below for usage details.
 * -----------------------------------------------------------------
 */

/***************************************************************************

                  FKINSOL Interface Package

 The FKINSOL Interface Package is a package of C functions which support the
 use of the KINSOL solver for the solution of nonlinear systems f(u) = 0, 
 in a mixed Fortran/C setting. While KINSOL is written in C, it is assumed 
 here that the user's calling program and user-supplied problem-defining 
 routines are written in Fortran. This package provides the necessary
 interface to KINSOL for both the serial and the parallel NVECTOR
 implementations.

 The user-callable functions, with the corresponding KINSOL functions,
 are as follows:

   FNVINITS and FNVINITP : initialize serial and parallel vector
                           computations, respectively
   FKINMALLOC :  interfaces to KINMalloc 
   FKINSPGMR : interfaces to KINSpgmr
   FKINSOL : interfaces to KINSol
   FKINFREE : interfaces to KINFree
   FNVFREES and FNVFREEP : finalize serial and parallel vector
                           computations, respectively

 The user-supplied functions, each with the corresponding interface function
 which calls it (and its type within KINSOL), are as follows:

   FKFUN : called by the interface function FKINfunc of type SysFn
   FKJTIMES : called by the interface function FKINJtimes of type
              KINSpgmrJacTimesVecFn
   FKPSOL : called by the interface function FKINPSol of type
            KINSpgmrPrecSolveFn
   FKPSET : called by the interface function FKINPSet of type
            KINSpgmrPrecSetupFn

 In contrast to the case of direct use of KINSOL, the names of all 
 user-supplied routines here are fixed, in order to maximize portability for
 the resulting mixed-language program.

 =========================================================================

                  Usage of the FKINSOL Interface Package

 The usage of FKINSOL requires calls to several interface functions, and 
 to a few user-supplied routines which define the problem to be solved.
 These function calls and user routines are summarized separately below.

 Some details are omitted, and the user is referred to the KINSOL manual
 for more complete documentation. Information on the arguments of any
 given user-callable interface routine, or of a given user-supplied
 function called by an interface function, can be found in the
 documentation on the corresponding function in the KINSOL package.

 The number labels on the instructions below end with "s" for instructions
 that apply to the serial version of CVODE only, and end with "p" for
 those that apply to the parallel version only.

 (1) User-supplied system routine: FKFUN

     The user must in all cases supply the following Fortran routine:

       SUBROUTINE FKFUN (UU, FVAL)
       DIMENSION UU(*), FVAL(*)

     It must set the FVAL array to f(u), the system function, as a
     function of the array UU = u. Here UU and FVAL are arrays representing
     vectors, which are distributed vectors in the parallel case.

 (2) Optional user-supplied Jacobian-vector product routine: FKJTIMES

     As an option, the user may supply a routine that computes the product
     of the system Jacobian and a given vector. This has the following form:

       SUBROUTINE FKJTIMES(V, Z, NEWU, UU, IER)
       DIMENSION V(*), Z(*), UU(*)

     This must set the array Z to the product J*V, where J is the Jacobian
     matrix J = dF/du, and V is a given array. Here UU is an array containing
     the current value of the unknown vector u. NEWU is an input integer 
     indicating whether UU has changed since FKJTIMES was last called 
     (1 = yes, 0 = no). If FKJTIMES computes and saves Jacobian data, then 
     no such computation is necessary when NEWU = 0. Here V, Z, and UU are 
     arrays of length NEQ, the problem size, or the local length of all 
     distributed vectors in the parallel case. FKJTIMES should return IER = 0 
     if successful, or a nonzero IER otherwise.

 (3) Initialization:  FNVINITS/FNVINITP and FKINMALLOC

 (3.1s) To initialize the serial machine environment, the user must make
        the following call:

          CALL FNVINITS (NEQ, IER)

        The arguments are:
          NEQ = size of vectors
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (3.1p) To initialize the parallel machine environment, the user must make 
        the following call:

          CALL FNVINITP (NLOCAL, NGLOBAL, IER)

        The arguments are:
          NLOCAL  = local size of vectors for this process
          NGLOBAL = the system size, and the global size of vectors
                    (the sum of all values of NLOCAL)
          IER     = return completion flag. Values are 0 = success,
                    -1 = failure.

 (3.2) To set various problem and solution parameters and allocate
       internal memory, make the following call:

         CALL FKINMALLOC(MSBPRE, FNORMTOL, SCSTEPTOL, CONSTRAINTS,
                         OPTIN, IOPT, ROPT, IER)

       The arguments are:
         MSBPRE      = maximum number of preconditioning solve calls
                       without calling the preconditioning setup routine;
                       0 indicates default.
         FNORMTOL    = tolerance on the norm of f(u) to accept convergence
         SCSTEPTOL   = tolerance on minimum scaled step size
         CONSTRAINTS = array of constraint values, on components of the
                       solution UU
         INOPT       = integer used as a flag to indicate whether possible
                       input values in IOPT are to be used for input:
                       0 = no, 1 = yes.
         IOPT        = array for integer optional inputs and outputs
                       (declare as INTEGER*4 or INTEGER*8 according to
                       C type long int)
         ROPT        = array of real optional inputs and outputs
         IER         = return completion flag. Values are 0 = success, and
                       -1 = failure.

       Note: See printed message for details in case of failure.

 (4) Specification of linear system solution method:

     The solution method in KINSOL involves the solution of linear systems 
     related to the Jacobian J = dF/du of the nonlinear system.

 (4.1) SPGMR treatment of the linear systems:

       For the Scaled Preconditioned GMRES solution of the linear systems,
       the user must make the call:

         CALL FKINSPGMR(MAXL, MAXLRST, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         MAXLRST  = maximum number of linear system restarts; 0 indicates
                    default (SPGMR only).
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

       If the user program includes the FKJTIMES routine for the evaluation
       of the Jacobian vector product, the following call must be made:

         CALL FKINSPGMRSETJAC(FLAG, IER)

       The argument FLAG = 0 specifies using the internal finite differences
       approximation to the Jacobian vector product, while FLAG = 1 specifies
       that FKJTIMES is provided.

       Usage of the user-supplied routines FKPSET and FKPSOL for the setup and
       solution of the preconditioned linear system is specified by calling:

         CALL FKINSPGMRSETPREC(FLAG, IER)

       where FLAG = 0 indicates no FKPSET or FKPSOL (default) and FLAG = 1
       specifies using FKPSET and FKPSOL. The user-supplied routines FKPSET
       and FKPSOL must be of the form:

         SUBROUTINE FKPSET (UU, USCALE, FVAL, FSCALE, VTEMP1, VTEMP2, IER)
         DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEMP1(*), VTEMP2(*)

       It must perform any evaluation of Jacobian-related data and
       preprocessing needed for the solution of the preconditioned linear
       systems by FKPSOL. The variables UU through FSCALE are for use in the
       preconditioning setup process. Typically, the system function FKFUN is
       called, so that FVAL will have been updated. UU is the current solution
       iterate. VTEMP1 and VTEMP2 are available for work space. If scaling is
       being used, USCALE and FSCALE are available for those operatins
       requiring scaling. NEQ is the (global) problem size.

       On return, set IER = 0 if FKPSET was successful, set IER = 1 if
       an error occurred.

         SUBROUTINE FKPSOL (UU, USCALE, FVAL, FSCALE, VTEM, FTEM, IER)
         DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)

       Typically this routine will use only UU, FVAL, VTEM and FTEM.
       It must solve the preconditioned linear system Pz = r, where
       r = VTEM is input, and store the solution z in VTEM as well. Here
       P is the right preconditioner. If scaling is being used, the
       routine supplied must also account for scaling on either coordinate
       or function value.

 (5) The solver: FKINSOL

     Solving the nonlinear system is accomplished by making the following
     call:

       CALL FKINSOL (UU, GLOBALSTRAT, USCALE, FSCALE, IER)

     The arguments are:
       UU          = array containing the initial guess on input, and the
                     solution on return
       GLOBALSTRAT = (INTEGER) a number defining the global strategy choice:
                     1 = InexactNewton, 2 = LineSearch
       USCALE      = array of scaling factors for the UU vector
       FSCALE      = array of scaling factors for the FVAL (function) vector
       IER         = INTEGER error flag as returned by KINSOL:
                     0 means success,
                     1 means initial guess satisfies f(u) = 0 (approx.),
                     2 means apparent stalling (small step),
                     a value < 0 means other error or failure.

     Note: See KINSOL documentation for detailed information.

 (6) Memory freeing: FKINFREE and FNVFREES/FNVFREEP

     To the free the internal memory created by the calls to FKINMALLOC
     and either FNVINITS or FNVINITP, make the following calls, in this
     order:

 (6.1s) CALL FKINFREE
        CALL FNVFREES

 (6.1p) CALL FKINFREE
        CALL FNVFREEP

 (7) Optional inputs and outputs: IOPT/ROPT

     The optional inputs available by way of IOPT and ROPT have the
     following names, locations, and descriptions. For further details
     see the KINSOL documentation. Note: A zero value results in the
     default.

       PRINTFL         = IOPT(1) = optional output print flag
       MXITER          = IOPT(2) = maximum Newton iterations
       PRECOND_NO_INIT = IOPT(3) = flag to suppress initial preconditioner
                                   setup call
       ETACHOICE       = IOPT(8) = choice of forcing term (1 = Choice 1,
                                   2 = Choice 2 and 3 = constant)
       NO_MIN_EPS      = IOPT(9) = flag to suppress minimum tolerance (eps)
       MXNEWTSTEP      = ROPT(1) = maximum size of Newton step
       RELFUNC         = ROPT(2) = relative error in computing f(u)
       ETACONST        = ROPT(5) and
       ETAGAMMA        = ROPT(6) and
       ETAALPHA        = ROPT(7) = constants in optional choices of forcing
                                   terms

     The optional outputs available by way of IOPT and ROPT have the
     following names, locations, and descriptions. For further details see
     the KINSOL documentation.
 
       NNI    = IOPT(4) = number of Newton iterations
       NFE    = IOPT(5) = number of f evaluations
       NBCF   = IOPT(6) = number of line search beta condition failures
       NBKTRK = IOPT(7) = number of line search backtracks
       FNORM  = ROPT(3) = final scaled norm of f(u)
       STEPL  = ROPT(4) = scaled last step length

     The following optional outputs are specific to the SPGMR module:

       NLI    = IOPT(11) = number of linear (Krylov) iterations
       NPE    = IOPT(12) = number of preconditioner evaluations
       NPS    = IOPT(13) = number of preconditioner solves
       NCFL   = IOPT(14) = number of linear convergence failures
       LSFLAG = IOPT(15) = last flag returned by linear solver

*******************************************************************************/

#ifndef _FKINSOL_H
#define _FKINSOL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * header files
 * -----------------------------------------------------------------
 */

#include "nvector.h"        /* definition of type N_Vector */
#include "sundialstypes.h"  /* definition of type realtype */

/*
 * -----------------------------------------------------------------
 * generic names are translated through the define statements below
 * -----------------------------------------------------------------
 */

#if defined(F77_FUNC)

#define FKIN_MALLOC       F77_FUNC(fkinmalloc, FKINMALLOC)
#define FKIN_SPGMR        F77_FUNC(fkinspgmr, FKINSPGMR)
#define FKIN_SPGMRSETJAC  F77_FUNC(fkinspgmrsetjac, FKINSPGMRSETJAC)
#define FKIN_SPGMRSETPREC F77_FUNC(fkinspgmrsetprec, FKINSPGMRSETPREC)
#define FKIN_SOL          F77_FUNC(fkinsol, FKINSOL)
#define FKIN_FREE         F77_FUNC(fkinfree, FKINFREE)
#define FK_FUN            F77_FUNC(fkfun, FKFUN)
#define FK_PSET           F77_FUNC(fkpset, FKPSET)
#define FK_PSOL           F77_FUNC(fkpsol, FKPSOL)
#define FK_JTIMES         F77_FUNC(fkjtimes, FKJTIMES)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC       fkinmalloc
#define FKIN_SPGMR        fkinspgmr
#define FKIN_SPGMRSETJAC  fkinspgmrsetjac
#define FKIN_SPGMRSETPREC fkinspgmrsetprec
#define FKIN_SOL          fkinsol
#define FKIN_FREE         fkinfree
#define FK_FUN            fkfun
#define FK_PSET           fkpset
#define FK_PSOL           fkpsol
#define FK_JTIMES         fkjtimes

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC       FKINMALLOC
#define FKIN_SPGMR        FKINSPGMR
#define FKIN_SPGMRSETJAC  FKINSPGMRSETJAC
#define FKIN_SPGMRSETPREC FKINSPGMRSETPREC
#define FKIN_SOL          FKINSOL
#define FKIN_FREE         FKINFREE
#define FK_FUN            FKFUN
#define FK_PSET           FKPSET
#define FK_PSOL           FKPSOL
#define FK_JTIMES         FKJTIMES

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC       fkinmalloc_
#define FKIN_SPGMR        fkinspgmr_
#define FKIN_SPGMRSETJAC  fkinspgmrsetjac_
#define FKIN_SPGMRSETPREC fkinspgmrsetprec_
#define FKIN_SOL          fkinsol_
#define FKIN_FREE         fkinfree_
#define FK_FUN            fkfun_
#define FK_PSET           fkpset_
#define FK_PSOL           fkpsol_
#define FK_JTIMES         fkjtimes_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC       FKINMALLOC_
#define FKIN_SPGMR        FKINSPGMR_
#define FKIN_SPGMRSETJAC  FKINSPGMRSETJAC_
#define FKIN_SPGMRSETPREC FKINSPGMRSETPREC_
#define FKIN_SOL          FKINSOL_
#define FKIN_FREE         FKINFREE_
#define FK_FUN            FKFUN_
#define FK_PSET           FKPSET_
#define FK_PSOL           FKPSOL_
#define FK_JTIMES         FKJTIMES_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC       fkinmalloc__
#define FKIN_SPGMR        fkinspgmr__
#define FKIN_SPGMRSETJAC  fkinspgmrsetjac__
#define FKIN_SPGMRSETPREC fkinspgmrsetprec__
#define FKIN_SOL          fkinsol__
#define FKIN_FREE         fkinfree__
#define FK_FUN            fkfun__
#define FK_PSET           fkpset__
#define FK_PSOL           fkpsol__
#define FK_JTIMES         fkjtimes__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC       FKINMALLOC__
#define FKIN_SPGMR        FKINSPGMR__
#define FKIN_SPGMRSETJAC  FKINSPGMRSETJAC__
#define FKIN_SPGMRSETPREC FKINSPGMRSETPREC__
#define FKIN_SOL          FKINSOL__
#define FKIN_FREE         FKINFREE__
#define FK_FUN            FKFUN__
#define FK_PSET           FKPSET__
#define FK_PSOL           FKPSOL__
#define FK_JTIMES         FKJTIMES__

#endif

/*
 * -----------------------------------------------------------------
 * Prototypes : exported functions
 * -----------------------------------------------------------------
 */

void FKIN_MALLOC(long int *msbpre, realtype *fnormtol, realtype *scsteptol,
		 realtype *constraints, int *optin, long int *iopt,
		 realtype *ropt, int *ier);
void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier);
void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier);
void FKIN_FREE(void);
void FKIN_SPGMRSETJAC(int *flag, int *ier);
void FKIN_SPGMRSETPREC(int *flag, int *ier);

/*
 * -----------------------------------------------------------------
 * Prototypes : functions called by the solver
 * -----------------------------------------------------------------
 */

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

/*
 * -----------------------------------------------------------------
 * declarations for global variables shared amongst various
 * routines
 * -----------------------------------------------------------------
 */

extern N_Vector F2C_vec;
extern realtype *data_F2C_vec;
extern void *KIN_mem;
extern long int *KIN_iopt;
extern realtype *KIN_ropt;
extern int KIN_ls;  /* Linear Solver: 1 = SPGMR */

#ifdef __cplusplus
}
#endif

#endif
