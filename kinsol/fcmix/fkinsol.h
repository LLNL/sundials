/*
 * -----------------------------------------------------------------
 * $Revision: 1.41 $
 * $Date: 2006-01-25 22:18:31 $
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

   FNVINITS and FNVINITP initialize serial and parallel vector
                         computations, respectively
   FKINMALLOC interfaces to KINMalloc
   FKINSETIIN, FKINSETRIN, FKINSETVIN interface to KINSet* functions
   FKINDENSE interfaces to KINDense
   FKINSPTFQMR interfaces to KINSptfqmr
   FKINSPGMR interfaces to KINSpgmr
   FKINSPBCG interfaces to KINSpbcg
   FKINSOL interfaces to KINSol and KINGet* functions
   FKINFREE interfaces to KINFree

 The user-supplied functions, each with the corresponding interface function
 which calls it (and its type within KINSOL), are as follows:

   FKFUN    : called by the interface function FKINfunc of type KINSysFn
   FKDJAC   : called by the interface function FKINDenseJac of type
              KINDenseJacFn
   FKBJAC   : called by the interface function FKINBandJac of type
              KINBandJacFn
   FKJTIMES : called by the interface function FKINJtimes of type
              KINSpilsJacTimesVecFn
   FKPSOL   : called by the interface function FKINPSol of type
              KINSpilsPrecSolveFn
   FKPSET   : called by the interface function FKINPSet of type
              KINSpilsPrecSetupFn

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

       SUBROUTINE FKFUN (UU, FVAL, IER)
       DIMENSION UU(*), FVAL(*)

     It must set the FVAL array to f(u), the system function, as a
     function of the array UU = u. Here UU and FVAL are arrays representing
     vectors, which are distributed vectors in the parallel case.
     IER is a return flag (currently not used).

 (2s) Optional user-supplied dense Jacobian approximation routine: FKDJAC
  
     As an option when using the DENSE linear solver, the user may supply a
     routine that computes a dense approximation of the system Jacobian 
     J = df/dy. If supplied, it must have the following form:
        
       SUBROUTINE FKDJAC(N, UU, FU, DJAC, WK1, WK2, IER)
       DIMENSION UU(*), FU(*), DJAC(N,*), WK1(*), WK2(*)

     This routine must compute the Jacobian and store it columnwise in DJAC.
     FKDJAC should return IER = 0 if successful, or a nonzero IER otherwise.

 (3s) Optional user-supplied band Jacobian approximation routine: FKBJAC
  
     As an option when using the BAND linear solver, the user may supply a
     routine that computes a band approximation of the system Jacobian 
     J = df/dy. If supplied, it must have the following form:
        
       SUBROUTINE FKBJAC(N, MU, ML, MDIM, UU, FU, BJAC, WK1, WK2, IER)
       DIMENSION UU(*), FU(*), BJAC(MDIM,*), WK1(*), WK2(*)

     This routine must load the MDIM by N array BJAC with the Jacobian matrix.
     FKBJAC should return IER = 0 if successful, or a nonzero IER otherwise.

 (4) Optional user-supplied Jacobian-vector product routine: FKJTIMES

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

 (5) Initialization:  FNVINITS/FNVINITP and FKINMALLOC

 (5.1s) To initialize the serial machine environment, the user must make
        the following call:

          CALL FNVINITS (3, NEQ, IER)

        The arguments are:
          NEQ = size of vectors
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (5.1p) To initialize the parallel machine environment, the user must make 
        the following call:

          CALL FNVINITP (3, NLOCAL, NGLOBAL, IER)

        The arguments are:
          NLOCAL  = local size of vectors for this process
          NGLOBAL = the system size, and the global size of vectors
                    (the sum of all values of NLOCAL)
          IER     = return completion flag. Values are 0 = success,
                    -1 = failure.

 (5.2) To allocate internal memory, make the following call:

         CALL FKINMALLOC(IOUT, ROUT, IER)

       The arguments are:
         IOUT        = array of length at least 15 for integer optional outputs
                       (declare as INTEGER*4 or INTEGER*8 according to
                       C type long int)
         ROUT        = array of length at least 2 for real optional outputs
         IER         = return completion flag. Values are 0 = success, and
                       -1 = failure.

       Note: See printed message for details in case of failure.

 (5.3) To set various integer optional inputs, make the folowing call:

          CALL FKINSETIIN(KEY, VALUE, IER)

       to set the optional input specified by the character key KEY to the 
       integer value VALUE.
       KEY is one of the following: PRNT_LEVEL, MAX_NITERS, ETA_FORM, 
       MAX_SETUPS, MAX_SP_SETUPS, NO_INIT_SETUP, NO_MIN_EPS, NO_RES_MON.

       To set various real optional inputs, make the folowing call:

         CALL FKINSETRIN(KEY, VALUE, IER)

      to set the optional input specified by the character key KEY to the
      real value VALUE.
      KEY is one of the following: FNORM_TOL, SSTEP_TOL, MAX_STEP, RERR_FUNC,
      ETA_CONST, ETA_PARAMS, RMON_CONST, RMON_PARAMS.
      Note that if KEY is ETA_PARAMS or RMON_PARAMS, then VALUE must be an
      array of dimension 2.

      To set the vector of constraints on the solution, make the following call:

        CALL FKINSETVIN(KEY, ARRAY, IER)

      where ARRAY is an array of reals and KEY is 'CONSTR_VEC'.

      FKINSETIIN, FKINSETRIN, and FKINSETVIN return IER=0 if successful and 
      IER<0 if an error occured.

 (6) Specification of linear system solution method:

     The solution method in KINSOL involves the solution of linear systems 
     related to the Jacobian J = dF/du of the nonlinear system.

 (6.1s) DENSE treatment of the linear systems (NVECTOR_SERIAL only):

       The user must make the following call:

         CALL FKINDENSE(NEQ, IER)

       In the above routine, the arguments are as follows:
         NEQ = problem size.
         IER = return completion flag.

       If the user program includes the FKDJAC routine for the evaluation
       of the dense approximation to the system Jacobian, the following call
       must be made:

         CALL FKINDENSESETJAC(FLAG, IER)

       with FLAG = 1 to specify that FKDJAC is provided.  (FLAG = 0 specifies
       using the internal finite difference approximation to the Jacobian.)

 (6.2s) BAND treatment of the linear systems (NVECTOR_SERIAL only):

       The user must make the following call:

         CALL FKINBAND(NEQ, MU, ML, IER)

       In the above routine, the arguments are as follows:
         NEQ = problem size.
         MU  = upper half-bandwidth
         ML  = lower half-bandwidth
         IER = return completion flag.

       If the user program includes the FKBJAC routine for the evaluation
       of the band approximation to the system Jacobian, the following call
       must be made:

         CALL FKINBANDSETJAC(FLAG, IER)

       with FLAG = 1 to specify that FKBJAC is provided.  (FLAG = 0 specifies
       using the internal finite difference approximation to the Jacobian.)

 (6.3) SPTFQMR treatment of the linear systems:

       For the Scaled Preconditioned TFQMR solution of the linear systems,
       the user must make the call:

         CALL FKINSPTFQMR(MAXL, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

       If the user program includes the FKJTIMES routine for the evaluation
       of the Jacobian-vector product, the following call must be made:

         CALL FKINSPTFQMRSETJAC(FLAG, IER)

       The argument FLAG = 0 specifies using the internal finite differences
       approximation to the Jacobian-vector product, while FLAG = 1 specifies
       that FKJTIMES is provided.

       Usage of the user-supplied routine FKPSOL for solution of the
       preconditioned linear system is specified by calling:

         CALL FKINSPTFQMRSETPREC(FLAG, IER)

       where FLAG = 0 indicates no FKPSOL or FKPSET (default) and FLAG = 1
       specifies using both FKPSOL and FKPSET.

       The user-supplied routine FKPSOL must be of the form:

         SUBROUTINE FKPSOL (UU, USCALE, FVAL, FSCALE, VTEM, FTEM, IER)
         DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)

       Typically this routine will use only UU, FVAL, VTEM and FTEM.
       It must solve the preconditioned linear system Pz = r, where
       r = VTEM is input, and store the solution z in VTEM as well. Here
       P is the right preconditioner. If scaling is being used, the
       routine supplied must also account for scaling on either coordinate
       or function value.

       The user-supplied routine FKPSET must be of the form:

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

 (6.4) SPBCG treatment of the linear systems:

       For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
       the user must make the call:

         CALL FKINSPBCG(MAXL, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

       If the user program includes the FKJTIMES routine for the evaluation
       of the Jacobian-vector product, the following call must be made:

         CALL FKINSPBCGSETJAC(FLAG, IER)

       The argument FLAG = 0 specifies using the internal finite differences
       approximation to the Jacobian-vector product, while FLAG = 1 specifies
       that FKJTIMES is provided.

       Usage of the user-supplied routine FKPSOL for solution of the
       preconditioned linear system is specified by calling:

         CALL FKINSPBCGSETPREC(FLAG, IER)

       where FLAG = 0 indicates no FKPSOL or FKPSET (default) and FLAG = 1
       specifies using both FKPSOL and FKPSET.

       The user-supplied routine FKPSOL must be of the form:

         SUBROUTINE FKPSOL (UU, USCALE, FVAL, FSCALE, VTEM, FTEM, IER)
         DIMENSION UU(*), USCALE(*), FVAL(*), FSCALE(*), VTEM(*), FTEM(*)

       Typically this routine will use only UU, FVAL, VTEM and FTEM.
       It must solve the preconditioned linear system Pz = r, where
       r = VTEM is input, and store the solution z in VTEM as well. Here
       P is the right preconditioner. If scaling is being used, the
       routine supplied must also account for scaling on either coordinate
       or function value.

       The user-supplied routine FKPSET must be of the form:

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

 (6.5) SPGMR treatment of the linear systems:

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
       of the Jacobian-vector product, the following call must be made:

         CALL FKINSPGMRSETJAC(FLAG, IER)

       The argument FLAG = 0 specifies using the internal finite differences
       approximation to the Jacobian-vector product, while FLAG = 1 specifies
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

 (7) The solver: FKINSOL

     Solving the nonlinear system is accomplished by making the following
     call:

       CALL FKINSOL (UU, GLOBALSTRAT, USCALE, FSCALE, IER)

     The arguments are:
       UU          = array containing the initial guess on input, and the
                     solution on return
       GLOBALSTRAT = (INTEGER) a number defining the global strategy choice:
                     0 = No globalization, 1 = LineSearch
       USCALE      = array of scaling factors for the UU vector
       FSCALE      = array of scaling factors for the FVAL (function) vector
       IER         = INTEGER error flag as returned by KINSOL:
                     0 means success,
                     1 means initial guess satisfies f(u) = 0 (approx.),
                     2 means apparent stalling (small step),
                     a value < 0 means other error or failure.

     Note: See KINSOL documentation for detailed information.

 (8) Memory freeing: FKINFREE

     To the free the internal memory created by the calls to FKINMALLOC
     and either FNVINITS or FNVINITP, make the following call:

       CALL FKINFREE

 (9) Optional outputs: IOUT/ROUT

     The optional outputs available by way of IOUT and ROUT have the
     following names, locations, and descriptions. For further details see
     the KINSOL documentation.
 
       LENRW  = IOUT(1) = real workspace size
       LENRW  = IOUT(2) = real workspace size
       NNI    = IOuT(3) = number of Newton iterations
       NFE    = IOUT(4) = number of f evaluations
       NBCF   = IOUT(5) = number of line search beta condition failures
       NBKTRK = IOUT(6) = number of line search backtracks

       FNORM  = ROUT(1) = final scaled norm of f(u)
       STEPL  = ROUT(2) = scaled last step length

     The following optional outputs are specific to the SPGMR/SPBCG/SPTFQMR
     module:

       LRW    = IOUT( 7) = real workspace size for the linear solver module
       LIW    = IOUT( 8) = integer workspace size for the linear solver module
       LSTF   = IOUT( 9) = last flag returned by linear solver
       NFE    = IOUT(10) = number of f evaluations for DQ Jacobian
       NJE    = IOUT(11) = number of Jacobian-vector product evaluations
       NPE    = IOUT(12) = number of preconditioner evaluations
       NPS    = IOUT(13) = number of preconditioner solves
       NLI    = IOUT(14) = number of linear (Krylov) iterations
       NCFL   = IOUT(15) = number of linear convergence failures

     The following optional outputs are specific to the DENSE/BAND module:

       LRW    = IOUT( 7) = real workspace size for the linear solver module
       LIW    = IOUT( 8) = integer workspace size for the linear solver module
       LSTF   = IOUT( 9) = last flag returned by linear solver
       NFE    = IOUT(10) = number of f evaluations for DQ Jacobian
       NJE    = IOUT(11) = number of Jacobian evaluations

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

#include "sundials_dense.h"   /* definition of DenseMat      */
#include "sundials_band.h"    /* definition of BandMat      */
#include "sundials_nvector.h" /* definition of type N_Vector */
#include "sundials_types.h"   /* definition of type realtype */

/*
 * -----------------------------------------------------------------
 * generic names are translated through the define statements below
 * -----------------------------------------------------------------
 */

#if defined(F77_FUNC)

#define FKIN_MALLOC         F77_FUNC(fkinmalloc, FKINMALLOC)
#define FKIN_SETIIN         F77_FUNC(fkinsetiin, FKINSETIIN)
#define FKIN_SETRIN         F77_FUNC(fkinsetrin, FKINSETRIN)
#define FKIN_SETVIN         F77_FUNC(fkinsetvin, FKINSETVIN)
#define FKIN_DENSE          F77_FUNC(fkindense, FKINDENSE)
#define FKIN_DENSESETJAC    F77_FUNC(fkindensesetjac, FKINDENSESETJAC)
#define FKIN_BAND           F77_FUNC(fkinband, FKINBAND)
#define FKIN_BANDSETJAC     F77_FUNC(fkinbandsetjac, FKINBANDSETJAC)
#define FKIN_SPTFQMR        F77_FUNC(fkinsptfqmr, FKINSPTFQMR)
#define FKIN_SPTFQMRSETJAC  F77_FUNC(fkinsptfqmrsetjac, FKINSPTFQMRSETJAC)
#define FKIN_SPTFQMRSETPREC F77_FUNC(fkinsptfqmrsetprec, FKINSPTFQMRSETPREC)
#define FKIN_SPBCG          F77_FUNC(fkinspbcg, FKINSPBCG)
#define FKIN_SPBCGSETJAC    F77_FUNC(fkinspbcgsetjac, FKINSPBCGSETJAC)
#define FKIN_SPBCGSETPREC   F77_FUNC(fkinspbcgsetprec, FKINSPBCGSETPREC)
#define FKIN_SPGMR          F77_FUNC(fkinspgmr, FKINSPGMR)
#define FKIN_SPGMRSETJAC    F77_FUNC(fkinspgmrsetjac, FKINSPGMRSETJAC)
#define FKIN_SPGMRSETPREC   F77_FUNC(fkinspgmrsetprec, FKINSPGMRSETPREC)
#define FKIN_SOL            F77_FUNC(fkinsol, FKINSOL)
#define FKIN_FREE           F77_FUNC(fkinfree, FKINFREE)
#define FK_FUN              F77_FUNC(fkfun, FKFUN)
#define FK_PSET             F77_FUNC(fkpset, FKPSET)
#define FK_PSOL             F77_FUNC(fkpsol, FKPSOL)
#define FK_JTIMES           F77_FUNC(fkjtimes, FKJTIMES)
#define FK_DJAC             F77_FUNC(fkdjac, FKDJAC)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC         fkinmalloc
#define FKIN_SETIIN         fkinsetiin
#define FKIN_SETRIN         fkinsetrin
#define FKIN_SETVIN         fkinsetvin
#define FKIN_DENSE          fkindense
#define FKIN_DENSESETJAC    fkindensesetjac
#define FKIN_BAND           fkinband
#define FKIN_BANDSETJAC     fkinbandsetjac
#define FKIN_SPTFQMR        fkinsptfqmr
#define FKIN_SPTFQMRSETJAC  fkinsptfqmrsetjac
#define FKIN_SPTFQMRSETPREC fkinsptfqmrsetprec
#define FKIN_SPBCG          fkinspbcg
#define FKIN_SPBCGSETJAC    fkinspbcgsetjac
#define FKIN_SPBCGSETPREC   fkinspbcgsetprec
#define FKIN_SPGMR          fkinspgmr
#define FKIN_SPGMRSETJAC    fkinspgmrsetjac
#define FKIN_SPGMRSETPREC   fkinspgmrsetprec
#define FKIN_SOL            fkinsol
#define FKIN_FREE           fkinfree
#define FK_FUN              fkfun
#define FK_PSET             fkpset
#define FK_PSOL             fkpsol
#define FK_JTIMES           fkjtimes
#define FK_DJAC             fkdjac

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC         FKINMALLOC
#define FKIN_SETIIN         FKINSETIIN
#define FKIN_SETRIN         FKINSETRIN
#define FKIN_SETVIN         FKINSETVIN
#define FKIN_DENSE          FKINDENSE
#define FKIN_DENSESETJAC    FKINDENSESETJAC
#define FKIN_BAND           FKINBAND
#define FKIN_BANDSETJAC     FKINBANDSETJAC
#define FKIN_SPTFQMR        FKINSPTFQMR
#define FKIN_SPTFQMRSETJAC  FKINSPTFQMRSETJAC
#define FKIN_SPTFQMRSETPREC FKINSPTFQMRSETPREC
#define FKIN_SPBCG          FKINSPBCG
#define FKIN_SPBCGSETJAC    FKINSPBCGSETJAC
#define FKIN_SPBCGSETPREC   FKINSPBCGSETPREC
#define FKIN_SPGMR          FKINSPGMR
#define FKIN_SPGMRSETJAC    FKINSPGMRSETJAC
#define FKIN_SPGMRSETPREC   FKINSPGMRSETPREC
#define FKIN_SOL            FKINSOL
#define FKIN_FREE           FKINFREE
#define FK_FUN              FKFUN
#define FK_PSET             FKPSET
#define FK_PSOL             FKPSOL
#define FK_JTIMES           FKJTIMES
#define FK_DJAC             FKDJAC

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC         fkinmalloc_
#define FKIN_SETIIN         fkinsetiin_
#define FKIN_SETRIN         fkinsetrin_
#define FKIN_SETVIN         fkinsetvin_
#define FKIN_DENSE          fkindense_
#define FKIN_DENSESETJAC    fkindensesetjac_
#define FKIN_BAND           fkinband_
#define FKIN_BANDSETJAC     fkinbandsetjac_
#define FKIN_SPTFQMR        fkinsptfqmr_
#define FKIN_SPTFQMRSETJAC  fkinsptfqmrsetjac_
#define FKIN_SPTFQMRSETPREC fkinsptfqmrsetprec_
#define FKIN_SPBCG          fkinspbcg_
#define FKIN_SPBCGSETJAC    fkinspbcgsetjac_
#define FKIN_SPBCGSETPREC   fkinspbcgsetprec_
#define FKIN_SPGMR          fkinspgmr_
#define FKIN_SPGMRSETJAC    fkinspgmrsetjac_
#define FKIN_SPGMRSETPREC   fkinspgmrsetprec_
#define FKIN_SOL            fkinsol_
#define FKIN_FREE           fkinfree_
#define FK_FUN              fkfun_
#define FK_PSET             fkpset_
#define FK_PSOL             fkpsol_
#define FK_JTIMES           fkjtimes_
#define FK_DJAC             fkdjac_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC         FKINMALLOC_
#define FKIN_SETIIN         FKINSETIIN_
#define FKIN_SETRIN         FKINSETRIN_
#define FKIN_SETVIN         FKINSETVIN_
#define FKIN_DENSE          FKINDENSE_
#define FKIN_DENSESETJAC    FKINDENSESETJAC_
#define FKIN_BAND           FKINBAND_
#define FKIN_BANDSETJAC     FKINBANDSETJAC_
#define FKIN_SPTFQMR        FKINSPTFQMR_
#define FKIN_SPTFQMRSETJAC  FKINSPTFQMRSETJAC_
#define FKIN_SPTFQMRSETPREC FKINSPTFQMRSETPREC_
#define FKIN_SPBCG          FKINSPBCG_
#define FKIN_SPBCGSETJAC    FKINSPBCGSETJAC_
#define FKIN_SPBCGSETPREC   FKINSPBCGSETPREC_
#define FKIN_SPGMR          FKINSPGMR_
#define FKIN_SPGMRSETJAC    FKINSPGMRSETJAC_
#define FKIN_SPGMRSETPREC   FKINSPGMRSETPREC_
#define FKIN_SOL            FKINSOL_
#define FKIN_FREE           FKINFREE_
#define FK_FUN              FKFUN_
#define FK_PSET             FKPSET_
#define FK_PSOL             FKPSOL_
#define FK_JTIMES           FKJTIMES_
#define FK_DJAC             FKDJAC_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FKIN_MALLOC         fkinmalloc__
#define FKIN_SETIIN         fkinsetiin__
#define FKIN_SETRIN         fkinsetrin__
#define FKIN_SETVIN         fkinsetvin__
#define FKIN_DENSE          fkindense__
#define FKIN_DENSESETJAC    fkindensesetjac__
#define FKIN_BAND           fkinband__
#define FKIN_BANDSETJAC     fkinbandsetjac__
#define FKIN_SPTFQMR        fkinsptfqmr__
#define FKIN_SPTFQMRSETJAC  fkinsptfqmrsetjac__
#define FKIN_SPTFQMRSETPREC fkinsptfqmrsetprec__
#define FKIN_SPBCG          fkinspbcg__
#define FKIN_SPBCGSETJAC    fkinspbcgsetjac__
#define FKIN_SPBCGSETPREC   fkinspbcgsetprec__
#define FKIN_SPGMR          fkinspgmr__
#define FKIN_SPGMRSETJAC    fkinspgmrsetjac__
#define FKIN_SPGMRSETPREC   fkinspgmrsetprec__
#define FKIN_SOL            fkinsol__
#define FKIN_FREE           fkinfree__
#define FK_FUN              fkfun__
#define FK_PSET             fkpset__
#define FK_PSOL             fkpsol__
#define FK_JTIMES           fkjtimes__
#define FK_DJAC             fkdjac__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FKIN_MALLOC         FKINMALLOC__
#define FKIN_SETIIN         FKINSETIIN__
#define FKIN_SETRIN         FKINSETRIN__
#define FKIN_SETVIN         FKINSETVIN__
#define FKIN_DENSE          FKINDENSE__
#define FKIN_DENSESETJAC    FKINDENSESETJAC__
#define FKIN_BAND           FKINBAND__
#define FKIN_BANDSETJAC     FKINBANDSETJAC__
#define FKIN_SPTFQMR        FKINSPTFQMR__
#define FKIN_SPTFQMRSETJAC  FKINSPTFQMRSETJAC__
#define FKIN_SPTFQMRSETPREC FKINSPTFQMRSETPREC__
#define FKIN_SPBCG          FKINSPBCG__
#define FKIN_SPBCGSETJAC    FKINSPBCGSETJAC__
#define FKIN_SPBCGSETPREC   FKINSPBCGSETPREC__
#define FKIN_SPGMR          FKINSPGMR__
#define FKIN_SPGMRSETJAC    FKINSPGMRSETJAC__
#define FKIN_SPGMRSETPREC   FKINSPGMRSETPREC__
#define FKIN_SOL            FKINSOL__
#define FKIN_FREE           FKINFREE__
#define FK_FUN              FKFUN__
#define FK_PSET             FKPSET__
#define FK_PSOL             FKPSOL__
#define FK_JTIMES           FKJTIMES__
#define FK_DJAC             FKDJAC__

#endif

/*
 * -----------------------------------------------------------------
 * Prototypes : exported functions
 * -----------------------------------------------------------------
 */

void FKIN_MALLOC(long int *iout, realtype *rout, int *ier);
void FKIN_SETIIN(char key_name[], long int *ival, int *ier, int key_len);
void FKIN_SETRIN(char key_name[], realtype *rval, int *ier, int key_len);
void FKIN_SETVIN(char key_name[], realtype *vval, int *ier, int key_len);
void FKIN_DENSE(long int *neq, int *ier);
void FKIN_DENSESETJAC(int *flag, int *ier);
void FKIN_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
void FKIN_BANDSETJAC(int *flag, int *ier);
void FKIN_SPTFQMR(int *maxl, int *ier);
void FKIN_SPBCG(int *maxl, int *ier);
void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier);
void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier);
void FKIN_FREE(void);
void FKIN_SPTFQMRSETJAC(int *flag, int *ier);
void FKIN_SPTFQMRSETPREC(int *flag, int *ier);
void FKIN_SPBCGSETJAC(int *flag, int *ier);
void FKIN_SPBCGSETPREC(int *flag, int *ier);
void FKIN_SPGMRSETJAC(int *flag, int *ier);
void FKIN_SPGMRSETPREC(int *flag, int *ier);

/*
 * -----------------------------------------------------------------
 * Prototypes : functions called by the solver
 * -----------------------------------------------------------------
 */

int FKINfunc(N_Vector uu, N_Vector fval, void *f_data);

int FKINDenseJac(long int N, DenseMat J, N_Vector uu, N_Vector fval,
                 void *jac_data, N_Vector vtemp1, N_Vector vtemp2);

int FKINBandJac(long int N, long int mupper, long int mlower,
                BandMat J, N_Vector uu, N_Vector fval, void *jac_data,
                N_Vector vtemp1, N_Vector vtemp2);

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

extern N_Vector F2C_KINSOL_vec;
extern void *KIN_kinmem;
extern long int *KIN_iout;
extern realtype *KIN_rout;
extern int KIN_ls;

/* Linear solver IDs */

enum { KIN_LS_SPGMR = 1, KIN_LS_SPBCG = 2, KIN_LS_SPTFQMR = 3, 
       KIN_LS_DENSE = 4, KIN_LS_BAND  = 5 };

#ifdef __cplusplus
}
#endif

#endif
