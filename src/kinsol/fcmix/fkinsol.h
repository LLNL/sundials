/*
 * -----------------------------------------------------------------
 * $Revision: 4905 $
 * $Date: 2016-09-14 16:04:36 -0700 (Wed, 14 Sep 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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
 interface to KINSOL for the serial and parallel NVECTOR
 implementations.

 The user-callable functions, with the corresponding KINSOL functions,
 are as follows:

   FNVINITS, FNVINITP, FNVINITOMP, FNVINITPTS
          initialize serial, distributed memory parallel, or threaded 
          vector computations
   FKINMALLOC interfaces to KINInit
   FKINCREATE interfaces to KINCreate
   FKININIT interfaces to KINInit
   FKINSETIIN, FKINSETRIN, FKINSETVIN interface to KINSet* functions
   FKINDENSE interfaces to KINDense
   FKINKLU interfaces to KINKLU
   FKINSUPERLUMT interfaces to KINSUPERLUMT
   FKINSPARSESETJAC interfaces to KINSlsSetSparseJacFn
   FKINSPTFQMR interfaces to KINSptfqmr
   FKINSPGMR interfaces to KINSpgmr
   FKINSPFGMR interfaces to KINSpfgmr
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
   FKINSPJAC: called by the interface function FKINSparseJac of type 
              KINSlsSparseJacFn
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
 that apply to the serial version of KINSOL only, and end with "p" for
 those that apply to the parallel version only.

 (1) User-supplied system routine: FKFUN

     The user must in all cases supply the following Fortran routine:

       SUBROUTINE FKFUN (UU, FVAL, IER)
       DIMENSION UU(*), FVAL(*)

     It must set the FVAL array to f(u), the system function, as a
     function of the array UU = u. Here UU and FVAL are arrays representing
     vectors, which are distributed vectors in the parallel case.
     IER is a return flag, which should be 0 if FKFUN was successful.
     Return IER > 0 if a recoverable error occurred (and KINSOL is to try
     to recover).  Return IER < 0 if an unrecoverable error occurred.
     

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

  (4.1s) User-supplied sparse Jacobian approximation routine: FKINSPJAC

     Required when using the KINKLU or KINSuperLUMT linear solvers, the 
     user must supply a routine that computes a compressed-sparse-column 
     [or compressed-sparse-row] approximation of the system Jacobian 
     J = dF(y)/dy.  If supplied, it must have the following form:

       SUBROUTINE FKINSPJAC(Y, FY, N, NNZ, JDATA, JRVALS, 
      &                     JCPTRS, WK1, WK2, IER)

     Typically this routine will use only N, NNZ, JDATA, JRVALS and 
     JCPTRS. It must load the N by N compressed sparse column [or compressed 
     sparse row] matrix with storage for NNZ nonzeros, stored in the arrays 
     JDATA (nonzero values), JRVALS (row [or column] indices for each nonzero), 
     JCOLPTRS (indices for start of each column [or row]), with the Jacobian 
     matrix at the current (y) in CSC [or CSR] form (see sundials_sparse.h for 
     more information).

     The arguments are:
         Y    -- array containing state variables [realtype, input]
         FY   -- array containing residual values [realtype, input]
         N    -- number of matrix rows/columns in Jacobian [int, input]
         NNZ  -- allocated length of nonzero storage [int, input]
        JDATA -- nonzero values in Jacobian
                 [realtype of length NNZ, output]
       JRVALS -- row [or column] indices for each nonzero in Jacobian
                  [int of length NNZ, output]
       JCPTRS -- pointers to each Jacobian column [or row] in preceding arrays
                 [int of length N+1, output]
         WK*  -- array containing temporary workspace of same size as Y 
                 [realtype, input]
         IER  -- return flag [int, output]:
                    0 if successful, 
                   >0 if a recoverable error occurred,
                   <0 if an unrecoverable error ocurred.

 (5) Initialization:  FNVINITS/FNVINITP/FNVINITOMP/FNVINITPTS and 
                      FKINCREATE and FKININIT 

 (5.1s) To initialize the serial machine environment, the user must make
        the following call:

          CALL FNVINITS (3, NEQ, IER)

        The arguments are:
          NEQ = size of vectors
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (5.1p) To initialize the distributed memory parallel machine environment, 
        the user must make the following call:

          CALL FNVINITP (3, NLOCAL, NGLOBAL, IER)

        The arguments are:
          NLOCAL  = local size of vectors for this process
          NGLOBAL = the system size, and the global size of vectors
                    (the sum of all values of NLOCAL)
          IER     = return completion flag. Values are 0 = success,
                    -1 = failure.

 (5.1omp) To initialize the openMP threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (3, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (5.1pts) To initialize the Pthreads threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (3, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (5.2) To create the internal memory structure, make the following call:

         CALL FKINCREATE(IER)

       The arguments are:
         IER         = return completion flag. Values are 0 = success, and
                       -1 = failure.

       Note: See printed message for details in case of failure.

 (5.3) To set various integer optional inputs, make the folowing call:

          CALL FKINSETIIN(KEY, VALUE, IER)

       to set the optional input specified by the character key KEY to the 
       integer value VALUE.
       KEY is one of the following: 'PRNT_LEVEL', 'MAX_NITERS', 'ETA_FORM', 'MAA',
       'MAX_SETUPS', 'MAX_SP_SETUPS', 'NO_INIT_SETUP', 'NO_MIN_EPS', 'NO_RES_MON'.

       To set various real optional inputs, make the folowing call:

         CALL FKINSETRIN(KEY, VALUE, IER)

      to set the optional input specified by the character key KEY to the
      real value VALUE.
      KEY is one of the following: 'FNORM_TOL', 'SSTEP_TOL', 'MAX_STEP',
      'RERR_FUNC', 'ETA_CONST', 'ETA_PARAMS', 'RMON_CONST', 'RMON_PARAMS'.
      Note that if KEY is 'ETA_PARAMS' or 'RMON_PARAMS', then VALUE must be an
      array of dimension 2.

      To set the vector of constraints on the solution, make the following call:

        CALL FKINSETVIN(KEY, ARRAY, IER)

      where ARRAY is an array of reals and KEY is 'CONSTR_VEC'.

      FKINSETIIN, FKINSETRIN, and FKINSETVIN return IER=0 if successful and 
      IER<0 if an error occured.

 (5.4) To allocate and initialize the internal memory structure, 
       make the following call:

         CALL FKININIT(IOUT, ROUT, IER)

       The arguments are:
         IOUT        = array of length at least 15 for integer optional outputs
                       (declare as INTEGER*4 or INTEGER*8 according to
                       C type long int)
         ROUT        = array of length at least 2 for real optional outputs
         IER         = return completion flag. Values are 0 = success, and
                       -1 = failure.

       Note: See printed message for details in case of failure.

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


 (6.3s) SPARSE treatment of the linear system using the KLU solver.

     The user must make the call

       CALL FKINKLU(NEQ, NNZ, SPARSETYPE, ORDERING, IER)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	SPARSETYPE = choice between CSC and CSR format
           (0 = CSC, 1 = CSR) [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the KLU package (0 = AMD, 1 = COLAMD) [int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     When using the KLU solver the user must provide the FKINSPJAC routine for the 
     evalution of the sparse approximation to the Jacobian. To indicate that this
     routine has been provided, after the call to FKINKLU, the following call must 
     be made    

       CALL FKINSPARSESETJAC(IER) 

     The int return flag IER=0 if successful, and nonzero otherwise.


     The KINSOL KLU solver will reuse much of the factorization information from one
     nonlinear iteration to the next.  If at any time the user wants to force a full
     refactorization or if the number of nonzeros in the Jacobian matrix changes, the
     user should make the call

       CALL FKINKLUREINIT(NEQ, NNZ, REINIT_TYPE)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be destroyed and 
          a new one will be allocated with NNZ nonzeros.  For a value of 2, 
	  only symbolic and numeric factorizations will be completed. 
 
     When using FKINKLU, the user is required to supply the FKINSPJAC 
     routine for the evaluation of the sparse approximation to the 
     Jacobian, as discussed above with the other user-supplied routines.
 
     Optional outputs specific to the KLU case are:
        LSTF    = IOUT(8)  from KINSlsGetLastFlag
        NJES    = IOUT(10) from KINSlsGetNumJacEvals
     See the KINSOL manual for descriptions.
 
 (6.4s) SPARSE treatment of the linear system using the SuperLUMT solver.

     The user must make the call

       CALL FKINSUPERLUMT(NTHREADS, NEQ, NNZ, ORDERING, IER)

     The arguments are:
        NTHREADS = desired number of threads to use [int; input]
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the SuperLU_MT package [int; input]
           0 = Natural
           1 = Minimum degree on A^T A
           2 = Minimum degree on A^T + A
           3 = COLAMD
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     At this time, there is no reinitialization capability for the SUNDIALS 
     interfaces to the SuperLUMT solver.

     When using FKINSUPERLUMT, the user is required to supply the FKINSPJAC 
     routine for the evaluation of the CSC approximation to the 
     Jacobian (note: the current SuperLU_MT interface in SUNDIALS does not 
     support CSR matrices). To indicate that this routine has been provided, 
     after the call to FKINSUPERLUMT, the following call must be made    

         CALL FKINSPARSESETJAC(IER) 

     The int return flag IER=0 if successful, and nonzero otherwise.
 
     Optional outputs specific to the SUPERLUMT case are:
        LSTF    = IOUT(8)  from KINSlsGetLastFlag
        NJES    = IOUT(10) from KINSlsGetNumJacEvals
     See the KINSOL manual for descriptions.
  
 (6.5) SPTFQMR treatment of the linear systems:

       For the Scaled Preconditioned TFQMR solution of the linear systems,
       the user must make the call:

         CALL FKINSPTFQMR(MAXL, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

 (6.6) SPBCG treatment of the linear systems:

       For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
       the user must make the call:

         CALL FKINSPBCG(MAXL, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

 (6.7) SPGMR and SPFGMR treatment of the linear systems:

       For the Scaled Preconditioned GMRES or Scaled Preconditioned Flexible 
       GMRES solution of the linear systems, the user must make one of the calls:

         CALL FKINSPGMR(MAXL, MAXLRST, IER)
         CALL FKINSPFGMR(MAXL, MAXLRST, IER)

       In the above routine, the arguments are as follows:
         MAXL     = maximum Krylov subspace dimension; 0 indicates default.
         MAXLRST  = maximum number of linear system restarts; 0 indicates
                    default (SPGMR and SPFGMR only).
         IER      = return completion flag.  Values are 0 = succes, and
                    -1 = failure.

       Note: See printed message for details in case of failure.

 (6.8) Specifying user-provided functions for the iterative linear solvers

       If the user program includes the FKJTIMES routine for the evaluation
       of the Jacobian-vector product, the following call must be made:

         CALL FKINSPILSSETJAC(FLAG, IER)

       The argument FLAG = 0 specifies using the internal finite differences
       approximation to the Jacobian-vector product, while FLAG = 1 specifies
       that FKJTIMES is provided.

       Usage of the user-supplied routines FKPSET and FKPSOL for the setup and
       solution of the preconditioned linear system is specified by calling:

         CALL FKINSPILSSETPREC(FLAG, IER)

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
                     0 = No globalization, 1 = LineSearch, 2 = Picard, 
                     3 = Fixed Point 
       USCALE      = array of scaling factors for the UU vector
       FSCALE      = array of scaling factors for the FVAL (function) vector
       IER         = INTEGER error flag as returned by KINSOL:
                     0 means success,
                     1 means initial guess satisfies f(u) = 0 (approx.),
                     2 means apparent stalling (small step),
                     a value < 0 means other error or failure.

     Note: See KINSOL documentation for detailed information.

 (8) Memory freeing: FKINFREE

     To the free the internal memory created by the calls to FKINCREATE and 
     FKININIT and any FNVINIT**, make the following call:

       CALL FKINFREE

 (9) Optional outputs: IOUT/ROUT

     The optional outputs available by way of IOUT and ROUT have the
     following names, locations, and descriptions. For further details see
     the KINSOL documentation.
 
       LENRW  = IOUT(1) = real workspace size
       LENRW  = IOUT(2) = real workspace size
       NNI    = IOUT(3) = number of Newton iterations
       NFE    = IOUT(4) = number of f evaluations
       NBCF   = IOUT(5) = number of line search beta condition failures
       NBKTRK = IOUT(6) = number of line search backtracks

       FNORM  = ROUT(1) = final scaled norm of f(u)
       STEPL  = ROUT(2) = scaled last step length

     The following optional outputs are specific to the SPGMR/SPFGMR/SPBCG/SPTFQMR
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

/*
 * -----------------------------------------------------------------
 * header files
 * -----------------------------------------------------------------
 */
#include <kinsol/kinsol.h>
#include <sundials/sundials_direct.h>  /* definition of type DlsMat   */
#include <sundials/sundials_sparse.h>  /* definition of type SlsMat   */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * generic names are translated through the define statements below
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_F77_FUNC)

#define FKIN_MALLOC         SUNDIALS_F77_FUNC(fkinmalloc, FKINMALLOC)
#define FKIN_CREATE         SUNDIALS_F77_FUNC(fkincreate, FKINCREATE)
#define FKIN_INIT           SUNDIALS_F77_FUNC(fkininit,   FKININIT)
#define FKIN_SETIIN         SUNDIALS_F77_FUNC(fkinsetiin, FKINSETIIN)
#define FKIN_SETRIN         SUNDIALS_F77_FUNC(fkinsetrin, FKINSETRIN)
#define FKIN_SETVIN         SUNDIALS_F77_FUNC(fkinsetvin, FKINSETVIN)
#define FKIN_DENSE          SUNDIALS_F77_FUNC(fkindense, FKINDENSE)
#define FKIN_DENSESETJAC    SUNDIALS_F77_FUNC(fkindensesetjac, FKINDENSESETJAC)
#define FKIN_BAND           SUNDIALS_F77_FUNC(fkinband, FKINBAND)
#define FKIN_BANDSETJAC     SUNDIALS_F77_FUNC(fkinbandsetjac, FKINBANDSETJAC)
#define FKIN_LAPACKDENSE       SUNDIALS_F77_FUNC(fkinlapackdense, FKINLAPACKDENSE)
#define FKIN_LAPACKDENSESETJAC SUNDIALS_F77_FUNC(fkinlapackdensesetjac, FKINLAPACKDENSESETJAC)
#define FKIN_LAPACKBAND        SUNDIALS_F77_FUNC(fkinlapackband, FKINLAPACKBAND)
#define FKIN_LAPACKBANDSETJAC  SUNDIALS_F77_FUNC(fkinlapackbandsetjac, FKINLAPACKBANDSETJAC)
#define FKIN_KLU            SUNDIALS_F77_FUNC(fkinklu, FKLUKLU)
#define FKIN_KLUREINIT      SUNDIALS_F77_FUNC(fkinklureinit, FKLUKLUREINIT)
#define FKIN_SUPERLUMT      SUNDIALS_F77_FUNC(fkinsuperlumt, FKLUSUPERLUMT)
#define FKIN_SPARSESETJAC   SUNDIALS_F77_FUNC(fkinsparsesetjac, FKINSPARSESETJAC)  
#define FKIN_SPTFQMR        SUNDIALS_F77_FUNC(fkinsptfqmr, FKINSPTFQMR)
#define FKIN_SPBCG          SUNDIALS_F77_FUNC(fkinspbcg, FKINSPBCG)
#define FKIN_SPGMR          SUNDIALS_F77_FUNC(fkinspgmr, FKINSPGMR)
#define FKIN_SPFGMR         SUNDIALS_F77_FUNC(fkinspfgmr, FKINSPFGMR)
#define FKIN_SPILSSETJAC    SUNDIALS_F77_FUNC(fkinspilssetjac, FKINSPILSSETJAC)
#define FKIN_SPILSSETPREC   SUNDIALS_F77_FUNC(fkinspilssetprec, FKINSPILSSETPREC)
#define FKIN_SOL            SUNDIALS_F77_FUNC(fkinsol, FKINSOL)
#define FKIN_FREE           SUNDIALS_F77_FUNC(fkinfree, FKINFREE)
#define FK_FUN              SUNDIALS_F77_FUNC(fkfun, FKFUN)
#define FK_PSET             SUNDIALS_F77_FUNC(fkpset, FKPSET)
#define FK_PSOL             SUNDIALS_F77_FUNC(fkpsol, FKPSOL)
#define FK_JTIMES           SUNDIALS_F77_FUNC(fkjtimes, FKJTIMES)
#define FK_DJAC             SUNDIALS_F77_FUNC(fkdjac, FKDJAC)
#define FK_BJAC             SUNDIALS_F77_FUNC(fkbjac, FKBJAC)

#else

#define FKIN_MALLOC         fkinmalloc_
#define FKIN_CREATE         fkincreate_
#define FKIN_INIT           fkininit_
#define FKIN_SETIIN         fkinsetiin_
#define FKIN_SETRIN         fkinsetrin_
#define FKIN_SETVIN         fkinsetvin_
#define FKIN_DENSE          fkindense_
#define FKIN_DENSESETJAC    fkindensesetjac_
#define FKIN_BAND           fkinband_
#define FKIN_BANDSETJAC     fkinbandsetjac_
#define FKIN_LAPACKDENSE       fkinlapackdense_
#define FKIN_LAPACKDENSESETJAC fkinlapackdensesetjac_
#define FKIN_LAPACKBAND        fkinlapackband_
#define FKIN_LAPACKBANDSETJAC  fkinlapackbandsetjac_
#define FKIN_KLU            fkinklu_
#define FKIN_KLUREINIT      fkinklureinit_
#define FKIN_SUPERLUMT      fkinsuperlumt_
#define FKIN_SPARSESETJAC   fkinsparsesetjac_
#define FKIN_SPTFQMR        fkinsptfqmr_
#define FKIN_SPBCG          fkinspbcg_
#define FKIN_SPGMR          fkinspgmr_
#define FKIN_SPFGMR         fkinspgmr_
#define FKIN_SPILSSETJAC    fkinspilssetjac_
#define FKIN_SPILSSETPREC   fkinspilssetprec_
#define FKIN_SOL            fkinsol_
#define FKIN_FREE           fkinfree_
#define FK_FUN              fkfun_
#define FK_PSET             fkpset_
#define FK_PSOL             fkpsol_
#define FK_JTIMES           fkjtimes_
#define FK_DJAC             fkdjac_
#define FK_BJAC             fkbjac_

#endif

/*
 * -----------------------------------------------------------------
 * Prototypes : exported functions
 * -----------------------------------------------------------------
 */

void FKIN_MALLOC(long int *iout, realtype *rout, int *ier);
void FKIN_CREATE(int *ier);
void FKIN_INIT(long int *iout, realtype *rout, int *ier);

void FKIN_SETIIN(char key_name[], long int *ival, int *ier);
void FKIN_SETRIN(char key_name[], realtype *rval, int *ier);
void FKIN_SETVIN(char key_name[], realtype *vval, int *ier);

void FKIN_DENSE(long int *neq, int *ier);
void FKIN_DENSESETJAC(int *flag, int *ier);

void FKIN_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
void FKIN_BANDSETJAC(int *flag, int *ier);

void FKIN_LAPACKDENSE(int *neq, int *ier);
void FKIN_LAPACKDENSESETJAC(int *flag, int *ier);
void FKIN_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
void FKIN_LAPACKBANDSETJAC(int *flag, int *ier);

void FKIN_KLU(int *neq, int *nnz, int *sparsetype, int *ordering, int *ier);
void FKIN_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier);
void FKIN_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier);
void FKIN_SPARSESETJAC(int *ier);

void FKIN_SPTFQMR(int *maxl, int *ier);
void FKIN_SPBCG(int *maxl, int *ier);
void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier);
void FKIN_SPFGMR(int *maxl, int *maxlrst, int *ier);

void FKIN_SPILSSETJAC(int *flag, int *ier);
void FKIN_SPILSSETPREC(int *flag, int *ier);

void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier);

void FKIN_FREE(void);

/*
 * -----------------------------------------------------------------
 * Prototypes : functions called by the solver
 * -----------------------------------------------------------------
 */

int FKINfunc(N_Vector uu, N_Vector fval, void *user_data);

int FKINDenseJac(long int N,
                 N_Vector uu, N_Vector fval,
                 DlsMat J, void *user_data, 
                 N_Vector vtemp1, N_Vector vtemp2);

int FKINBandJac(long int N, long int mupper, long int mlower,
                N_Vector uu, N_Vector fval, 
                DlsMat J, void *user_data,
                N_Vector vtemp1, N_Vector vtemp2);

int FKINLapackDenseJac(long int N,
                       N_Vector uu, N_Vector fval,
                       DlsMat J, void *user_data, 
                       N_Vector vtemp1, N_Vector vtemp2);

int FKINLapackBandJac(long int N, long int mupper, long int mlower,
                      N_Vector uu, N_Vector fval, 
                      DlsMat J, void *user_data,
                      N_Vector vtemp1, N_Vector vtemp2);

int FKINSparseJac(N_Vector y, N_Vector fval, SlsMat J,
		  void *user_data, N_Vector vtemp1, N_Vector vtemp2);

int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *user_data,
             N_Vector vtemp1, N_Vector vtemp2);

int FKINPSol(N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale, 
             N_Vector vv, void *user_data,
             N_Vector vtemp);

int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *user_data);

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

  enum { KIN_LS_SPGMR = 1, KIN_LS_SPFGMR = 2, KIN_LS_SPBCG = 3, KIN_LS_SPTFQMR = 4, 
	 KIN_LS_DENSE = 5, KIN_LS_BAND  = 6,
	 KIN_LS_LAPACKDENSE = 7, KIN_LS_LAPACKBAND = 8,
         KIN_LS_KLU = 9, KIN_LS_SUPERLUMT = 10 };

#ifdef __cplusplus
}
#endif

#endif
