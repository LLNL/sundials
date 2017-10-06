/*
 * -----------------------------------------------------------------
 * $Revision: 4795 $
 * $Date: 2016-07-01 15:10:25 -0700 (Fri, 01 Jul 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
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
 * This is the header file for FIDA, the Fortran interface to
 * the IDA package.
 * -----------------------------------------------------------------
 */

/*
 * =============================================================================
 *
 *                  FIDA Interface Package
 *
 * The FIDA Interface Package is a package of C functions which support
 * the use of the IDA solver, for the solution of DAE systems, in a
 * mixed Fortran/C setting.  While IDA is written in C, it is assumed
 * here that the user's calling program and user-supplied problem-defining
 * routines are written in Fortran.  This package provides the necessary
 * interface to IDA for both the serial and the parallel NVECTOR
 * implementations.
 *
 * The user-callable functions, with the corresponding IDA functions,
 * are as follows:
 *
 *   FNVINITS* and FNVINITP*  interface to N_VNew_Serial and
 *                            N_VNew_Parallel, respectively
 *   FNVINITOMP               N_VNew_OpenMP
 *   FNVINITPTS               N_VNew_Pthreads
 *
 *   FIDAMALLOC  interfaces to IDACreate and IDAInit
 *
 *   FIDAREINIT  interfaces to IDAReInit
 *
 *   FIDASETIIN, FIDASETRIN, FIDASETVIN interface to IDASet*
 *
 *   FIDATOLREINIT  interfaces to IDASetTolerances
 *
 *   FIDACALCIC  interfaces to IDACalcIC
 *
 *   FIDAEWTSET  interfaces to IDAWFtolerances
 *
 *   FIDADENSE        interfaces to IDADense
 *   FIDADENSESETJAC  interfaces to IDADenseSetJacFn
 *
 *   FIDABAND        interfaces to IDABand
 *   FIDABANDSETJAC  interfaces to IDABandSetJacFn
 *
 *   FIDAKLU         interfaces to IDAKLU
 *   FIDAKLUReinit   interfaces to IDAKLUReinit
 *   FIDASUPERLUMT   interfaces to IDASuperLUMT
 *   FIDASPARSESETJAC    interfaces to IDASlsSetSparseJacFn
 *
 *   FIDASPTFQMR/FIDASPTFQMRREINIT  interface to IDASptfqmr and IDASptfqmrSet*
 *   FIDASPBCG/FIDASPBCGREINIT      interface to IDASpbcg and IDASpbcgSet*
 *   FIDASPGMR/FIDASPGMRREINIT      interface to IDASpgmr and IDASpgmrSet*
 *   FIDASPILSSETJAC                interfaces to IDASpilsSetJacFn
 *   FIDASPILSSETPREC               interfaces to IDASpilsSetPreconditioner
 *
 *   FIDASOLVE  interfaces to IDASolve, IDAGet*, and IDA*Get*
 *
 *   FIDAGETDKY  interfaces to IDAGetDky
 *
 *   FIDAGETERRWEIGHTS  interfaces to IDAGetErrWeights
 *
 *   FIDAGETESTLOCALERR  interfaces to IDAGetEstLocalErrors
 *
 *   FIDAFREE  interfaces to IDAFree
 *
 * The user-supplied functions, each listed with the corresponding interface
 * function which calls it (and its type within IDA), are as follows:
 *   FIDARESFUN is called by the interface function FIDAresfn of type IDAResFn
 *   FIDADJAC   is called by the interface fn. FIDADenseJac of type IDADenseJacFn
 *   FIDABJAC   is called by the interface fn. FIDABandJac of type IDABandJacFn
 *   FIDAPSOL   is called by the interface fn. FIDAPSol of type IDASpilsPrecSolveFn
 *   FIDAPSET   is called by the interface fn. FIDAPSet of type IDASpilsPrecSetupFn
 *   FIDAJTIMES is called by interface fn. FIDAJtimes of type IDASpilsJacTimesVecFn
 *   FIDASPJAC  is called by interface fn. FIDASparseJac of type IDASlsSparseJacFn
 *   FIDAEWT    is called by interface fn. FIDAEwtSet of type IDAEwtFn
 * In contrast to the case of direct use of IDA, the names of all user-supplied
 * routines here are fixed, in order to maximize portability for the resulting
 * mixed-language program.
 *
 * Important note on portability:
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions in the
 * header file fida.h.
 *
 * =============================================================================
 *
 *                  Usage of the FIDA Interface Package
 *
 * The usage of FIDA requires calls to a few different interface
 * functions, depending on the method options selected, and one or more
 * user-supplied routines which define the problem to be solved.  These
 * function calls and user routines are summarized separately below.
 *
 * Some details are omitted, and the user is referred to the user documents
 * on IDA for more complete documentation.  Information on the
 * arguments of any given user-callable interface routine, or of a given
 * user-supplied function called by an interface function, can be found in
 * the documentation on the corresponding function in the IDA package.
 *
 * The number labels on the instructions below end with s for instructions
 * that apply to the serial version of IDA only, and end with p for
 * those that apply to the parallel version only.
 *
 * -----------------------------------------------------------------------------
 *
 * (1) User-supplied residual routine: FIDARESFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FIDARESFUN(T, Y, YP, R, IPAR, RPAR, IER)
 *       DIMENSION Y(*), YP(*), R(*), IPAR(*), RPAR(*)
 * It must set the R array to F(t,y,y'), the residual of the DAE 
 * system, as a function of T = t, the array Y = y, and the array YP = y'.
 * Here Y, YP and R are distributed vectors. 
 * IPAR and RPAR are arrays of integer and real user data, respectively,
 * as passed to FIDAMALLOC.
 *
 * (2s) Optional user-supplied dense Jacobian approximation routine: FIDADJAC
 * As an option when using the DENSE linear solver, the user may supply a
 * routine that computes a dense approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FIDADJAC(NEQ, T, Y, YP, R, DJAC, CJ, EWT, H,
 *      1                    IPAR, RPAR, WK1, WK2, WK3, IER)
 *       DIMENSION Y(*), YP(*), R(*), EWT(*), DJAC(NEQ,*),
 *      1          IPAR(*), RPAR(*), WK1(*), WK2(*), WK3(*)
 * This routine must compute the Jacobian and store it columnwise in DJAC.
 * IPAR and RPAR are user (integer and real) arrays passed to FIDAMALLOC.
 *
 * (3s) Optional user-supplied band Jacobian approximation routine: FIDABJAC
 * As an option when using the BAND linear solver, the user may supply a
 * routine that computes a band approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FIDABJAC(NEQ, MU, ML, MDIM, T, Y, YP, R, CJ,
 *      1                    BJAC, EWT, H, IPAR, RPAR, WK1, WK2, WK3, IER)
 *       DIMENSION Y(*), YP(*), R(*), EWT(*), BJAC(MDIM,*),
 *      1          IPAR(*), RPAR(*), WK1(*), WK2(*), WK3(*)
 * This routine must load the MDIM by N array BJAC with the Jacobian matrix at the
 * current (t,y,y') in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
 * with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.
 * IPAR and RPAR are user (integer and real) arrays passed to FIDAMALLOC.
 *
 * (4) Optional user-supplied Jacobian-vector product routine: FIDAJTIMES
 * As an option when using the SPGMR/SPBCG/SPTFQMR linear solver, the user may
 * supply a routine that computes the product of the system Jacobian J = df/dy
 * and a given vector v.  If supplied, it must have the following form:
 *       SUBROUTINE FIDAJTIMES(T, Y, YP, R, V, FJV, CJ, EWT, H, 
 *      1                      IPAR, RPAR, WK1, WK2, IER)
 *       DIMENSION V(*), FJV(*), Y(*), YP(*), R(*), EWT(*),
 *      1          IPAR(*), RPAR(*), WK1(*), WK2(*)
 * This routine must compute the product vector Jv, where the vector v is stored
 * in V, and store the product in FJV.  On return, set IER = 0 if FIDAJTIMES was
 * successful, and nonzero otherwise.
 * IPAR and RPAR are user (integer and real) arrays passed to FIDAMALLOC.
 *
 * (5s) User-supplied sparse Jacobian approximation routine: FIDASPJAC
 *
 * Required when using the IDAKLU or IDASuperLUMT linear solvers, the 
 * user must supply a routine that computes a compressed-sparse-column [or 
 * compressed-sparse-row] approximation of the system Jacobian 
 * J = dF/dy' + c_j*dF/dy.  If supplied, it must have the following form:
 *
 *       SUBROUTINE FIDASPJAC(T, CJ, Y, YP, R, N, NNZ, JDATA, JRVALS, 
 *      &                    JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)
 *
 * It must load the N by N compressed sparse column [row] matrix 
 * with storage for NNZ nonzeros, stored in the arrays JDATA (nonzero
 * values), JRVALS (row [column] indices for each nonzero), JCOLPTRS (indices 
 * for start of each column [row]), with the Jacobian matrix at the current
 * (t,y) in CSC [CSR] form (see sundials_sparse.h for more information).
 *
 * The arguments are:
 *         T    -- current time [realtype, input]
 *         CJ   -- Scalar in the system Jacobian proportional 
 *                 to inverse step size [realtype, input]
 *         Y    -- array containing state variables [realtype, input]
 *         YP   -- array containing state derivatives [realtype, input]
 *         R    -- array containing system residual F(T, Y, YP) [realtype, input]
 *         N    -- number of matrix rows/columns in Jacobian [int, input]
 *         NNZ  -- allocated length of nonzero storage [int, input]
 *        JDATA -- nonzero values in Jacobian
 *                 [realtype of length NNZ, output]
 *       JRVALS -- row [column] indices for each nonzero in Jacobian
 *                  [int of length NNZ, output]
 *       JCPTRS -- pointers to each Jacobian column [row] in preceding arrays
 *                 [int of length N+1, output]
 *         H    -- current step size [realtype, input]
 *         IPAR -- array containing integer user data that was passed to
 *                 FIDAMALLOC [long int, input]
 *         RPAR -- array containing real user data that was passed to
 *                 FIDAMALLOC [realtype, input]
 *         WK*  -- array containing temporary workspace of same size as Y 
 *                 [realtype, input]
 *         IER  -- return flag [int, output]:
 *                    0 if successful, 
 *                   >0 if a recoverable error occurred,
 *                   <0 if an unrecoverable error ocurred.
 *
 * (6) Optional user-supplied error weight vector routine: FIDAEWT
 * As an option to providing the relative and absolute tolerances, the user
 * may supply a routine that computes the weights used in the WRMS norms.
 * If supplied, it must have the following form:
 *       SUBROUTINE FIDAEWT(Y, EWT, IPAR, RPAR, IER)
 *       DIMENSION Y(*), EWT(*)
 * It must store the error weights in EWT, given the current solution vector Y.
 * On return, set IER = 0 if successful, and nonzero otherwise.
 * IPAR and RPAR are user (integer and real) arrays passed to FIDAMALLOC.
 *
 * -----------------------------------------------------------------------------
 *
 * (7) Initialization:  FNVINITS/FNVINITP/FNVINITOMP/FNVINITPTS, 
 *                      FIDAMALLOC, FIDAREINIT,
 *                      FIDATOLREINIT, and FIDACALCIC
 *
 * (7.1s) To initialize the serial machine environment, the user must make
 * the following call:
 *        CALL FNVINITS(KEY, NEQ, IER)
 * The arguments are:
 * KEY = 2 for IDA
 * NEQ = size of vectors
 * IER = return completion flag. Values are 0 = success, -1 = failure.
 *
 * (7.1p) To initialize the distributed memory parallel machine environment, 
 * the user must make one of the following calls:
 *        CALL FNVINITP(KEY, NLOCAL, NGLOBAL, IER)
 *                     -or-
 *        CALL FNVINITP(COMM, KEY, NLOCAL, NGLOBAL, IER)
 * The arguments are:
 * COMM = MPI communicator (e.g., MPI_COMM_WORLD)
 * KEY = 2 for IDA
 * NLOCAL  = local size of vectors on this processor
 * NGLOBAL = the system size, and the global size of vectors (the sum 
 *           of all values of NLOCAL)
 * IER     = return completion flag. Values are 0 = success, -1 = failure.
 * NOTE: The COMM argument passed to the FNVINITP routine is only supported if
 * the MPI implementation used to build SUNDIALS includes the MPI_Comm_f2c
 * function from the MPI-2 specification.  To check if the function is supported
 * look for the line "#define SUNDIALS_MPI_COMM_F2C 1" in the sundials_config.h
 * header file.
 *
 (7.1omp) To initialize the openMP threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (KEY, NEQ, NUM_THREADS, IER)

        The arguments are:
	  KEY = 2 for IDA
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (7.1pts) To initialize the Pthreads threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (KEY, NEQ, NUM_THREADS, IER)

        The arguments are:
	  KEY = 2 for IDA
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.
 *
 * (7.2) To set various problem and solution parameters and allocate
 * internal memory, make the following call:
 *       CALL FIDAMALLOC(T0, Y0, YP0, IATOL, RTOL, ATOL, 
 *      1                IOUT, ROUT, IPAR, RPAR, IER)
 * The arguments are:
 * T0    = initial value of t
 * Y0    = array of initial conditions, y(t0)
 * YP0   = value of y'(t0)
 * IATOL = type for absolute tolerance ATOL: 1 = scalar, 2 = array.
 *         If IATOL = 3, then the user must supply a routine FIDAEWT to compute
 *         the error weight vector.
 * RTOL  = relative tolerance (scalar)
 * ATOL  = absolute tolerance (scalar or array)
 * IOUT  = array of length at least 21 for integer optional outputs
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * ROUT  = array of length at least 6 for real optional outputs
 * IPAR  = array with user integer data
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * RPAR  = array with user real data
 * IER   = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
 *         See printed message for details in case of failure.
 *
 * The user data arrays IPAR and RPAR are passed unmodified to all subsequent
 * calls to user-provided routines. Modifications to either array inside a
 * user-provided routine will be propagated. Using these two arrays, the user
 * can dispense with Common blocks to pass data betwen user-provided routines.
 * 
 * The optional outputs are:
 *           LENRW   = IOUT( 1) -> IDAGetWorkSpace
 *           LENIW   = IOUT( 2) -> IDAGetWorkSpace
 *           NST     = IOUT( 3) -> IDAGetNumSteps
 *           NRE     = IOUT( 4) -> IDAGetNumResEvals
 *           NETF    = IOUT( 5) -> IDAGetNumErrTestFails
 *           NCFN    = IOUT( 6) -> IDAGetNumNonlinSolvConvFails
 *           NNI     = IOUT( 7) -> IDAGetNumNonlinSolvIters
 *           NSETUPS = IOUT( 8) -> IDAGetNumLinSolvSetups
 *           KLAST   = IOUT( 9) -> IDAGetLastOrder
 *           KCUR    = IOUT(10) -> IDAGetCurrentOrder
 *           NBCKTRK = IOUT(11) -> IDAGetNumBacktrackOps
 *           NGE     = IOUT(12) -> IDAGetNumGEvals
 *
 *           HINUSED = ROUT( 1) -> IDAGetActualInitStep
 *           HLAST   = ROUT( 2) -> IDAGetLastStep
 *           HCUR    = ROUT( 3) -> IDAGetCurrentStep
 *           TCUR    = ROUT( 4) -> IDAGetCurrentTime
 *           TOLSFAC = ROUT( 5) -> IDAGetTolScaleFactor
 *           UNITRND = ROUT( 6) -> UNIT_ROUNDOFF
 *
 *
 * If the user program includes the FIDAEWT routine for the evaluation of the 
 * error weights, the following call must be made
 *       CALL FIDAEWTSET(FLAG, IER)
 * with FLAG = 1 to specify that FIDAEWT is provided.
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * (7.3) To set various integer optional inputs, make the folowing call:
 *       CALL FIDASETIIN(KEY, VALUE, IER)
 * to set the optional input specified by the character key KEY to the 
 * integer value VAL.
 * KEY is one of the following: MAX_ORD, MAX_NSTEPS, MAX_ERRFAIL, MAX_NITERS, 
 * MAX_CONVFAIL, SUPPRESS_ALG, MAX_NSTEPS_IC, MAX_NITERS_IC, MAX_NJE_IC, LS_OFF_IC.
 *
 * To set various real optional inputs, make the folowing call:
 *       CALL FIDASETRIN(KEY, VALUE, IER)
 * to set the optional input specified by the character key KEY to the
 * real value VAL.
 * KEY is one of the following: INIT_STEP, MAX_STEP, MIIN_STEP, STOP_TIME,
 * NLCONV_COEF.
 *
 * To set the vector of variable IDs or the vector of constraints, make
 * the following call:
 *       CALL FIDASETVIN(KEY, ARRAY, IER)
 * where ARRAY is an array of reals and KEY is 'ID_VEC' or 'CONSTR_VEC'.
 *
 * FIDASETIIN, FIDASETRIN, and FIDASETVIN return IER=0 if successful and 
 * IER<0 if an error occured.
 *
 * (7.4) To re-initialize the FIDA solver for the solution of a new problem
 * of the same size as one already solved, make the following call:
 *       CALL FIDAREINIT(T0, Y0, YP0, IATOL, RTOL, ATOL, ID, CONSTR, IER)
 * The arguments have the same names and meanings as those of FIDAMALLOC.
 * FIDAREINIT performs the same initializations as FIDAMALLOC, but does no memory 
 * allocation for IDA data structures, using instead the existing internal memory
 * created by the previous FIDAMALLOC call.  The call to specify the linear system
 * solution method may or may not be needed.  See below.
 *
 * (7.5) To modify the tolerance parameters, make the following call:
 *       CALL FIDATOLREINIT(IATOL, RTOL, ATOL, IER)
 * The arguments have the same names and meanings as those of FIDAMALLOC.
 * FIDATOLREINIT simple calls IDASetTolerances with the given arguments.
 *
 * (7.6) To compute consistent initial conditions for an index-one DAE system,
 * make the following call:
 *       CALL FIDACALCIC(ICOPT, TOUT, IER)
 * The arguments are:
 * ICOPT = specifies the option: 1 = IDA_YP_YDP_INIT, 2 = IDA_Y_INIT.
 *         (See user guide for additional details.)
 * TOUT  = the first value of t at which a solution will be requested
 *         (from FIDASOLVE).
 * IER   = return completion flag.
 * 
 * -----------------------------------------------------------------------------
 *
 * (8) Specification of linear system solution method.
 * FIDA presently includes four choices for the treatment of these systems,
 * and the user of FIDA must call a routine with a specific name to make the
 * desired choice.
 * 
 * (8.1s) DENSE treatment of the linear system.
 * The user must make the call
 *       CALL FIDADENSE(NEQ, IER)
 * The arguments are:
 * NEQ = size of vectors
 * IER = error return flag: 0 = success , negative value = an error occured
 * 
 * If the user program includes the FIDADJAC routine for the evaluation of the 
 * dense approximation to the Jacobian, the following call must be made
 *       CALL FIDADENSESETJAC(FLAG, IER)
 * with FLAG = 1 to specify that FIDADJAC is provided.  (FLAG = 0 specifies
 * using the internal finite differences approximation to the Jacobian.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 * 
 * Optional outputs specific to the DENSE case are:
 *        LENRWLS = IOUT(13) -> IDADenseGetWorkSpace
 *        LENIWLS = IOUT(14) -> IDADenseGetWorkSpace
 *        LSTF    = IOUT(15) -> IDADenseGetLastFlag
 *        NRELS   = IOUT(16) -> IDADenseGetNumResEvals
 *        NJE     = IOUT(17) -> IDADenseGetNumJacEvals
 *
 * (8.2s) BAND treatment of the linear system
 * The user must make the call
 *       CALL FIDABAND(NEQ, MU, ML, IER)
 * The arguments are:
 * NEQ = size of vectors
 * MU  = upper bandwidth
 * ML  = lower bandwidth
 * IER = error return flag: 0 = success , negative value = an error occured
 * 
 * If the user program includes the FIDABJAC routine for the evaluation of the 
 * band approximation to the Jacobian, the following call must be made
 *       CALL FIDABANDSETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDABJAC is provided.  (FLAG = 0 specifies
 * using the internal finite differences approximation to the Jacobian.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * Optional outputs specific to the BAND case are:
 *        LENRWLS = IOUT(13) -> IDABandGetWorkSpace
 *        LENIWLS = IOUT(14) -> IDABandGetWorkSpace
 *        LSTF    = IOUT(15) -> IDABandGetLastFlag
 *        NRELS   = IOUT(16) -> IDABandGetNumResEvals
 *        NJE     = IOUT(17) -> IDABandGetNumJacEvals
 *
 *  (8.3s) SPARSE treatment of the linear system using the KLU solver.
 *
 *     The user must make the call
 *
 *       CALL FIDAKLU(NEQ, NNZ, SPARSETYPE, ORDERING, IER)
 *
 *     The arguments are:
 *        NEQ = the problem size [int; input]
 *        NNZ = the maximum number of nonzeros [int; input]
 *        SPARSETYPE = choice between CSC and CSR format
 *           (0 = CSC, 1 = CSR) [int; input]
 *        ORDERING = the matrix ordering desired, possible values
 *           come from the KLU package (0 = AMD, 1 = COLAMD) [int; input]
 *        IER = error return flag [int, output]: 
 *	         0 = success, 
 *		 negative = error.
 * 
 *     When using the KLU solver the user must provide the FIDASPJAC routine for the 
 *     evalution of the sparse approximation to the Jacobian. To indicate that this
 *     routine has been provided, after the call to FIDAKLU, the following call must 
 *     be made    
 *
 *       CALL FIDASPARSESETJAC(IER) 
 *
 *     The int return flag IER=0 if successful, and nonzero otherwise.
 *
 *  
 *     The IDA KLU solver will reuse much of the factorization information from one
 *     nonlinear iteration to the next.  If at any time the user wants to force a full
 *     refactorization or if the number of nonzeros in the Jacobian matrix changes, the
 *     user should make the call
 *
 *       CALL FIDAKLUREINIT(NEQ, NNZ, REINIT_TYPE)
 *
 *     The arguments are:
 *        NEQ = the problem size [int; input]
 *        NNZ = the maximum number of nonzeros [int; input]
 *	REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be destroyed and 
 *          a new one will be allocated with NNZ nonzeros.  For a value of 2, 
 *	  only symbolic and numeric factorizations will be completed. 
 * 
 *     When using FIDAKLU, the user is required to supply the FIDASPJAC 
 *     routine for the evaluation of the sparse approximation to the 
 *     Jacobian, as discussed above with the other user-supplied routines.
 * 
 *     Optional outputs specific to the KLU case are:
 *        LSTF    = IOUT(16) from IDASlsGetLastFlag
 *        NJES    = IOUT(18) from IDASlsGetNumJacEvals
 *     See the IDA manual for descriptions.
 * 
 * (8.4s) SPARSE treatment of the linear system using the SuperLUMT solver.
 *
 *     The user must make the call
 *
 *       CALL FIDASUPERLUMT(NTHREADS, NEQ, NNZ, ORDERING, IER)
 *
 *     The arguments are:
 *        NTHREADS = desired number of threads to use [int; input]
 *        NEQ = the problem size [int; input]
 *        NNZ = the maximum number of nonzeros [int; input]
 *	ORDERING = the matrix ordering desired, possible values
 *	   come from the SuperLU_MT package [int; input]
 *           0 = Natural
 *           1 = Minimum degree on A^T A
 *           2 = Minimum degree on A^T + A
 *           3 = COLAMD
 *	IER = error return flag [int, output]: 
 *	         0 = success, 
 *		 negative = error.
 *	 
 *     At this time, there is no reinitialization capability for the SUNDIALS 
 *     interfaces to the SuperLUMT solver.
 *
 *     When using FIDASUPERLUMT, the user is required to supply the FIDASPJAC 
 *     routine for the evaluation of the CSC approximation to the Jacobian (note: the 
 *     current SuperLU_MT interface in SUNDIALS does not support CSR matrices). To 
 *     indicate that this routine has been provided, after the call to FIDASUPERLUMT,
 *     the following call must be made    
 *
 *       CALL FIDASPARSESETJAC(IER) 
 *
 *     The int return flag IER=0 if successful, and nonzero otherwise.
 * 
 *     Optional outputs specific to the SUPERLUMT case are:
 *        LSTF    = IOUT(14) from IDASlsGetLastFlag
 *        NJES    = IOUT(16) from IDASlsGetNumJacEvals
 *     See the IDA manual for descriptions.
 * 
 * (8.5) SPGMR treatment of the linear systems.
 * For the Scaled Preconditioned GMRES solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPGMR(MAXL, IGSTYPE, MAXRS, EPLIFAC, DQINCFAC, IER)
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * IGSTYPE  = specifies the type of Gram-Schmidt orthogonalization to be used:
 *            1 = MODIFIED_GS, 2 = CLASSICAL_GS
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to the matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * Optional outputs specific to the SPGMR case are:
 *        LENRWLS = IOUT(13) -> IDASpgmrGetWorkSpace
 *        LENIWLS = IOUT(14) -> IDASpgmrGetWorkSpace
 *        LSTF    = IOUT(15) -> IDASpgmrGetLastFlag
 *        NRELS   = IOUT(16) -> IDASpgmrGetResEvals
 *        NJE     = IOUT(17) -> IDASpgmrGetJtimesEvals
 *        NPE     = IOUT(18) -> IDASpgmrGetPrecEvals
 *        NPS     = IOUT(19) -> IDASpgmrGetPrecSolves
 *        NLI     = IOUT(20) -> IDASpgmrGetLinIters
 *        NLCF    = IOUT(21) -> IDASpgmrGetConvFails
 *
 * If a sequence of problems of the same size is being solved using the
 * SPGMR linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPGMRREINIT routine is needed if any of IGSTYPE, MAXRS, EPLIFAC, or
 * DQINCFAC is being changed.  In that case, call FIDASPGMRREINIT as follows:
 *       CALL FIDASPGMRREINIT (IGSTYPE, MAXRS, EPLIFAC, DQINCFAC, IER)              
 * The arguments have the same meanings as for FIDASPGMR.  If MAXL is being
 * changed, then call FIDASPGMR instead.
 *
 * (8.6) SPBCG treatment of the linear systems.
 * For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPBCG(MAXL, EPLIFAC, DQINCFAC, IER)              
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * Optional outputs specific to the SPBCG case are:
 *        LENRWLS = IOUT(13) -> IDASpbcgGetWorkSpace
 *        LENIWLS = IOUT(14) -> IDASpbcgGetWorkSpace
 *        LSTF    = IOUT(15) -> IDASpbcgGetLastFlag
 *        NRELS   = IOUT(16) -> IDASpbcgGetResEvals
 *        NJE     = IOUT(17) -> IDASpbcgGetJtimesEvals
 *        NPE     = IOUT(18) -> IDASpbcgGetPrecEvals
 *        NPS     = IOUT(19) -> IDASpbcgGetPrecSolves
 *        NLI     = IOUT(20) -> IDASpbcgGetLinIters
 *        NLCF    = IOUT(21) -> IDASpbcgGetConvFails
 *
 *      If a sequence of problems of the same size is being solved using the
 * SPBCG linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPBCGREINIT routine is needed if MAXL, EPLIFAC, or DQINCFAC is
 * being changed.  In that case, call FIDASPBCGREINIT as follows:
 *       CALL FIDASPBCGREINIT(MAXL, EPLIFAC, DQINCFAC, IER)
 * The arguments have the same meanings as for FIDASPBCG.
 *
 * (8.7) SPTFQMR treatment of the linear systems.
 * For the Scaled Preconditioned TFQMR solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPTFQMR(MAXL, EPLIFAC, DQINCFAC, IER)              
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * Optional outputs specific to the SPTFQMR case are:
 *        LENRWLS = IOUT(13) -> IDASptfqmrGetWorkSpace
 *        LENIWLS = IOUT(14) -> IDASptfqmrGetWorkSpace
 *        LSTF    = IOUT(15) -> IDASptfqmrGetLastFlag
 *        NRELS   = IOUT(16) -> IDASptfqmrGetResEvals
 *        NJE     = IOUT(17) -> IDASptfqmrGetJtimesEvals
 *        NPE     = IOUT(18) -> IDASptfqmrGetPrecEvals
 *        NPS     = IOUT(19) -> IDASptfqmrGetPrecSolves
 *        NLI     = IOUT(20) -> IDASptfqmrGetLinIters
 *        NLCF    = IOUT(21) -> IDASptfqmrGetConvFails
 *
 *      If a sequence of problems of the same size is being solved using the
 * SPTFQMR linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPTFQMRREINIT routine is needed if MAXL, EPLIFAC, or DQINCFAC is
 * being changed.  In that case, call FIDASPTFQMRREINIT as follows:
 *       CALL FIDASPTFQMRREINIT (MAXL, EPLIFAC, DQINCFAC, IER)
 * The arguments have the same meanings as for FIDASPTFQMR.
 *
 * (8.8) Using user-provided functions for the iterative linear solvers
 * 
 * If the user program includes the FIDAJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FIDASPILSSETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDAJTIMES is provided.  (FLAG = 0 specifies
 * using and internal finite difference approximation to this product.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * Usage of the user-supplied routines FIDAPSOL and FIDAPSET for solution of the 
 * preconditioner linear system requires the following call:
 *       CALL FIDASPILSSETPREC(FLAG, IER)
 * with FLAG = 1. The return flag IER is 0 if successful, nonzero otherwise.
 * The user-supplied routine FIDAPSOL must have the form:
 *       SUBROUTINE FIDAPSOL(T, Y, YP, R, RV, ZV, CJ, DELTA, EWT, 
 *      1                    IPAR, RPAR, WRK, IER)
 *       DIMENSION Y(*), YP(*), R(*), RV(*), ZV(*), 
 *      1          IPAR(*), RPAR(*), EWT(*), WRK(*)
 * This routine must solve the preconditioner linear system Pz = r, where r = RV
 * is input, and store the solution z in ZV.
 *
 * The user-supplied routine FIDAPSET must be of the form:
 *       SUBROUTINE FIDAPSET(T, Y, YP, R, CJ, EWT, H, IPAR, RPAR, 
 *      1                    WK1, WK2, WK3, IER)
 *       DIMENSION Y(*), YP(*), R(*), EWT(*), IPAR(*), RPAR(*), 
 *      1          WK1(*), WK2(*), WK3(*)
 * This routine must perform any evaluation of Jacobian-related data and
 * preprocessing needed for the solution of the preconditioner linear systems
 * by FIDAPSOL.  On return, set IER = 0 if FIDAPSET was successful, set IER
 * positive if a recoverable error occurred, and set IER negative if a
 * non-recoverable error occurred.
 * IPAR and RPAR are user (integer and real) arrays passed to FIDAMALLOC.
 *
 * -----------------------------------------------------------------------------
 *
 * (9) The solver: FIDASOLVE
 * To solve the DAE system, make the following call:
 *       CALL FIDASOLVE(TOUT, TRET, Y, YP, ITASK, IER)
 * The arguments are:
 * TOUT  = next value of t at which a solution is desired (input)
 * TRET  = value of t reached by the solver on output
 * Y     = array containing the computed solution on output
 * YP    = array containing current value of y'
 * ITASK = task indicator: 1 = normal mode (overshoot TOUT and interpolate)
 *         2 = one-step mode (return after each internal step taken)
 *         3 = normal tstop mode (like 1, but integration never proceeds past 
 *             TSTOP, which must be specified through a call to FIDASETRIN
 *             using the key 'STOP_TIME'
 *         4 = one step tstop (like 2, but integration never goes past TSTOP)
 * IER   = completion flag: 0 = success, 1 = tstop return, 2 = root return, 
 *         values -1 ... -10 are various failure modes (see IDA manual).
 * The current values of the optional outputs are available in IOUT and ROUT.
 *
 * -----------------------------------------------------------------------------
 *
 * (10) Getting current solution derivative: FIDAGETDKY
 * To obtain interpolated values of y and y' for any value of t in the last
 * internal step taken by IDA, make the following call:
 *       CALL FIDAGETDKY(T, K, DKY, IER)
 * The arguments are:
 * T   = value of t at which solution is desired, in [TCUR - HU,TCUR].
 * K   = order of derivative requested.
 * DKY = array containing computed K-th derivative of the solution y.
 * IER = return flag: = 0 for success, < 0 for illegal argument.
 *
 * -----------------------------------------------------------------------------
 *
 * (11) Memory freeing: FIDAFREE
 * To the free the internal memory created by the calls to FIDAMALLOC and
 * FNVINITS or FNVINITP, depending on the version (serial/parallel), make
 * the following call:
 *       CALL FIDAFREE
 *
 * =============================================================================
 */

#ifndef _FIDA_H
#define _FIDA_H

#include <ida/ida.h>                   /* definition of type IDAResFn */
#include <sundials/sundials_direct.h>  /* definition of type DlsMat  */
#include <sundials/sundials_sparse.h>  /* definition of type SlsMat  */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)

#define FIDA_MALLOC         SUNDIALS_F77_FUNC(fidamalloc, FIDAMALLOC)
#define FIDA_REINIT         SUNDIALS_F77_FUNC(fidareinit, FIDAREINIT)
#define FIDA_SETIIN         SUNDIALS_F77_FUNC(fidasetiin, FIDASETIIN)
#define FIDA_SETRIN         SUNDIALS_F77_FUNC(fidasetrin, FIDASETRIN)
#define FIDA_SETVIN         SUNDIALS_F77_FUNC(fidasetvin, FIDASETVIN)
#define FIDA_TOLREINIT      SUNDIALS_F77_FUNC(fidatolreinit, FIDATOLREINIT)
#define FIDA_SOLVE          SUNDIALS_F77_FUNC(fidasolve, FIDASOLVE)
#define FIDA_FREE           SUNDIALS_F77_FUNC(fidafree, FIDAFREE)
#define FIDA_CALCIC         SUNDIALS_F77_FUNC(fidacalcic, FIDACALCIC)
#define FIDA_BAND           SUNDIALS_F77_FUNC(fidaband, FIDABAND)
#define FIDA_BANDSETJAC     SUNDIALS_F77_FUNC(fidabandsetjac, FIDABANDSETJAC)
#define FIDA_DENSE          SUNDIALS_F77_FUNC(fidadense, FIDADENSE)
#define FIDA_DENSESETJAC    SUNDIALS_F77_FUNC(fidadensesetjac, FIDADENSESETJAC)
#define FIDA_LAPACKBAND        SUNDIALS_F77_FUNC(fidalapackband, FIDALAPACKBAND)
#define FIDA_LAPACKBANDSETJAC  SUNDIALS_F77_FUNC(fidalapackbandsetjac, FIDALAPACKBANDSETJAC)
#define FIDA_LAPACKDENSE       SUNDIALS_F77_FUNC(fidalapackdense, FIDALAPACKDENSE)
#define FIDA_LAPACKDENSESETJAC SUNDIALS_F77_FUNC(fidalapackdensesetjac, FIDALAPACKDENSESETJAC)
#define FIDA_KLU            SUNDIALS_F77_FUNC(fidaklu, FIDAKLU)
#define FIDA_KLUREINIT      SUNDIALS_F77_FUNC(fidaklureinit, FIDAKLUREINIT)
#define FIDA_SUPERLUMT      SUNDIALS_F77_FUNC(fidasuperlumt, FIDASUPERLUMT)
#define FIDA_SPARSESETJAC   SUNDIALS_F77_FUNC(fidasparsesetjac, FIDASPARSESETJAC)
#define FIDA_SPTFQMR        SUNDIALS_F77_FUNC(fidasptfqmr, FIDASPTFQMR)
#define FIDA_SPBCG          SUNDIALS_F77_FUNC(fidaspbcg, FIDASPBCG)
#define FIDA_SPGMR          SUNDIALS_F77_FUNC(fidaspgmr, FIDASPGMR)
#define FIDA_SPTFQMRREINIT  SUNDIALS_F77_FUNC(fidasptfqmrreinit, FIDASPTFQMRREINIT)
#define FIDA_SPBCGREINIT    SUNDIALS_F77_FUNC(fidaspbcgreinit, FIDASPBCGREINIT)
#define FIDA_SPGMRREINIT    SUNDIALS_F77_FUNC(fidaspgmrreinit, FIDASPGMRREINIT)
#define FIDA_SPILSSETJAC    SUNDIALS_F77_FUNC(fidaspilssetjac, FIDASPILSSETJAC)
#define FIDA_SPILSSETPREC   SUNDIALS_F77_FUNC(fidaspilssetprec, FIDASPILSSETPREC)
#define FIDA_RESFUN         SUNDIALS_F77_FUNC(fidaresfun, FIDARESFUN)
#define FIDA_DJAC           SUNDIALS_F77_FUNC(fidadjac, FIDADJAC)
#define FIDA_BJAC           SUNDIALS_F77_FUNC(fidabjac, FIDABJAC)
#define FIDA_PSET           SUNDIALS_F77_FUNC(fidapset, FIDAPSET)
#define FIDA_PSOL           SUNDIALS_F77_FUNC(fidapsol, FIDAPSOL)
#define FIDA_JTIMES         SUNDIALS_F77_FUNC(fidajtimes, FIDAJTIMES)
#define FIDA_EWT            SUNDIALS_F77_FUNC(fidaewt, FIDAEWT)
#define FIDA_EWTSET         SUNDIALS_F77_FUNC(fidaewtset, FIDAEWTSET)
#define FIDA_GETDKY         SUNDIALS_F77_FUNC(fidagetdky, FIDAGETDKY)
#define FIDA_GETERRWEIGHTS  SUNDIALS_F77_FUNC(fidageterrweights, FIDAGETERRWEIGHTS)
#define FIDA_GETESTLOCALERR SUNDIALS_F77_FUNC(fidagetestlocalerr, FIDAGETESTLOCALERR)

#else

#define FIDA_MALLOC         fidamalloc_
#define FIDA_REINIT         fidareinit_
#define FIDA_SETIIN         fidasetiin_
#define FIDA_SETRIN         fidasetrin_
#define FIDA_SETVIN         fidasetvin_
#define FIDA_TOLREINIT      fidatolreinit_
#define FIDA_SOLVE          fidasolve_
#define FIDA_FREE           fidafree_
#define FIDA_CALCIC         fidacalcic_
#define FIDA_BAND           fidaband_
#define FIDA_BANDSETJAC     fidabandsetjac_
#define FIDA_DENSE          fidadense_
#define FIDA_DENSESETJAC    fidadensesetjac_
#define FIDA_LAPACKBAND        fidalapackband_
#define FIDA_LAPACKBANDSETJAC  fidalapackbandsetjac_
#define FIDA_LAPACKDENSE       fidalapackdense_
#define FIDA_LAPACKDENSESETJAC fidalapackdensesetjac_
#define FIDA_KLU            fidaklu_
#define FIDA_KLUREINIT      fidaklureinit_
#define FIDA_SUPERLUMT      fidasuperlumt_
#define FIDA_SPARSESETJAC   fidasparsesetjac_
#define FIDA_SPTFQMR        fidasptfqmr_
#define FIDA_SPBCG          fidaspbcg_
#define FIDA_SPGMR          fidaspgmr_
#define FIDA_SPTFQMRREINIT  fidasptfqmrreinit_
#define FIDA_SPBCGREINIT    fidaspbcgreinit_
#define FIDA_SPGMRREINIT    fidaspgmrreinit_
#define FIDA_SPILSSETJAC    fidaspilssetjac_
#define FIDA_SPILSSETPREC   fidaspilssetprec_
#define FIDA_RESFUN         fidaresfun_
#define FIDA_DJAC           fidadjac_
#define FIDA_BJAC           fidabjac_
#define FIDA_PSET           fidapset_
#define FIDA_PSOL           fidapsol_
#define FIDA_JTIMES         fidajtimes_
#define FIDA_EWT            fidaewt_
#define FIDA_GETDKY         fidagetdky_
#define FIDA_GETERRWEIGHTS  fidageterrweights_
#define FIDA_GETESTLOCALERR fidagetestlocalerr_

#endif

/* Type for user data */

typedef struct {
  realtype *rpar;
  long int *ipar;
} *FIDAUserData;

/* Prototypes of exported functions */

void FIDA_MALLOC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout,
                 long int *ipar, realtype *rpar,
                 int *ier);
void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 int *ier);

void FIDA_SETIIN(char key_name[], long int *ival, int *ier);

void FIDA_SETRIN(char key_name[], realtype *rval, int *ier);

void FIDA_SETVIN(char key_name[], realtype *vval, int *ier);

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier);
void FIDA_CALCIC(int *icopt, realtype *tout1, int *ier);

void FIDA_DENSE(long int *neq, int *ier);
void FIDA_DENSESETJAC(int *flag, int *ier);
void FIDA_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
void FIDA_BANDSETJAC(int *flag, int *ier);

void FIDA_LAPACKDENSE(int *neq, int *ier);
void FIDA_LAPACKDENSESETJAC(int *flag, int *ier);
void FIDA_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
void FIDA_LAPACKBANDSETJAC(int *flag, int *ier);

void FIDA_KLU(int *neq, int *nnz, int *sparsetype, int *ordering, int *ier);
void FIDA_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier);
void FIDA_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier);
void FIDA_SPARSESETJAC(int *ier);

void FIDA_SPTFQMR(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPBCG(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPGMR(int *maxl, int *gstype, int *maxrs, realtype *eplifac,
                realtype *dqincfac, int *ier);
void FIDA_SPTFQMRREINIT(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPBCGREINIT(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPGMRREINIT(int *gstype, int *maxrs, realtype *eplifac,
                      realtype *dqincfac, int *ier);
void FIDA_SPILSSETJAC(int *flag, int *ier);
void FIDA_SPILSSETPREC(int *flag, int *ier);

void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
                realtype *ypret, int *itask, int *ier);
void FIDA_FREE(void);
void FIDA_EWTSET(int *flag, int *ier);
void FIDA_GETDKY(realtype *t, int *k, realtype *dky, int *ier);
void FIDA_GETERRWEIGHTS(realtype *eweight, int *ier);
void FIDA_GETESTLOCALERR(realtype *ele, int *ier);

/* Prototypes: Functions Called by the IDA Solver */

int FIDAresfn(realtype t, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data);

int FIDADenseJac(long int N, realtype t, realtype c_j, 
                 N_Vector yy, N_Vector yp, N_Vector rr,
                 DlsMat Jac, void *user_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDABandJac(long int N, long int mupper, long int mlower,
                realtype t, realtype c_j, 
                N_Vector yy, N_Vector yp, N_Vector rr,
                DlsMat Jac, void *user_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDALapackDenseJac(long int N, realtype t, realtype c_j, 
                       N_Vector yy, N_Vector yp, N_Vector rr,
                       DlsMat Jac, void *user_data,
                       N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDALapackBandJac(long int N, long int mupper, long int mlower,
                      realtype t, realtype c_j, 
                      N_Vector yy, N_Vector yp, N_Vector rr,
                      DlsMat J, void *user_data,
                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDASparseJac(realtype t, realtype cj, N_Vector y, N_Vector yp,
		  N_Vector fval, SlsMat J, void *user_data, 
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDAJtimes(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
               N_Vector v, N_Vector Jv,
               realtype c_j, void *user_data,
               N_Vector vtemp1, N_Vector vtemp2);

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             realtype c_j, void *user_data,
             N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             N_Vector rvec, N_Vector zvec,
             realtype c_j, realtype delta, void *user_data,
             N_Vector vtemp1);

int FIDAEwtSet(N_Vector yy, N_Vector ewt, void *user_data);

/* Declarations for global variables shared amongst various routines */

extern N_Vector F2C_IDA_vec;    /* defined in FNVECTOR module */

extern N_Vector F2C_IDA_ypvec;  /* defined in fida.c */
extern N_Vector F2C_IDA_ewtvec; /* defined in fida.c */
extern void *IDA_idamem;        /* defined in fida.c */
extern long int *IDA_iout;      /* defined in fida.c */
extern realtype *IDA_rout;      /* defined in fida.c */  
extern int IDA_ls;              /* defined in fida.c */
extern int IDA_nrtfn;           /* defined in fida.c */

/* Linear solver IDs */

enum { IDA_LS_DENSE = 1, IDA_LS_BAND = 2, 
       IDA_LS_LAPACKDENSE = 3, IDA_LS_LAPACKBAND = 4, 
       IDA_LS_KLU = 5, IDA_LS_SUPERLUMT = 6, 
       IDA_LS_SPGMR = 7, IDA_LS_SPBCG = 8, IDA_LS_SPTFQMR = 9 };

#ifdef __cplusplus
}
#endif

#endif
