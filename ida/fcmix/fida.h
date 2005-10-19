/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2005-10-19 20:23:11 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
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
 *
 *   FIDAMALLOC  interfaces to IDACreate and IDAMalloc
 *
 *   FIDAREINIT  interfaces to IDAReInit
 *
 *   FIDASETIIN, FIDASETRIN, FIDASETVIN interface to IDASet*
 *
 *   FIDATOLREINIT  interfaces to IDASetTolerances
 *
 *   FIDACALCIC  interfaces to IDACalcIC
 *
 *   FIDAEWTSET  interfaces to IDASetEwtFn
 *
 *   FIDADENSE        interfaces to IDADense
 *   FIDADENSESETJAC  interfaces to IDADenseSetJacFn
 *
 *   FIDABAND        interfaces to IDABand
 *   FIDABANDSETJAC  interfaces to IDABandSetJacFn
 *
 *   FIDASPTFQMR/FIDASPTFQMRREINIT  interface to IDASptfqmr and IDASptfqmrSet*
 *   FIDASPTFQMRSETJAC              interfaces to IDASptfqmrSetJacFn
 *   FIDASPTFQMRSETPREC             interfaces to IDASptfqmrSetPreconditioner
 *
 *   FIDASPBCG/FIDASPBCGREINIT  interface to IDASpbcg and IDASpbcgSet*
 *   FIDASPBCGSETJAC            interfaces to IDASpbcgSetJacFn
 *   FIDASPBCGSETPREC           interfaces to IDASpbcgSetPreconditioner
 *
 *   FIDASPGMR/FIDASPGMRREINIT  interface to IDASpgmr and IDASpgmrSet*
 *   FIDASPGMRSETJAC            interfaces to IDASpgmrSetJacFn
 *   FIDASPGMRSETPREC           interfaces to IDASpgmrSetPreconditioner
 *
 *   FIDASOLVE  interfaces to IDASolve, IDAGet*, and IDA*Get*
 *
 *   FIDAGETSOL  interfaces to IDAGetSolution
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
 * (1) User-supplied residual routine: FIDARESFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FIDARESFUN (T, Y, YP, R, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*)
 * It must set the R array to F(t,y,y'), the residual of the DAE 
 * system, as a function of T = t, the array Y = y, and the array YP = y'.
 * Here Y, YP and R are distributed vectors.
 *
 * (2s) Optional user-supplied dense Jacobian approximation routine: FIDADJAC
 * As an option when using the DENSE linear solver, the user may supply a
 * routine that computes a dense approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FIDADJAC (NEQ, T, Y, YP, R, DJAC, CJ, EWT, H,
 *      1                     WK1, WK2, WK3, IER)
 *       INTEGER NEQ, IER
 *       DIMENSION T, Y(*), YP(*), R(*), EWT(*), DJAC(NEQ,*), CJ, H,
 *      1          WK1(*), WK2(*), WK3(*)
 * This routine must compute the Jacobian and store it columnwise in DJAC.
 *
 * (3s) Optional user-supplied band Jacobian approximation routine: FIDABJAC
 * As an option when using the BAND linear solver, the user may supply a
 * routine that computes a band approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FIDABJAC (NEQ, MU, ML, MDIM, T, Y, YP, R, CJ,
 *      1                     BJAC, EWT, H, WK1, WK2, WK3, IER)
 *       INTEGER NEQ, MU, ML, MDIM, IER
 *       DIMENSION T, Y(*), YP(*), R(*), CJ, EWT(*), H, BJAC(MDIM,*), WK1(*),
 *      1          WK2(*), WK3(*)
 * This routine must load the MDIM by N array BJAC with the Jacobian matrix at the
 * current (t,y,y') in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
 * with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.
 *
 * (4) Optional user-supplied Jacobian-vector product routine: FIDAJTIMES
 * As an option when using the SPGMR/SPBCG/SPTFQMR linear solver, the user may
 * supply a routine that computes the product of the system Jacobian J = df/dy
 * and a given vector v.  If supplied, it must have the following form:
 *       SUBROUTINE FIDAJTIMES (T, Y, YP, R, V, FJV, CJ, EWT, H, WK1, WK2, IER)
 *       INTEGER IER
 *       DIMENSION T, V(*), FJV(*), Y(*), YP(*), R(*), CJ, H, EWT(*),
 *      1          WK1(*), WK2(*)
 * This routine must compute the product vector Jv, where the vector v is stored
 * in V, and store the product in FJV.  On return, set IER = 0 if FIDAJTIMES was
 * successful, and nonzero otherwise.
 *
 * (5) Optional user-supplied error weight vector routine: FIDAEWT
 * As an option to providing the relative and absolute tolerances, the user
 * may supply a routine that computes the weights used in the WRMS norms.
 * If supplied, it must have the following form:
 *       SUBROUTINE FIDAEWT (Y, EWT, IER)
 *       INTEGER IER
 *       DIMENSION Y(*), EWT(*)
 * It must store the error weights in EWT, given the current solution vector Y.
 * On return, set IER = 0 if successful, and nonzero otherwise.
 *
 * (6) Initialization:  FNVINITS / FNVINITP , FIDAMALLOC, FIDAREINIT,
 *                      FIDATOLREINIT, and FIDACALCIC
 *
 * (6.1s) To initialize the serial machine environment, the user must make
 * the following call:
 *        CALL FNVINITS (KEY, NEQ, IER)
 * The arguments are:
 * KEY = 2 for IDA
 * NEQ = size of vectors
 * IER = return completion flag. Values are 0 = success, -1 = failure.
 *
 * (6.1p) To initialize the parallel machine environment, the user must make 
 * one of the following calls:
 *        CALL FNVINITP (KEY, NLOCAL, NGLOBAL, IER)
 *                     -or-
 *        CALL FNVINITP (COMM, KEY, NLOCAL, NGLOBAL, IER)
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
 * (6.2) To set various problem and solution parameters and allocate
 * internal memory, make the following call:
 *       CALL FIDAMALLOC(T0, Y0, YP0, IATOL, RTOL, ATOL, IOUT, ROUT, IER)
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
 * IER   = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
 *         See printed message for details in case of failure.
 *
 * If the user program includes the FIDAEWT routine for the evaluation of the 
 * error weights, the following call must be made
 *       CALL FIDAEWTSET (FLAG, IER)
 * with FLAG = 1 to specify that FIDAEWT is provided.
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * (6.3) To set various integer optional inputs, make the folowing call:
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
 * (6.4) To re-initialize the FIDA solver for the solution of a new problem
 * of the same size as one already solved, make the following call:
 *       CALL FIDAREINIT (T0, Y0, YP0, IATOL, RTOL, ATOL, ID, CONSTR, IER)
 * The arguments have the same names and meanings as those of FIDAMALLOC.
 * FIDAREINIT performs the same initializations as FIDAMALLOC, but does no memory 
 * allocation for IDA data structures, using instead the existing internal memory
 * created by the previous FIDAMALLOC call.  The call to specify the linear system
 * solution method may or may not be needed.  See below.
 *
 * (6.5) To modify the tolerance parameters, make the following call:
 *       CALL FIDATOLREINIT (IATOL, RTOL, ATOL, IER)
 * The arguments have the same names and meanings as those of FIDAMALLOC.
 * FIDATOLREINIT simple calls IDASetTolerances with the given arguments.
 *
 * (6.6) To compute consistent initial conditions for an index-one DAE system,
 * make the following call:
 *       CALL FIDACALCIC(T0, Y0, YP0, ICOPT, TOUT, IER)
 * The arguments are:
 * T0    = initial value of t
 * Y0    = initial condition vector y(t0)
 * YP0   = initial condition vector y'(t0)
 * ICOPT = specifies the option: 1 = IDA_YP_YDP_INIT, 2 = IDA_Y_INIT.
 *         (See user guide for additional details.)
 * TOUT  = the first value of t at which a solution will be requested
 *         (from FIDASOLVE).
 * IER   = return completion flag.
 * 
 * (7) Specification of linear system solution method.
 * FIDA presently includes four choices for the treatment of these systems,
 * and the user of FIDA must call a routine with a specific name to make the
 * desired choice.
 * 
 * (7.1s) DENSE treatment of the linear system.
 * The user must make the call
 *       CALL FIDADENSE (NEQ, IER)
 * The arguments are:
 * NEQ = size of vectors
 * IER = error return flag: 0 = success , negative value = an error occured
 * 
 * If the user program includes the FIDADJAC routine for the evaluation of the 
 * dense approximation to the Jacobian, the following call must be made
 *       CALL FIDADENSESETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDADJAC is provided.  (FLAG = 0 specifies
 * using the internal finite differences approximation to the Jacobian.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 * 
 * Optional outputs specific to the DENSE case are:
 *        LENRWD = IOUT(13) -> IDADenseGetWorkSpace
 *        LENIWD = IOUT(14) -> IDADenseGetWorkSpace
 *        LSTF   = IOUT(15) -> IDADenseGetLastFlag
 *        NRED   = IOUT(16) -> IDADenseGetNumResEvals
 *        NJED   = IOUT(17) -> IDADenseGetNumJacEvals
 *
 * (7.2s) BAND treatment of the linear system
 * The user must make the call
 *       CALL FIDABAND (NEQ, MU, ML, IER)
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
 *        LENRWB = IOUT(13) -> IDABandGetWorkSpace
 *        LENIWB = IOUT(14) -> IDABandGetWorkSpace
 *        LSTF   = IOUT(15) -> IDABandGetLastFlag
 *        NREB   = IOUT(16) -> IDABandGetNumResEvals
 *        NJEB   = IOUT(17) -> IDABandGetNumJacEvals
 *
 * (7.3) SPGMR treatment of the linear systems.
 * For the Scaled Preconditioned GMRES solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPGMR (MAXL, IGSTYPE, MAXRS, EPLIFAC, DQINCFAC, IER)
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * IGSTYPE  = specifies the type of Gram-Schmidt orthogonalization to be used:
 *            1 = MODIFIED_GS, 2 = CLASSICAL_GS
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to the matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * If the user program includes the FIDAJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FIDASPGMRSETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDAJTIMES is provided.  (FLAG = 0 specifies
 * using and internal finite difference approximation to this product.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * Usage of the user-supplied routines FIDAPSOL and FIDAPSET for solution of the 
 * preconditioner linear system requires the following call:
 *       CALL FIDASPGMRSETPREC (FLAG, IER)
 * with FLAG = 1. The return flag IER is 0 if successful, nonzero otherwise.
 * The user-supplied routine FIDAPSOL must have the form:
 *       SUBROUTINE FIDAPSOL (T, Y, YP, R, RV, ZV, CJ, DELTA, EWT, WRK, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), RV(*), ZV(*), CJ, DELTA, EWT(*), WRK(*)
 * This routine must solve the preconditioner linear system Pz = r, where r = RV
 * is input, and store the solution z in ZV.
 *
 * The user-supplied routine FIDAPSET must be of the form:
 *       SUBROUTINE FIDAPSET(T, Y, YP, R, CJ, EWT, H, WK1, WK2, WK3, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), CJ, EWT(*), H, WK1(*), WK2(*), WK3(*)
 * This routine must perform any evaluation of Jacobian-related data and
 * preprocessing needed for the solution of the preconditioner linear systems
 * by FIDAPSOL.  On return, set IER = 0 if FIDAPSET was successful, set IER
 * positive if a recoverable error occurred, and set IER negative if a
 * non-recoverable error occurred.
 *
 * Optional outputs specific to the SPGMR case are:
 *        LENRWG = IOUT(13) -> IDASpgmrGetWorkSpace
 *        LENIWG = IOUT(14) -> IDASpgmrGetWorkSpace
 *        LSTF   = IOUT(15) -> IDASpgmrGetLastFlag
 *        NRE    = IOUT(16) -> IDASpgmrGetResEvals
 *        NJE    = IOUT(17) -> IDASpgmrGetJtimesEvals
 *        NPE    = IOUT(18) -> IDASpgmrGetPrecEvals
 *        NPS    = IOUT(19) -> IDASpgmrGetPrecSolves
 *        NLI    = IOUT(20) -> IDASpgmrGetLinIters
 *        NLCF   = IOUT(21) -> IDASpgmrGetConvFails
 *
 * If a sequence of problems of the same size is being solved using the
 * SPGMR linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPGMRREINIT routine is needed if any of IGSTYPE, MAXRS, EPLIFAC, or
 * DQINCFAC is being changed.  In that case, call FIDASPGMRREINIT as follows:
 *       CALL FIDASPGMRREINIT (IGSTYPE, MAXRS, EPLIFAC, DQINCFAC, IER)              
 * The arguments have the same meanings as for FIDASPGMR.  If MAXL is being
 * changed, then call FIDASPGMR instead.
 *
 * (7.4) SPBCG treatment of the linear systems.
 * For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPBCG (MAXL, EPLIFAC, DQINCFAC, IER)              
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * If the user program includes the FIDAJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FIDASPBCGSETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDAJTIMES is provided.  (FLAG = 0 specifies
 * using and internal finite difference approximation to this product.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * Usage of the user-supplied routines FIDAPSOL and FIDAPSET for solution of the 
 * preconditioner linear system requires the following call:
 *       CALL FIDASPBCGSETPREC (FLAG, IER)
 * with FLAG = 1. The return flag IER is 0 if successful, nonzero otherwise.
 * The user-supplied routine FIDAPSOL must have the form:
 *       SUBROUTINE FIDAPSOL (T, Y, YP, R, RV, ZV, CJ, DELTA, EWT, WRK, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), RV(*), ZV(*), CJ, DELTA, EWT(*), WRK(*)
 * This routine must solve the preconditioner linear system Pz = r, where r = RV
 * is input, and store the solution z in ZV.
 *
 * The user-supplied routine FIDAPSET must be of the form:
 *       SUBROUTINE FIDAPSET (T, Y, YP, R, CJ, EWT, H, WK1, WK2, WK3, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), CJ, EWT(*), H, WK1(*), WK2(*), WK3(*)
 * This routine must perform any evaluation of Jacobian-related data and
 * preprocessing needed for the solution of the preconditioner linear systems
 * by FIDAPSOL.  On return, set IER = 0 if FIDAPSET was successful, set IER
 * positive if a recoverable error occurred, and set IER negative if a
 * non-recoverable error occurred.
 *
 * Optional outputs specific to the SPBCG case are:
 *        LENRWC = IOUT(13) -> IDASpbcgGetWorkSpace
 *        LENIWC = IOUT(14) -> IDASpbcgGetWorkSpace
 *        LSTF   = IOUT(15) -> IDASpbcgGetLastFlag
 *        NRE    = IOUT(16) -> IDASpbcgGetResEvals
 *        NJE    = IOUT(17) -> IDASpbcgGetJtimesEvals
 *        NPE    = IOUT(18) -> IDASpbcgGetPrecEvals
 *        NPS    = IOUT(19) -> IDASpbcgGetPrecSolves
 *        NLI    = IOUT(20) -> IDASpbcgGetLinIters
 *        NLCF   = IOUT(21) -> IDASpbcgGetConvFails
 *
 *      If a sequence of problems of the same size is being solved using the
 * SPBCG linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPBCGREINIT routine is needed if MAXL, EPLIFAC, or DQINCFAC is
 * being changed.  In that case, call FIDASPBCGREINIT as follows:
 *       CALL FIDASPBCGREINIT (MAXL, EPLIFAC, DQINCFAC, IER)
 * The arguments have the same meanings as for FIDASPBCG.
 *
 * (7.5) SPTFQMR treatment of the linear systems.
 * For the Scaled Preconditioned TFQMR solution of the linear systems,
 * the user must make the following call:
 *       CALL FIDASPTFQMR (MAXL, EPLIFAC, DQINCFAC, IER)              
 * The arguments are:
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * EPLIFAC  = factor in the linear iteration convergence test constant
 * DQINCFAC = factor in the increments to y used in the difference quotient
 *            approximations to matrix-vector products Jv
 * IER      = error return flag: 0 = success; negative value = an error occured
 *
 * If the user program includes the FIDAJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FIDASPTFQMRSETJAC (FLAG, IER)
 * with FLAG = 1 to specify that FIDAJTIMES is provided.  (FLAG = 0 specifies
 * using and internal finite difference approximation to this product.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * Usage of the user-supplied routines FIDAPSOL and FIDAPSET for solution of the 
 * preconditioner linear system requires the following call:
 *       CALL FIDASPTFQMRSETPREC (FLAG, IER)
 * with FLAG = 1. The return flag IER is 0 if successful, nonzero otherwise.
 * The user-supplied routine FIDAPSOL must have the form:
 *       SUBROUTINE FIDAPSOL (T, Y, YP, R, RV, ZV, CJ, DELTA, EWT, WRK, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), RV(*), ZV(*), CJ, DELTA, EWT(*), WRK(*)
 * This routine must solve the preconditioner linear system Pz = r, where r = RV
 * is input, and store the solution z in ZV.
 *
 * The user-supplied routine FIDAPSET must be of the form:
 *       SUBROUTINE FIDAPSET (T, Y, YP, R, CJ, EWT, H, WK1, WK2, WK3, IER)
 *       INTEGER IER
 *       DIMENSION T, Y(*), YP(*), R(*), CJ, EWT(*), H, WK1(*), WK2(*), WK3(*)
 * This routine must perform any evaluation of Jacobian-related data and
 * preprocessing needed for the solution of the preconditioner linear systems
 * by FIDAPSOL.  On return, set IER = 0 if FIDAPSET was successful, set IER
 * positive if a recoverable error occurred, and set IER negative if a
 * non-recoverable error occurred.
 *
 * Optional outputs specific to the SPTFQMR case are:
 *        LENRWC = IOUT(13) -> IDASptfqmrGetWorkSpace
 *        LENIWC = IOUT(14) -> IDASptfqmrGetWorkSpace
 *        LSTF   = IOUT(15) -> IDASptfqmrGetLastFlag
 *        NRE    = IOUT(16) -> IDASptfqmrGetResEvals
 *        NJE    = IOUT(17) -> IDASptfqmrGetJtimesEvals
 *        NPE    = IOUT(18) -> IDASptfqmrGetPrecEvals
 *        NPS    = IOUT(19) -> IDASptfqmrGetPrecSolves
 *        NLI    = IOUT(20) -> IDASptfqmrGetLinIters
 *        NLCF   = IOUT(21) -> IDASptfqmrGetConvFails
 *
 *      If a sequence of problems of the same size is being solved using the
 * SPTFQMR linear solver, then following the call to FIDAREINIT, a call to the
 * FIDASPTFQMRREINIT routine is needed if MAXL, EPLIFAC, or DQINCFAC is
 * being changed.  In that case, call FIDASPTFQMRREINIT as follows:
 *       CALL FIDASPTFQMRREINIT (MAXL, EPLIFAC, DQINCFAC, IER)
 * The arguments have the same meanings as for FIDASPTFQMR.
 *
 * (8) The solver: FIDASOLVE
 * To solve the DAE system, make the following call:
 *       CALL FIDASOLVE (TOUT, TRET, Y, YP, ITASK, IER)
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
 * (9) Getting current solution: FIDAGETSOL
 * To obtain interpolated values of y and y' for any value of t in the last
 * internal step taken by IDA, make the following call:
 *       CALL FIDAGETSOL (T, YRET, YPRET, IER)
 * The arguments are:
 * T   = value of t at which solution is desired, in [TCUR-HU,TCUR].
 * Y   = array containing interpolated y
 * YP  = array containing the derivative of the computed solution, y'(tret)
 * IER = return flag: = 0 for success, < 0 for illegal argument.
 *
 * (10) Memory freeing: FIDAFREE
 * To the free the internal memory created by the calls to FIDAMALLOC and
 * FNVINITS or FNVINITP, depending on the version (serial/parallel), make
 * the following call:
 *       CALL FIDAFREE
 *
 * =============================================================================
 */

#ifndef _FIDA_H
#define _FIDA_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "band.h"           /* definition of type BandMat  */
#include "ida.h"            /* definition of type IDAResFn */
#include "dense.h"          /* definition of type DenseMat */
#include "nvector.h"        /* definition of type N_Vector */
#include "sundialstypes.h"  /* definition of type realtype */

#if defined(F77_FUNC)

#define FIDA_MALLOC         F77_FUNC(fidamalloc, FIDAMALLOC)
#define FIDA_REINIT         F77_FUNC(fidareinit, FIDAREINIT)
#define FIDA_SETIIN         F77_FUNC(fidasetiin, FIDASETIIN)
#define FIDA_SETRIN         F77_FUNC(fidasetrin, FIDASETRIN)
#define FIDA_SETVIN         F77_FUNC(fidasetvin, FIDASETVIN)
#define FIDA_TOLREINIT      F77_FUNC(fidatolreinit, FIDATOLREINIT)
#define FIDA_SOLVE          F77_FUNC(fidasolve, FIDASOLVE)
#define FIDA_FREE           F77_FUNC(fidafree, FIDAFREE)
#define FIDA_CALCIC         F77_FUNC(fidacalcic, FIDACALCIC)
#define FIDA_SPTFQMR        F77_FUNC(fidasptfqmr, FIDASPTFQMR)
#define FIDA_SPBCG          F77_FUNC(fidaspbcg, FIDASPBCG)
#define FIDA_SPGMR          F77_FUNC(fidaspgmr, FIDASPGMR)
#define FIDA_SPTFQMRREINIT  F77_FUNC(fidasptfqmrreinit, FIDASPTFQMRREINIT)
#define FIDA_SPBCGREINIT    F77_FUNC(fidaspbcgreinit, FIDASPBCGREINIT)
#define FIDA_SPGMRREINIT    F77_FUNC(fidaspgmrreinit, FIDASPGMRREINIT)
#define FIDA_BAND           F77_FUNC(fidaband, FIDABAND)
#define FIDA_BJAC           F77_FUNC(fidabjac, FIDABJAC)
#define FIDA_BANDSETJAC     F77_FUNC(fidabandsetjac, FIDABANDSETJAC)
#define FIDA_DENSE          F77_FUNC(fidadense, FIDADENSE)
#define FIDA_DENSESETJAC    F77_FUNC(fidadensesetjac, FIDADENSESETJAC)
#define FIDA_DJAC           F77_FUNC(fidadjac, FIDADJAC)
#define FIDA_SPTFQMRSETJAC  F77_FUNC(fidasptfqmrsetjac, FIDASPTFQMRSETJAC)
#define FIDA_SPBCGSETJAC    F77_FUNC(fidaspbcgsetjac, FIDASPBCGSETJAC)
#define FIDA_SPGMRSETJAC    F77_FUNC(fidaspgmrsetjac, FIDASPGMRSETJAC)
#define FIDA_JTIMES         F77_FUNC(fidajtimes, FIDAJTIMES)
#define FIDA_SPTFQMRSETPREC F77_FUNC(fidasptfqmrsetprec, FIDASPTFQMRSETPREC)
#define FIDA_SPBCGSETPREC   F77_FUNC(fidaspbcgsetprec, FIDASPBCGSETPREC)
#define FIDA_SPGMRSETPREC   F77_FUNC(fidaspgmrsetprec, FIDASPGMRSETPREC)
#define FIDA_PSET           F77_FUNC(fidapset, FIDAPSET)
#define FIDA_PSOL           F77_FUNC(fidapsol, FIDAPSOL)
#define FIDA_RESFUN         F77_FUNC(fidaresfun, FIDARESFUN)
#define FIDA_EWT            F77_FUNC(fidaewt, FIDAEWT)
#define FIDA_GETSOL         F77_FUNC(fidagetsol, FIDAGETSOL)
#define FIDA_GETERRWEIGHTS  F77_FUNC(fidageterrweights, FIDAGETERRWEIGHTS)
#define FIDA_GETESTLOCALERR F77_FUNC(fidagetestlocalerr, FIDAGETESTLOCALERR)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FIDA_MALLOC         fidamalloc
#define FIDA_REINIT         fidareinit
#define FIDA_SETIIN         fidasetiin
#define FIDA_SETRIN         fidasetrin
#define FIDA_SETVIN         fidasetvin
#define FIDA_TOLREINIT      fidatolreinit
#define FIDA_SOLVE          fidasolve
#define FIDA_FREE           fidafree
#define FIDA_CALCIC         fidacalcic
#define FIDA_SPTFQMR        fidasptfqmr
#define FIDA_SPBCG          fidaspbcg
#define FIDA_SPGMR          fidaspgmr
#define FIDA_SPTFQMRREINIT  fidasptfqmrreinit
#define FIDA_SPBCGREINIT    fidaspbcgreinit
#define FIDA_SPGMRREINIT    fidaspgmrreinit
#define FIDA_BAND           fidaband
#define FIDA_BJAC           fidabjac
#define FIDA_BANDSETJAC     fidabandsetjac
#define FIDA_DENSE          fidadense
#define FIDA_DENSESETJAC    fidadensesetjac
#define FIDA_DJAC           fidadjac
#define FIDA_SPTFQMRSETJAC  fidasptfqmrsetjac
#define FIDA_SPBCGSETJAC    fidaspbcgsetjac
#define FIDA_SPGMRSETJAC    fidaspgmrsetjac
#define FIDA_JTIMES         fidajtimes
#define FIDA_SPTFQMRSETPREC fidasptfqmrsetprec
#define FIDA_SPBCGSETPREC   fidaspbcgsetprec
#define FIDA_SPGMRSETPREC   fidaspgmrsetprec
#define FIDA_PSET           fidapset
#define FIDA_PSOL           fidapsol
#define FIDA_RESFUN         fidaresfun
#define FIDA_EWT            fidaewt
#define FIDA_GETSOL         fidagetsol
#define FIDA_GETERRWEIGHTS  fidageterrweights
#define FIDA_GETESTLOCALERR fidagetestlocalerr

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FIDA_MALLOC         FIDAMALLOC
#define FIDA_REINIT         FIDAREINIT
#define FIDA_SETIIN         FIDASETIIN
#define FIDA_SETRIN         FIDASETRIN
#define FIDA_SETVIN         FIDASETVIN
#define FIDA_TOLREINIT      FIDATOLREINIT
#define FIDA_SOLVE          FIDASOLVE
#define FIDA_FREE           FIDAFREE
#define FIDA_CALCIC         FIDACALCIC
#define FIDA_SPTFQMR        FIDASPTFQMR
#define FIDA_SPBCG          FIDASPBCG
#define FIDA_SPGMR          FIDASPGMR
#define FIDA_SPTFQMRREINIT  FIDASPTFQMRREINIT
#define FIDA_SPBCGREINIT    FIDASPBCGREINIT
#define FIDA_SPGMRREINIT    FIDASPGMRREINIT
#define FIDA_BAND           FIDABAND
#define FIDA_BJAC           FIDABJAC
#define FIDA_BANDSETJAC     FIDABANDSETJAC
#define FIDA_DENSE          FIDADENSE
#define FIDA_DENSESETJAC    FIDADENSESETJAC
#define FIDA_DJAC           FIDADJAC
#define FIDA_SPTFQMRSETJAC  FIDASPTFQMRSETJAC
#define FIDA_SPBCGSETJAC    FIDASPBCGSETJAC
#define FIDA_SPGMRSETJAC    FIDASPGMRSETJAC
#define FIDA_JTIMES         FIDAJTIMES
#define FIDA_SPTFQMRSETPREC FIDASPTFQMRSETPREC
#define FIDA_SPBCGSETPREC   FIDASPBCGSETPREC
#define FIDA_SPGMRSETPREC   FIDASPGMRSETPREC
#define FIDA_PSET           FIDAPSET
#define FIDA_PSOL           FIDAPSOL
#define FIDA_RESFUN         FIDARESFUN
#define FIDA_EWT            FIDAEWT
#define FIDA_GETSOL         FIDAGETSOL
#define FIDA_GETERRWEIGHTS  FIDAGETERRWEIGHTS
#define FIDA_GETESTLOCALERR FIDAGETESTLOCALERR

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FIDA_MALLOC         fidamalloc_
#define FIDA_REINIT         fidareinit_
#define FIDA_SETIIN         fidasetiin_
#define FIDA_SETRIN         fidasetrin_
#define FIDA_SETVIN         fidasetvin_
#define FIDA_TOLREINIT      fidatolreinit_
#define FIDA_SOLVE          fidasolve_
#define FIDA_FREE           fidafree_
#define FIDA_CALCIC         fidacalcic_
#define FIDA_SPTFQMR        fidasptfqmr_
#define FIDA_SPBCG          fidaspbcg_
#define FIDA_SPGMR          fidaspgmr_
#define FIDA_SPTFQMRREINIT  fidasptfqmrreinit_
#define FIDA_SPBCGREINIT    fidaspbcgreinit_
#define FIDA_SPGMRREINIT    fidaspgmrreinit_
#define FIDA_BAND           fidaband_
#define FIDA_BJAC           fidabjac_
#define FIDA_BANDSETJAC     fidabandsetjac_
#define FIDA_DENSE          fidadense_
#define FIDA_DENSESETJAC    fidadensesetjac_
#define FIDA_DJAC           fidadjac_
#define FIDA_SPTFQMRSETJAC  fidasptfqmrsetjac_
#define FIDA_SPBCGSETJAC    fidaspbcgsetjac_
#define FIDA_SPGMRSETJAC    fidaspgmrsetjac_
#define FIDA_JTIMES         fidajtimes_
#define FIDA_SPTFQMRSETPREC fidasptfqmrsetprec_
#define FIDA_SPBCGSETPREC   fidaspbcgsetprec_
#define FIDA_SPGMRSETPREC   fidaspgmrsetprec_
#define FIDA_PSET           fidapset_
#define FIDA_PSOL           fidapsol_
#define FIDA_RESFUN         fidaresfun_
#define FIDA_EWT            fidaewt_
#define FIDA_GETSOL         fidagetsol_
#define FIDA_GETERRWEIGHTS  fidageterrweights_
#define FIDA_GETESTLOCALERR fidagetestlocalerr_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FIDA_MALLOC         FIDAMALLOC_
#define FIDA_REINIT         FIDAREINIT_
#define FIDA_SETIIN         FIDASETIIN_
#define FIDA_SETRIN         FIDASETRIN_
#define FIDA_SETVIN         FIDASETVIN_
#define FIDA_TOLREINIT      FIDATOLREINIT_
#define FIDA_SOLVE          FIDASOLVE_
#define FIDA_FREE           FIDAFREE_
#define FIDA_CALCIC         FIDACALCIC_
#define FIDA_SPTFQMR        FIDASPTFQMR_
#define FIDA_SPBCG          FIDASPBCG_
#define FIDA_SPGMR          FIDASPGMR_
#define FIDA_SPTFQMRREINIT  FIDASPTFQMRREINIT_
#define FIDA_SPBCGREINIT    FIDASPBCGREINIT_
#define FIDA_SPGMRREINIT    FIDASPGMRREINIT_
#define FIDA_BAND           FIDABAND_
#define FIDA_BJAC           FIDABJAC_
#define FIDA_BANDSETJAC     FIDABANDSETJAC_
#define FIDA_DENSE          FIDADENSE_
#define FIDA_DENSESETJAC    FIDADENSESETJAC_
#define FIDA_DJAC           FIDADJAC_
#define FIDA_SPTFQMRSETJAC  FIDASPTFQMRSETJAC_
#define FIDA_SPBCGSETJAC    FIDASPBCGSETJAC_
#define FIDA_SPGMRSETJAC    FIDASPGMRSETJAC_
#define FIDA_JTIMES         FIDAJTIMES_
#define FIDA_SPTFQMRSETPREC FIDASPTFQMRSETPREC_
#define FIDA_SPBCGSETPREC   FIDASPBCGSETPREC_
#define FIDA_SPGMRSETPREC   FIDASPGMRSETPREC_
#define FIDA_PSET           FIDAPSET_
#define FIDA_PSOL           FIDAPSOL_
#define FIDA_RESFUN         FIDARESFUN_
#define FIDA_EWT            FIDAEWT_
#define FIDA_GETSOL         FIDAGETSOL_
#define FIDA_GETERRWEIGHTS  FIDAGETERRWEIGHTS_
#define FIDA_GETESTLOCALERR FIDAGETESTLOCALERR_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FIDA_MALLOC         fidamalloc__
#define FIDA_REINIT         fidareinit__
#define FIDA_SETIIN         fidasetiin__
#define FIDA_SETRIN         fidasetrin__
#define FIDA_SETVIN         fidasetvin__
#define FIDA_TOLREINIT      fidatolreinit__
#define FIDA_SOLVE          fidasolve__
#define FIDA_FREE           fidafree__
#define FIDA_CALCIC         fidacalcic__
#define FIDA_SPTFQMR        fidasptfqmr__
#define FIDA_SPBCG          fidaspbcg__
#define FIDA_SPGMR          fidaspgmr__
#define FIDA_SPTFQMRREINIT  fidasptfqmrreinit__
#define FIDA_SPBCGREINIT    fidaspbcgreinit__
#define FIDA_SPGMRREINIT    fidaspgmrreinit__
#define FIDA_BAND           fidaband__
#define FIDA_BJAC           fidabjac__
#define FIDA_BANDSETJAC     fidabandsetjac__
#define FIDA_DENSE          fidadense__
#define FIDA_DENSESETJAC    fidadensesetjac__
#define FIDA_DJAC           fidadjac__
#define FIDA_SPTFQMRSETJAC  fidasptfqmrsetjac__
#define FIDA_SPBCGSETJAC    fidaspbcgsetjac__
#define FIDA_SPGMRSETJAC    fidaspgmrsetjac__
#define FIDA_JTIMES         fidajtimes__
#define FIDA_SPTFQMRSETPREC fidasptfqmrsetprec__
#define FIDA_SPBCGSETPREC   fidaspbcgsetprec__
#define FIDA_SPGMRSETPREC   fidaspgmrsetprec__
#define FIDA_PSET           fidapset__
#define FIDA_PSOL           fidapsol__
#define FIDA_RESFUN         fidaresfun__
#define FIDA_EWT            fidaewt__
#define FIDA_GETSOL         fidagetsol__
#define FIDA_GETERRWEIGHTS  fidageterrweights__
#define FIDA_GETESTLOCALERR fidagetestlocalerr__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FIDA_MALLOC         FIDAMALLOC__
#define FIDA_REINIT         FIDAREINIT__
#define FIDA_SETIIN         FIDASTIIN__
#define FIDA_SETRIN         FIDASETRIN__
#define FIDA_SETVIN         FIDASETVIN__
#define FIDA_TOLREINIT      FIDATOLREINIT__
#define FIDA_SOLVE          FIDASOLVE__
#define FIDA_FREE           FIDAFREE__
#define FIDA_CALCIC         FIDACALCIC__
#define FIDA_SPTFQMR        FIDASPTFQMR__
#define FIDA_SPBCG          FIDASPBCG__
#define FIDA_SPGMR          FIDASPGMR__
#define FIDA_SPTFQMRREINIT  FIDASPTFQMRREINIT__
#define FIDA_SPBCGREINIT    FIDASPBCGREINIT__
#define FIDA_SPGMRREINIT    FIDASPGMRREINIT__
#define FIDA_BAND           FIDABAND__
#define FIDA_BJAC           FIDABJAC__
#define FIDA_BANDSETJAC     FIDABANDSETJAC__
#define FIDA_DENSE          FIDADENSE__
#define FIDA_DENSESETJAC    FIDADENSESETJAC__
#define FIDA_DJAC           FIDADJAC__
#define FIDA_SPTFQMRSETJAC  FIDASPTFQMRSETJAC__
#define FIDA_SPBCGSETJAC    FIDASPBCGSETJAC__
#define FIDA_SPGMRSETJAC    FIDASPGMRSETJAC__
#define FIDA_JTIMES         FIDAJTIMES__
#define FIDA_SPTFQMRSETPREC FIDASPTFQMRSETPREC__
#define FIDA_SPBCGSETPREC   FIDASPBCGSETPREC__
#define FIDA_SPGMRSETPREC   FIDASPGMRSETPREC__
#define FIDA_PSET           FIDAPSET__
#define FIDA_PSOL           FIDAPSOL__
#define FIDA_RESFUN         FIDARESFUN__
#define FIDA_EWT            FIDAEWT__
#define FIDA_GETSOL         FIDAGETSOL__
#define FIDA_GETERRWEIGHTS  FIDAGETERRWEIGHTS__
#define FIDA_GETESTLOCALERR FIDAGETESTLOCALERR__

#endif

/* Prototypes of exported functions */

void FIDA_MALLOC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout,
                 int *ier);
void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 int *ier);

void FIDA_SETIIN(char key_name[], long int *ival, int *ier, int key_len);

void FIDA_SETRIN(char key_name[], realtype *rval, int *ier, int key_len);

void FIDA_SETVIN(char key_name[], realtype *vval, int *ier, int key_len);

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier);
void FIDA_CALCIC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *icopt, realtype *tout1, int *ier);
void FIDA_SPTFQMR(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPBCG(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPGMR(int *maxl, int *gstype, int *maxrs, realtype *eplifac,
                realtype *dqincfac, int *ier);
void FIDA_DENSE(long int *neq, int *ier);
void FIDA_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
void FIDA_SPTFQMRREINIT(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPBCGREINIT(int *maxl, realtype *eplifac, realtype *dqincfac, int *ier);
void FIDA_SPGMRREINIT(int *gstype, int *maxrs, realtype *eplifac,
                      realtype *dqincfac, int *ier);
void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
                realtype *ypret, int *itask, int *ier);
void FIDA_FREE(void);
void FIDA_BANDSETJAC(int *flag, int *ier);
void FIDA_DENSESETJAC(int *flag, int *ier);
void FIDA_SPGMRSETJAC(int *flag, int *ier);
void FIDA_SPBCGSETJAC(int *flag, int *ier);
void FIDA_SPTFQMRSETJAC(int *flag, int *ier);
void FIDA_SPGMRSETPREC(int *flag, int *ier);
void FIDA_SPBCGSETPREC(int *flag, int *ier);
void FIDA_SPTFQMRSETPREC(int *flag, int *ier);
void FIDA_EWTSET(int *flag, int *ier);
void FIDA_GETSOL(realtype *t, realtype *yret, realtype *ypret, int *ier);
void FIDA_GETERRWEIGHTS(realtype *eweight, int *ier);
void FIDA_GETESTLOCALERR(realtype *ele, int *ier);

/* Prototypes: Functions Called by the IDA Solver */

int FIDAresfn(realtype t, N_Vector yy, N_Vector yp, N_Vector rr, void *res_data);

int FIDADenseJac(long int N, realtype t,
                 N_Vector yy, N_Vector yp, N_Vector rr,
                 realtype c_j, void *jac_data, DenseMat Jac,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDABandJac(long int N, long int mupper, long int mlower,
                BandMat J, realtype t,
                N_Vector yy, N_Vector yp, N_Vector rr,
                realtype c_j, void *jac_data, BandMat Jac,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDAJtimes(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
               N_Vector v, N_Vector Jv,
               realtype c_j, void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2);

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             realtype c_j, void *prec_data,
             N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             N_Vector rvec, N_Vector zvec,
             realtype c_j, realtype delta, void *prec_data,
             N_Vector vtemp1);

int FIDAEwtSet(N_Vector yy, N_Vector ewt, void *e_data);

/* Declarations for global variables shared amongst various routines */

extern N_Vector F2C_IDA_vec;

extern N_Vector F2C_IDA_ypvec, F2C_IDA_ewtvec;
extern void *IDA_idamem;
extern long int *IDA_iout;
extern realtype *IDA_rout;
extern int IDA_ls;
extern int IDA_nrtfn;

/* Linear solver IDs */

enum { IDA_LS_DENSE = 1, IDA_LS_BAND = 2, IDA_LS_SPGMR = 3, 
       IDA_LS_SPBCG = 4, IDA_LS_SPTFQMR = 5 };

#ifdef __cplusplus
}
#endif

#endif
