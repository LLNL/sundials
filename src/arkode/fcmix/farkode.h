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
 * This is the header file for FARKODE, the Fortran interface to
 * the ARKODE package.
 *--------------------------------------------------------------*/

/*===============================================================
                FARKODE Interface Package

 The FARKODE Interface Package is a package of C functions which
 support the use of the ARKODE solver, for the solution of ODE
 systems
         M(t) dy/dt = fe(t,y) + fi(t,y),
 in a mixed Fortran/C setting.  While ARKODE is written in C, it
 is assumed here that the user's calling program and user-supplied
 problem-defining routines are written in Fortran.  This package
 provides the necessary interface to ARKODE for any acceptable
 NVECTOR implementation.

 A summary of the user-callable functions, with the corresponding
 ARKODE functions, are as follows:

   Fortran                    ARKODE
   ---------------------      --------------------------------
   FNVINITS                   N_VNew_Serial
   FNVINITP                   N_VNew_Parallel
   FNVINITOMP                 N_VNew_OpenMP
   FNVINITPTS                 N_VNew_Pthreads

   FSUNBANDMATINIT            SUNBandMatrix
   FSUNDENSEMATINIT           SUNDenseMatrix
   FSUNSPARSEMATINIT          SUNSparseMatrix

   FSUNBANDMASSMATINIT        SUNBandMatrix
   FSUNDENSEMASSMATINIT       SUNDenseMatrix
   FSUNSPARSEMASSMATINIT      SUNSparseMatrix

   FSUNBANDLINSOLINIT         SUNBandLinearSolver
   FSUNDENSELINSOLINIT        SUNDenseLinearSolver
   FSUNKLUINIT                SUNKLU
   FSUNKLUREINIT              SUNKLUReinit
   FSUNLAPACKBANDINIT         SUNLapackBand
   FSUNLAPACKDENSEINIT        SUNLapackDense
   FSUNPCGINIT                SUNPCG
   FSUNSPBCGSINIT             SUNSPBCGS
   FSUNSPFGMRINIT             SUNSPFGMR
   FSUNSPGMRINIT              SUNSPGMR
   FSUNSPTFQMRINIT            SUNSPTFQMR
   FSUNSUPERLUMTINIT          SUNSuperLUMT

   FSUNMASSBANDLINSOLINIT     SUNBandLinearSolver
   FSUNMASSDENSELINSOLINIT    SUNDenseLinearSolver
   FSUNMASSKLUINIT            SUNKLU
   FSUNMASSKLUREINIT          SUNKLUReinit
   FSUNMASSLAPACKBANDINIT     SUNLapackBand
   FSUNMASSLAPACKDENSEINIT    SUNLapackDense
   FSUNMASSPCGINIT            SUNPCG
   FSUNMASSSPBCGSINIT         SUNSPBCGS
   FSUNMASSSPFGMRINIT         SUNSPFGMR
   FSUNMASSSPGMRINIT          SUNSPGMR
   FSUNMASSSPTFQMRINIT        SUNSPTFQMR
   FSUNMASSSUPERLUMTINIT      SUNSuperLUMT

   FARKMALLOC                 ARKodeCreate, ARKodeSetUserData,
                                 and ARKodeInit
   FARKREINIT                 ARKReInit
   FARKRESIZE                 ARKResize
   FARKSETIIN                 ARKodeSet* (integer arguments)
   FARKSETRIN                 ARKodeSet* (real arguments)
   FARKSETADAPTIVITYMETHOD    ARKodeSetAdaptivityMethod
   FARKSETDEFAULTS            ARKodeSetDefaults
   FARKSETERKTABLE            ARKodeSetERKTable
   FARKSETIRKTABLE            ARKodeSetIRKTable
   FARKSETARKTABLES           ARKodeSetARKTables
   FARKSETRESTOLERANCE        ARKodeResStolerance, ARKodeResVtolerance
   FARKEWTSET                 ARKodeWFtolerances
   FARKADAPTSET               ARKodeSetAdaptivityFn
   FARKEXPSTABSET             ARKodeSetStabilityFn
   FARKSETDIAGNOSTICS         ARKodeSetDiagnostics
   FARKSTOPDIAGNOSTICS        (none)

   FARKDLSINIT                ARKDlsSetLinearSolver
   FARKDENSESETJAC            ARKDlsSetJacFn
   FARKBANDSETJAC             ARKDlsSetJacFn
   FARKSPARSESETJAC           ARKDlsSetJacFn

   FARKDLSMASSINIT            ARKDlsSetMassLinearSolver
   FARKDENSESETMASS           ARKDlsSetMassFn
   FARKBANDSETMASS            ARKDlsSetMassFn
   FARKSPARSESETMASS          ARKDlsSetMassFn

   FARKSPILSINIT              ARKSpilsSetLinearSolver
   FARKSPILSSETEPSLIN         ARKSpilsSetEpsLin

   FARKSPILSMASSINIT          ARKSpilsSetMassLinearSolver
   FARKSPILSSETMASSEPSLIN     ARKSpilsSetMassEpsLin

   FARKSPILSSETJAC            ARKSpilsSetJacTimes
   FARKSPILSSETPREC           ARKSpilsSetPreconditioner

   FARKSPILSSETMASS           ARKSpilsSetMassTimes
   FARKSPILSSETMASSPREC       ARKSpilsSetMassPreconditioner

   FARKODE                    ARKode, ARKodeGet*, and ARK*Get*
   FARKDKY                    ARKodeGetDky

   FARKGETERRWEIGHTS          ARKodeGetErrWeights
   FARKGETRESWEIGHTS          ARKodeGetResWeights
   FARKGETESTLOCALERR         ARKodeGetEstLocalErrors

   FARKFREE                   ARKodeFree
   ---------------------      --------------------------------


 The user-supplied functions, each listed with the corresponding interface
 function which calls it (and its type within ARKODE), are as follows:

   Fortran:           Interface Fcn:           ARKODE Type:
   -------------      ------------------       -----------------------
   FARKIFUN           FARKfi                   ARKRhsFn
   FARKEFUN           FARKfe                   ARKRhsFn
   FARKDJAC           FARKDenseJac             ARKDlsJacFn
   FARKDMASS          FARKDenseMass            ARKDlsMassFn
   FARKBJAC           FARKBandJac              ARKDlsJacFn
   FARKBMASS          FARKBandMass             ARKDlsMassFn
   FARKSPJAC          FARKSparseJac            ARKDlsJacFn
   FARKSPMASS         FARKSparseMass           ARKDlsMassFn
   FARKPSET           FARKPSet                 ARKSpilsPrecSetupFn
   FARKPSOL           FARKPSol                 ARKSpilsPrecSolveFn
   FARKMASSPSET       FARKMassPSet             ARKSpilsMassPrecSetupFn
   FARKMASSPSOL       FARKMassPSol             ARKSpilsMassPrecSolveFn
   FARKJTSETUP        FARKJTSetup              ARKSpilsJacTimesSetupFn
   FARKJTIMES         FARKJtimes               ARKSpilsJacTimesVecFn
   FARKMTSETUP        FARKMTSetup              ARKSpilsMassTimesSetupFn
   FARKMTIMES         FARKMtimes               ARKSpilsMassTimesVecFn
   FARKEWT            FARKEwt                  ARKEwtFn
   FARKADAPT          FARKAdapt                ARKAdaptFn
   FARKEXPSTAB        FARKExpStab              ARKExpStabFn
   -------------      ------------------       -----------------------

 In contrast to the case of direct use of ARKODE, and of most Fortran ODE
 solvers, the names of all user-supplied routines here are fixed, in
 order to maximize portability for the resulting mixed-language program.

 Important note on portability:  In this package, the names of the 
 interface functions, and the names of the Fortran user routines called by
 them, appear as dummy names which are mapped to actual values by a series 
 of definitions, in this and other header files.

 =============================================================================

                  Usage of the FARKODE Interface Package

 The usage of FARKODE requires calls to a variety of interface
 functions, depending on the method options selected, and one or more
 user-supplied routines which define the problem to be solved.  These
 function calls and user routines are summarized separately below.

 Some details are omitted, and the user is referred to the ARKODE user 
 documentation for more complete information.  Information on the arguments 
 of any given user-callable interface routine, or of a given user-supplied 
 function called by an interface function, can be found in the 
 documentation on the corresponding function in the ARKODE package.

 The number labels on the instructions below end with s for instructions
 that are specific to use with the serial/OpenMP/PThreads NVector package; 
 similarly those that end with p are specific to use with the N_VParallel package.

 -----------------------------------------------------------------------------

                               Data Types

 Throughout this documentation, we will refer to data types according to 
 their usage in SUNDIALS.  The equivalent types to these may vary, 
 depending on your computer architecture and on how SUNDIALS was compiled.  
 A Fortran user should take care that all arguments passed through this 
 Fortran/C interface are declared of the appropriate type.

 Integers: SUNDIALS uses 'int', 'long int' and 'sunindextype' types.  At 
 compilation, SUNDIALS allows the configuration of the 'index' type, that 
 accepts values of 32-bit signed and 64-bit signed.  This choice dictates 
 the size of a SUNDIALS 'sunindextype' variable.
   int      -- equivalent to an INTEGER or INTEGER*4 in Fortran
   long int -- equivalent to an INTEGER*8 in Fortran (Linux/UNIX/OSX), or 
               equivalent to an INTEGER in Windows
   sunindextype -- this will depend on the SUNDIALS configuration:
               32-bit -- equivalent to an INTEGER or INTEGER*4 in Fortran
               64-bit -- equivalent to an INTEGER*8 in Fortran
	      
 Real numbers:  At compilation, SUNDIALS allows the configuration option 
 '--with-precision', that accepts values of 'single', 'double' or 
 'extended' (the default is 'double').  This choice dictates the size of a 
 SUNDIALS 'realtype' variable.  The corresponding Fortran types for these 
 'realtype' sizes are:
   single   -- equivalent to a REAL or REAL*4 in Fortran
   double   -- equivalent to a DOUBLE PRECISION or REAL*8 in Fortran
   extended -- equivalent to a REAL*16 in Fortran

 -----------------------------------------------------------------------------

 (1) User-supplied right-hand side routines: FARKIFUN and FARKEFUN

     The user must in all cases supply at least one of the following 
     Fortran routines:

       SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)

     Sets the YDOT array to fi(T,Y), the implicit portion of the right-hand
     side of the ODE system, as function of time T and the state variable
     array Y.

       SUBROUTINE FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)

     Sets the YDOT array to fe(t,y), the explicit portion of the right-hand
     side of the ODE system, as function of time T and the state variable
     array Y.

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       YDOT -- array containing state derivatives [realtype, output]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.

 (2s) Optional user-supplied dense Jacobian approximation routine: FARKDJAC

     As an option when using the Dense or LapackDense linear solvers, the 
     user may supply a routine that computes a dense approximation of the 
     system Jacobian J = dfi(t,y)/dy.  If supplied, it must have the 
     following form:

       SUBROUTINE FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2,
      &                    WK3, IER)

     Typically this routine will use only NEQ, T, Y, and DJAC. It must 
     compute the Jacobian and store it column-wise in DJAC.

     The arguments are:
       NEQ  -- number of rows in the matrix [long int, input]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       DJAC -- 2D array containing the jacobian entries [realtype of size
               (NEQ,NEQ), output]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WK*  -- array containing temporary workspace of same size as Y
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.

 (2s) Optional user-supplied band Jacobian approximation routine: FARKBJAC

     As an option when using the Band or LapackBand linear solvers, the 
     user may supply a routine that computes a band approximation of the 
     system Jacobian J = dfi(t,y)/dy. If supplied, it must have the 
     following form:

       SUBROUTINE FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H,
      &                    IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. It
     must load the MDIM by N array BJAC with the Jacobian matrix at the
     current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element
     J(i,j)  with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.

     The arguments are:
       NEQ  -- number of rows in the matrix [long int, input]
       MU   -- upper half-bandwidth of the matrix [long int, input]
       ML   -- lower half-bandwidth of the matrix [long int, input]
       MDIM -- leading dimension of BJAC array [long int, input]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       BJAC -- 2D array containing the jacobian entries [realtype of size
               (MDIM,NEQ), output]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WK*  -- array containing temporary workspace of same size as Y
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.
 
 (2s) User-supplied sparse Jacobian approximation routine: FARKSPJAC

     Required when using the KLU or SuperLUMT linear solvers, the 
     user must supply a routine that computes a compressed-sparse-column 
     [or compressed-sparse-row] approximation of the system Jacobian 
     J = dfi(t,y)/dy.  This routine must have the following form:

       SUBROUTINE FARKSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS,
      &                     JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only M, N, NNZ, JDATA, JRVALS and 
     JCPTRS. It must load the N by N compressed sparse column [or 
     compressed sparse row] matrix with storage for NNZ nonzeros, stored 
     in the arrays JDATA (nonzero values), JRVALS (row [or column] indices 
     for each nonzero), JCOLPTRS (indices for start of each column 
     [or row]), with the Jacobian matrix at the current (t,y) in CSC [or 
     CSR] format (see sunmatrix_sparse.h for more information).

     The arguments are:
         T    -- current time [realtype, input]
         Y    -- array containing state variables [realtype, input]
         FY   -- array containing state derivatives [realtype, input]
         N    -- number of matrix rows/columns in Jacobian [long int, input]
         NNZ  -- allocated length of nonzero storage [long int, input]
        JDATA -- nonzero values in Jacobian
                 [realtype of length NNZ, output]
       JRVALS -- row [or column] indices for each nonzero in Jacobian
                 [long int of length NNZ, output]
       JCPTRS -- pointers to each Jacobian column [or row] in preceding 
                 arrays [long int of length N+1, output]
         H    -- current step size [realtype, input]
         IPAR -- array containing integer user data that was passed to
                 FARKMALLOC [long int, input]
         RPAR -- array containing real user data that was passed to
                 FARKMALLOC [realtype, input]
         WK*  -- array containing temporary workspace of same size as Y
                 [realtype, input]
         IER  -- return flag [int, output]:
                    0 if successful,
                   >0 if a recoverable error occurred,
                   <0 if an unrecoverable error ocurred.
 
     NOTE: this may ONLY be used if SUNDIALS has been configured with 
     sunindextype set to 64-bit integers.
 
 (2) Optional user-supplied Jacobian-vector product setup routine: 
     FARKJTSETUP

     As an option when using the ARKSpils iterative linear solver 
     interface, the user may supply a routine that computes the product of 
     the system Jacobian J = dfi(t,y)/dy and a given vector v, as well as 
     a routine to set up any user data structures in preparation for the 
     matrix-vector product.  If a 'setup' routine is supplied, it must 
     have the following form:

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
 
 (2) Optional user-supplied Jacobian-vector product routine: FARKJTIMES

     As an option when using the ARKSpils linear solver interface, the 
     user may supply a routine that computes the product of the system 
     Jacobian J = dfi(t,y)/dy and a given vector v.  If supplied, it 
     must have the following form:

       SUBROUTINE FARKJTIMES(V, JV, T, Y, FY, H, IPAR, RPAR, WORK, IER)

     Typically this routine will use only NEQ, T, Y, V, and JV.  It must
     compute the product vector J*v where v is stored in the vector V 
     and the result J*v is stored in JV.  

     The arguments are:
       V    -- array containing vector to multiply [realtype, input]
       JV   -- array containing product vector [realtype, output]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WORK -- array containing temporary workspace of same size as Y
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                  nonzero if an error.

 (3) Optional user-supplied preconditioner setup/solve routines: FARKPSET 
   and FARKPSOL

     As an option when using the ARKSPILS linear solver interface, the 
     user may supply routines to setup and apply the preconditioner.  
     If supplied, these must have the following form:

       SUBROUTINE FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,IER)

     This routine must set up the preconditioner P to be used in the 
     subsequent call to FARKPSOL.  The preconditioner (or the product of 
     the left and right preconditioners if using both) should be an 
     approximation to the matrix  A(T) = M(T) - GAMMA*J  
     (M(T) = mass matrix, J = Jacobian),

     The arguments are:
       T = current time [realtype, input]
       Y = current state variable array [realtype, input]
       FY = current state variable derivative array [realtype, input]
       JOK = flag indicating whether Jacobian-related data needs to be 
           recomputed [int, input]:
                  0 = recompute, 
                  1 = reuse with the current value of GAMMA
       JCUR = return flag to denote if Jacobian data was recomputed
           [realtype, output], 1=yes, 0=no
       GAMMA = Jacobian scaling factor [realtype, input]
       H = current time step [realtype, input]
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 >0 = recoverable failure
                 <0 = non-recoverable failure

     The user-supplied routine FARKPSOL must have the form:

       SUBROUTINE FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,IER)

     Typically this routine will use only T, Y, GAMMA, R, LR, and Z.  It
     must solve the preconditioner linear system Pz = r.  The preconditioner
     (or the product of the left and right preconditioners if both are 
     nontrivial) should be an approximation to the matrix  M(T) - GAMMA*J  
     (M(T) = mass matrix, J = Jacobian).

     The arguments are:
       T = current time [realtype, input]
       Y = current state variable array [realtype, input]
       FY = current state variable derivative array [realtype, input]
       R = right-hand side array [realtype, input]
       Z = solution array [realtype, output]
       GAMMA = Jacobian scaling factor [realtype, input]
       DELTA = desired residual tolerance [realtype, input]
       LR = flag denoting to solve the right or left preconditioner system
                  1 = left preconditioner
                  2 = right preconditioner
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 >0 = recoverable failure
                 <0 = non-recoverable failure

 (4s) Optional user-supplied dense mass matrix routine: FARKDMASS

     Required when using the Dense or LapackDense mass matrix linear 
     solvers, the user must supply a routine that computes the system mass 
     matrix M.  This routine must have the following form: 

       SUBROUTINE FARKDMASS(NEQ, T, DMASS, IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only NEQ, T, and DMASS. It must compute
     the mass matrix and store it column-wise in DMASS.

     The arguments are:
       NEQ   -- number of rows in the matrix [long int, input]
       T     -- current time [realtype, input]
       DMASS -- 2D array containing the mass matrix entries [realtype 
                of size (NEQ,NEQ), output]
       IPAR  -- array containing integer user data that was passed to
                FARKMALLOC [long int, input]
       RPAR  -- array containing real user data that was passed to
                FARKMALLOC [realtype, input]
       WK*   -- array containing temporary workspace of same size as Y 
                [realtype, input]
       IER   -- return flag [int, output]:
                   0 if successful, 
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.
 
 (4s) Optional user-supplied band mass matrix routine: FARKBMASS

     Required when using the Band or LapackBand mass matrix linear solvers, 
     the user must supply a routine that computes the system mass matrix 
     M.  This routine must have the following form: 

       SUBROUTINE FARKBMASS(NEQ, MU, ML, MDIM, T, BMASS, IPAR, 
      &                     RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only NEQ, MU, ML, T, and BMASS. It 
     must load the MDIM by N array BMASS with the system mass matrix at 
     the current (t) in band form.  Store in BMASS(k,j) the mass matrix 
     element M(i,j)  with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and 
     j = 1 ... N.

     The arguments are:
       NEQ   -- number of rows in the matrix [long int, input]
       MU    -- upper half-bandwidth of the matrix [long int, input]
       ML    -- lower half-bandwidth of the matrix [long int, input]
       MDIM  -- leading dimension of BMASS array [long int, input]
       T     -- current time [realtype, input]
       BMASS -- 2D array containing the mass matrix entries [realtype 
                of size (MDIM,NEQ), output]
       IPAR  -- array containing integer user data that was passed to
                FARKMALLOC [long int, input]
       RPAR  -- array containing real user data that was passed to
                FARKMALLOC [realtype, input]
       WK*   -- array containing temporary workspace of same size as Y 
                [realtype, input]
       IER   -- return flag [int, output]:
                   0 if successful, 
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.
 
 (4s) User-supplied sparse mass matrix routine: FARKSPMASS

     Required when using the KLU or SuperLUMT mass matrix linear solvers, 
     the user must supply a routine that computes a 
     compressed-sparse-column [or compressed-sparse-row] version of the 
     (possibly time-dependent) system mass matrix M(t).  If supplied, 
     it must have the following form:

       SUBROUTINE FARKSPMASS(T, N, NNZ, MDATA, MRVALS, MCPTRS,
      &                      IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only M, N, NNZ, MDATA, MRVALS and 
     MCPTRS. It must load the N by N compressed sparse column [or 
     compressed sparse row] matrix with storage for NNZ nonzeros, stored 
     in the arrays MDATA (nonzero values), MRVALS (row [or column] indices 
     for each nonzero), MCOLPTRS (indices for start of each column 
     [or row]), with the system mass matrix at the current (t) in CSC [or 
     CSR] format (see sundials_sparse.h for more information).

     The arguments are:
         T    -- current time [realtype, input]
         N    -- number of rows/columns in mass matrix [long int, input]
         NNZ  -- allocated length of nonzero storage [long int, input]
        MDATA -- nonzero values in mass matrix
                 [realtype of length NNZ, output]
       MRVALS -- row [or column] indices for each nonzero in mass matrix
                 [long int of length NNZ, output]
       MCPTRS -- pointers to each mass matrix column [or row] in preceding 
                 arrays [long int of length N+1, output]
         IPAR -- array containing integer user data that was passed to
                 FARKMALLOC [long int, input]
         RPAR -- array containing real user data that was passed to
                 FARKMALLOC [realtype, input]
         WK*  -- array containing temporary workspace of same size as Y
                 [realtype, input]
         IER  -- return flag [int, output]:
                    0 if successful,
                   >0 if a recoverable error occurred,
                   <0 if an unrecoverable error ocurred.
 
 (4) User-supplied mass-matrix-vector product routine: FARKMTSETUP

     When using the ARKSpils iterative mass matrix linear solver interface, 
     the user is required to supply a routine that computes the product of
     the system Jacobian J = dfi(t,y)/dy and a given vector v.  As an 
     option, the user may also supply a routine to set up any user data
     structures in preparation for the mass matrix-vector product.  If this 
     is supplied, it must have the following form:

       SUBROUTINE FARKMTSETUP(T, IPAR, RPAR, IER)

     It must perform any relevant preparations for subsequent calls to the 
     user-provided FARKMTIMES routine (see below).  

     The arguments are:
       T    -- current time [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.

 (4) User-supplied mass-matrix-vector product routine: FARKMTIMES

     Required when using the SP* linear solvers, the user should supply a
     routine that computes the product of the system mass matrix M and a
     given vector v.  It must have the following form:

       SUBROUTINE FARKMTIMES(V, MV, T, IPAR, RPAR, IER)

     Typically this routine will use only NEQ, T, V, and MV.  It must
     compute the product vector M*v where v is stored in the vector V
     and the result M*v is stored in MV.

     The arguments are:
       V    -- array containing vector to multiply [realtype, input]
       MV   -- array containing product vector [realtype, output]
       T    -- current time [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.
 
 (5) Optional user-supplied mass matrix preconditioner setup/solve 
   routines: FARKMASSPSET and FARKMASSPSOL

     As an option when using the SP* mass matrix linear solver, the user 
     may supply routines to setup and apply the preconditioner.  If 
     supplied, these must have the following form:

       SUBROUTINE FARKMASSPSET(T, IPAR, RPAR, IER)

     This routine must set up the preconditioner P to be used in the 
     subsequent call to FARKMASSPSOL.  The preconditioner (or the product 
     of the left and right preconditioners if using both) should be an 
     approximation to the system mass matrix M.

     The arguments are:
       T = current time [realtype, input]
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

     The user-supplied routine FARKMASSPSOL must have the form:

       SUBROUTINE FARKMASSPSOL(T, R, Z, DELTA, LR, IPAR, RPAR, IER)

       DIMENSION R(*), Z(*), IPAR(*), RPAR(*)

     Typically this routine will use only T, R, LR, and Z.  It
     must solve the preconditioner linear system Pz = r.  The 
     preconditioner (or the product of the left and right 
     preconditioners if both are nontrivial) should be an 
     approximation to the system mass matrix M.

     The arguments are:
       T = current time [realtype, input]
       R = right-hand side array [realtype, input]
       Z = solution array [realtype, output]
       DELTA = desired residual tolerance [realtype, input]
       LR = flag denoting to solve the right or left preconditioner system
                  1 = left preconditioner
		  2 = right preconditioner
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

 (6) Optional user-supplied error weight vector routine: FARKEWT
 
     As an option to providing the relative and absolute tolerances, the 
     user may supply a routine that computes the weights used in the WRMS 
     norms.  If supplied, it must have the following form:

       SUBROUTINE FARKEWT(Y, EWT, IPAR, RPAR, IER)

     It must store the error weights in EWT, given the current solution Y. 

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       EWT  -- array containing the error weight vector [realtype, output]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                  nonzero if an error.

 (7) Optional user-supplied error weight vector routine: FARKADAPT
 
     As an option to providing the time step adaptivity, the user
     may supply a routine that computes the new time step for ARKode to 
     use, based on information about the three most recent errors and 
     previous time steps taken by the solver. If supplied, it must have 
     the following form:

       SUBROUTINE FARKADAPT(Y, T, H1, H2, H3, E1, E2, E3, Q, P, HNEW, IPAR, RPAR, IER)

     It must output the new time step in HNEW. 

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       T    -- current time [realtype, input]
       H1   -- current step size [realtype, input]
       H2   -- previous step size [realtype, input]
       H3   -- previous-previous step size [realtype, input]
       E1   -- estimated temporal error in current step [realtype, input]
       E2   -- estimated temporal error in previous step [realtype, input]
       E3   -- estimated temporal error in previous-previous step 
               [realtype, input]
       Q    -- global order of accuracy for RK method [int, input]
       P    -- global order of accuracy for RK embedding [int, input]
       HNEW -- predicted next step size [realtype, output]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.

 (8) Optional user-supplied explicitly stable time step routine: FARKEXPSTAB
 
     As an option, the user may provide a routine to return the maximum 
     stable time step size for the explicit ODE RHS function.  If supplied, 
     it must have the following form:

       SUBROUTINE FARKEXPSTAB(Y, T, HSTAB, IPAR, RPAR, IER)

     It must store the explicitly stable step size in the ouptut HSTAB.

     The arguments are:
       Y     -- array containing state variables [realtype, input]
       T     -- current time [realtype, input]
       HSTAB -- explicitly-stable step size [realtype, output]
       IPAR  -- array containing integer user data that was passed to
                FARKMALLOC [long int, input]
       RPAR  -- array containing real user data that was passed to
                FARKMALLOC [realtype, input]
       IER   -- return flag [int, output]:
                  0 if successful,
                  nonzero if an error.

 -----------------------------------------------------------------------------

 (9) Initialization:  FNVINITS / FNVINITP / FNVINITOMP / FNVINITPTS, 
                      FSUNBANDMATINIT / FSUNDENSEMATINIT / FSUNSPARSEMATINIT,
                      FSUNBANDMASSMATINIT / FSUNDENSEMASSMATINIT /
                         FSUNSPARSEMASSMATINIT,
                      FSUNBANDLINSOLINIT / FSUNDENSELINSOLINIT / 
                         FSUNKLUINIT / FSUNKLUREINIT / FSUNKLUSETORDERING / 
                         FSUNLAPACKBANDINIT / FSUNLAPACKDENSEINIT / 
                         FSUNPCGINIT / FSUNSPBCGSINIT / FSUNSPFGMRINIT / 
                         FSUNSPGMRINIT / FSUNSPTFQMRINIT / FSUNSUPERLUMTINIT /
                         FSUNSUPERLUMTSETORDERING,
                      FSUNMASSBANDLINSOLINIT / FSUNMASSDENSELINSOLINIT /
                         FSUNMASSKLUINIT / FSUNMASSKLUREINIT / 
                         FSUNMASSKLUSETORDERING / FSUNMASSLAPACKBANDINIT / 
                         FSUNMASSLAPACKDENSEINIT / FSUNMASSPCGINIT / 
                         FSUNMASSSPBCGSINIT / FSUNMASSSPFGMRINIT / 
                         FSUNMASSSPGMRINIT / FSUNMASSSPTFQMRINIT / 
                         FSUNMASSSUPERLUMTINIT / FSUNMASSSUPERLUMTSETORDERING,
                      FARKMALLOC,
                      FARKDLSINIT / FARKSPILSINIT,
                      FARKREINIT, FARKRESIZE

     NOTE: the initialization order is important!  It *must* proceed as 
     shown: vector, matrix (if used), linear solver (if used), ARKode, 
     ARKDls/ARKSpils, reinit/resize.
 
 (9.1s) To initialize the a vector specification for storing the solution 
     data, the user must make one of the following calls:

       (serial)   
          CALL FNVINITS(4, NEQ, IER)
       (MPI parallel)
          CALL FNVINITP(COMM, 4, NLOCAL, NGLOBAL, IER)
       (OpenMP threaded)
          CALL FNVINITOMP(4, NEQ, NUM_THREADS, IER)
       (PThreads threaded)
          CALL FNVINITPTS(4, NEQ, NUM_THREADS, IER)

     In each of these, one argument is an int containing the ARKODE solver 
     ID (4). 

     The other arguments are:
        NEQ = size of vectors [long int, input]
        COMM = the MPI communicator [int, input]
        NLOCAL = local size of vectors on this processor 
           [long int, input]
        NGLOBAL = the system size, and the global size of vectors (the sum 
           of all values of NLOCAL) [long int, input]
        NUM_THREADS = number of threads
        IER = return completion flag [int, output]:
                  0 = success, 
                 -1 = failure.

 (9.2) To initialize a band/dense/sparse matrix structure for 
     storing the system Jacobian and for use within a direct linear solver,
     the user must make one of the following calls:
 
          CALL FSUNBANDMATINIT(4, N, MU, ML, SMU, IER)
          CALL FSUNDENSEMATINIT(4, M, N, IER)
          CALL FSUNSPARSEMATINIT(4, M, N, NNZ, SPARSETYPE, IER)

     In each of these, one argument is an int containing the ARKODE solver 
     ID (4). 

     The other arguments are:

        M = the number of rows of the matrix [long int, input]
        N = the number of columns of the matrix [long int, input]
        MU = the number of upper bands (diagonal not included) in a banded 
           matrix [long int, input]
        ML = the number of lower bands (diagonal not included) in a banded 
           matrix [long int, input]
        SMU = the number of upper bands to store (diagonal not included) 
           for factorization of a banded matrix [long int, input]
        NNZ = the storage size (upper bound on the number of nonzeros) for 
           a sparse matrix [long int, input]
        SPARSETYPE = integer denoting use of CSC (0) vs CSR (1) storage 
           for a sparse matrix [int, input]
        IER = return completion flag [int, output]:
                  0 = success, 
                 -1 = failure.

 (9.3) To initialize a band/dense/sparse matrix structure for 
     storing the mass matrix and for use within a direct mass matrix linear 
     solver, the user must make one of the following calls:
  
          CALL FSUNBANDMASSMATINIT(N, MU, ML, SMU, IER)
          CALL FSUNDENSEMASSMATINIT(M, N, IER)
          CALL FSUNSPARSEMASSMATINIT(M, N, NNZ, SPARSETYPE, IER)

     The arguments have the same meanings as with the Jacobian matrix 
     initialization routines above.

 (9.4) To initialize a linear solver structure for solving linear systems 
     arising from implicit or IMEX treatment of the IVP, the user must make 
     one of the following calls:

          CALL FSUNBANDLINSOLINIT(4, IER)
          CALL FSUNDENSELINSOLINIT(4, IER)
          CALL FSUNKLUINIT(4, IER)
          CALL FSUNLAPACKBANDINIT(4, IER)
          CALL FSUNLAPACKDENSEINIT(4, IER)
          CALL FSUNPCGINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPBCGSINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPFGMRINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPGMRINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSPTFQMRINIT(4, PRETYPE, MAXL, IER)
          CALL FSUNSUPERLUMTINIT(4, NUM_THREADS, IER)

     Or once these have been initialized, their solver parameters may be
     modified via calls to the functions

          CALL FSUNKLUSETORDERING(4, ORD_CHOICE, IER)
          CALL FSUNSUPERLUMTSETORDERING(4, ORD_CHOICE, IER)

          CALL FSUNPCGSETPRECTYPE(4, PRETYPE, IER)
          CALL FSUNPCGSETMAXL(4, MAXL, IER)
          CALL FSUNSPBCGSSETPRECTYPE(4, PRETYPE, IER)
          CALL FSUNSPBCGSSETMAXL(4, MAXL, IER)
          CALL FSUNSPFGMRSETGSTYPE(4, GSTYPE, IER)
          CALL FSUNSPFGMRSETPRECTYPE(4, PRETYPE, IER)
          CALL FSUNSPGMRSETGSTYPE(4, GSTYPE, IER)
          CALL FSUNSPGMRSETPRECTYPE(4, PRETYPE, IER)
          CALL FSUNSPTFQMRSETPRECTYPE(4, PRETYPE, IER)
          CALL FSUNSPTFQMRSETMAXL(4, MAXL, IER)

     In all of the above, one argument is an int containing the ARKODE solver 
     ID (4). 

     The other arguments are:

        NNZ = the storage size (upper bound on the number of nonzeros) for 
           a sparse matrix [long int, input]
        ORD_CHOICE = integer denoting ordering choice (see 
           SUNKLUSetOrdering and SUNSuperLUMTSetOrdering documentation 
           for details) [int, input]
        PRETYPE = type of preconditioning to perform (0=none, 1=left, 
           2=right, 3=both) [int, input]
        MAXL = maximum Krylov subspace dimension [int, input]
        GSTYPE = choice of Gram-Schmidt orthogonalization algorithm 
           (0=modified, 1=classical) [int, input]
        IER = return completion flag [int, output]:
                  0 = success, 
                 -1 = failure.

 (9.5) To initialize a linear solver structure for solving linear systems 
     arising from use of a non-identity mass matrix, the user must make 
     one of the following calls:

          CALL FSUNMASSBANDLINSOLINIT(IER)
          CALL FSUNMASSDENSELINSOLINIT(IER)
          CALL FSUNMASSKLUINIT(IER)
          CALL FSUNMASSLAPACKBANDINIT(IER)
          CALL FSUNMASSLAPACKDENSEINIT(IER)
          CALL FSUNMASSPCGINIT(PRETYPE, MAXL, IER)
          CALL FSUNMASSSPBCGSINIT(PRETYPE, MAXL, IER)
          CALL FSUNMASSSPFGMRINIT(PRETYPE, MAXL, IER)
          CALL FSUNMASSSPGMRINIT(PRETYPE, MAXL, IER)
          CALL FSUNMASSSPTFQMRINIT(PRETYPE, MAXL, IER)
          CALL FSUNMASSSUPERLUMTINIT(NUM_THREADS, IER)
          
     Or once these have been initialized, their solver parameters may be
     modified via calls to the functions

          CALL FSUNMASSKLUSETORDERING(ORD_CHOICE, IER)
          CALL FSUNMASSSUPERLUMTSETORDERING(ORD_CHOICE, IER)

          CALL FSUNMASSPCGSETPRECTYPE(PRETYPE, IER)
          CALL FSUNMASSPCGSETMAXL(MAXL, IER)
          CALL FSUNMASSSPBCGSSETPRECTYPE(PRETYPE, IER)
          CALL FSUNMASSSPBCGSSETMAXL(MAXL, IER)
          CALL FSUNMASSSPFGMRSETGSTYPE(GSTYPE, IER)
          CALL FSUNMASSSPFGMRSETPRECTYPE(PRETYPE, IER)
          CALL FSUNMASSSPGMRSETGSTYPE(GSTYPE, IER)
          CALL FSUNMASSSPGMRSETPRECTYPE(PRETYPE, IER)
          CALL FSUNMASSSPTFQMRSETPRECTYPE(PRETYPE, IER)
          CALL FSUNMASSSPTFQMRSETMAXL(MAXL, IER)

     The arguments have the same meanings as with the Jacobian matrix 
     initialization routines above.

 (9.6) To set various problem and solution parameters and allocate
     internal memory, make the following call:

       CALL FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT, ROUT,
      &                IPAR, RPAR, IER)

     The arguments are:
        T0 = initial value of t [realtype, input]
        Y0 = array of initial conditions [realtype, input]
        IMEX = flag denoting basic integration method [int, input]: 
                  0 = implicit, 
                  1 = explicit, 
                  2 = imex
        IATOL = type for absolute tolerance ATOL [int, input]:
                  1 = scalar,
                  2 = array,
                  3 = user-supplied function; the user must supply a routine
                      FARKEWT to compute the error weight vector.
        RTOL = scalar relative tolerance [realtype, input]
        ATOL = scalar or array absolute tolerance [realtype, input]
        IOUT = array of length at least 29 for integer optional outputs
               [long int, output]
        ROUT = array of length 6 for real optional outputs [realtype, output]
        IPAR = array of user integer data [long int, input/output]
        RPAR = array with user real data [realtype, input/output]
        IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

     The user data arrays IPAR and RPAR are passed unmodified to all 
     subsequent calls to user-provided routines. Modifications to either 
     array inside a user-provided routine will be propagated. Using these 
     two arrays, the user can dispense with COMMON blocks to pass data 
     betwen user-provided routines. 
 
     The optional outputs are:
           LENRW   = IOUT( 1) from ARKodeGetWorkSpace
           LENIW   = IOUT( 2) from ARKodeGetWorkSpace
           NST     = IOUT( 3) from ARKodeGetNumSteps
           NST_STB = IOUT( 4) from ARKodeGetNumExpSteps
           NST_ACC = IOUT( 5) from ARKodeGetNumAccSteps
           NST_ATT = IOUT( 6) from ARKodeGetNumStepAttempts
           NFE     = IOUT( 7) from ARKodeGetNumRhsEvals (fe(t,y) calls)
           NFI     = IOUT( 8) from ARKodeGetNumRhsEvals (fi(t,y) calls)
           NSETUPS = IOUT( 9) from ARKodeGetNumLinSolvSetups
           NETF    = IOUT(10) from ARKodeGetNumErrTestFails
           NNI     = IOUT(11) from ARKodeGetNumNonlinSolvIters
           NCFN    = IOUT(12) from ARKodeGetNumNonlinSolvConvFails
           NGE     = IOUT(13) from ARKodeGetNumGEvals

           H0U     = ROUT( 1) from ARKodeGetActualInitStep
           HU      = ROUT( 2) from ARKodeGetLastStep
           HCUR    = ROUT( 3) from ARKodeGetCurrentStep
           TCUR    = ROUT( 4) from ARKodeGetCurrentTime
           TOLSF   = ROUT( 5) from ARKodeGetTolScaleFactor
           UROUND  = ROUT( 6) from UNIT_ROUNDOFF
     See the ARKODE manual for details.

 (9.7) If a direct linear solver was created in step (7.4) then it must be 
     attached to ARKode.  If the user called any one of FSUNBANDLINSOLINIT, 
     FSUNDENSELINSOLINIT, FSUNKLUINIT, FSUNLAPACKBANDINIT, 
     FSUNLAPACKDENSEINIT, or FSUNSUPERLUMTINIT, then this must be 
     attached to the ARKDLS interface using the command:

       CALL FARKDLSINIT(IER)

     The arguments are:
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

 (9.7) If an iterative linear solver was created in step (7.4) then it must 
     be attached to ARKode.  If the user called any one of FSUNPCGINIT, 
     FSUNSPBCGSINIT, FSUNSPFGMRINIT, FSUNSPGMRINIT, or FSUNSPTFQMRINIT, 
     then this must be attached to the ARKSPILS interface using the command:

       CALL FARKSPILSINIT(IER)

     The arguments are:
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

 (9.8) If a mass matrix linear solver was created in step (7.5) then it 
     must be attached to ARKode.  If the user called any one of 
     FSUNMASSBANDLINSOLINIT, FSUNMASSDENSELINSOLINIT, 
     FSUNMASSKLUINIT, FSUNMASSLAPACKBANDINIT, FSUNMASSLAPACKDENSEINIT, 
     or FSUNMASSSUPERLUMTINIT, then this must be attached to the 
     ARKDLS interface using the command:

       CALL FARKDLSMASSINIT(TIME_DEP, IER)

     The arguments are:
	TIME_DEP = flag indicating whether the mass matrix is 
           time-dependent (1) or not (0) [int, output]:
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

 (9.9) If a mass matrix linear solver was created in step (7.5) then it 
     must be attached to ARKode.  If the user called any one of 
     FSUNMASSPCGINIT, FSUNMASSSPBCGSINIT, FSUNMASSSPFGMRINIT, 
     FSUNMASSSPGMRINIT, or FSUNMASSSPTFQMRINIT, then this must be attached
     to the ARKDLS interface using the command:

       CALL FARKSPILSMASSINIT(TIME_DEP, IER)

     The arguments are:
	TIME_DEP = flag indicating whether the mass matrix is 
           time-dependent (1) or not (0) [int, output]:
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

 (9.10) If the user program includes the FARKEWT routine for the evaluation 
     of the error weights, the following call must be made

       CALL FARKEWTSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKEWT is provided and 
     should be used; FLAG = 0 resets to the default EWT formulation.  The 
     int return flag IER is 0 if successful, and nonzero otherwise.

 (9.11) If the user program includes the FARKADAPT routine for performing 
     step adaptivity, the following call must be made

       CALL FARKADAPTSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKADAPT is provided 
     and should be used; FLAG = 0 resets to the default adaptivity 
     formulation. The int return flag IER is 0 if successful, and nonzero 
     otherwise.

 (9.12) If the user program includes the FARKEXPSTAB routine for 
     calculation of the maximum explicitly stable step size, the following 
     call must be made

       CALL FARKEXPSTABSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKEXPSTAB is provided 
     and should be used; FLAG = 0 resets to the default explicit stability 
     formulation.  The int return flag IER is 0 if successful, and nonzero 
     otherwise.

 (9.13) If the user program includes the FARKBJAC routine for the 
     evaluation of the band approximation to the Jacobian, then following 
     the call to FARKDLSINIT, the following call must be made 

       CALL FARKBANDSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKBJAC is provided and should be 
     used; FLAG=0 specifies a reset to the internal finite difference 
     Jacobian approximation.  The int return flag IER=0 if successful, 
     nonzero otherwise.
 
     If the user program includes the FARKDJAC routine for the evaluation 
     of the dense approximation to the Jacobian, then after the call to 
     FARKDLSINIT, the following call must be made 

       CALL FARKDENSESETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKDJAC is provided and should be 
     used; FLAG=0 specifies a reset to the internal finite difference 
     Jacobian approximation.  The int return flag IER=0 if successful, and 
     nonzero otherwise.
 
     When using a sparse matrix and linear solver the user must provide the
     FARKSPJAC routine for the evaluation of the sparse approximation to 
     the Jacobian.  To indicate that this routine has been provided, after 
     the call to FARKDLSINIT, the following call must be made 

       CALL FARKSPARSESETJAC(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.

 (9.14) If the user program includes the FARKJTSETUP and FARKJTIMES 
     routines for setup of a Jacobian-times-vector product (for use with 
     the ARKSpils interface), then after creating the ARKSpils interface, 
     the following call must be made:

       CALL FARKSPILSSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKJTSETUP and FARKJTIMES are 
     provided and should be used; FLAG=0 specifies a reset to the internal 
     finite difference approximation to this product).  The int return 
     flag IER=0 if successful, and nonzero otherwise.
 
 (9.16) If the user program includes the FARKPSET and FARKPSOL routines 
     for supplying a preconditioner to an iterative linear solver, then 
     after creating the ARKSpils interface, the following call must be made

       CALL FARKSPILSSETPREC(FLAG, IER)

     with the int FLAG=1.  If FLAG=0 then preconditioning with these 
     routines will be disabled. The return flag IER=0 if successful, 
     nonzero otherwise.

 (9.17) If the user wishes to use one of ARKode's built-in preconditioning 
     modules, FARKBP or FARKBBD, then that should be initialized after 
     creating the ARKSpils interface using one of the calls

       CALL FARKBPINIT(NEQ, MU, ML, IER)
       CALL FARKBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)

     Detailed explanation of the inputs to these functions, as well as any 
     requirements of user-supplied functions on which these preconditioning 
     modules rely, may be found in the header files for each module, 
     farkbp.h or farkbbd.h, respectively.

 (9.18) When using a band / dense / sparse mass matrix and 
     corresponding linear solver the user must provide the FARKBMASS / 
     FARKDMASS / FARKSPMASS routine for the evaluation of the 
     approximation to the mass matrix.  To indicate that the 
     appropriate routine has been provided, after the call to 
     FARKDLSMASSINIT, one of the following calls must be made 

       CALL FARKBANDSETMASS(IER)
       CALL FARKDENSESETMASS(IER)
       CALL FARKSPARSESETMASS(IER)

     The int return flag IER=0 if successful, nonzero otherwise.


 (9.19) When using the ARKSPILS interface for the mass matrix solvers, the 
     user must supply a mass matrix-times-vector routine FARKMTIMES, and 
     an associated 'setup' routine FARKMTSETUP.  After creating the 
     ARKSpilsMass interface, the user must make the following call to 
     signal that FARKMTSETUP and FARKMTIMES have been provided

       CALL FARKSPILSSETMASS(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.

 (9.20) If the user program includes the FARKMASSPSET and FARKMASSPSOL 
     routines for supplying a preconditioner to an iterative mass matrix 
     linear solver, then after creating the ARKSpilsMass interface, the 
     following call must be made

       CALL FARKSPILSSETMASSPREC(FLAG, IER)

     with the int FLAG=1.  If FLAG=0 then preconditioning with these 
     routines will be disabled. The return flag IER=0 if successful, 
     nonzero otherwise.

 (9.21) To re-initialize the ARKODE solver for the solution of a new problem
     of the same size as one already solved, make the following call:

       CALL FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)

     The arguments have the same names and meanings as those of FARKMALLOC. 
     FARKREINIT performs the same initializations as FARKMALLOC, but does 
     no memory allocation, using instead the existing internal memory 
     created by the previous FARKMALLOC call.  The subsequent calls to 
     attach the linear system or mass matrix system solvers are only needed 
     if those objects have been re-created.
 
 (9.22) To re-initialize the ARKODE solver for the solution of a new problem
     of a different size as one already solved, but with the same dynamical 
     time scale and method choice, make the following call:

       CALL FARKRESIZE(T0, Y0, HSCALE, ITOL, RTOL, ATOL, IER)

     The arguments are:
        T0 = initial value of t [realtype, input]
	Y0 = array of initial conditions [realtype, input]
	HSCALE = desired step size scale factor [realtype, input]
	          1.0 is the default
		  any value <= 0.0 results in the default.
        ITOL = flag denoting that a new relative tolerance and vector of
	   absolute tolerances are supplied in the RTOL and ATOL arguments
	   [int, input]:
                  0 = retain the current relative tolerance and current
		      scalar-valued or user-supplied function
                  1 = RTOL contains the new scalar-valued relative tolerance
                      and ATOL contains a new array of absolute tolerances
	RTOL = scalar-valued relative tolerance [realtype, input]
	ATOL = array of absolute tolerances [realtype, input]
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

     FARKRESIZE performs the opposite set of of operations as FARKREINIT: 
     it does not reinitialize any of the time-step heuristics, but it does 
     perform memory reallocation.  

     Since the problem has changed size, it is likely that the linear 
     solver object(s) must also be reconstructed.  The previous solver 
     object(s) may be freed and new ones created by calling the respective 
     FSUN***INIT routine again with modified arguments.  This should be 
     performed _after_ FARKRESIZE has been called.  Following 
     reconstruction of the linear solver object(s), the ARKSPILS or ARKDLS 
     interface must also be reconstructed, and any options (e.g. setting
     the Jacobian-times-vector product routine) must be specified again.

 (9.23) The SUNKLU solver will reuse much of the factorization information 
     from one solve to the next.  If at any time the user wants to force a 
     full refactorization or if the number of nonzeros in the Jacobian 
     matrix changes, the user should make the call

        CALL FSUNKLUREINIT(4, NNZ, REINIT_TYPE, IER)

     Similarly, if a user wants to force a full refactorization of the mass
     matrix, or if the number of nonzeros in the mass matrix changes, the 
     user should make the call

        CALL FSUNMASSKLUREINIT(NNZ, REINIT_TYPE, IER)

     The arguments are:
        NNZ = the maximum number of nonzeros [int; input]
        REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be 
          destroyed and a new one will be allocated with NNZ nonzeros.  
          For a value of 2, only symbolic and numeric factorizations will 
          be completed. 
 
 (9.24) To set various integer optional inputs, make the folowing call:

       CALL FARKSETIIN(KEY, VALUE, IER)

     to set the integer value VALUE to the optional input specified by the 
     quoted character string KEY. VALUE must be a Fortran integer of size 
     commensurate with a C "long int".  KEY must be one of the following: 
     ORDER, DENSE_ORDER, LINEAR, NONLINEAR, FIXEDPOINT, NEWTON, EXPLICIT, 
     IMPLICIT, IMEX, IRK_TABLE_NUM, ERK_TABLE_NUM, ARK_TABLE_NUM (pass in 
     an int array of length 2, implicit method first), MAX_NSTEPS, 
     HNIL_WARNS, PREDICT_METHOD, MAX_ERRFAIL, MAX_CONVFAIL, MAX_NITERS, 
     ADAPT_SMALL_NEF or LSETUP_MSBP.  The int return flag IER is 0 if 
     successful, and nonzero otherwise. 

 (9.25) To set various real optional inputs, make the following call:

       CALL FARKSETRIN(KEY, VALUE, IER)

     to set the realtype value VALUE to the optional input specified by the 
     quoted character string KEY.  VALUE must be a Fortran real-valued 
     number of size commensurate with the SUNDIALS "realtype".  KEY must 
     one of the following: INIT_STEP, MAX_STEP, MIN_STEP, STOP_TIME, 
     NLCONV_COEF, ADAPT_CFL, ADAPT_SAFETY, ADAPT_BIAS, ADAPT_GROWTH, 
     ADAPT_BOUNDS (pass in a realtype array of length 2), ADAPT_ETAMX1, 
     ADAPT_ETAMXF, ADAPT_ETACF, NONLIN_CRDOWN, NONLIN_RDIV, LSETUP_DGMAX, 
     or FIXED_STEP.  The int return flag IER is 0 if successful, and nonzero 
     otherwise.

 (9.26) To set the time step adaptivity method (and its associated 
     parameters), make the following call: 

       CALL FARKSETADAPTIVITYMETHOD(IMETHOD, IDEFAULT, IPQ, PARAMS, IER)

     The arguments are:
       IMETHOD  = the adaptivity method to use [integer, input]
       IDEFAULT = flag to use (1) or not (0) the default adaptivity 
                  parameters [integer, input]
       IPQ      = flag to use the embedding order p (0) or the method order 
                  q (1) for error-based step adaptivity [integer, input]
       PARAMS   = if IDEFAULT=0, this should be a realtype array of length 
                  2 containing the custom adaptivity parameters to use in 
		  the method [realtype, input].
       IER      = integer error flag (0 = success, 1 = failure) [integer, output]

 (9.27) To reset all optional inputs to their default values, make the 
     following call:

       CALL FARKSETDEFAULTS(IER)

 (9.28) To set a custom explicit Runge-Kutta table, make the following call:

       CALL FARKSETERKTABLE(S, Q, P, C, A, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       A = array of length S*S containing the ERK coefficients (stored in
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

     To set a custom diagonally-implicit Runge-Kutta table, make the 
     following call:

       CALL FARKSETIRKTABLE(S, Q, P, C, A, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       A = array of length S*S containing the DIRK coefficients (stored in
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

     To set a custom additive Runge-Kutta table, make the following call:

       CALL FARKSETARKTABLES(S, Q, P, CI, CE, AI, AE, BI, BE, B2I, B2E, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       CI = array of length S containing the implicit stage times
           [realtype, input]
       CE = array of length S containing the explicit stage times
           [realtype, input]
       AI = array of length S*S containing the DIRK coefficients (stored in
           row-major, "C", order) [realtype, input]
       AE = array of length S*S containing the ERK coefficients (stored in
           row-major, "C", order) [realtype, input]
       BI = array of length S containing the implicit solution coefficients
           [realtype, input]
       BE = array of length S containing the explicit solution coefficients
           [realtype, input]
       B2I = array of length S containing the implicit embedding coefficients
           [realtype, input]
       B2E = array of length S containing the explicit embedding coefficients
           [realtype, input]

 (9.29) When using a non-identity mass matrix, to set an absolute residual 
     tolerance (scalar or vector), call:

       CALL FARKSETRESTOLERANCE(IATOL, ATOL, IER)

     The arguments are:
       IATOL = type for absolute tolerance ATOL [int, input]:
                 1 = scalar,
                 2 = array
	ATOL = scalar or array absolute residual tolerance [realtype, input]
	IER  = return completion flag [int, output]:
                 0 = SUCCESS,
                -1 = failure (see printed message for failure details).

 (9.30) To set a solver diagnostics output file, make the folowing call:

       CALL FARKSETDIAGNOSTICS(FNAME, FLEN, IER)

     The desired diagnostics filename should be supplied by the 
     quoted character string FNAME.  The integer argument FLEN should 
     contain the length (in characters) of FNAME (for portability).  The 
     int return flag IER is 0 if successful (able to open file), and 
     nonzero otherwise.

 (9.31) To close the solver diagnostics output file, make the folowing call:

       CALL FARKSTOPDIAGNOSTICS(IER)

     The int return flag IER is 0 if successful (able to close file), and
     nonzero otherwise.


 -----------------------------------------------------------------------------

 (10) Optional outputs from DLS and SPILS linear solvers (stored in the 
     IOUT array that was passed to FARKMALLOC)

     Optional outputs specific to the ARKDLS interface:
        LENRWLS  = IOUT(14) from ARKDlsGetWorkSpace (realtype space)
        LENIWLS  = IOUT(15) from ARKDlsGetWorkSpace (integer space)
        LSTF     = IOUT(16) from ARKDlsGetLastFlag
        NFELS    = IOUT(17) from ARKDlsGetNumRhsEvals
        NJED     = IOUT(18) from ARKDlsGetNumJacEvals

     Optional outputs specific to the ARKDLSMASS interface:
        LENRWMS  = IOUT(23) from ARKDlsGetMassWorkSpace (realtype space)
        LENIWMS  = IOUT(24) from ARKDlsGetMassWorkSpace (integer space)
        LSTMF    = IOUT(25) from ARKDlsGetLastMassFlag
        NMSETUP  = IOUT(26) from ARKDlsGetNumMassSetups
        NMSOLVES = IOUT(27) from ARKDlsGetNumMassSolves
        NMMULTS  = IOUT(28) from ARKDlsGetNumMassMult

     Optional outputs specific to the ARKSPILS interface:
        LENRWLS  = IOUT(14) from ARKSpilsGetWorkSpace
        LENIWLS  = IOUT(15) from ARKSpilsGetWorkSpace
        LSTF     = IOUT(16) from ARKSpilsGetLastFlag
        NFELS    = IOUT(17) from ARKSpilsGetNumRhsEvals
        NJTV     = IOUT(18) from ARKSpilsGetNumJtimesEvals
        NPE      = IOUT(19) from ARKSpilsGetNumPrecEvals
        NPS      = IOUT(20) from ARKSpilsGetNumPrecSolves
        NLI      = IOUT(21) from ARKSpilsGetNumLinIters
        NCFL     = IOUT(22) from ARKSpilsGetNumConvFails
 
     Optional outputs specific to the ARKSPILSMASS interface:
        LENRWMS  = IOUT(23) from ARKSpilsGetMassWorkSpace
        LENIWMS  = IOUT(24) from ARKSpilsGetMassWorkSpace
        LSTMF    = IOUT(25) from ARKSpilsGetLastMassFlag
        NMPE     = IOUT(26) from ARKSpilsGetNumMassPrecEvals
        NMPS     = IOUT(27) from ARKSpilsGetNumMassPrecSolves
        NMLI     = IOUT(28) from ARKSpilsGetNumMassIters
        NMCFL    = IOUT(29) from ARKSpilsGetNumMassConvFails
 
     See the ARKODE manual for more detailed descriptions of any of the 
     above.

 -----------------------------------------------------------------------------

 (11) The integrator: FARKODE

     Carrying out the integration is accomplished by making calls as follows:

       CALL FARKODE(TOUT, T, Y, ITASK, IER)

     The arguments are:
       TOUT = next value of t at which a solution is desired [realtype, input]
       T = value of t reached by the solver [realtype, output]
       Y = array containing state variables on output [realtype, output]
       ITASK = task indicator [int, input]:
                  1 = normal mode (overshoot TOUT and interpolate)
                  2 = one-step mode (return after each internal step taken)
                  3 = normal tstop mode (like 1, but integration never 
                      proceeds past TSTOP, which must be specified through a 
                      call to FARKSETRIN using the key 'STOP_TIME')
                  4 = one step tstop (like 2, but integration never goes 
                      past TSTOP)
       IER = completion flag [int, output]: 
                  0 = success, 
                  1 = tstop return, 
                  2 = root return, 
                  values -1 ... -10 are failure modes (see ARKODE manual).
     The current values of the optional outputs are immediately available in
     the IOUT and ROUT arrays.
 
 -----------------------------------------------------------------------------

 (12) Computing solution derivatives: FARKDKY

     To obtain a derivative of the solution, of order up to the method order,
     make the following call:

       CALL FARKDKY(T, K, DKY, IER)

     The arguments are:
       T = time at which solution derivative is desired, within the interval
           [TCUR-HU,TCUR], [realtype, input].
       K = derivative order (0 .le. K .le. QU) [int, input]
       DKY = array containing computed K-th derivative of y [realtype, output]
       IER = return flag [int, output]: 0=success, <0 = illegal argument.
 
 -----------------------------------------------------------------------------

 (13) Get the current error weight vector: FARKGETERRWEIGHTS

     To obtain the current error weight vector, make the following call:

       CALL FARKGETERRWEIGHTS(EWT, IER)

     The arguments are:
       EWT = array containing the error weight vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (14) Get the current residual weight vector: FARKGETRESWEIGHTS

     To obtain the current residual weight vector, make the following call:

       CALL FARKGETRESWEIGHTS(RWT, IER)

     The arguments are:
       RWT = array containing the residual weight vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (15) Get an estimate of the local error: FARKGETESTLOCALERR

     To obtain the current error estimate vector, make the following call:

       CALL FARKGETESTLOCALERR(ELE, IER)

     The arguments are:
       ELE = array with the estimated local error vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (16) Memory freeing: FARKFREE 

     To free the internal memory created by the calls to FARKMALLOC, 
     FARKDLSINIT/FARKSPILSINIT, the generic linear solver and matrix modules, 
     and FNVINIT*, make the call

       CALL FARKFREE()
 
===============================================================*/

#ifndef _FARKODE_H
#define _FARKODE_H

/* header files  */
#include <arkode/arkode.h>
#include <sundials/sundials_linearsolver.h>  /* definition of type SUNLinearSolver */
#include <sundials/sundials_matrix.h>        /* definition of type SUNMatrix */
#include <sundials/sundials_nvector.h>       /* definition of type N_Vector */
#include <sundials/sundials_types.h>         /* definition of type realtype */

/*=============================================================*/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_IMP_FUN             SUNDIALS_F77_FUNC(farkifun,                FARKIFUN)
#define FARK_EXP_FUN             SUNDIALS_F77_FUNC(farkefun,                FARKEFUN)
#define FARK_MALLOC              SUNDIALS_F77_FUNC(farkmalloc,              FARKMALLOC)
#define FARK_REINIT              SUNDIALS_F77_FUNC(farkreinit,              FARKREINIT)
#define FARK_RESIZE              SUNDIALS_F77_FUNC(farkresize,              FARKRESIZE)
#define FARK_SETDEFAULTS         SUNDIALS_F77_FUNC(farksetdefaults,         FARKSETDEFAULTS)
#define FARK_SETIIN              SUNDIALS_F77_FUNC(farksetiin,              FARKSETIIN)
#define FARK_SETRIN              SUNDIALS_F77_FUNC(farksetrin,              FARKSETRIN)
#define FARK_SETADAPTMETHOD      SUNDIALS_F77_FUNC(farksetadaptivitymethod, FARKSETADAPTIVITYMETHOD)
#define FARK_SETERKTABLE         SUNDIALS_F77_FUNC(farkseterktable,         FARKSETERKTABLE)
#define FARK_SETIRKTABLE         SUNDIALS_F77_FUNC(farksetirktable,         FARKSETIRKTABLE)
#define FARK_SETARKTABLES        SUNDIALS_F77_FUNC(farksetarktables,        FARKSETARKTABLES)
#define FARK_SETRESTOLERANCE     SUNDIALS_F77_FUNC(farksetrestolerance,     FARKSETRESTOLERANCE)
#define FARK_SETDIAGNOSTICS      SUNDIALS_F77_FUNC(farksetdiagnostics,      FARKSETDIAGNOSTICS)
#define FARK_STOPDIAGNOSTICS     SUNDIALS_F77_FUNC(farkstopdiagnostics,     FARKSTOPDIAGNOSTICS)
#define FARK_DLSINIT             SUNDIALS_F77_FUNC(farkdlsinit,             FARKDLSINIT)
#define FARK_DLSMASSINIT         SUNDIALS_F77_FUNC(farkdlsmassinit,         FARKDLSMASSINIT)
#define FARK_SPILSINIT           SUNDIALS_F77_FUNC(farkspilsinit,           FARKSPILSINIT)
#define FARK_SPILSSETEPSLIN      SUNDIALS_F77_FUNC(farkspilssetepslin,      FARKSPILSSETEPSLIN)
#define FARK_SPILSMASSINIT       SUNDIALS_F77_FUNC(farkspilsmassinit,       FARKSPILSMASSINIT)
#define FARK_SPILSSETMASSEPSLIN  SUNDIALS_F77_FUNC(farkspilssetmassepslin,  FARKSPILSSETMASSEPSLIN)
#define FARK_ARKODE              SUNDIALS_F77_FUNC(farkode,                 FARKODE)
#define FARK_DKY                 SUNDIALS_F77_FUNC(farkdky,                 FARKDKY)
#define FARK_GETERRWEIGHTS       SUNDIALS_F77_FUNC(farkgeterrweights,       FARKGETERRWEIGHTS)
#define FARK_GETRESWEIGHTS       SUNDIALS_F77_FUNC(farkgetresweights,       FARKGETRESWEIGHTS)
#define FARK_GETESTLOCALERR      SUNDIALS_F77_FUNC(farkgetestlocalerr,      FARKGETESTLOCALERR)
#define FARK_FREE                SUNDIALS_F77_FUNC(farkfree,                FARKFREE)
#define FARK_WRITEPARAMETERS     SUNDIALS_F77_FUNC(farkwriteparameters,     FARKWRITEPARAMETERS)

#define FARK_DENSESETJAC         SUNDIALS_F77_FUNC(farkdensesetjac,         FARKDENSESETJAC)
#define FARK_DJAC                SUNDIALS_F77_FUNC(farkdjac,                FARKDJAC)

#define FARK_BANDSETJAC          SUNDIALS_F77_FUNC(farkbandsetjac,          FARKBANDSETJAC)
#define FARK_BJAC                SUNDIALS_F77_FUNC(farkbjac,                FARKBJAC)

#define FARK_SPARSESETJAC        SUNDIALS_F77_FUNC(farksparsesetjac,        FARKSPARSESETJAC)
#define FARK_SPJAC               SUNDIALS_F77_FUNC(farkspjac,               FARKSPJAC)

#define FARK_DENSESETMASS        SUNDIALS_F77_FUNC(farkdensesetmass,        FARKDENSESETMASS)
#define FARK_DMASS               SUNDIALS_F77_FUNC(farkdmass,               FARKDMASS)

#define FARK_BANDSETMASS         SUNDIALS_F77_FUNC(farkbandsetmass,         FARKBANDSETMASS)
#define FARK_BMASS               SUNDIALS_F77_FUNC(farkbmass,               FARKBMASS)

#define FARK_SPARSESETMASS       SUNDIALS_F77_FUNC(farksparsesetmass,       FARKSPARSESETMASS)
#define FARK_SPMASS              SUNDIALS_F77_FUNC(farkspmass,              FARKSPMASS)

#define FARK_SPILSSETJAC         SUNDIALS_F77_FUNC(farkspilssetjac,         FARKSPILSSETJAC)
#define FARK_JTSETUP             SUNDIALS_F77_FUNC(farkjtsetup,             FARKJTSETUP)
#define FARK_JTIMES              SUNDIALS_F77_FUNC(farkjtimes,              FARKJTIMES)

#define FARK_SPILSSETPREC        SUNDIALS_F77_FUNC(farkspilssetprec,        FARKSPILSSETPREC)
#define FARK_PSOL                SUNDIALS_F77_FUNC(farkpsol,                FARKPSOL)
#define FARK_PSET                SUNDIALS_F77_FUNC(farkpset,                FARKPSET)

#define FARK_SPILSSETMASS        SUNDIALS_F77_FUNC(farkspilssetmass,        FARKSPILSSETMASS)
#define FARK_MTSETUP             SUNDIALS_F77_FUNC(farkmtsetup,             FARKMTSETUP)
#define FARK_MTIMES              SUNDIALS_F77_FUNC(farkmtimes,              FARKMTIMES)

#define FARK_SPILSSETMASSPREC    SUNDIALS_F77_FUNC(farkspilssetmassprec,    FARKSPILSSETMASSPREC)
#define FARK_MASSPSOL            SUNDIALS_F77_FUNC(farkmasspsol,            FARKMASSPSOL)
#define FARK_MASSPSET            SUNDIALS_F77_FUNC(farkmasspset,            FARKMASSPSET)

#define FARK_EWTSET              SUNDIALS_F77_FUNC(farkewtset,              FARKEWTSET)
#define FARK_EWT                 SUNDIALS_F77_FUNC(farkewt,                 FARKEWT)

#define FARK_ADAPTSET            SUNDIALS_F77_FUNC(farkadaptset,            FARKADAPTSET)
#define FARK_ADAPT               SUNDIALS_F77_FUNC(farkadapt,               FARKADAPT)

#define FARK_EXPSTABSET          SUNDIALS_F77_FUNC(farkexpstabset,          FARKEXPSTABSET)
#define FARK_EXPSTAB             SUNDIALS_F77_FUNC(farkexpstab,             FARKEXPSTAB)

#else

#define FARK_IMP_FUN             farkifun_
#define FARK_EXP_FUN             farkefun_
#define FARK_MALLOC              farkmalloc_
#define FARK_REINIT              farkreinit_
#define FARK_RESIZE              farkresize_
#define FARK_SETDEFAULTS         farksetdefaults_
#define FARK_SETIIN              farksetiin_
#define FARK_SETRIN              farksetrin_
#define FARK_SETADAPTMETHOD      farksetadaptivitymethod_
#define FARK_SETERKTABLE         farkseterktable_
#define FARK_SETIRKTABLE         farksetirktable_
#define FARK_SETARKTABLES        farksetarktables_
#define FARK_SETRESTOLERANCE     farksetrestolerance_
#define FARK_SETDIAGNOSTICS      farksetdiagnostics_
#define FARK_STOPDIAGNOSTICS     farkstopdiagnostics_
#define FARK_DLSINIT             farkdlsinit_
#define FARK_DLSMASSINIT         farkdlsmassinit_
#define FARK_SPILSINIT           farkspilsinit_
#define FARK_SPILSSETEPSLIN      farkspilssetepslin_
#define FARK_SPILSMASSINIT       farkspilsmassinit_
#define FARK_SPILSSETMASSEPSLIN  farkspilssetmassepslin_
#define FARK_ARKODE              farkode_
#define FARK_DKY                 farkdky_
#define FARK_GETERRWEIGHTS       farkgeterrweights_
#define FARK_GETRESWEIGHTS       farkgetresweights_
#define FARK_GETESTLOCALERR      farkgetestlocalerr_
#define FARK_FREE                farkfree_
#define FARK_WRITEPARAMETERS     farkwriteparameters_

#define FARK_DENSESETJAC         farkdensesetjac_
#define FARK_DJAC                farkdjac_

#define FARK_BANDSETJAC          farkbandsetjac_
#define FARK_BJAC                farkbjac_

#define FARK_SPARSESETJAC        farksparsesetjac_
#define FARK_SPJAC               farkspjac_

#define FARK_DENSESETMASS        farkdensesetmass_
#define FARK_DMASS               farkdmass_

#define FARK_BANDSETMASS         farkbandsetmass_
#define FARK_BMASS               farkbmass_

#define FARK_SPARSESETMASS       farksparsesetmass_
#define FARK_SPMASS              farkspmass_


#define FARK_SPILSSETJAC         farkspilssetjac_
#define FARK_JTSETUP             farkjtsetup_
#define FARK_JTIMES              farkjtimes_

#define FARK_SPILSSETPREC        farkspilssetprec_
#define FARK_PSOL                farkpsol_
#define FARK_PSET                farkpset_

#define FARK_SPILSSETMASS        farkspilssetmass_
#define FARK_MTSETUP             farkmtsetup_
#define FARK_MTIMES              farkmtimes_

#define FARK_SPILSSETMASSPREC    farkspilssetmassprec_
#define FARK_MASSPSOL            farkmasspsol_
#define FARK_MASSPSET            farkmasspset_

#define FARK_EWTSET              farkewtset_
#define FARK_EWT                 farkewt_

#define FARK_ADAPTSET            farkadaptset_
#define FARK_ADAPT               farkadapt_

#define FARK_EXPSTABSET          farkexpstabset_
#define FARK_EXPSTAB             farkexpstab_

#endif

  /* Type for user data */
  typedef struct {
    realtype *rpar;
    long int *ipar;
  } *FARKUserData;

  /* Prototypes of exported functions */
  void FARK_MALLOC(realtype *t0, realtype *y0, int *imex, 
		   int *iatol, realtype *rtol, realtype *atol, 
		   long int *iout, realtype *rout, 
		   long int *ipar, realtype *rpar, int *ier);

  void FARK_REINIT(realtype *t0, realtype *y0, int *imex,
		   int *iatol, realtype *rtol, realtype *atol,
		   int *ier);

  void FARK_RESIZE(realtype *t0, realtype *y0, realtype *hscale,
		   int *itol, realtype *rtol, realtype *atol, int *ier);

  void FARK_SETDEFAULTS(int *ier);
  void FARK_SETIIN(char key_name[], long int *ival, int *ier);
  void FARK_SETRIN(char key_name[], realtype *rval, int *ier);

  void FARK_SETADAPTMETHOD(int *imethod, int *idefault, int *ipq, 
			   realtype *params, int *ier);

  void FARK_SETERKTABLE(int *s, int *q, int *p, realtype *c, realtype *A, 
			realtype *b, realtype *b2, int *ier);
  void FARK_SETIRKTABLE(int *s, int *q, int *p, realtype *c,
			realtype *A, realtype *b, realtype *b2, int *ier);
  void FARK_SETARKTABLES(int *s, int *q, int *p, realtype *ci, 
                         realtype *ce, realtype *Ai, realtype *Ae, 
                         realtype *bi, realtype *be, realtype *b2i, 
                         realtype *b2e, int *ier);

  void FARK_SETRESTOLERANCE(int *itol, realtype *atol, int *ier);
  void FARK_SETDIAGNOSTICS(char fname[], int *flen, int *ier);
  void FARK_STOPDIAGNOSTICS(int *ier);

  void FARK_DLSINIT(int *ier);
  void FARK_DLSMASSINIT(int *time_dep, int *ier);

  void FARK_SPILSINIT(int *ier);
  void FARK_SPILSSETEPSLIN(realtype *eplifac, int *ier);

  void FARK_SPILSMASSINIT(int *time_dep, int *ier);
  void FARK_SPILSSETMASSEPSLIN(realtype *eplifac, int *ier);

  void FARK_ARKODE(realtype *tout, realtype *t, realtype *y, 
                   int *itask, int *ier);
  void FARK_DKY(realtype *t, int *k, realtype *dky, int *ier);

  void FARK_GETERRWEIGHTS(realtype *eweight, int *ier);
  void FARK_GETRESWEIGHTS(realtype *rweight, int *ier);
  void FARK_GETESTLOCALERR(realtype *ele, int *ier);

  void FARK_FREE(void);

  void FARK_WRITEPARAMETERS(int *ier);

  void FARK_DENSESETJAC(int *flag, int *ier);
  void FARK_BANDSETJAC(int *flag, int *ier);
  void FARK_SPARSESETJAC(int *ier);

  void FARK_DENSESETMASS(int *ier);
  void FARK_BANDSETMASS(int *ier);
  void FARK_SPARSESETMASS(int *ier);


  void FARK_SPILSSETJAC(int *flag, int *ier);
  void FARK_SPILSSETPREC(int *flag, int *ier);
  void FARK_SPILSSETMASS(int *ier);
  void FARK_SPILSSETMASSPREC(int *flag, int *ier);

  void FARK_EWTSET(int *flag, int *ier);
  void FARK_ADAPTSET(int *flag, int *ier);
  void FARK_EXPSTABSET(int *flag, int *ier);


  /* Prototypes: Functions Called by the ARKODE Solver */
  int FARKfe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  int FARKfi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  
  int FARKDenseJac(realtype t, N_Vector y, N_Vector fy, 
		   SUNMatrix J, void *user_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKBandJac(realtype t, N_Vector y, N_Vector fy,
		  SUNMatrix J, void *user_data,
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKSparseJac(realtype t, N_Vector y, N_Vector fy, 
		    SUNMatrix J, void *user_data, N_Vector vtemp1, 
		    N_Vector vtemp2, N_Vector vtemp3);

  
  int FARKDenseMass(realtype t, SUNMatrix M, void *user_data, 
                    N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKBandMass(realtype t, SUNMatrix M, void *user_data, 
                    N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKSparseMass(realtype t, SUNMatrix M, void *user_data, 
                    N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  

  int FARKPSet(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
	       booleantype *jcurPtr, realtype gamma, void *user_data);
  
  int FARKMassPSet(realtype tn, void *user_data);
  
  int FARKPSol(realtype tn, N_Vector y, N_Vector fy, N_Vector r, 
               N_Vector z, realtype gamma, realtype delta, int lr, 
	       void *user_data);
  
  int FARKMassPSol(realtype tn, N_Vector r, N_Vector z, realtype delta, 
		   int lr, void *user_data);
  
  int FARKJTSetup(realtype t, N_Vector y, N_Vector fy, void *user_data);
  
  int FARKJtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
		 N_Vector fy, void *user_data, N_Vector work);
  
  int FARKMTSetup(realtype t, void *user_data);
  
  int FARKMtimes(N_Vector v, N_Vector Mv, realtype t, void *user_data);

  int FARKEwt(N_Vector y, N_Vector ewt, void *user_data);

  int FARKAdapt(N_Vector y, realtype t, realtype h1, realtype h2,
		realtype h3, realtype e1, realtype e2, realtype e3,
		int q, int p, realtype *hnew, void *user_data);

  int FARKExpStab(N_Vector y, realtype t, realtype *hstab, void *user_data);

  /* Declarations for global variables shared amongst various routines */
  extern N_Vector F2C_ARKODE_vec;             /* defined in FNVECTOR module */
  extern SUNMatrix F2C_ARKODE_matrix;         /* defined in FSUNMATRIX module */
  extern SUNMatrix F2C_ARKODE_mass_matrix;  
  extern SUNLinearSolver F2C_ARKODE_linsol;   /* defined in FSUNLINSOL module */
  extern SUNLinearSolver F2C_ARKODE_mass_sol; 

  extern void *ARK_arkodemem;     /* defined in farkode.c */
  extern long int *ARK_iout;      /* defined in farkode.c */
  extern realtype *ARK_rout;      /* defined in farkode.c */
  extern int ARK_nrtfn;           /* defined in farkode.c */
  extern int ARK_ls;              /* defined in farkode.c */
  extern int ARK_mass_ls;         /* defined in farkode.c */

  /* Linear solver IDs */
  enum { ARK_LS_ITERATIVE   = 0, 
         ARK_LS_DIRECT      = 1, 
	 ARK_LS_CUSTOM      = 2 };

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
