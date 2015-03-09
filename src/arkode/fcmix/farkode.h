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
 provides the necessary interface to ARKODE any acceptable NVECTOR
 implementation.
 
 A summary of the user-callable functions, with the corresponding 
 ARKODE functions, are as follows:
 
   Fortran                    ARKODE
   ---------------------      --------------------------------
   FNVINITS                   N_VNew_Serial
   FNVINITP                   N_VNew_Parallel
   FNVINITOMP                 N_VNew_OpenMP
   FNVINITPTS                 N_VNew_Pthreads

   FARKMALLOC                 ARKodeCreate, ARKodeSetUserData, 
                                 and ARKodeInit
   FARKREINIT                 ARKReInit
   FARKRESIZE                 ARKResize
   FARKSETIIN                 ARKodeSet* (integer arguments)
   FARKSETRIN                 ARKodeSet* (real arguments)
   FARKSETADAPTIVITYMETHOD    ARKodeSetAdaptivityMethod
   FARKEWTSET                 ARKodeWFtolerances
   FARKADAPTSET               ARKodeSetAdaptivityFn
   FARKEXPSTABSET             ARKodeSetStabilityFn
   FARKSETDIAGNOSTICS         ARKodeSetDiagnostics
   FARKSTOPDIAGNOSTICS        (none)

   FARKDENSE                  ARKDense
   FARKDENSESETJAC            ARKDlsSetDenseJacFn
   FARKMASSDENSE              ARKMassDense
   FARKDENSESETMASS           ARKDlsSetDenseMassFn

   FARKBAND                   ARKBand
   FARKBANDSETJAC             ARKDlsSetBandJacFn
   FARKMASSBAND               ARKMassBand
   FARKBANDSETMASS            ARKDlsSetBandMassFn

   FARKLAPACKDENSE            ARKLapackDense
   FARKLAPACKDENSESETJAC      ARKDlsSetDenseJacFn
   FARKMASSLAPACKDENSE        ARKMassLapackDense
   FARKLAPACKDENSESETMASS     ARKDlsSetDenseMassFn

   FARKLAPACKBAND             ARKLapackBand
   FARKLAPACKBANDSETJAC       ARKDlsSetBandJacFn
   FARKMASSLAPACKBAND         ARKMassLapackBand
   FARKLAPACKBANDSETMASS      ARKDlsSetBandMassFn

   FARKKLU                    ARKKLU
   FARKKLUREINIT              ARKKLUReinit
   FARKSUPERLUMT              ARKSuperLUMT
   FARKSPARSESETJAC           ARKSlsSetSparseJacFn
   FARKMASSKLU                ARKMassKLU
   FARKMASSKLUREINIT          ARKMassKLUReinit
   FARKMASSSUPERLUMT          ARKMassSuperLUMT
   FARKSPARSESETMASS          ARKSlsSetSparseMassFn

   FARKSPGMR                  ARKSpgmr and ARKSpilsSet*
   FARKSPGMRREINIT            ARKSpilsSet*
   FARKMASSSPGMR              ARKMassSpgmr and ARKSpilsSetMass*
   FARKMASSSPGMRREINIT        ARKSpilsSetMass*

   FARKSPBCG                  ARKSpbcg and ARKSpilsSet*
   FARKSPBCGREINIT            ARKSpilsSet*
   FARKMASSSPBCG              ARKMassSpbcg and ARKSpilsSetMass*
   FARKMASSSPBCGREINIT        ARKSpilsSetMass*

   FARKSPTFQMR                ARKSptfqmr and ARKSpilsSet*
   FARKSPTFQMRREINIT          ARKSpilsSet*
   FARKMASSSPTFQMR            ARKMassSptfqmr and ARKSpilsSetMass*
   FARKMASSSPTFQMRREINIT      ARKSpilsSetMass*

   FARKSPFGMR                 ARKSpfgmr and ARKSpilsSet*
   FARKSPFGMRREINIT           ARKSpilsSet*
   FARKMASSSPFGMR             ARKwMASSMassSpfgmr and ARKSpilsSetMASS*
   FARKMASSSPFGMRREINIT       ARKSpilsSetMASS*

   FARKPCG                    ARKPcg and ARKSpilsSet*
   FARKPCGREINIT              ARKSpilsSet*
   FARKMASSPCG                ARKMASSPcg and ARKSpilsSetMASS*
   FARKMASSPCGREINIT          ARKSpilsSetMASS*

   FARKSPILSSETJAC            ARKSpilsSetJacTimesVecFn
   FARKSPILSSETPREC           ARKSpilsSetPreconditioner
   FARKSPILSSETMASS           ARKSpilsSetMassTimesVecFn
   FARKSPILSSETMASSPREC       ARKSpilsSetMassPreconditioner
 
   FARKODE                    ARKode, ARKodeGet*, and ARK*Get*
   FARKDKY                    ARKodeGetDky

   FARKGETERRWEIGHTS          ARKodeGetErrWeights
   FARKGETESTLOCALERR         ARKodeGetEstLocalErrors

   FARKFREE                   ARKodeFree
   ---------------------      --------------------------------

 
 The user-supplied functions, each listed with the corresponding interface
 function which calls it (and its type within ARKODE), are as follows:

   Fortran:           Interface Fcn:           ARKODE Type:
   -------------      ------------------       -----------------------
   FARKIFUN           FARKfi                   ARKRhsFn
   FARKEFUN           FARKfe                   ARKRhsFn
   FARKDJAC           FARKDenseJac             ARKDlsDenseJacFn
   FARKDMASS          FARKDenseMass            ARKDlsDenseMassFn
   FARKBJAC           FARKBandJac              ARKDlsBandJacFn
   FARKBMASS          FARKBandMass             ARKDlsBandMassFn
   FARKSPJAC          FARKSparseJac            ARKSlsSparseJacFn
   FARKSPMASS         FARKSparseMass           ARKSlsSparseMassFn
   FARKPSOL           FARKPSol                 ARKSpilsPrecSolveFn
   FARKPSET           FARKPSet                 ARKSpilsPrecSetupFn
   FARKMASSPSOL       FARKMassPSol             ARKSpilsMassPrecSolveFn
   FARKMASSPSET       FARKMassPSet             ARKSpilsMassPrecSetupFn
   FARKJTIMES         FARKJtimes               ARKSpilsJacTimesVecFn
   FARKMTIMES         FARKMtimes               ARKSpilsMassTimesVecFn
   FARKEWT            FARKEwt                  ARKEwtFn
   FARKADAPT          FARKAdapt                ARKAdaptFn
   FARKEXPSTAB        FARKExpStab              ARKExpStabFn
   -------------      ------------------       -----------------------

 In contrast to the case of direct use of ARKODE, and of most Fortran ODE
 solvers, the names of all user-supplied routines here are fixed, in
 order to maximize portability for the resulting mixed-language program.
 
 Important note on portability:  In this package, the names of the interface 
 functions, and the names of the Fortran user routines called by them, appear 
 as dummy names which are mapped to actual values by a series of definitions, 
 in this and other header files.
 
 =============================================================================
 
                  Usage of the FARKODE Interface Package
 
 The usage of FARKODE requires calls to five or more interface
 functions, depending on the method options selected, and one or more
 user-supplied routines which define the problem to be solved.  These
 function calls and user routines are summarized separately below.
 
 Some details are omitted, and the user is referred to the ARKODE user 
 documentation for more complete information.  Information on the arguments 
 of any given user-callable interface routine, or of a given user-supplied 
 function called by an interface function, can be found in the documentation 
 on the corresponding function in the ARKODE package.
 
 The number labels on the instructions below end with s for instructions
 that are specific to use with the N_VSerial package; similarly those that 
 end with p are specific to use with the N_VParallel package.

 -----------------------------------------------------------------------------

                               Data Types

 Throughout this documentation, we will refer to data types according to their
 usage in SUNDIALS.  The equivalent types to these may vary, depending on your
 computer architecture and on how SUNDIALS was compiled.  A Fortran user 
 should take care that all arguments passed through this Fortran/C interface 
 are declared of the appropriate type.
 
 Integers: SUNDIALS uses both 'int' and 'long int' types:
   int      -- equivalent to an INTEGER or INTEGER*4 in Fortran
   long int -- this will depend on the computer architecture:
                 32-bit -- equivalent to an INTEGER or INTEGER*4 in Fortran
                 64-bit -- equivalent to an INTEGER*8 in Fortran
	      
 Real numbers:  At compilation, SUNDIALS allows the configuration option 
 '--with-precision', that accepts values of 'single', 'double' or 'extended' 
 (the default is 'double').  This choice dictates the size of a SUNDIALS 
 'realtype' variable.  The corresponding Fortran types for these 'realtype' 
 sizes are:
   single   -- equivalent to a REAL or REAL*4 in Fortran
   double   -- equivalent to a DOUBLE PRECISION or REAL*8 in Fortran
   extended -- equivalent to a REAL*16 in Fortran

 -----------------------------------------------------------------------------

 (1) User-supplied right-hand side routines: FARKIFUN and FARKEFUN
     The user must in all cases supply at least one of the following Fortran 
     routines:

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

     As an option when using the DENSE linear solver, the user may supply a
     routine that computes a dense approximation of the system Jacobian 
     J = dfi(t,y)/dy.  If supplied, it must have the following form:

       SUBROUTINE FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, 
      &                    WK3, IER)

     Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
     the Jacobian and store it column-wise in DJAC.

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
 
 (2s) Optional user-supplied dense mass matrix routine: FARKDMASS

     Required when using the DENSE or LAPACKDENSE mass matrix linear solvers,
     the user must supply a routine that computes the system mass matrix M.  
     This routine must have the following form: 

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
 
 (3s) Optional user-supplied band Jacobian approximation routine: FARKBJAC

     As an option when using the BAND linear solver, the user may supply a
     routine that computes a band approximation of the system Jacobian 
     J = dfi(t,y)/dy. If supplied, it must have the following form:

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
 
 (4s) Optional user-supplied band Mass matrix routine: FARKBMASS

     Required when using the BAND or LAPACKBAND mass matrix linear solvers, 
     the user must supply a routine that computes the system mass matrix M.  
     This routine must have the following form: 

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
 
 (5s) User-supplied sparse Jacobian approximation routine: FARKSPJAC

     Required when using the ARKKLU or ARKSuperLUMT linear solvers, the 
     user must supply a routine that computes a compressed-sparse-column 
     approximation of the system Jacobian J = dfi(t,y)/dy.  This routine
     must have the following form:

       SUBROUTINE FARKSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS, 
      &                     JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only M, N, NNZ, JDATA, JRVALS and 
     JCPTRS. It must load the N by N compressed sparse column matrix 
     with storage for NNZ nonzeros, stored in the arrays JDATA (nonzero
     values), JRVALS (row indices for each nonzero), JCOLPTRS (indices 
     for start of each column), with the Jacobian matrix at the current
     (t,y) in CSC form (see sundials_sparse.h for more information).

     The arguments are:
         T    -- current time [realtype, input]
         Y    -- array containing state variables [realtype, input]
         FY   -- array containing state derivatives [realtype, input]
         N    -- number of matrix rows/columns in Jacobian [int, input]
         NNZ  -- allocated length of nonzero storage [int, input]
        JDATA -- nonzero values in Jacobian
                 [realtype of length NNZ, output]
       JRVALS -- row indices for each nonzero in Jacobian
                  [int of length NNZ, output]
       JCPTRS -- pointers to each Jacobian column in preceding arrays
                 [int of length N+1, output]
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
 
 (6s) User-supplied sparse Mass matrix routine: FARKSPMASS

     Required when using the ARKMassKLU or ARKMassSuperLUMT mass matrix 
     linear solvers, the user must supply a routine that computes a 
     compressed-sparse-column version of the (possibly time-dependent) 
     system mass matrix M(t).  If supplied, it must have the following 
     form:

       SUBROUTINE FARKSPMASS(T, N, NNZ, MDATA, MRVALS, MCPTRS, 
      &                      IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only M, N, NNZ, MDATA, MRVALS and 
     MCPTRS. It must load the N by N compressed sparse column matrix 
     with storage for NNZ nonzeros, stored in the arrays MDATA (nonzero
     values), MRVALS (row indices for each nonzero), MCOLPTRS (indices 
     for start of each column), with the system mass matrix at the current
     (t) in CSC form (see sundials_sparse.h for more information).

     The arguments are:
         T    -- current time [realtype, input]
         N    -- number of rows/columns in mass matrix [int, input]
         NNZ  -- allocated length of nonzero storage [int, input]
        MDATA -- nonzero values in mass matrix
                 [realtype of length NNZ, output]
       MRVALS -- row indices for each nonzero in mass matrix
                  [int of length NNZ, output]
       MCPTRS -- pointers to each mass matrix column in preceding arrays
                 [int of length N+1, output]
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
 
 (7) Optional user-supplied Jacobian-vector product routine: FARKJTIMES

     As an option when using the SP* linear solver, the user may supply a
     routine that computes the product of the system Jacobian 
     J = dfi(t,y)/dy and a given vector v.  If supplied, it must have the 
     following form:

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
 
 (8) User-supplied mass-matrix-vector product routine: FARKMTIMES

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
 
 (9) Optional user-supplied error weight vector routine: FARKEWT
 
     As an option to providing the relative and absolute tolerances, the user
     may supply a routine that computes the weights used in the WRMS norms.
     If supplied, it must have the following form:

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

 (10) Optional user-supplied error weight vector routine: FARKADAPT
 
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
       E3   -- estimated temporal error in previous-previous step [realtype, input]
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

 (11) Optional user-supplied explicitly stable time step routine: FARKEXPSTAB
 
     As an option, the user may provide a routine to return the maximum stable 
     time step size for the explicit ODE RHS function.  If supplied, it must 
     have the following form:

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

 (12) Initialization:  FNVINITS/FNVINITP/FNVINITOMP/FNVINITPTS, 
                       FARKMALLOC, FARKREINIT, FARKRESIZE
 
 (12.1s) To initialize the serial vector specification, the user must make the
     following call:

        CALL FNVINITS(4, NEQ, IER)

     where the first argument is an int containing the ARKODE solver ID (4). 
     The other arguments are:
        NEQ = size of vectors [long int, input]
	IER = return completion flag [int, output]:
	          0 = success, 
		 -1 = failure.
 
 (12.1p) To initialize the distributed memory parallel machine environment, 
        the user must make the following call:

        CALL FNVINITP(COMM, 4, NLOCAL, NGLOBAL, IER)

     The arguments are:
        COMM = the MPI communicator [int, input]
        NLOCAL = local size of vectors on this processor [long int, input]
        NGLOBAL = the system size, and the global size of vectors (the sum 
           of all values of NLOCAL) [long int, input]
        IER = return completion flag [int, ouptut]. 
                  0 = success, 
                 -1 = failure.

     Note: If MPI was initialized by the user, the communicator must be set to
     MPI_COMM_WORLD.  If not, this routine initializes MPI and sets the
     communicator equal to MPI_COMM_WORLD.
 
 (12.1omp) To initialize the openMP threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP(3, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (12.1pts) To initialize the Pthreads threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITPTS(3, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (12.2) To set various problem and solution parameters and allocate
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
	IOUT = array of length 22 for integer optional outputs
	   [long int, output]
	ROUT = array of length 6 for real optional outputs [realtype, output]
	IPAR = array of user integer data [long int, input/output]
	RPAR = array with user real data [realtype, input/output]
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

     The user data arrays IPAR and RPAR are passed unmodified to all 
     subsequent calls to user-provided routines. Modifications to either array
     inside a user-provided routine will be propagated. Using these two 
     arrays, the user can dispense with COMMON blocks to pass data betwen 
     user-provided routines. 
 
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

 (12.3) If the user program includes the FARKEWT routine for the evaluation 
     of the error weights, the following call must be made

       CALL FARKEWTSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKEWT is provided and 
     should be used; FLAG = 0 resets to the default EWT formulation.  The int 
     return flag IER is 0 if successful, and nonzero otherwise.

 (12.4) If the user program includes the FARKADAPT routine for performing 
     step adaptivity, the following call must be made

       CALL FARKADAPTSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKADAPT is provided and 
     should be used; FLAG = 0 resets to the default adaptivity formulation. 
     The int return flag IER is 0 if successful, and nonzero otherwise.

 (12.5) If the user program includes the FARKEXPSTAB routine for calculation of
     the maximum explicitly stable step size, the following call must be made

       CALL FARKEXPSTABSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKEXPSTAB is provided and
     should be used; FLAG = 0 resets to the default explicit stability formulation.
     The int return flag IER is 0 if successful, and nonzero otherwise.

 (12.6) To re-initialize the ARKODE solver for the solution of a new problem
     of the same size as one already solved, make the following call:

       CALL FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)

     The arguments have the same names and meanings as those of FARKMALLOC.  
     FARKREINIT performs the same initializations as FARKMALLOC, but does no 
     memory allocation, using instead the existing internal memory created by 
     the previous FARKMALLOC call.  The subsequent call to specify the linear 
     system solution method may or may not be needed; see paragraph (9) below.
 
 (12.7) To re-initialize the ARKODE solver for the solution of a new problem
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

     FARKRESIZE performs the opposite set of of operations as FARKREINIT: it 
     does not reinitialize any of the time-step heuristics, but it does 
     perform memory reallocation.  The subsequent call to specify the linear 
     system solution method is required, since its internal memory 
     structures will no longer be the correct size.; see paragraph (9) below.
 
 (12.8) To set various integer optional inputs, make the folowing call:

       CALL FARKSETIIN(KEY, VALUE, IER)

     to set the integer value VALUE to the optional input specified by the 
     quoted character string KEY. VALUE must be a Fortran integer of size 
     commensurate with a C "long int".  KEY must be one of the following: 
     ORDER, DENSE_ORDER, LINEAR, NONLINEAR, FIXEDPOINT, NEWTON, EXPLICIT, 
     IMPLICIT, IMEX, IRK_TABLE_NUM, ERK_TABLE_NUM, ARK_TABLE_NUM (pass in an 
     int array of length 2, implicit method first), MAX_NSTEPS, HNIL_WARNS, 
     PREDICT_METHOD, MAX_ERRFAIL, MAX_CONVFAIL, MAX_NITERS, ADAPT_SMALL_NEF 
     or LSETUP_MSBP.  The int return flag IER is 0 if successful, and 
     nonzero otherwise.

     To set various real optional inputs, make the following call:

       CALL FARKSETRIN(KEY, VALUE, IER)

     to set the realtype value VALUE to the optional input specified by the 
     quoted character string KEY.  VALUE must be a Fortran real-valued number
     of size commensurate with the SUNDIALS "realtype".  KEY must one of the 
     following: INIT_STEP, MAX_STEP, MIN_STEP, STOP_TIME, NLCONV_COEF, 
     ADAPT_CFL, ADAPT_SAFETY, ADAPT_BIAS, ADAPT_GROWTH, ADAPT_BOUNDS (pass in 
     a realtype array of length 2), ADAPT_ETAMX1, ADAPT_ETAMXF, ADAPT_ETACF, 
     NONLIN_CRDOWN, NONLIN_RDIV, LSETUP_DGMAX, or FIXED_STEP.

     To set the time step adaptivity method (and its associated parameters), 
     make the following call: 

       CALL FARKSETADAPTIVITYMETHOD(IMETHOD, IDEFAULT, IPQ, PARAMS, IER)
       
     The arguments are:
       IMETHOD  = the adaptivity method to use
       IDEFAULT = flag to use (1) or not (0) the default adaptivity parameters
       IPQ      = flag to use the embedding order p (0) or the method order 
                  q (1) for error-based step adaptivity
       PARAMS   = if IDEFAULT=0, this should be a realtype array of length 2 
                  containing the custom adaptivity parameters to use in the 
		  method.
       IMETHOD  = integer error flag (0 = success, 1 = failure)

     To reset all optional inputs to their default values, make the following 
     call:

       CALL FARKSETDEFAULTS(IER)

     To set a custom explicit Runge-Kutta table, make the following call:

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

     To set a custom diagonally-implicit Runge-Kutta table, make the following
     call:

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

       CALL FARKSETARKTABLES(S, Q, P, C, AI, AE, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       AI = array of length S*S containing the DIRK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       AE = array of length S*S containing the ERK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

 (12.9) To set a solver diagnostics output file, make the folowing call:

       CALL FARKSETDIAGNOSTICS(FNAME, FLEN, IER)

     The desired diagnostics filename should be supplied by the 
     quoted character string FNAME.  The integer argument FLEN should contain
     the length (in characters) of FNAME (for portability).  The int return 
     flag IER is 0 if successful (able to open file), and nonzero otherwise.

 (12.10) To close the solver diagnostics output file, make the folowing call:

       CALL FARKSTOPDIAGNOSTICS(IER)

     The int return flag IER is 0 if successful (able to close file), and
     nonzero otherwise.


 -----------------------------------------------------------------------------

 (13) Specification of linear system solution method.

     In the case of using either an implicit or ImEx method and a Newton 
     iteration, the solution of each Runge-Kutta stage may involve the 
     solution of linear systems related to the Jacobian J = dfi(t,y)/dy of the 
     implicit portion of the ODE system. ARKode presently includes a variety 
     of choices for the treatment of these systems, and the user of FARKODE 
     must call a routine with a specific name to make the desired choice.  
     Following any call to FARKMALLOC or FARKRESIZE, one of these solver 
     specification routines must be called again to set up the linear solver
     memory.

 (13.1s) DENSE treatment of the linear system.

     The user must make the call

       CALL FARKDENSE(NEQ, IER)

     The arguments are:
        NEQ = the problem size [long int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     If the user program includes the FARKDJAC routine for the evaluation of 
     the dense approximation to the Jacobian, then after the call to 
     FARKDENSE, the following call must be made 

       CALL FARKDENSESETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKDJAC is provided and should be 
     used; FLAG=0 specifies a reset to the internal finite difference Jacobian
     approximation.  The int return flag IER=0 if successful, and nonzero 
     otherwise.
 
     Optional outputs specific to the DENSE case are:
        LENRWLS = IOUT(14) from ARKDlsGetWorkSpace (realtype space)
        LENIWLS = IOUT(15) from ARKDlsGetWorkSpace (integer space)
        LSTF    = IOUT(16) from ARKDlsGetLastFlag
        NFELS   = IOUT(17) from ARKDlsGetNumRhsEvals
        NJED    = IOUT(18) from ARKDlsGetNumJacEvals
     See the ARKODE manual for descriptions.
 
 (13.2s) DENSE treatment of the mass matrix linear system.

     The user must make the call

       CALL FARKMASSDENSE(NEQ, IER)

     The arguments are:
        NEQ = the problem size [long int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     When using MASSDENSE, the user program must provide the FARKDMASS routine 
     for the evaluation of the dense mass matrix.  To indicate that this routine
     has been provided, following the call to FARKMASSDENSE, the following call 
     must be made 

       CALL FARKDENSESETMASS(IER)

     The int return flag IER=0 if successful, nonzero otherwise.
 
     Optional outputs specific to the MASSDENSE case are:
        LENRWMS = IOUT(23) from ARKDlsGetMassWorkSpace (realtype space)
        LENIWMS = IOUT(24) from ARKDlsGetMassWorkSpace (integer space)
        LSTMF   = IOUT(25) from ARKDlsGetLastMassFlag
        NME     = IOUT(26) from ARKDlsGetNumMassEvals
     See the ARKODE manual for descriptions.
 
 (13.3s) BAND treatment of the linear system

     The user must make the call

       CALL FARKBAND(NEQ, MU, ML, IER)

     The arguments are:
        NEQ = problem size [long int, input]
	MU = upper bandwidth [long int, input]
	ML = lower bandwidth [long int, input]
        IER = return flag [int, output]; 0 if successful, nonzero otherwise.
 
     If the user program includes the FARKBJAC routine for the evaluation of 
     the band approximation to the Jacobian, then following the call to 
     FARKBAND, the following call must be made 

       CALL FARKBANDSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKBJAC is provided and should be 
     used; FLAG=0 specifies a reset to the internal finite difference Jacobian 
     approximation.  The int return flag IER=0 if successful, nonzero otherwise.
 
     Optional outputs specific to the BAND case are:
        LENRWLS = IOUT(14) from ARKDlsGetWorkSpace (realtype space)
        LENIWLS = IOUT(15) from ARKDlsGetWorkSpace (integer space)
        LSTF    = IOUT(16) from ARKDlsGetLastFlag
        NFELS   = IOUT(17) from ARKDlsGetNumRhsEvals
        NJED    = IOUT(18) from ARKDlsGetNumJacEvals
     See the ARKODE manual for descriptions.

 (13.4s) BAND treatment of the mass matrix linear system

     The user must make the call

       CALL FARKMASSBAND(NEQ, MU, ML, IER)

     The arguments are:
        NEQ = problem size [long int, input]
	MU = upper bandwidth [long int, input]
	ML = lower bandwidth [long int, input]
        IER = return flag [int, output]; 0 if successful, nonzero otherwise.
 
     When using MASSBAND, the user program must provide the FARKBMASS routine 
     for the evaluation of the band mass matrix.  To indicate that this routine
     has been provided, following the call to FARKMASSBAND, the following call 
     must be made 

       CALL FARKBANDSETMASS(IER)

     The int return flag IER=0 if successful, nonzero otherwise.
 
     Optional outputs specific to the MASSBAND case are:
        LENRWMS = IOUT(23) from ARKDlsGetMassWorkSpace (realtype space)
        LENIWMS = IOUT(24) from ARKDlsGetMassWorkSpace (integer space)
        LSTMF   = IOUT(25) from ARKDlsGetLastMassFlag
        NME     = IOUT(26) from ARKDlsGetNumMassEvals
     See the ARKODE manual for descriptions.

 (13.5s) LAPACK dense treatment of the linear system

     The user must make the call

       CALL FARKLAPACKDENSE(NEQ, IER)

     The arguments match those for FARKDENSE, except that NEQ is now a 
     normal int (and not a long int).

     Following the call to FARKLAPACKDENSE, the user may optionally call

       CALL FARKLAPACKDENSESETJAC(FLAG, IER)
       
     with the int FLAG=1 if the user provides the function FARKDJAC (and 
     wishes to use it); FLAG=0 implies a reset to the default finite 
     difference Jacobian approximation.

     The optional outputs when using FARKLAPACKDENSE match those from 
     FARKDENSE.

 (13.6s) LAPACK dense treatment of the mass matrix linear system

     The user must make the call

       CALL FARKMASSLAPACKDENSE(NEQ, IER)

     The arguments match those for FARKMASSDENSE, except that NEQ is 
     now a normal int (and not a long int).

     Following the call to FARKMASSLAPACKDENSE, the user must call 

       CALL FARKLAPACKDENSESETMASS(IER)

     The int return flag IER=0 if successful, nonzero otherwise.
 
     The optional outputs when using FARKMASSLAPACKDENSE match those 
     from FARKMASSDENSE.

 (13.7s) LAPACK band treatment of the linear system

     The user must make the call

       CALL FARKLAPACKBAND(NEQ, MUPPER, MLOWER, IER)

     The arguments match those for FARKBAND, except that now all arguments 
     have type 'int'.

     Following the call to FARKLAPACKBAND, the user may optionally call

       CALL FARKLAPACKBANDSETJAC(FLAG, IER)

     with the int FLAG=1 if the user provides the function FARKBJAC (and 
     wishes to use it); FLAG=0 implies a reset to the default finite 
     difference Jacobian approximation.

     The optional outputs when using FARKLAPACKBAND match those from FARKBAND.

 (13.8s) LAPACK band treatment of the mass matrix linear system

     The user must make the call

       CALL FARKMASSLAPACKBAND(NEQ, MUPPER, MLOWER, IER)

     The arguments match those for FARKMASSBAND, except that now all 
     arguments have type 'int'.

     Following the call to FARKMASSLAPACKBAND, the user must call 

       CALL FARKLAPACKBANDSETMASS(IER)

     The int return flag IER=0 if successful, nonzero otherwise.
 
     The optional outputs when using FARKMASSLAPACKBAND match those from 
     FARKMASSBAND.

 (13.9s) SPARSE treatment of the linear system using the KLU solver.

     The user must make the call

       CALL FARKKLU(NEQ, NNZ, ORDERING, IER)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the KLU package (0 = AMD, 1 = COLAMD) [int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     When using the KLU solver the user must provide the FARKSPJAC routine for the 
     evaluation of the sparse approximation to the Jacobian.  To indicate that this
     routine has been provided, after the call to FARKKLU, the following call must 
     be made 

       CALL FARKSPARSESETJAC(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.
 
 
     The ARKODE KLU solver will reuse much of the factorization information from one
     nonlinear iteration to the next.  If at any time the user wants to force a full
     refactorization or if the number of nonzeros in the Jacobian matrix changes, the
     user should make the call

       CALL FARKKLUREINIT(NEQ, NNZ, REINIT_TYPE)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be destroyed and 
          a new one will be allocated with NNZ nonzeros.  For a value of 2, 
	  only symbolic and numeric factorizations will be completed. 
 
     Optional outputs specific to the KLU case are:
        LSTF    = IOUT(16) from ARKSlsGetLastFlag
        NJES    = IOUT(18) from ARKSlsGetNumJacEvals
     See the ARKODE manual for descriptions.
 
 (13.10s) SPARSE treatment of the mass matrix linear system using the KLU solver.

     The user must make the call

       CALL FARKMASSKLU(NEQ, NNZ, ORDERING, IER)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of mass matrix nonzeros [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the KLU package (0 = AMD, 1 = COLAMD) [int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     When using the MASSKLU solver the user must provide the FARKSPMASS routine for 
     the evaluation of the sparse mass matrix.  To indicate that this routine has 
     been provided, after the call to FARKMASSKLU, the following call must be made 

       CALL FARKSPARSESETMASS(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.
 

     The ARKODE KLU solver will reuse much of the factorization information from one
     mass matrix solve to the next.  If at any time the user wants to force a full
     refactorization or if the number of nonzeros in the mass matrix changes, the
     user should make the call

       CALL FARKMASSKLUREINIT(NEQ, NNZ, REINIT_TYPE)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be destroyed and 
          a new one will be allocated with NNZ nonzeros.  For a value of 2, 
	  only symbolic and numeric factorizations will be completed. 
 
     Optional outputs specific to the MASSKLU case are:
        LSTMF = IOUT(25) from ARKSlsGetLastMassFlag
        NME = IOUT(26) from ARKSlsGetNumMassEvals
     See the ARKODE manual for descriptions.
 
 (13.11s) SPARSE treatment of the linear system using the SuperLUMT solver.

     The user must make the call

       CALL FARKSUPERLUMT(NTHREADS, NEQ, NNZ, ORDERING, IER)

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

     When using the SUPERLUMT solver the user must provide the FARKSPJAC routine
     for the evaluation of the sparse approximation to the Jacobian.  To indicate 
     that this routine has been provided, after the call to FARKKLU, the following 
     call must be made 

       CALL FARKSPARSESETJAC(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.
 
     Optional outputs specific to the SUPERLUMT case are:
        LSTF    = IOUT(16) from ARKSlsGetLastFlag
        NJES    = IOUT(18) from ARKSlsGetNumJacEvals
     See the ARKODE manual for descriptions.
 
 (13.12s) SPARSE treatment of the mass matrix linear system using the SuperLUMT solver.

     The user must make the call

       CALL FARKMASSSUPERLUMT(NTHREADS, NEQ, NNZ, ORDERING, IER)

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

     When using the MASSSUPERLUMT solver the user must provide the FARKSPJAC 
     routine for the evaluation of the mass matrix.  To indicate that this 
     routine has been provided, after the call to FARKKLU, the following call
     must be made 

       CALL FARKSPARSESETMASS(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.
 
     Optional outputs specific to the MASSSUPERLUMT case are:
        LSTMF = IOUT(25) from ARKSlsGetLastMassFlag
        NME = IOUT(26) from ARKSlsGetNumMassEvals
     See the ARKODE manual for descriptions.
 
 (13.13) SPGMR treatment of the linear systems.

     For the Scaled Preconditioned GMRES solution of the linear systems,
     the user must make the following call:

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
	DELT = linear convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     Optional outputs specific to the SPGMR case are:
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
 
     If a sequence of problems of the same size is being solved using the
     SPGMR linear solver, then following the call to FARKREINIT, a call to 
     the FARKSPGMRREINIT routine is needed if any of IPRETYPE, IGSTYPE, DELT 
     is being changed.  In that case, call FARKSPGMRREINIT as follows:

       CALL FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKSPGMR.  If MAXL is being
     changed, then the user should call FARKSPGMR instead.  

     Note: if the problem has been resized using FARKRESIZE, then FARKSPGMR 
     must be called again. 
 
 (13.14) SPGMR treatment of the mass matrix linear systems.

     For the Scaled Preconditioned GMRES solution of the mass matrix 
     linear systems, the user must make the following call:

       CALL FARKMASSSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)              

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
	DELT = linear convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     Optional outputs specific to the MASSSPGMR case are:
        LENRWMS = IOUT(23) from ARKSpilsGetMassWorkSpace
        LENIWMS = IOUT(24) from ARKSpilsGetMassWorkSpace
        LSTMF   = IOUT(25) from ARKSpilsGetLastMassFlag
        NMPE    = IOUT(26) from ARKSpilsGetNumMassPrecEvals
        NMPS    = IOUT(27) from ARKSpilsGetNumMassPrecSolves
        NMLI    = IOUT(28) from ARKSpilsGetNumMassIters
        NMCFL   = IOUT(29) from ARKSpilsGetNumMassConvFails
     See the ARKODE manual for descriptions.
 
     If a sequence of problems of the same size is being solved using the
     MASSSPGMR linear solver, then following the call to FARKREINIT, a call
     to the FARKMASSSPGMRREINIT routine is needed if any of IPRETYPE, 
     IGSTYPE, DELT is being changed.  In that case, call FARKMASSSPGMRREINIT 
     as follows:

       CALL FARKMASSSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKMASSSPGMR.  If MAXL is 
     being changed, then the user should call FARKMASSSPGMR instead.  

     Note: if the problem has been resized using FARKRESIZE, then 
     FARKMASSSPGMR must be called again. 
 
 (13.15) SPBCG treatment of the linear systems.

     For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
     the user must make the following call:

       CALL FARKSPBCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the SPBCG case are identical to 
     those from SPGMR.
 
     If a sequence of problems of the same size is being solved using the
     SPBCG linear solver, then following the call to FARKREINIT, a call to the
     FARKSPBCGREINIT routine is needed if any of its arguments is being 
     changed.  The call is:

       CALL FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKSPBCG.

     Note: if the problem has been resized using FARKRESIZE, then FARKSPBCG 
     must be called again. 

 (13.16) SPBCG treatment of the mass matrix linear systems.

     For the Scaled Preconditioned Bi-CGSTAB solution of the mass matrix 
     linear systems, the user must make the following call:

       CALL FARKMASSSPBCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the MASSSPBCG case are identical to 
     those from MASSSPGMR.
 
     If a sequence of problems of the same size is being solved using the
     MASSSPBCG linear solver, then following the call to FARKREINIT, a call 
     to the FARKMASSSPBCGREINIT routine is needed if any of its arguments 
     is being changed.  The call is:

       CALL FARKMASSSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKMASSSPBCG.

     Note: if the problem has been resized using FARKRESIZE, then 
     FARKMASSSPBCG must be called again. 

 (13.17) SPTFQMR treatment of the linear systems.

     For the Scaled Preconditioned TFQMR solution of the linear systems, the
     user must make the following call:

       CALL FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]
	      0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the SPTFQMR case are identical to 
     those from SPGMR.
 
     If a sequence of problems of the same size is being solved using the
     SPTFQMR linear solver, then following the call to FARKREINIT, a call to 
     the FARKSPTFQMRREINIT routine is needed if any of its arguments is
     being changed.  The call is:

       CALL FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKSPTFQMR.

     Note: if the problem has been resized using FARKRESIZE, then FARKSPTFQMR 
     must be called again. 

 (13.18) SPTFQMR treatment of the mass matrix linear systems.

     For the Scaled Preconditioned TFQMR solution of the mass matrix 
     linear systems, the user must make the following call:

       CALL FARKMASSSPTFQMR(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]
	      0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the MASSSPTFQMR case are identical to 
     those from MASSSPGMR.
 
     If a sequence of problems of the same size is being solved using the
     MASSSPTFQMR linear solver, then following the call to FARKREINIT, a 
     call to the FARKMASSSPTFQMRREINIT routine is needed if any of its 
     arguments is being changed.  The call is:

       CALL FARKMASSSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKMASSSPTFQMR.

     Note: if the problem has been resized using FARKRESIZE, then 
     FARKMASSSPTFQMR must be called again. 

 (13.19) SPFGMR treatment of the linear systems.

     For the Scaled Preconditioned Flexible GMRES solution of the linear 
     systems, the user must make the following call:

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
	DELT = linear convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     Optional outputs specific to the SPFGMR case are identical to 
     those from SPGMR.
 
     If a sequence of problems of the same size is being solved using the
     SPFGMR linear solver, then following the call to FARKREINIT, a call to 
     the FARKSPFGMRREINIT routine is needed if any of IPRETYPE, IGSTYPE, DELT 
     is being changed.  In that case, call FARKSPFGMRREINIT as follows:

       CALL FARKSPFGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKSPFGMR.  If MAXL is being
     changed, then the user should call FARKSPFGMR instead.  

     Note: if the problem has been resized using FARKRESIZE, then FARKSPFGMR 
     must be called again. 
 
 (13.20) SPFGMR treatment of the mass matrix linear systems.

     For the Scaled Preconditioned Flexible GMRES solution of the mass 
     matrix linear systems, the user must make the following call:

       CALL FARKMASSSPFGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)              

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
	DELT = linear convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     Optional outputs specific to the MASSSPFGMR case are identical to 
     those from MASSSPGMR.
 
     If a sequence of problems of the same size is being solved using the
     MASSSPFGMR linear solver, then following the call to FARKREINIT, a 
     call to the FARKMASSSPFGMRREINIT routine is needed if any of IPRETYPE, 
     IGSTYPE, DELT is being changed.  In that case, call 
     FARKMASSSPFGMRREINIT as follows:

       CALL FARKMASSSPFGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKMASSSPFGMR.  If MAXL 
     is being changed, then the user should call FARKMASSSPFGMR instead.  

     Note: if the problem has been resized using FARKRESIZE, then 
     FARKMASSSPFGMR must be called again. 
 
 (13.21) PCG treatment of the linear systems.

     For the Preconditioned Conjugate Gradient solution of the linear systems,
     the user must make the following call:

       CALL FARKPCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the PCG case are identical to 
     those from SPGMR.
 
     If a sequence of problems of the same size is being solved using the
     PCG linear solver, then following the call to FARKREINIT, a call to the
     FARKPCGREINIT routine is needed if any of its arguments is being 
     changed.  The call is:

       CALL FARKPCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKPCG.

     Note: if the problem has been resized using FARKRESIZE, then FARKPCG
     must be called again. 

 (13.22) PCG treatment of the mass matrix linear systems.

     For the Preconditioned Conjugate Gradient solution of the mass matrix 
     linear systems, the user must make the following call:

       CALL FARKMASSPCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the MASSPCG case are identical to 
     those from MASSSPGMR.
 
     If a sequence of problems of the same size is being solved using the
     MASSPCG linear solver, then following the call to FARKREINIT, a call 
     to the FARKMASSPCGREINIT routine is needed if any of its arguments 
     is being changed.  The call is:

       CALL FARKMASSPCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKMASSPCG.

     Note: if the problem has been resized using FARKRESIZE, then 
     FARKMASSPCG must be called again. 

 (13.23) Usage of user-supplied routines for the system Krylov solvers

     If the user program includes the FARKJTIMES routine for the evaluation of
     the Jacobian vector product, then after specifying the linear solver 
     choice (e.g. FARKSPGMR), the following call must be made

       CALL FARKSPILSSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKJTIMES is provided and should
     be used; FLAG=0 specifies a reset to the internal finite difference 
     approximation to this product).  The int return flag IER=0 if 
     successful, and nonzero otherwise.
 
     Usage of the user-supplied routines FARKPSOL and FARKPSET for solution of
     the preconditioner linear system similarly requires the following call after 
     specifying the linear solver module:

       CALL FARKSPILSSETPREC(FLAG, IER)

     with the int FLAG=1.  If FLAG=0 then preconditioning with these routines 
     will be disabled. The return flag IER=0 if successful, nonzero otherwise.

     NOTE: following any call to FARKRESIZE, either of the above routines must 
     again be called following re-specification of the linear solver module.

     The user-supplied routine FARKPSET must have the form:

       SUBROUTINE FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,V1,V2,V3,IER)

     This routine must set up the preconditioner P to be used in the 
     subsequent call to FARKPSOL.  The preconditioner (or the product of the 
     left and right preconditioners if using both) should be an approximation 
     to the matrix  M(T) - GAMMA*J  (M(T) = mass matrix, J = Jacobian),

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
       V* -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

     The user-supplied routine FARKPSOL must have the form:

       SUBROUTINE FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)

       DIMENSION Y(*), FY(*), VT(*), R(*), Z(*), IPAR(*), RPAR(*)

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
       VT -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

 (13.24) Usage of user-supplied routines for the mass matrix Krylov solvers

     When using any of the MASSSPGMR, MASSSPBCG, MASSSPTFQMR, MASSSPFGMR or 
     MASSPCG mass matrix linear solvers, the user must supply the 
     FARKMTIMES routine for the evaluation of the mass matrix vector product.  
     After specifying the linear solver choice (e.g. FARKMASSSPGMR), the 
     following call must be made

       CALL FARKSPILSSETMASS(IER)

     The int return flag IER=0 if successful, and nonzero otherwise.

 
     Usage of the user-supplied routines FARKMASSPSOL and FARKMASSPSET for 
     solution of the mass matrix preconditioner linear system similarly 
     requires the following call after specifying the mass matrix linear 
     solver module:

       CALL FARKSPILSSETMASSPREC(FLAG, IER)

     with the int FLAG=1. If FLAG=0 then preconditioning with these routines 
     will be disabled. The return flag IER=0 if successful, nonzero otherwise.

     NOTE: following any call to FARKRESIZE, the above routine must again 
     be called following re-specification of the linear solver module.

     The user-supplied routine FARKMASSPSET must have the form:

       SUBROUTINE FARKMASSPSET(T,IPAR,RPAR,V1,V2,V3,IER)

     This routine must set up the preconditioner P to be used in the 
     subsequent call to FARKMASSPSOL.  The preconditioner (or the product 
     of the left and right preconditioners if using both) should be an 
     approximation to the system mass matrix M.

     The arguments are:
       T = current time [realtype, input]
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       V* -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

     The user-supplied routine FARKMASSPSOL must have the form:

       SUBROUTINE FARKMASSPSOL(T,R,Z,DELTA,LR,IPAR,RPAR,VT,IER)

       DIMENSION VT(*), R(*), Z(*), IPAR(*), RPAR(*)

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
       VT -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

 -----------------------------------------------------------------------------

 (14) The integrator: FARKODE

     Carrying out the integration is accomplished by making calls as follows:

       CALL FARKODE(TOUT, T, Y, ITASK, IER)

     The arguments are:
       TOUT = next value of t at which a solution is desired [realtype, input]
       T = value of t reached by the solver [realtype, output]
       Y = array containing state variables on output (realtype, output]
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

 (15) Computing solution derivatives: FARKDKY

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

 (16) Get the current weight vector: FARKGETERRWEIGHTS

     To obtain the current weight vector, make the following call:

       CALL FARKGETERRWEIGHTS(EWT, IER)

     The arguments are:
       EWT = array containing the error weight vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (17) Get an estimate of the local error: FARKGETESTLOCALERR

     To obtain the current error estimate vector, make the following call:

       CALL FARKGETESTLOCALERR(ELE, IER)

     The arguments are:
       ELE = array with the estimated local error vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (18) Memory freeing: FARKFREE 

     To free the internal memory created by the calls to FARKMALLOC and
     FNVINITS, FNVINITOMP, FNVINITPTS or FNVINITP, make the call

       CALL FARKFREE()
 
===============================================================*/

#ifndef _FARKODE_H
#define _FARKODE_H

/* header files  */
#include <arkode/arkode.h>
#include <sundials/sundials_direct.h>  /* definition of type DlsMat   */
#include <sundials/sundials_sparse.h>  /* definition of type SlsMat   */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_MALLOC              SUNDIALS_F77_FUNC(farkmalloc,              FARKMALLOC)
#define FARK_REINIT              SUNDIALS_F77_FUNC(farkreinit,              FARKREINIT)
#define FARK_RESIZE              SUNDIALS_F77_FUNC(farkresize,              FARKRESIZE)
#define FARK_SETDEFAULTS         SUNDIALS_F77_FUNC(farksetdefaults,         FARKSETDEFAULTS)
#define FARK_SETIIN              SUNDIALS_F77_FUNC(farksetiin,              FARKSETIIN)
#define FARK_SETRIN              SUNDIALS_F77_FUNC(farksetrin,              FARKSETRIN)
#define FARK_SETERKTABLE         SUNDIALS_F77_FUNC(farkseterktable,         FARKSETERKTABLE)
#define FARK_SETIRKTABLE         SUNDIALS_F77_FUNC(farksetirktable,         FARKSETIRKTABLE)
#define FARK_SETARKTABLES        SUNDIALS_F77_FUNC(farksetarktables,        FARKSETARKTABLES)
#define FARK_SETDIAGNOSTICS      SUNDIALS_F77_FUNC(farksetdiagnostics,      FARKSETDIAGNOSTICS)
#define FARK_STOPDIAGNOSTICS     SUNDIALS_F77_FUNC(farkstopdiagnostics,     FARKSTOPDIAGNOSTICS)
#define FARK_DENSE               SUNDIALS_F77_FUNC(farkdense,               FARKDENSE)
#define FARK_DENSESETJAC         SUNDIALS_F77_FUNC(farkdensesetjac,         FARKDENSESETJAC)
#define FARK_MASSDENSE           SUNDIALS_F77_FUNC(farkmassdense,           FARKMASSDENSE)
#define FARK_DENSESETMASS        SUNDIALS_F77_FUNC(farkdensesetmass,        FARKDENSESETMASS)
#define FARK_BAND                SUNDIALS_F77_FUNC(farkband,                FARKBAND)
#define FARK_BANDSETJAC          SUNDIALS_F77_FUNC(farkbandsetjac,          FARKBANDSETJAC)
#define FARK_MASSBAND            SUNDIALS_F77_FUNC(farkmassband,            FARKMASSBAND)
#define FARK_BANDSETMASS         SUNDIALS_F77_FUNC(farkbandsetmass,         FARKBANDSETMASS)
#define FARK_LAPACKDENSE         SUNDIALS_F77_FUNC(farklapackdense,         FARKLAPACKDENSE)
#define FARK_LAPACKDENSESETJAC   SUNDIALS_F77_FUNC(farklapackdensesetjac,   FARKLAPACKDENSESETJAC)
#define FARK_MASSLAPACKDENSE     SUNDIALS_F77_FUNC(farkmasslapackdense,     FARKMASSLAPACKDENSE)
#define FARK_LAPACKDENSESETMASS  SUNDIALS_F77_FUNC(farklapackdensesetmass,  FARKLAPACKDENSESETMASS)
#define FARK_LAPACKBAND          SUNDIALS_F77_FUNC(farklapackband,          FARKLAPACKBAND)
#define FARK_LAPACKBANDSETJAC    SUNDIALS_F77_FUNC(farklapackbandsetjac,    FARKLAPACKBANDSETJAC)
#define FARK_MASSLAPACKBAND      SUNDIALS_F77_FUNC(farkmasslapackband,      FARKMASSLAPACKBAND)
#define FARK_LAPACKBANDSETMASS   SUNDIALS_F77_FUNC(farklapackbandsetmass,   FARKLAPACKBANDSETMASS)
#define FARK_KLU                 SUNDIALS_F77_FUNC(farkklu,                 FARKKLU)
#define FARK_KLUREINIT           SUNDIALS_F77_FUNC(farkklureinit,           FARKKLUREINIT)
#define FARK_SUPERLUMT           SUNDIALS_F77_FUNC(farksuperlumt,           FARKSUPERLUMT)
#define FARK_SPARSESETJAC        SUNDIALS_F77_FUNC(farksparsesetjac,        FARKSPARSESETJAC)
#define FARK_MASSKLU             SUNDIALS_F77_FUNC(farkmassklu,             FARKMASSKLU)
#define FARK_MASSKLUREINIT       SUNDIALS_F77_FUNC(farkmassklureinit,       FARKMASSKLUREINIT)
#define FARK_MASSSUPERLUMT       SUNDIALS_F77_FUNC(farkmasssuperlumt,       FARKMASSSUPERLUMT)
#define FARK_SPARSESETMASS       SUNDIALS_F77_FUNC(farksparsesetmass,       FARKSPARSESETMASS)
#define FARK_SPTFQMR             SUNDIALS_F77_FUNC(farksptfqmr,             FARKSPTFQMR)
#define FARK_SPTFQMRREINIT       SUNDIALS_F77_FUNC(farksptfqmrreinit,       FARKSPTFQMRREINIT)
#define FARK_MASSSPTFQMR         SUNDIALS_F77_FUNC(farkmasssptfqmr,         FARKMASSSPTFQMR)
#define FARK_MASSSPTFQMRREINIT   SUNDIALS_F77_FUNC(farkmasssptfqmrreinit,   FARKMASSSPTFQMRREINIT)
#define FARK_SPBCG               SUNDIALS_F77_FUNC(farkspbcg,               FARKSPBCG)
#define FARK_SPBCGREINIT         SUNDIALS_F77_FUNC(farkspbcgreinit,         FARKSPBCGREINIT)
#define FARK_MASSSPBCG           SUNDIALS_F77_FUNC(farkmassspbcg,           FARKMASSSPBCG)
#define FARK_MASSSPBCGREINIT     SUNDIALS_F77_FUNC(farkmassspbcgreinit,     FARKMASSSPBCGREINIT)
#define FARK_SPGMR               SUNDIALS_F77_FUNC(farkspgmr,               FARKSPGMR)
#define FARK_SPGMRREINIT         SUNDIALS_F77_FUNC(farkspgmrreinit,         FARKSPGMRREINIT)
#define FARK_MASSSPGMR           SUNDIALS_F77_FUNC(farkmassspgmr,           FARKMASSSPGMR)
#define FARK_MASSSPGMRREINIT     SUNDIALS_F77_FUNC(farkmassspgmrreinit,     FARKMASSSPGMRREINIT)
#define FARK_SPFGMR              SUNDIALS_F77_FUNC(farkspfgmr,              FARKSPFGMR)
#define FARK_SPFGMRREINIT        SUNDIALS_F77_FUNC(farkspfgmrreinit,        FARKSPFGMRREINIT)
#define FARK_MASSSPFGMR          SUNDIALS_F77_FUNC(farkmassspfgmr,          FARKMASSSPFGMR)
#define FARK_MASSSPFGMRREINIT    SUNDIALS_F77_FUNC(farkmassspfgmrreinit,    FARKMASSSPFGMRREINIT)
#define FARK_PCG                 SUNDIALS_F77_FUNC(farkpcg,                 FARKPCG)
#define FARK_PCGREINIT           SUNDIALS_F77_FUNC(farkpcgreinit,           FARKPCGREINIT)
#define FARK_MASSPCG             SUNDIALS_F77_FUNC(farkmasspcg,             FARKMASSPCG)
#define FARK_MASSPCGREINIT       SUNDIALS_F77_FUNC(farkmasspcgreinit,       FARKMASSPCGREINIT)
#define FARK_SPILSSETJAC         SUNDIALS_F77_FUNC(farkspilssetjac,         FARKSPILSSETJAC)
#define FARK_SPILSSETPREC        SUNDIALS_F77_FUNC(farkspilssetprec,        FARKSPILSSETPREC)
#define FARK_SPILSSETMASS        SUNDIALS_F77_FUNC(farkspilssetmass,        FARKSPILSSETMASS)
#define FARK_SPILSSETMASSPREC    SUNDIALS_F77_FUNC(farkspilssetmassprec,    FARKSPILSSETMASSPREC)
#define FARK_ARKODE              SUNDIALS_F77_FUNC(farkode,                 FARKODE)
#define FARK_DKY                 SUNDIALS_F77_FUNC(farkdky,                 FARKDKY)
#define FARK_FREE                SUNDIALS_F77_FUNC(farkfree,                FARKFREE)
#define FARK_IMP_FUN             SUNDIALS_F77_FUNC(farkifun,                FARKIFUN)
#define FARK_EXP_FUN             SUNDIALS_F77_FUNC(farkefun,                FARKEFUN)
#define FARK_DJAC                SUNDIALS_F77_FUNC(farkdjac,                FARKDJAC)
#define FARK_DMASS               SUNDIALS_F77_FUNC(farkdmass,               FARKDMASS)
#define FARK_BJAC                SUNDIALS_F77_FUNC(farkbjac,                FARKBJAC)
#define FARK_BMASS               SUNDIALS_F77_FUNC(farkbmass,               FARKBMASS)
#define FARK_SPJAC               SUNDIALS_F77_FUNC(farkspjac,               FARKSPJAC)
#define FARK_SPMASS              SUNDIALS_F77_FUNC(farkspmass,              FARKSPMASS)
#define FARK_PSOL                SUNDIALS_F77_FUNC(farkpsol,                FARKPSOL)
#define FARK_PSET                SUNDIALS_F77_FUNC(farkpset,                FARKPSET)
#define FARK_MASSPSOL            SUNDIALS_F77_FUNC(farkmasspsol,            FARKMASSPSOL)
#define FARK_MASSPSET            SUNDIALS_F77_FUNC(farkmasspset,            FARKMASSPSET)
#define FARK_JTIMES              SUNDIALS_F77_FUNC(farkjtimes,              FARKJTIMES)
#define FARK_MTIMES              SUNDIALS_F77_FUNC(farkmtimes,              FARKMTIMES)
#define FARK_EWT                 SUNDIALS_F77_FUNC(farkewt,                 FARKEWT)
#define FARK_EWTSET              SUNDIALS_F77_FUNC(farkewtset,              FARKEWTSET)
#define FARK_GETERRWEIGHTS       SUNDIALS_F77_FUNC(farkgeterrweights,       FARKGETERRWEIGHTS)
#define FARK_GETESTLOCALERR      SUNDIALS_F77_FUNC(farkgetestlocalerr,      FARKGETESTLOCALERR)
#define FARK_WRITEPARAMETERS     SUNDIALS_F77_FUNC(farkwriteparameters,     FARKWRITEPARAMETERS)
#define FARK_ADAPT               SUNDIALS_F77_FUNC(farkadapt,               FARKADAPT)
#define FARK_ADAPTSET            SUNDIALS_F77_FUNC(farkadaptset,            FARKADAPTSET)
#define FARK_EXPSTAB             SUNDIALS_F77_FUNC(farkexpstab,             FARKEXPSTAB)
#define FARK_EXPSTABSET          SUNDIALS_F77_FUNC(farkexpstabset,          FARKEXPSTABSET)
#define FARK_SETADAPTMETHOD      SUNDIALS_F77_FUNC(farksetadaptivitymethod, FARKSETADAPTIVITYMETHOD)

#else

#define FARK_MALLOC              farkmalloc_
#define FARK_REINIT              farkreinit_
#define FARK_RESIZE              farkresize_
#define FARK_SETDEFAULTS         farksetdefaults_
#define FARK_SETIIN              farksetiin_
#define FARK_SETRIN              farksetrin_
#define FARK_SETERKTABLE         farkseterktable_
#define FARK_SETIRKTABLE         farksetirktable_
#define FARK_SETARKTABLES        farksetarktables_
#define FARK_SETDIAGNOSTICS      farksetdiagnostics_
#define FARK_STOPDIAGNOSTICS     farkstopdiagnostics_
#define FARK_DENSE               farkdense_
#define FARK_DENSESETJAC         farkdensesetjac_
#define FARK_MASSDENSE           farkmassdense_
#define FARK_DENSESETMASS        farkdensesetmass_
#define FARK_BAND                farkband_
#define FARK_BANDSETJAC          farkbandsetjac_
#define FARK_MASSBAND            farkmassband_
#define FARK_BANDSETMASS         farkbandsetmass_
#define FARK_LAPACKDENSE         farklapackdense_
#define FARK_LAPACKDENSESETJAC   farklapackdensesetjac_
#define FARK_MASSLAPACKDENSE     farkmasslapackdense_
#define FARK_LAPACKDENSESETMASS  farklapackdensesetmass_
#define FARK_LAPACKBAND          farklapackband_
#define FARK_LAPACKBANDSETJAC    farklapackbandsetjac_
#define FARK_MASSLAPACKBAND      farkmasslapackband_
#define FARK_LAPACKBANDSETMASS   farklapackbandsetmass_
#define FARK_KLU                 farkklu_
#define FARK_KLUREINIT           farkklureinit_
#define FARK_SUPERLUMT           farksuperlumt_
#define FARK_SPARSESETJAC        farksparsesetjac_
#define FARK_MASSKLU             farkmassklu_
#define FARK_MASSKLUREINIT       farkmassklureinit_
#define FARK_MASSSUPERLUMT       farkmasssuperlumt_
#define FARK_SPARSESETMASS       farksparsesetmass_
#define FARK_SPTFQMR             farksptfqmr_
#define FARK_SPTFQMRREINIT       farksptfqmrreinit_
#define FARK_MASSSPTFQMR         farkmasssptfqmr_
#define FARK_MASSSPTFQMRREINIT   farkmasssptfqmrreinit_
#define FARK_SPBCG               farkspbcg_
#define FARK_SPBCGREINIT         farkspbcgreinit_
#define FARK_MASSSPBCG           farkmassspbcg_
#define FARK_MASSSPBCGREINIT     farkmassspbcgreinit_
#define FARK_SPGMR               farkspgmr_
#define FARK_SPGMRREINIT         farkspgmrreinit_
#define FARK_MASSSPGMR           farkmassspgmr_
#define FARK_MASSSPGMRREINIT     farkmassspgmrreinit_
#define FARK_SPFGMR              farkspfgmr_
#define FARK_SPFGMRREINIT        farkspfgmrreinit_
#define FARK_MASSSPFGMR          farkmassspfgmr_
#define FARK_MASSSPFGMRREINIT    farkmassspfgmrreinit_
#define FARK_PCG                 farkpcg_
#define FARK_PCGREINIT           farkpcgreinit_
#define FARK_MASSPCG             farkmasspcg_
#define FARK_MASSPCGREINIT       farkmasspcgreinit_
#define FARK_SPILSSETJAC         farkspilssetjac_
#define FARK_SPILSSETPREC        farkspilssetprec_
#define FARK_SPILSSETMASS        farkspilssetmass_
#define FARK_SPILSSETMASSPREC    farkmassspilssetprec_
#define FARK_ARKODE              farkode_
#define FARK_DKY                 farkdky_
#define FARK_FREE                farkfree_
#define FARK_IMP_FUN             farkifun_
#define FARK_EXP_FUN             farkefun_
#define FARK_DJAC                farkdjac_
#define FARK_DMASS               farkdmass_
#define FARK_BJAC                farkbjac_
#define FARK_BMASS               farkbmass_
#define FARK_SPJAC               farkspjac_
#define FARK_SPMASS              farkspmass_
#define FARK_PSOL                farkpsol_
#define FARK_PSET                farkpset_
#define FARK_MASSPSOL            farkmasspsol_
#define FARK_MASSPSET            farkmasspset_
#define FARK_JTIMES              farkjtimes_
#define FARK_MTIMES              farkmtimes_
#define FARK_EWT                 farkewt_
#define FARK_EWTSET              farkewtset_
#define FARK_GETERRWEIGHTS       farkgeterrweights_
#define FARK_GETESTLOCALERR      farkgetestlocalerr_
#define FARK_WRITEPARAMETERS     farkwriteparameters_
#define FARK_ADAPT               farkadapt_
#define FARK_ADAPTSET            farkadaptset_
#define FARK_EXPSTAB             farkexpstab_
#define FARK_EXPSTABSET          farkexpstabset_
#define FARK_SETADAPTMETHOD      farksetadaptivitymethod_

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

  void FARK_SETIIN(char key_name[], long int *ival, int *ier);
  void FARK_SETRIN(char key_name[], realtype *rval, int *ier);
  void FARK_SETDEFAULTS(int *ier);
  void FARK_SETERKTABLE(int *s, int *q, int *p, realtype *c, realtype *A, 
			realtype *b, realtype *b2, int *ier);
  void FARK_SETIRKTABLE(int *s, int *q, int *p, realtype *c, 
			realtype *A, realtype *b, realtype *b2, int *ier);
  void FARK_SETARKTABLES(int *s, int *q, int *p, realtype *c, realtype *Ai, 
			 realtype *Ae, realtype *b, realtype *b2, int *ier);
  void FARK_SETDIAGNOSTICS(char fname[], int *flen, int *ier);
  void FARK_STOPDIAGNOSTICS(int *ier);

  void FARK_EWTSET(int *flag, int *ier);
  void FARK_ADAPTSET(int *flag, int *ier);
  void FARK_EXPSTABSET(int *flag, int *ier);

  void FARK_DENSE(long int *neq, int *ier);
  void FARK_DENSESETJAC(int *flag, int *ier);
  void FARK_MASSDENSE(long int *neq, int *ier);
  void FARK_DENSESETMASS(int *ier);

  void FARK_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
  void FARK_BANDSETJAC(int *flag, int *ier);
  void FARK_MASSBAND(long int *neq, long int *mupper, long int *mlower, int *ier);
  void FARK_BANDSETMASS(int *ier);

  void FARK_LAPACKDENSE(int *neq, int *ier);
  void FARK_LAPACKDENSESETJAC(int *flag, int *ier);
  void FARK_MASSLAPACKDENSE(int *neq, int *ier);
  void FARK_LAPACKDENSESETMASS(int *ier);

  void FARK_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
  void FARK_LAPACKBANDSETJAC(int *flag, int *ier);
  void FARK_MASSLAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
  void FARK_LAPACKBANDSETMASS(int *ier);

  void FARK_KLU(int *neq, int *nnz, int *ordering, int *ier);
  void FARK_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier);
  void FARK_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier);
  void FARK_SPARSESETJAC(int *ier);
  void FARK_MASSKLU(int *neq, int *nnz, int *ordering, int *ier);
  void FARK_MASSKLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier);
  void FARK_MASSSUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier);
  void FARK_SPARSESETMASS(int *ier);

  void FARK_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FARK_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);
  void FARK_MASSSPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);

  void FARK_SPFGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FARK_SPFGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);
  void FARK_MASSSPFGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPFGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);

  void FARK_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_SPBCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPBCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPBCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FARK_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_SPTFQMRREINIT(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSSPTFQMRREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FARK_PCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_PCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSPCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_MASSPCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FARK_SPILSSETJAC(int *flag, int *ier);
  void FARK_SPILSSETPREC(int *flag, int *ier);
  void FARK_SPILSSETMASS(int *ier);
  void FARK_SPILSSETMASSPREC(int *flag, int *ier);
  
  void FARK_ARKODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier);

  void FARK_DKY(realtype *t, int *k, realtype *dky, int *ier);

  void FARK_GETERRWEIGHTS(realtype *eweight, int *ier);
  void FARK_GETESTLOCALERR(realtype *ele, int *ier);

  void FARK_FREE(void);

  void FARK_WRITEPARAMETERS(int *ier);

  void FARK_SETADAPTMETHOD(int *imethod, int *idefault, int *ipq, 
			   realtype *params, int *ier);


  /* Prototypes: Functions Called by the ARKODE Solver */
  int FARKfi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  int FARKfe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  
  int FARKDenseJac(long int N, realtype t, 
		   N_Vector y, N_Vector fy, 
		   DlsMat J, void *user_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKDenseMass(long int N, realtype t, DlsMat M, void *user_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKBandJac(long int N, long int mupper, long int mlower,
		  realtype t, N_Vector y, N_Vector fy,
		  DlsMat J, void *user_data,
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKBandMass(long int N, long int mupper, long int mlower,
		   realtype t, DlsMat M, void *user_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKLapackDenseJac(long int N, realtype t, 
			 N_Vector y, N_Vector fy, 
			 DlsMat J, void *user_data,
			 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKLapackDenseMass(long int N, realtype t, DlsMat M, void *user_data,
			  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

  int FARKLapackBandJac(long int N, long int mupper, long int mlower,
			realtype t, N_Vector y, N_Vector fy,
			DlsMat J, void *user_data,
			N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKLapackBandMass(long int N, long int mupper, long int mlower,
			 realtype t, DlsMat M, void *user_data,
			 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKSparseJac(realtype t, N_Vector y, N_Vector fy, 
		    SlsMat J, void *user_data, N_Vector vtemp1, 
		    N_Vector vtemp2, N_Vector vtemp3);

  int FARKSparseMass(realtype t, SlsMat Mass, void *user_data, 
		     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKPSet(realtype tn, N_Vector y,N_Vector fy, booleantype jok,
	       booleantype *jcurPtr, realtype gamma, void *user_data,
	       N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKMassPSet(realtype tn, void *user_data, N_Vector vtemp1, 
		   N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKPSol(realtype tn, N_Vector y, N_Vector fy, 
	       N_Vector r, N_Vector z,
	       realtype gamma, realtype delta,
	       int lr, void *user_data, N_Vector vtemp);
  
  int FARKMassPSol(realtype tn, N_Vector r, N_Vector z, realtype delta, 
		   int lr, void *user_data, N_Vector vtemp);
  
  int FARKJtimes(N_Vector v, N_Vector Jv, realtype t, 
		 N_Vector y, N_Vector fy,
		 void *user_data, N_Vector work);
  
  int FARKMtimes(N_Vector v, N_Vector Mv, realtype t, void *user_data);
  
  int FARKEwt(N_Vector y, N_Vector ewt, void *user_data);

  int FARKAdapt(N_Vector y, realtype t, realtype h1, realtype h2, 
		realtype h3, realtype e1, realtype e2, realtype e3, 
		int q, int p, realtype *hnew, void *user_data);

  int FARKExpStab(N_Vector y, realtype t, realtype *hstab, void *user_data);

  /* Declarations for global variables shared amongst various routines */
  extern N_Vector F2C_ARKODE_vec;     /* defined in FNVECTOR module */

  extern void *ARK_arkodemem;     /* defined in farkode.c */
  extern long int *ARK_iout;      /* defined in farkode.c */
  extern realtype *ARK_rout;      /* defined in farkode.c */
  extern int ARK_nrtfn;           /* defined in farkode.c */
  extern int ARK_ls;              /* defined in farkode.c */
  extern int ARK_mass_ls;         /* defined in farkode.c */

  /* Linear solver IDs */
  enum { ARK_LS_DENSE       = 1, 
	 ARK_LS_BAND        = 2, 
         ARK_LS_LAPACKDENSE = 3, 
	 ARK_LS_LAPACKBAND  = 4,
	 ARK_LS_KLU         = 5, 
	 ARK_LS_SUPERLUMT   = 6, 
	 ARK_LS_SPGMR       = 7, 
	 ARK_LS_SPBCG       = 8, 
	 ARK_LS_SPTFQMR     = 9, 
	 ARK_LS_SPFGMR      = 10, 
	 ARK_LS_PCG         = 12 };

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
