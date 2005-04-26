C     ----------------------------------------------------------------
C     $Revision: 1.20 $
C     $Date: 2005-04-26 23:43:21 $
C     ----------------------------------------------------------------
C     Diagonal ODE example.  Stiff case, with diagonal preconditioner.
C     Uses FCVODE interfaces and FCVBBD interfaces.
C     Solves problem twice -- with left and right preconditioning.
C     ----------------------------------------------------------------
C     
C     Include MPI-Fortran header file for MPI_COMM_WORLD, MPI types.
      
      INCLUDE "mpif.h"
C     
      INTEGER NOUT, LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNPE
      INTEGER LNLI, LNPS, LNCFL, MYPE, IER, NPES, METH, ITMETH
      INTEGER IATOL, INOPT, ITASK, IPRE, IGS, IOUT
      INTEGER*4 IOPT(40)
      INTEGER*4 NEQ, NLOCAL, I, MUDQ, MLDQ, MU, ML, NETF
      INTEGER*4 NST, NFE, NPSET, NPE, NPS, NNI, NLI, NCFN, NCFL
      INTEGER*4 LENRPW, LENIPW, NGE
      DOUBLE PRECISION ALPHA, TOUT, ERMAX, AVDIM
      DOUBLE PRECISION ATOL, ERRI, RTOL, GERMAX, DTOUT, Y, ROPT, T
      DIMENSION Y(1024), ROPT(40)
C     
      DATA ATOL/1.0D-10/, RTOL/1.0D-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LNPE/18/, LNLI/19/, LNPS/20/, LNCFL/21/
C
      COMMON /PCOM/ ALPHA, NLOCAL, MYPE
C     
C     Get NPES and MYPE.  Requires initialization of MPI.
      CALL MPI_INIT(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,5) IER
 5       FORMAT(///' MPI_ERROR: MPI_INIT returned IER = ', I5)
         STOP
      ENDIF
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPES, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,6) IER
 6       FORMAT(///' MPI_ERROR: MPI_COMM_SIZE returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYPE, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,7) IER
 7       FORMAT(///' MPI_ERROR: MPI_COMM_RANK returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      
C     
C     Set input arguments.
      NLOCAL = 10
      NEQ = NPES * NLOCAL
      T = 0.0D0
      METH = 2
      ITMETH = 2
      IATOL = 1
      INOPT = 0
      ITASK = 1
      IPRE = 1
      IGS = 1
C     Set parameter alpha
      ALPHA  = 10.0D0
C     
      DO I = 1, NLOCAL
         Y(I) = 1.0D0
      ENDDO
C     
      IF (MYPE .EQ. 0) THEN
         WRITE(6,15) NEQ, ALPHA, RTOL, ATOL, NPES
 15      FORMAT('Diagonal test problem:'//' NEQ = ', I3, /
     &          ' parameter alpha = ', F8.3/
     &          ' ydot_i = -alpha*i * y_i (i = 1,...,NEQ)'/
     &          ' RTOL, ATOL = ', 2E10.1/
     &          ' Method is BDF/NEWTON/SPGMR'/
     &          ' Preconditioner is band-block-diagonal, using CVBBDPRE'
     &          /' Number of processors = ', I3/)
      ENDIF
C     
      CALL FNVINITP(MPI_COMM_WORLD, 1, NLOCAL, NEQ, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
 20      FORMAT(///' SUNDIALS_ERROR: FNVINITP returned IER = ', I5)
         CALL MPI_FINALIZE(IER)
         STOP
      ENDIF
C     
      CALL FCVMALLOC(T, Y, METH, ITMETH, IATOL, RTOL, ATOL,
     &               INOPT, IOPT, ROPT, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
 30      FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C     
      MUDQ = 0
      MLDQ = 0
      MU = 0
      ML = 0
      CALL FCVBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, 0.0D0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,35) IER
 35      FORMAT(///' SUNDIALS_ERROR: FCVBBDINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      CALL FCVBBDSPGMR(IPRE, IGS, 0, 0.0D0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,36) IER
 36      FORMAT(///' SUNDIALS_ERROR: FCVBBDSPGMR returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IF (MYPE .EQ. 0) WRITE(6,38)
 38   FORMAT(/'Preconditioning on left'/)
C     
C     Looping point for cases IPRE = 1 and 2.
C     
 40   CONTINUE
C     
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 60 IOUT = 1, NOUT
C     
         CALL FCVODE(TOUT, T, Y, ITASK, IER)
C     
         IF (MYPE .EQ. 0) WRITE(6,45) T, IOPT(LNST), IOPT(LNFE)
 45      FORMAT(' t = ', E10.2, 5X, 'no. steps = ', I5,
     &          '   no. f-s = ', I5)
C     
         IF (IER .NE. 0) THEN
            WRITE(6,50) IER, IOPT(26)
 50         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     &             '                 Linear Solver returned IER = ', I5)
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
            STOP
         ENDIF
C     
         TOUT = TOUT + DTOUT
 60   CONTINUE
C     
C     Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 65 I = 1, NLOCAL
         ERRI  = Y(I) - EXP(-ALPHA * (MYPE * NLOCAL + I) * T)
         ERMAX = MAX(ERMAX, ABS(ERRI))
 65   CONTINUE
C     Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE(ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     &                0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,70) IER
 70      FORMAT(///' MPI_ERROR: MPI_REDUCE returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      IF (MYPE .EQ. 0) WRITE(6,75) GERMAX
 75   FORMAT(/'Max. absolute error is', E10.2/)
C     
C     Print final statistics.
      IF (MYPE .EQ. 0) THEN
         NST = IOPT(LNST)
         NFE = IOPT(LNFE)
         NPSET = IOPT(LNSETUP)
         NPE = IOPT(LNPE)
         NPS = IOPT(LNPS)
         NNI = IOPT(LNNI)
         NLI = IOPT(LNLI)
         AVDIM = DBLE(NLI) / DBLE(NNI)
         NCFN = IOPT(LNCF)
         NCFL = IOPT(LNCFL)
         NETF = IOPT(LNETF)
         WRITE(6,80) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &               NCFL, NETF
 80      FORMAT(/'Final statistics:'//
     &          ' number of steps        = ', I5, 4X,
     &          ' number of f evals.     = ', I5/
     &          ' number of prec. setups = ', I5/
     &          ' number of prec. evals. = ', I5, 4X,
     &          ' number of prec. solves = ', I5/
     &          ' number of nonl. iters. = ', I5, 4X,
     &          ' number of lin. iters.  = ', I5/
     &          ' average Krylov subspace dimension (NLI/NNI) = ', F8.4/
     &          ' number of conv. failures.. nonlinear = ', I3,
     &          '  linear = ', I3/
     &          ' number of error test failures = ', I3/)
         CALL FCVBBDOPT(LENRPW, LENIPW, NGE)
         WRITE(6,82) LENRPW, LENIPW, NGE
 82      FORMAT('In CVBBDPRE:'//
     &          ' real/int local workspace = ', 2I5/
     &          ' number of g evals. = ', I5)
      ENDIF
C     
C     If IPRE = 1, re-initialize T, Y, and the solver, and loop for case IPRE = 2.
C     Otherwise jump to final block.
      IF (IPRE .EQ. 2) GO TO 99
C
      T = 0.0D0
      DO I = 1, NLOCAL
         Y(I) = 1.0D0
      ENDDO         
C
      CALL FCVREINIT(T, Y, IATOL, RTOL, ATOL,INOPT, IOPT, ROPT, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,91) IER
 91      FORMAT(///' SUNDIALS_ERROR: FCVREINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IPRE = 2
C
      CALL FCVBBDREINIT(NLOCAL, MUDQ, MLDQ, 0.0D0, IER) 
      IF (IER .NE. 0) THEN
         WRITE(6,92) IER
 92      FORMAT(///' SUNDIALS_ERROR: FCVBBDREINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      CALL FCVSPGMRREINIT(IPRE, IGS, 0.0D0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,93) IER
 93      FORMAT(///' SUNDIALS_ERROR: FCVSPGMRREINIT returned IER = ',I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IF (MYPE .EQ. 0) WRITE(6,95)
 95   FORMAT(//60('-')///'Preconditioning on right'/)
      GO TO 40
C     
C     Free the memory and finalize MPI.
 99   CALL FCVBBDFREE
      CALL FCVFREE
      CALL MPI_FINALIZE(IER)
C     
      STOP
      END
C
      SUBROUTINE FCVFUN(T, Y, YDOT)
C     Routine for right-hand side function f
C
      IMPLICIT NONE
C
      INTEGER MYPE
      INTEGER*4 I, NLOCAL
      DOUBLE PRECISION Y, YDOT, ALPHA, T
      DIMENSION Y(*), YDOT(*)
C
      COMMON /PCOM/ ALPHA, NLOCAL, MYPE
C     
      DO I = 1, NLOCAL
         YDOT(I) = -ALPHA * (MYPE * NLOCAL + I) * Y(I)
      ENDDO
C     
      RETURN
      END
C
      SUBROUTINE FCVGLOCFN(NLOC, T, YLOC, GLOC)
C     Routine to define local approximate function g, here the same as f. 
      IMPLICIT NONE
C
      INTEGER*4 NLOC
      DOUBLE PRECISION YLOC, GLOC, T
      DIMENSION YLOC(*), GLOC(*)
C     
      CALL FCVFUN(T, YLOC, GLOC)
C     
      RETURN
      END
      
      SUBROUTINE FCVCOMMFN(NLOC, T, YLOC)
C     Routine to perform communication required for evaluation of g.
      RETURN
      END
