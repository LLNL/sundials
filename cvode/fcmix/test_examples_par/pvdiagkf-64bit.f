C ------------------------------------------------------------------
C $Revision: 1.2 $
C $Date: 2004-05-17 18:50:18 $
C ------------------------------------------------------------------
C Diagonal ODE example. Stiff case, with BDF/SPGMR, diagonal
C preconditioner. Solved with preconditioning on left, then with
C preconditioning on right.
C ------------------------------------------------------------------
C
C Include MPI-Fortran header file for MPI_COMM_WORLD, MPI types.
C
      INCLUDE "mpif.h"
C
      INTEGER LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNPE, LNLI, LNPS
      INTEGER LNCFL, NOUT, MYPE, NPES, IER, METH, ITMETH, IATOL, INOPT
      INTEGER ITASK, IPRE, IGS, IOUT
      INTEGER*8 IOPT(40)
      INTEGER*8 NEQ, NLOCAL, I, NST, NFE, NPSET, NPE, NPS, NNI, NLI
      INTEGER*8 NCFL, NETF, NCFN
      DOUBLE PRECISION Y(1024), ROPT(40)
      DOUBLE PRECISION ATOL, DTOUT, T, ALPHA, RTOL, TOUT, ERMAX, ERRI
      DOUBLE PRECISION GERMAX, AVDIM
      DATA ATOL/1.0E-10/, RTOL/1.0E-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LNPE/18/, LNLI/19/, LNPS/20/, LNCFL/21/
C
      COMMON /PCOM/ NLOCAL, ALPHA, MYPE
C
C Get NPES and MYPE.  Requires initialization of MPI.
      CALL MPI_INIT (IER)
      IF (IER .NE. 0) THEN
        WRITE(6,5) IER
 5      FORMAT(///' MPI_ERROR: MPI_INIT returned IER =',I5)
        STOP
        ENDIF
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD, NPES, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,6) IER
 6      FORMAT(///' MPI_ERROR: MPI_COMM_SIZE returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
        ENDIF
      CALL MPI_COMM_RANK (MPI_COMM_WORLD, MYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,7) IER
 7      FORMAT(///' MPI_ERROR: MPI_COMM_RANK returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
        ENDIF
C
C Set input arguments.
      NLOCAL = 10
      NEQ = NPES*NLOCAL
      T = 0.0D0
      METH = 2
      ITMETH = 2
      IATOL = 1
      INOPT = 0
      ITASK = 1
      IPRE = 1
      IGS = 1
C Set parameter alpha.
      ALPHA  = 10.0D0
C
C Do remaining initializations for first case: IPRE = 1 (prec. on left).
C
      DO 10 I = 1,NLOCAL
  10    Y(I) = 1.0D0
C
      IF (MYPE .EQ. 0) THEN
        WRITE(6,11) NEQ, ALPHA
  11    FORMAT('Diagonal test problem, size NEQ =',I5,
     1         '  parameter alpha =',F8.3)
        WRITE(6,12)
  12    FORMAT('  ydot_i = -alpha*i * y_i (i = 1,...,NEQ)'/)
        WRITE(6,13)RTOL,ATOL
  13    FORMAT('RTOL, ATOL =',2E10.1/)
        WRITE(6,14)
  14    FORMAT('Method is BDF/NEWTON/SPGMR'/
     1         'Diagonal preconditioner uses approximate Jacobian'/)
        WRITE(6,15)NPES
  15    FORMAT('Number of processors =',I3)
        WRITE(6,16)
  16    FORMAT(///'Case 1: preconditioning on left')
        ENDIF
C
      CALL FNVINITP (NLOCAL, NEQ, IER)
C
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
  20    FORMAT(///' SUNDIALS_ERROR: FMENVINITP returned IER =',I5)
        CALL MPI_FINALIZE(IER)
        STOP
        ENDIF
C
      CALL FCVMALLOC (T, Y, METH, ITMETH, IATOL, RTOL, ATOL,
     1                INOPT, IOPT, ROPT, IER)
C
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
  30    FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
        ENDIF
C
      CALL FCVSPGMR (IPRE, IGS, 0, 0.0D0, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,35) IER
  35    FORMAT(///' SUNDIALS_ERROR: FCVSPGMR returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
        ENDIF
C
        CALL FCVSPGMRSETPSOL (1, IER)

C Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 70 IOUT = 1,NOUT
C
        CALL FCVODE (TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T,IOPT(LNST),IOPT(LNFE)
  40    FORMAT(/' t =',E10.2,5X,'no. steps =',I5,'   no. f-s =',I5)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER
  60      FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER =',I5)
          CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
          STOP
          ENDIF
C
        TOUT = TOUT + DTOUT
  70    CONTINUE
C
C Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 75 I = 1,NLOCAL
        ERRI  = Y(I) - EXP(-ALPHA*(MYPE*NLOCAL + I)*T)
  75    ERMAX = MAX(ERMAX,ABS(ERRI))
C Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE (ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     1                 0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,80) IER
  80    FORMAT(///' MPI_ERROR: MPI_REDUCE returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
        ENDIF
      IF (MYPE .EQ. 0) WRITE(6,85) GERMAX
  85  FORMAT(//'Max. absolute error is',E10.2)
C
C Print final statistics.
      NST = IOPT(LNST)
      NFE = IOPT(LNFE)
      NPSET = IOPT(LNSETUP)
      NPE = IOPT(LNPE)
      NPS = IOPT(LNPS)
      NNI = IOPT(LNNI)
      NLI = IOPT(LNLI)
      AVDIM = DBLE(NLI)/DBLE(NNI)
      NCFN = IOPT(LNCF)
      NCFL = IOPT(LNCFL)
      NETF = IOPT(LNETF)
      IF (MYPE .EQ. 0)
     1  WRITE (6,90) NST,NFE,NPSET,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL,NETF
  90  FORMAT(//' Final statistics..'/
     1 ' number of steps        =',I5,5X,'number of f evals.     =',I5/
     2 ' number of prec. setups =',I5/
     3 ' number of prec. evals. =',I5,5X,'number of prec. solves =',I5/
     4 ' number of nonl. iters. =',I5,5X,'number of lin. iters.  =',I5/
     5 ' average Krylov subspace dimension (NLI/NNI)  =',F8.4/
     6 ' number of conv. failures..  nonlinear =',I3,'  linear =',I3/
     7 ' number of error test failures =',I3)
C
C Re-initialize to run second case: IPRE = 2 (prec. on right).
      IPRE = 2
      T = 0.0D0
      DO 110 I = 1,NLOCAL
 110    Y(I) = 1.0D0
C
      IF (MYPE .EQ. 0)  WRITE(6,111) 
 111  FORMAT(///'Case 2: preconditioning on right')
C
      CALL FCVREINIT (T, Y, IATOL, RTOL, ATOL,
     1                INOPT, IOPT, ROPT, IER)
C
      IF (IER .NE. 0) THEN
        WRITE(6,130) IER
 130    FORMAT(///' SUNDIALS_ERROR: FCVREINIT returned IER =',I5)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
      ENDIF
C
      CALL FCVSPGMRREINIT (IPRE, IGS, 0.0D0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,140) IER
 140     FORMAT(///' SUNDIALS_ERROR: FCVSPGMRREINIT returned IER =',I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
C Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 170 IOUT = 1,NOUT
C
        CALL FCVODE (TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T,IOPT(LNST),IOPT(LNFE)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER
          CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
          STOP
          ENDIF
C
        TOUT = TOUT + DTOUT
 170    CONTINUE
C
C Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 175 I = 1,NLOCAL
        ERRI  = Y(I) - EXP(-ALPHA*(MYPE*NLOCAL + I)*T)
 175    ERMAX = MAX(ERMAX,ABS(ERRI))
C Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE (ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     1                 0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,80) IER
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
        STOP
      ENDIF
      IF (MYPE .EQ. 0) WRITE(6,85) GERMAX
C
C Print final statistics.
      NST = IOPT(LNST)
      NFE = IOPT(LNFE)
      NPSET = IOPT(LNSETUP)
      NPE = IOPT(LNPE)
      NPS = IOPT(LNPS)
      NNI = IOPT(LNNI)
      NLI = IOPT(LNLI)
      AVDIM = REAL(NLI)/REAL(NNI)
      NCFN = IOPT(LNCF)
      NCFL = IOPT(LNCFL)
      NETF = IOPT(LNETF)
      IF (MYPE .EQ. 0)
     1  WRITE (6,90) NST,NFE,NPSET,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL,NETF
C
C Free the memory and finalize MPI.
      CALL FCVFREE
      CALL FNVFREEP
      CALL MPI_FINALIZE(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,195) IER
 195    FORMAT(///' MPI_ERROR: MPI_FINALIZE returned IER =',I5)
        STOP
        ENDIF
C
      STOP
      END

      SUBROUTINE FCVFUN (T, Y, YDOT)
C Routine for right-hand side function f
      IMPLICIT NONE
      INTEGER MYPE
      INTEGER*8 I, NLOCAL
      DOUBLE PRECISION YDOT, ALPHA, Y, T
      DIMENSION Y(*), YDOT(*)
      COMMON /PCOM/ NLOCAL, ALPHA, MYPE
C
      DO 10 I = 1,NLOCAL
  10    YDOT(I) = -ALPHA*(MYPE*NLOCAL + I)*Y(I)
C
      RETURN
      END

      SUBROUTINE FCVPSOL(T, Y, FY, VT, GAMMA, EWT, DELTA,
     1                  R, LR, Z, IER)
C Routine to solve preconditioner linear system
C This routine uses a diagonal preconditioner P = I - gamma*J,
C where J is a diagonal approximation to the true Jacobian, given by:
C J = diag(0, 0, 0, -4*alpha, ..., -N*alpha).
C The vector r is copied to z, and the inverse of P (restricted to the
C local vector segment) is applied to the vector z.
      IMPLICIT NONE
      INTEGER IER, MYPE
      INTEGER*8 I, NLOCAL, ISTART, IBASE, LR
      DOUBLE PRECISION PSUBI, GAMMA, ALPHA, Y, T, FY, VT, EWT, DELTA
      DOUBLE PRECISION Z, R
      DIMENSION Y(*), R(*), Z(*)
      COMMON /PCOM/ NLOCAL, ALPHA, MYPE
C
      DO 10 I = 1,NLOCAL
 10     Z(I) = R(I)
C
      IBASE = MYPE*NLOCAL
      ISTART = MAX(1, 4 - IBASE)
      DO 20 I = ISTART,NLOCAL
        PSUBI = 1.0 + GAMMA*ALPHA*(IBASE + I)
        Z(I) = Z(I)/PSUBI
 20     CONTINUE
C
      RETURN
      END
