C File: sdiagkbf.f
C Sensitivity version of diagonal ODE example.  
C Stiff case, with diagonal preconditioner.
C Uses SENSFPVODE interfaces and SENSFPVBBD interfaces.
C Version of 25 August 2000
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C Include MPI-Fortran header file for MPI_COMM_WORLD, MPI_MAX,
C                                     MPI_DOUBLE_PRECISION
      INCLUDE "mpif.h"
C
      INTEGER*4 IOPT(40)
      DIMENSION Y(1024), ROPT(40)
      DIMENSION P(10), PBAR(10)
      DATA ATOL/1.0E-10/, RTOL/1.0E-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LQ/11/, LH/5/, LNPE/14/, LNLI/15/, LNPS/16/, LNCFL/17/
C
      COMMON /PCOM/ MYPE
C
C Get NPES and MYPE.  Requires initialization of MPI.
      CALL MPI_INIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,5) IER
 5      FORMAT(///' MPI_INIT returned IER =',I5)
        STOP
        ENDIF
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPES, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,6) IER
 6      FORMAT(///' MPI_COMM_SIZE returned IER =',I5)
        STOP
        ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,7) IER
 7      FORMAT(///' MPI_COMM_RANK returned IER =',I5)
        STOP
        ENDIF
C
C Set input arguments.
      NLOCAL = 10
      NY = NPES*NLOCAL
C
      NS = 1
      NTOTAL = (1+NS)*NY
C
      T = 0.0D0
      METH = 2
      ITMETH = 2
      IATOL = 1
      INOPT = 0
      ITASK = 0
      IPRE = 1
      IGS = 0
C Set parameter alpha
      ALPHA  = 10.0D0

C Set real parameter array and scaling factor(s)
      P(1) = ALPHA
      PBAR(1) = ALPHA
C
      RHOMAX = 0.0
C
      DO 10 I = 1, NLOCAL
  10    Y(I) = 1.0D0
C
C Initialize Sensitivity Variables
C
      DO 9 J = 1, NS
      DO 8 I = 1, NLOCAL
        Y(J*NLOCAL + I) = 0.0D0
 8    CONTINUE
 9    CONTINUE
C
      IF (MYPE .EQ. 0) THEN
        WRITE(6,11) NY, ALPHA
  11    FORMAT('Diagonal test problem, size NY =',I5,
     1         '  parameter alpha =',F8.3)
        WRITE(6,12)
  12    FORMAT('  ydot_i = -alpha*i * y_i (i = 1,...,NY)'/)
        WRITE(6,13)RTOL,ATOL
  13    FORMAT('RTOL, ATOL =',2E10.1/)
        WRITE(6,14)
  14    FORMAT('Method is BDF/NEWTON/SPGMR'/
     1         'Preconditioner is band-block-diagonal, using PVBBDPRE'/)
        WRITE(6,15)NPES
  15    FORMAT('Number of processors =',I3)
        WRITE(6,16)NS
  16    FORMAT('Number of sensitivity vectors =',i3)
        WRITE(6,17)NTOTAL
  17    FORMAT('Total number of equations =',i3)
        ENDIF
C
      CALL FPVINITMPI (NLOCAL, NY, IER)
C
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
  20    FORMAT(///' FPVINITMPI returned IER =',I5)
        STOP
        ENDIF
C
      CALL SFPVMALLOC (NY, NS, NTOTAL, T, Y, METH, ITMETH, IATOL,
     1                RTOL, ATOL, INOPT, IOPT, ROPT, IER, 
     2                P, PBAR, RHOMAX)
C
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
  30    FORMAT(///' SFPVMALLOC returned IER =',I5)
        STOP
        ENDIF
C
      MUDQ = 0
      MLDQ = 0
      MU = 0
      ML = 0
      CALL SFPVBBDIN (MUDQ,MLDQ, MU,ML, 0.0D0, IPRE, IGS, 0, 0.0D0, IER)
        IF (IER .NE. 0) THEN
          WRITE(6,35) IER
  35      FORMAT(///' SFPVBBDIN returned IER =',I5)
          STOP
          ENDIF
C
C Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 70 IOUT = 1, NOUT
C
        CALL FCVODE (TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T,IOPT(LNST),IOPT(LNFE)
  40    FORMAT(/' t =',E10.2,5X,'no. steps =',I5,'   no. f-s =',I5)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER
  60      FORMAT(///' FCVODE returned IER =',I5)
          STOP
          ENDIF
C
        TOUT = TOUT + DTOUT
  70    CONTINUE
C
C Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 75 I = 1, NLOCAL
        ERRI  = Y(I) - EXP(-P(1)*(MYPE*NLOCAL + I)*T)
        ERMAX = MAX(ERMAX,ABS(ERRI))
  75    CONTINUE
C
C Get max. absolute error in each local scaled sensitivity vector.
      DO 78 J = 1, NS
      DO 77 I = 1, NLOCAL
        ERRI  = Y(NS*NLOCAL+I) -
     1        (-P(1)*(MYPE*NLOCAL+I)*T*EXP(-P(1)*(MYPE*NLOCAL+I)*T))
        ERMAX = MAX(ERMAX,ABS(ERRI))        
  77  CONTINUE
  78  CONTINUE
C
C Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE (ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     1                 0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,80) IER
  80    FORMAT(///' MPI_REDUCE returned IER =',I5)
        STOP
        ENDIF
      IF (MYPE .EQ. 0) WRITE(6,85) GERMAX
  85  FORMAT(//'Max. absolute error is',E10.2)
C
C Print final statistics.
      IF (MYPE .EQ. 0) THEN
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
      WRITE (6,90) NST,NFE,NPSET,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL,NETF
  90  FORMAT(//'Final statistics..'/
     1 ' number of steps        =',I5,5X,'number of f evals.     =',I5/
     2 ' number of prec. setups =',I5/
     3 ' number of prec. evals. =',I5,5X,'number of prec. solves =',I5/
     4 ' number of nonl. iters. =',I5,5X,'number of lin. iters.  =',I5/
     5 ' average Krylov subspace dimension (NLI/NNI)  =',F8.4/
     6 ' number of conv. failures..  nonlinear =',I3,'  linear =',I3/
     7 ' number of error test failures =',I3)
      CALL FPVBBDOPT (LENRPW, LENIPW, NGE)
      WRITE (6,92) LENRPW, LENIPW, NGE
  92  FORMAT(/'In PVBBDPRE: real/integer local work space sizes =',2I5/
     1 ' number of g evals.     =',I5)
      ENDIF
C
C Free the memory and finalize MPI.
      CALL FPVBBDF
      CALL SFCVFREE
      CALL FPVFREEMPI
      CALL MPI_FINALIZE(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,95) IER
 95     FORMAT(///' MPI_FINALIZE returned IER =',I5)
        STOP
        ENDIF
C
      STOP
      END

      SUBROUTINE PVFUN (NLOC, T, Y, YDOT, P)
C Routine for right-hand side function f, with parameters
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION Y(*), YDOT(*)
      DIMENSION P(*)
      COMMON /PCOM/ MYPE
C
      DO 10 I = 1,NLOC
  10    YDOT(I) = -P(1)*(MYPE*NLOC + I)*Y(I)
C
      RETURN
      END

      SUBROUTINE PVLOCFN (NLOC, T, YLOC, GLOC, P)
C Routine to define local approximate function g, here the same as f. 
C Include parameters in the call sequence
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION YLOC(*), GLOC(*)
      DIMENSION P(*)
C
      CALL PVFUN(NLOC, T, YLOC, GLOC, P)
C
      RETURN
      END

      SUBROUTINE PVCOMMF (NLOC, T, YLOC)
C Routine to perform communication required for evaluation of g.
      RETURN
      END
