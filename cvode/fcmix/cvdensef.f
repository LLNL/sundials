C File: cvdensef.f
C ODE example with BDF/DENSE (serial NVECTOR).
C Version of 28 March 2002
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     
      INTEGER*4 IOPT(40)
      DIMENSION Y(3), ATOL(3), ROPT(40)
      DATA Y1/1.0/, Y2/0.0/, Y3/0.0/
      DATA ATOL1/1.0E-8/, ATOL2/1.0E-14/, ATOL3/1.0E-6/, RTOL/1.0E-4/
      DATA NOUT/12/
      DATA T1/0.4/, TMULT/10.0/
      DATA LNST/4/, LNFE/5/, LNNI/7/, LNCF/8/, LNETF/9/, LNCFL/19/
C     
C     
C     Set input arguments.
      NEQ = 3
      T = 0.0D0
      METH = 2
      ITMETH = 2
      IATOL = 2
      INOPT = 0
      ITASK = 0
      IPRE = 1
      IGS = 0
C     
      Y(1) = Y1
      Y(2) = Y2
      Y(3) = Y3
      
      ATOL(1) = ATOL1
      ATOL(2) = ATOL2
      ATOL(3) = ATOL3
      
      WRITE(6, 10) NEQ
 10   FORMAT('BDF/DENSE test problem, size NEQ =',I3)

C     
      CALL FMENVINITS(NEQ, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
 20      FORMAT(///' FMENVINITS returned IER =',I5)
         STOP
      ENDIF

C     
      CALL FCVMALLOC(NEQ, T, Y, METH, ITMETH, IATOL, RTOL, ATOL,
     &     INOPT, IOPT, ROPT, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
 30      FORMAT(///' FCVMALLOC returned IER =',I5)
         STOP
      ENDIF

C     
      CALL FCVDENSE0(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,35) IER
 35      FORMAT(///' FCVDENSE0 returned IER =',I5)
         STOP
      ENDIF

C     
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = T1
      DO 70 IOUT = 1, NOUT
C     
         CALL FCVODE(TOUT, T, Y, ITASK, IER)
C     
         WRITE(6,40) T, Y(1), Y(2), Y(3)
 40      FORMAT(/' t =', E10.2, 5X, 'Y =', E14.6, E14.6, E14.6)
C     
         IF (IER .NE. 0) THEN
            WRITE(6,60) IER
 60         FORMAT(///' FCVODE returned IER =',I5)
            STOP
         ENDIF
C     
         TOUT = TOUT * TMULT
 70   CONTINUE
C     
C     Print final statistics.
      NST = IOPT(LNST)
      NFE = IOPT(LNFE)
      NNI = IOPT(LNNI)
      NCFN = IOPT(LNCF)
      NCFL = IOPT(LNCFL)
      NETF = IOPT(LNETF)

      WRITE (6,90) NST,NFE,NNI,NCFN,NCFL,NETF

 90   FORMAT(//' Final statistics..'/
     &     ' number of steps        =',I5/
     &     ' number of f evals.     =',I5/
     &     ' number of nonl. iters. =',I5/
     &     ' number of conv. failures..  nonlinear =',I3,
     &                                   '  linear =',I3/
     &     ' number of error test failures =',I3)
C     
C     Free the memory.
      CALL FCVFREE
      CALL FMENVFREES
C     
      STOP
      END
      
C===============================================================================
C===============================================================================

      SUBROUTINE CVFUN (NEQ, T, Y, YDOT)
C     Routine for right-hand side function f
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION Y(*), YDOT(*)
      
      YDOT(1) = -0.04 * Y(1) + 1.0E4 * Y(2) * Y(3)
      YDOT(3) = 3.0E7 * Y(2) * Y(2)
      YDOT(2) = - YDOT(1) - YDOT(3)

      RETURN
      END

