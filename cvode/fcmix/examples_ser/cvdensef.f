C File cvdensef.f
C
C FCVODE Example Problem:  Robertson kinetics, dense user Jacobian.
C Version of 27 January 2004
C
C The following is a simple example problem, with the coding
C needed for its solution by CVODE.  The problem is from chemical
C kinetics, and consists of the following three rate equations..
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C
C The following coding solves this problem with CVODE, using the
C Fortran/C interface routine package.  This solution uses the BDF
C method and user-supplied Jacobian routine, and prints results at
C t = .4, 4., ..., 4.e10.  It uses ITOL = 2 and ATOL much smaller for
C y2 than y1 or y3 because y2 has much smaller values.
C At the end of the run, various counters of interest are printed.
C
      DOUBLE PRECISION RTOL, T, T0, TOUT
      DOUBLE PRECISION Y(3), ATOL(3), ROPT(40)
      INTEGER IOPT(40)
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LNJE/16/
C
      NEQ = 3
      T0 = 0.0D0
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      METH = 2
      ITMETH = 2
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-8
      ATOL(2) = 1.D-14
      ATOL(3) = 1.D-6
      INOPT = 0
      TOUT = 0.4D0
      ITASK = 0
C
      WRITE(6,10) NEQ
 10   FORMAT('Dense example problem: Robertson kinetics, NEQ = ',I2//)
C
      CALL FNVINITS (NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER =',I5)
        STOP
      ENDIF
C
      CALL FCVMALLOC (T0, Y, METH, ITMETH, ITOL, RTOL, ATOL,
     1                INOPT, IOPT, ROPT, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER =',I5)
        CALL FNVFREES
        STOP
      ENDIF
C
      CALL FCVDENSE (NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FCVDENSE returned IER =',I5)
        CALL FNVFREES
        CALL FCVFREE
        STOP
      ENDIF
C
      CALL FCVDENSESETJAC(1, IER)
C
      DO 70 IOUT = 1,12
C
        CALL FCVODE (TOUT, T, Y, ITASK, IER)
C
        WRITE(6,50)T,Y(1),Y(2),Y(3)
 50     FORMAT('At t =',D12.4,'   y =',3D14.6)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER
 60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER =',I5)
          CALL FNVFREES
          CALL FCVFREE
          STOP
        ENDIF
C
        TOUT = TOUT*10.0D0
 70   CONTINUE
C
      CALL FCVDKY (T, 1, Y, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,80) IER
 80     FORMAT(///' SUNDIALS_ERROR: FCVDKY returned IER =',I5)
        CALL FNVFREES
        CALL FCVFREE
        STOP
      ENDIF
      WRITE(6,85)y(1),y(2),y(3)
 85   FORMAT(/'Final value of ydot =',3D11.3)
C
      WRITE(6,90) IOPT(LNST), IOPT(LNFE), IOPT(LNJE), IOPT(LNSETUP),
     1            IOPT(LNNI), IOPT(LNCF), IOPT(LNETF)
 90   FORMAT(/'No. steps =',I4,'   No. f-s =',I4,
     1       '   No. J-s =',I4,'   No. LU-s =',I4/
     2       'No. nonlinear iterations =',I4/
     3       'No. nonlinear convergence failures =',I4/
     4       'No. error test failures =',I4/)
C
      CALL FCVFREE
      CALL FNVFREES
C
      STOP
      END

      SUBROUTINE FCVFUN (T, Y, YDOT)
C Fortran routine for right-hand side function.
      DOUBLE PRECISION T, Y(*), YDOT(*)
      YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
      YDOT(3) = 3.D7*Y(2)*Y(2)
      YDOT(2) = -YDOT(1) - YDOT(3)
      RETURN
      END

      SUBROUTINE FCVDJAC(N, T, Y, FY, JAC, V1, V2, V3)
C Fortran routine for dense user-supplied Jacobian.
      DOUBLE PRECISION T, Y(*), JAC(N,*), Y1, Y2, Y3
      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)
      JAC(1,1) = -0.04D0
      JAC(1,2) = 1.D4*Y3
      JAC(1,3) = 1.D4*Y2
      JAC(2,1) =  0.04D0
      JAC(2,2) = -1.D4*Y3 - 6.D7*Y2
      JAC(2,3) = -1.D4*Y2
      JAC(3,2) = 6.D7*Y2
      RETURN
      END
