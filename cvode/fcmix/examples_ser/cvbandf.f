C File cvbandf.f
C
C FCVODE Example Problem: Advection-diffusion, banded user Jacobian.
C Version of 30 March 2003
C
C The following is a simple example problem with a banded Jacobian.
C The problem is the semi-discrete form of the advection-diffusion
C equation in 2D:
C   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
C on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
C interval 0 <= t <= 1.  Homogeneous Dirichlet boundary conditions
C are posed, and the initial condition is
C   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy) .
C The PDE is discretized on a uniform MX+2 by MY+2 grid with
C central differencing, and with boundary values eliminated,
C leaving an ODE system of size NEQ = MX*MY.
C This program solves this problem with CVODE, using the Fortran/C 
C interface routine package.  This solution uses the BDF method,
C a user-supplied banded Jacobian routine, and scalar relative and
C absolute tolerances.  It prints results at t = .1, .2, ..., 1.0.
C At the end of the run, various counters of interest are printed.
C
      DOUBLE PRECISION RTOL, ATOL, T0, T, TOUT, DTOUT, UNORM 
      DOUBLE PRECISION U(50), ROPT(40)
      INTEGER*4 IOPT(40)
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LNJE/16/
C
      MX = 10
      MY = 5
      NEQ = MX*MY
      T0 = 0.0D0
      CALL INITBX (MX, MY, U)
      METH = 2
      ITMETH = 2
      IATOL = 1
      RTOL = 0.0D0
      ATOL = 1.D-5
      INOPT = 0
      MU = MY
      ML = MY
      DTOUT = 0.1D0
      ITASK = 0
C
      WRITE(6,10)NEQ
 10   FORMAT('Band example problem: Advection-diffusion, NEQ = ',I2//)
C
      CALL FMENVINITS(NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' FMENVINITS returned IER =',I5)
        STOP
      ENDIF
C
      CALL FCVMALLOC (T0, U, METH, ITMETH, IATOL, RTOL, ATOL,
     1                INOPT, IOPT, ROPT, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' FCVMALLOC returned IER =',I5)
        STOP
        ENDIF
C
      CALL FCVBAND1 (NEQ, MU, ML, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' FCVBAND1 returned IER =',I5)
        STOP
      ENDIF
C
      CALL MAXNORM (NEQ, U, UNORM)
      WRITE(6,45)T,UNORM
 45   FORMAT(' At t =',F6.2,'   max.norm(u) =',D14.6)
C
      TOUT = DTOUT
      DO 70 IOUT = 1,10
C
        CALL FCVODE (TOUT, T, U, ITASK, IER)
C
        CALL MAXNORM (NEQ, U, UNORM)
        WRITE(6,50)T,UNORM, IOPT(LNST)
 50     FORMAT(' At t =',F6.2,'   max.norm(u) =',D14.6,'   NST =',I4)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER
 60       FORMAT(///' FCVODE returned IER =',I5)
          STOP
          ENDIF
C
 70     TOUT = TOUT + DTOUT
C
      WRITE(6,80) IOPT(LNST), IOPT(LNFE), IOPT(LNJE), IOPT(LNSETUP),
     1            IOPT(LNNI), IOPT(LNCF), IOPT(LNETF)
 80   FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
     1       '   No. J-s =',I4,'   No. LU-s =',I4/
     2       ' No. nonlinear iterations =',I4/
     3       ' No. nonlinear convergence failures =',I4/
     4       ' No. error test failures =',I4/)
C
      CALL FCVFREE
      CALL FMENVFREES
C
      STOP
      END

      SUBROUTINE INITBX (MESHX, MESHY, U0)
C Load Common with problem constants and U0 with initial values
      DOUBLE PRECISION U0(MESHY,MESHX)
      DOUBLE PRECISION DX, DY, HDCOEF, HACOEF, VDCOEF
      DOUBLE PRECISION XMAX, YMAX, X, Y
      DATA XMAX/2.0D0/, YMAX/1.0D0/
      COMMON /PAR/ DX, DY, HDCOEF, HACOEF, VDCOEF, MX, MY
C
C Load constants in Common.
      MX = MESHX
      MY = MESHY
      DX = XMAX/(MX + 1)
      DY = YMAX/(MY + 1)
      HDCOEF = 1.0D0/(DX*DX)
      HACOEF = 0.5D0/(2.0D0*DX)
      VDCOEF = 1.0D0/(DY*DY)
C
C Loop over grid and load initial values.
      DO 20 I = 1,MX
        X = I*DX
        DO 10 J = 1,MY
          Y = J*DY
          U0(J,I) = X*(XMAX - X)*Y*(YMAX - Y)*EXP(5.0D0*X*Y)
 10       CONTINUE
 20     CONTINUE
C
      RETURN
      END

      SUBROUTINE MAXNORM (N, U, UNORM)
C Compute max-norm of array U
      DOUBLE PRECISION U(*), UNORM, TEMP
      TEMP = 0.0D0
      DO 10 I = 1,N
 10     TEMP = MAX(ABS(U(I)),TEMP)
      UNORM = TEMP
      RETURN
      END

      SUBROUTINE CVFUN (T, U, UDOT)
C Right-hand side routine
      DOUBLE PRECISION T, U(*), UDOT(*)
      DOUBLE PRECISION UIJ, UDN, UUP, ULT, URT, HDIFF, HADV, VDIFF
      DOUBLE PRECISION DX, DY, HDCOEF, HACOEF, VDCOEF
      COMMON /PAR/ DX, DY, HDCOEF, HACOEF, VDCOEF, MX, MY
C
C Loop over all grid points.
      DO 20 I = 1,MX
        IOFF = (I-1)*MY
        DO 10 J = 1,MY
C
C Extract u at x_i, y_j and four neighboring points.
          IJ = J+IOFF
          UIJ = U(IJ)
          UDN = 0.0D0
          IF (J .NE. 1)  UDN = U(IJ-1)
          UUP = 0.0D0
          IF (J .NE. MY) UUP = U(IJ+1)
          ULT = 0.0D0
          IF (I .NE. 1)  ULT = U(IJ-MY)
          URT = 0.0D0
          IF (I .NE. MX) URT = U(IJ+MY)
C
C Set diffusion and advection terms and load into UDOT.
          HDIFF = HDCOEF*(ULT - 2.0D0*UIJ + URT)
          HADV = HACOEF*(URT - ULT)
          VDIFF = VDCOEF*(UUP - 2.0D0*UIJ + UDN)
          UDOT(IJ) = HDIFF + HADV + VDIFF
 10       CONTINUE
 20     CONTINUE
C
      RETURN
      END

      SUBROUTINE CVBJAC(N, MU, ML, MDIM, T, U, FU,
     1                 BJAC, V1, V2, V3)
C Load banded Jacobian
      DOUBLE PRECISION BJAC(MDIM,*)
      DOUBLE PRECISION DX, DY, HDCOEF, HACOEF, VDCOEF
      COMMON /PAR/ DX, DY, HDCOEF, HACOEF, VDCOEF, MX, MY
C
      MU1 = MU + 1
      MU2 = MU + 2
      MBAND = MU + 1 + ML
C
C Loop over all grid points.
      DO 20 I = 1,MX
        IOFF = (I-1)*MY
        DO 10 J = 1,MY
          K = J + IOFF
C
C Set Jacobian elements in column k of Jb.
          BJAC(MU1,K) = -2.0D0*(VDCOEF+HDCOEF)
          IF (I .NE. 1)  BJAC(1,K) = HDCOEF + HACOEF
          IF (I .NE. MX) BJAC(MBAND,K) = HDCOEF - HACOEF
          IF (J .NE. 1)  BJAC(MU,K) = VDCOEF
          IF (J .NE. MY) BJAC(MU2,K) = VDCOEF
C
 10       CONTINUE
 20     CONTINUE
C
      RETURN
      END
