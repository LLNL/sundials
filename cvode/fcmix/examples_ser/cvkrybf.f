C     ----------------------------------------------------------------
C     $Revision: 1.4 $
C     $Date: 2004-09-08 21:06:17 $
C     ----------------------------------------------------------------
C     FCVODE Example Problem: 2D kinetics-transport, 
C     precond. Krylov solver. 
C     
C     An ODE system is generated from the following 2-species diurnal
C     kinetics advection-diffusion PDE system in 2 space dimensions:
C     
C     dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
C                           + Ri(c1,c2,t)      for i = 1,2,   where
C     R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
C     R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
C     Kv(z) = Kv0*exp(z/5) ,
C     Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
C     vary diurnally.
C
C     The problem is posed on the square
C     0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
C     with homogeneous Neumann boundary conditions, and for time t in
C     0 .le. t .le. 86400 sec (1 day).
C
C     The PDE system is treated by central differences on a uniform
C     10 x 10 mesh, with simple polynomial initial profiles.
C     The problem is solved with CVODE, with the BDF/GMRES method and
C     using the FCVBP banded preconditioner module
C     
C     The second and third dimensions of Y here must match the values of
C     MESHX and MESHY, for consistency with the output statements below.
C     ----------------------------------------------------------------
C
      DOUBLE PRECISION ATOL, AVDIM, DELT, FLOOR, RTOL, T, TOUT, TWOHR
      DOUBLE PRECISION Y(2,10,10), ROPT(40)
      INTEGER IOPT(40)
      DATA TWOHR/7200.0D0/, RTOL/1.0D-5/, FLOOR/100.0D0/,
     1     JPRETYPE/1/, IGSTYPE/1/, MAXL/0/, DELT/0.0D0/
      DATA LNST/4/, LNFE/5/, LNSETUP/6/, LNNI/7/, LNCF/8/, LNETF/9/,
     1     LQ/11/, LH/5/, LNPE/18/, LNLI/19/, LNPS/20/, LNCFL/21/
C     
C     Set mesh sizes
      MESHX = 10
      MESHZ = 10
C     Load Common and initial values in Subroutine INITKX
      CALL INITKX (MESHX, MESHZ, Y)
C     Set other input arguments.
      NEQ = 2*MESHX*MESHZ
      T = 0.0D0
      METH = 2
      ITMETH = 2
      IATOL = 1
      ATOL = RTOL*FLOOR
      INOPT = 0
      ITASK = 1
C
      WRITE(6,10)NEQ
 10   FORMAT('Krylov example problem: Kinetics-transport, NEQ = ',I4/)
C     
C     Initialize vector specification
      CALL FNVINITS (NEQ, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
 20      FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER =',I5)
         STOP
      ENDIF
C     
C     Initialize CVODE
      CALL FCVMALLOC (T, Y, METH, ITMETH, IATOL, RTOL, ATOL,
     1     INOPT, IOPT, ROPT, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
 30      FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER =',I5)
         CALL FNVFREES
         STOP
      ENDIF
C     
C     Initialize band preconditioner
      MU = 2
      ML = 2
      CALL FCVBPINIT(NEQ, MU, ML, IER) 
      IF (IER .NE. 0) THEN
         WRITE(6,40) IER
 40      FORMAT(///' SUNDIALS_ERROR: FCVBPINIT returned IER =',I5)
         CALL FNVFREES
         CALL FCVFREE
         STOP
      ENDIF
C     
C     Initialize SPGMR solver with band preconditioner
      CALL FCVBPSPGMR(JPRETYPE, IGSTYPE, MAXL, DELT, IER) 
      IF (IER .NE. 0) THEN
         WRITE(6,45) IER
 45      FORMAT(///' SUNDIALS_ERROR: FCVBPSPGMR returned IER =',I5)
         CALL FNVFREES
         CALL FCVFREE
         STOP
      ENDIF
C     
C     Loop over output points, call FCVODE, print sample solution values.
      TOUT = TWOHR
      DO 70 IOUT = 1,12
C
         CALL FCVODE (TOUT, T, Y, ITASK, IER)
C     
         WRITE(6,50) T,IOPT(LNST),IOPT(LQ),ROPT(LH)
 50      FORMAT(/' t =',D10.2,5X,'no. steps =',I5,
     1        '   order =',I3,'   stepsize =',D10.2)
         WRITE(6,55) Y(1,1,1), Y(1,5,5), Y(1,10,10),
     1        Y(2,1,1), Y(2,5,5), Y(2,10,10)
 55      FORMAT('  c1 (bot.left/middle/top rt.) =',3D12.3/
     1        '  c2 (bot.left/middle/top rt.) =',3D12.3)
C     
         IF (IER .NE. 0) THEN
            WRITE(6,60) IER
 60         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER =',I5)
            CALL FCVBPFREE
            CALL FNVFREES
            CALL FCVFREE
            STOP
         ENDIF
C     
         TOUT = TOUT + TWOHR
 70   CONTINUE
      
C     Print final statistics.
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
      WRITE (6,80) NST,NFE,NPSET,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL,NETF
 80   FORMAT(//'Final statistics..'/
     &     ' number of steps        =',I5,4X,
     &     ' number of f evals.     =',I5/
     &     ' number of prec. setups =',I5/
     &     ' number of prec. evals. =',I5,4X,
     &     ' number of prec. solves =',I5/
     &     ' number of nonl. iters. =',I5,4X,
     &     ' number of lin. iters.  =',I5/
     &     ' average Krylov subspace dimension (NLI/NNI)  =',F8.4/
     &     ' number of conv. failures..  nonlinear =',I3,
     &     ' linear =',I3/
     &     ' number of error test failures =',I3)
C     
      CALL FCVBPFREE
      CALL FCVFREE
      CALL FNVFREES
C     
      STOP
      END

      SUBROUTINE INITKX (MESHX, MESHZ, Y0)
C     Routine to set problem constants and initial values
      DIMENSION Y0(2, MESHX, MESHZ)
      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
      DOUBLE PRECISION CX, CZ, DKH, DKV0, DX, HALFDA, PI, VEL, 
     1     X, Y0, Z
      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/,
     1     PI/3.1415926535898D0/
C     
C     Load Common block of problem parameters.
      MX = MESHX
      MZ = MESHZ
      MM = MX*MZ
      Q1 = 1.63D-16
      Q2 = 4.66D-16
      A3 = 22.62D0
      A4 = 7.601D0
      OM = PI/HALFDA
      C3 = 3.7D16
      DX = 20.0D0/(MX - 1.0D0)
      DZ = 20.0D0/(MZ - 1.0D0)
      HDCO = DKH/DX**2
      HACO = VEL/(2.0D0*DX)
      VDCO = (1.0D0/DZ**2)*DKV0
C     
C     Set initial profiles.
      DO 20 JZ = 1,MZ
         Z = 30.0D0 + (JZ - 1.0D0)*DZ
         CZ = (0.1D0*(Z - 40.0D0))**2
         CZ = 1.0D0 - CZ + 0.5D0*CZ**2
         DO 10 JX = 1,MX
            X = (JX - 1.0D0)*DX
            CX = (0.1D0*(X - 10.0D0))**2
            CX = 1.0D0 - CX + 0.5D0*CX**2
            Y0(1,JX,JZ) = 1.0D6*CX*CZ
            Y0(2,JX,JZ) = 1.0D12*CX*CZ
 10      CONTINUE
 20   CONTINUE
C     
      RETURN
      END
      
      SUBROUTINE FCVFUN (T, Y, YDOT)
C     Routine for right-hand side function f
      DOUBLE PRECISION T, Y(2,*), YDOT(2,*)
      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
      DOUBLE PRECISION C1, C2, C1DN, C2DN, C1UP, C2UP, C1LT, C2LT,
     1     C1RT, C2RT, CZDN, CZUP, HORD1, HORD2, HORAD1, HORAD2,
     2     QQ1, QQ2, QQ3, QQ4, RKIN1, RKIN2, S, VERTD1, VERTD2, ZDN, ZUP
C     
C     Set diurnal rate coefficients.
      S = SIN(OM*T)
      IF (S .GT. 0.0D0) THEN
         Q3 = EXP(-A3/S)
         Q4 = EXP(-A4/S)
      ELSE
         Q3 = 0.0D0
         Q4 = 0.0D0
      ENDIF
C     
C     Loop over all grid points.
      DO 20 JZ = 1,MZ
         ZDN = 30.0D0 + (JZ - 1.5D0)*DZ
         ZUP = ZDN + DZ
         CZDN = VDCO*EXP(0.2D0*ZDN)
         CZUP = VDCO*EXP(0.2D0*ZUP)
         IBLOK0 = (JZ-1)*MX
         IDN = -MX
         IF (JZ .EQ. 1) IDN = MX
         IUP = MX
         IF (JZ .EQ. MZ) IUP = -MX
         DO 10 JX = 1,MX
            IBLOK = IBLOK0 + JX
            C1 = Y(1,IBLOK)
            C2 = Y(2,IBLOK)
C     Set kinetic rate terms.
            QQ1 = Q1*C1*C3
            QQ2 = Q2*C1*C2
            QQ3 = Q3*C3
            QQ4 = Q4*C2
            RKIN1 = -QQ1 - QQ2 + 2.0D0*QQ3 + QQ4
            RKIN2 = QQ1 - QQ2 - QQ4
C     Set vertical diffusion terms.
            C1DN = Y(1,IBLOK+IDN)
            C2DN = Y(2,IBLOK+IDN)
            C1UP = Y(1,IBLOK+IUP)
            C2UP = Y(2,IBLOK+IUP)
            VERTD1 = CZUP*(C1UP - C1) - CZDN*(C1 - C1DN)
            VERTD2 = CZUP*(C2UP - C2) - CZDN*(C2 - C2DN)
C     Set horizontal diffusion and advection terms.
            ILEFT = -1
            IF (JX .EQ. 1) ILEFT = 1
            IRIGHT = 1
            IF (JX .EQ. MX) IRIGHT = -1
            C1LT = Y(1,IBLOK+ILEFT)
            C2LT = Y(2,IBLOK+ILEFT)
            C1RT = Y(1,IBLOK+IRIGHT)
            C2RT = Y(2,IBLOK+IRIGHT)
            HORD1 = HDCO*(C1RT - 2.0D0*C1 + C1LT)
            HORD2 = HDCO*(C2RT - 2.0D0*C2 + C2LT)
            HORAD1 = HACO*(C1RT - C1LT)
            HORAD2 = HACO*(C2RT - C2LT)
C     Load all terms into YDOT.
            YDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1
            YDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END

