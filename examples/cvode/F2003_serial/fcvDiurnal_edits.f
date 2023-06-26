C     --------------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2023, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
C     --------------------------------------------------------------------
C     FCVODE Example Problem: 2D kinetics-transport, precond. Krylov
C     solver.
C
C     An ODE system is generated from the following 2-species diurnal
C     kinetics advection-diffusion PDE system in 2 space dimensions:
C
C     dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
C                           + Ri(c1,c2,t)      for i = 1,2,   where
C     R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
C     R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
C     Kv(y) = Kv0*exp(y/5) ,
C     Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
C     vary diurnally.
C
C     The problem is posed on the square
C     0 .le. x .le. 20,    30 .le. y .le. 50   (all in km),
C     with homogeneous Neumann boundary conditions, and for time t
C     in 0 .le. t .le. 86400 sec (1 day).
C
C     The PDE system is treated by central differences on a uniform
C     10 x 10 mesh, with simple polynomial initial profiles.
C     The problem is solved with CVODE, with the BDF/GMRES method and
C     the block-diagonal part of the Jacobian as a left
C     preconditioner.
C
C     The preconditioner system is block-diagonal with 2x2 blocks, so
C     it is solved using the analytical inverse.
C
C     The second and third dimensions of U here must match the values
C     of MX and MY, for consistency with the output statements below.
C     --------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY
      PARAMETER (MX=10, MY=10)
      INTEGER*4 LENIPAR, LENRPAR
      PARAMETER (LENIPAR=6+2*MX*MY, LENRPAR=12+8*MX*MY)
C
      INTEGER*4 METH,IATOL,ITASK,IER,LNCFL,LNPS
      INTEGER*4 LNST,LNFE,LNSETUP,LNNI,LNCF,LQ,LH,LNPE,LNLI,LNETF
      INTEGER*4 JOUT,JPRETYPE,IGSTYPE,MAXL
C The following declaration specification should match C type long int.
      INTEGER*8 NEQ, IOUT(25), IPAR(LENIPAR)
      INTEGER*4 NST,NFE,NPSET,NPE,NPS,NNI,NETF
      INTEGER*4 NLI,NCFN,NCFL
      DOUBLE PRECISION ATOL,AVDIM,T,TOUT,TWOHR,RTOL,FLOOR,DELT
      DOUBLE PRECISION U(2,MX,MY),ROUT(10),RPAR(LENRPAR)
C
      DATA TWOHR/7200.0D0/, RTOL/1.0D-5/, FLOOR/100.0D0/,
     &     JPRETYPE/1/, IGSTYPE/1/, MAXL/0/, DELT/0.0D0/
      DATA LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/,
     &     LQ/9/, LNPE/20/, LNLI/22/, LNPS/21/, LNCFL/23/
      DATA LH/2/
C
C     Load problem constants into IPAR, RPAR, and set initial values
      CALL INITKX(MX, MY, U, IPAR, RPAR)
C
C     Set other input arguments.
      NEQ = 2*MX*MY
      T = 0.0D0
      METH = 2
      IATOL = 1
      ATOL = RTOL * FLOOR
      ITASK = 1
C
      WRITE(6,10) NEQ
 10   FORMAT('Krylov example problem:'//
     &       ' Kinetics-transport, NEQ = ', I4/)
C
C     Initialize vector specification
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
        STOP
      ENDIF
C
C     Initialize SPGMR linear solver module
      call FSUNSPGMRINIT(1, JPRETYPE, MAXL, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRINIT IER = ', I5)
        STOP
      ENDIF
      call FSUNSPGMRSETGSTYPE(1, IGSTYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,27) IER
 27     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETGSTYPE IER = ', I5)
        STOP
      ENDIF
C
C     Initialize CVODE
      CALL FCVMALLOC(T, U, METH, IATOL, RTOL, ATOL,
     &     IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
        STOP
      ENDIF
C
C     attach linear solver module to CVLs interface
      CALL FCVLSINIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FCVLSINIT returned IER = ',I5)
        CALL FCVFREE
        STOP
      ENDIF
C
C     attach preconditioner to CVLs interface
      CALL FCVLSSETPREC(1, IER)
C
C Loop over output points, call FCVODE, print sample solution values.
      TOUT = TWOHR
      DO JOUT = 1, 12
C
         CALL FCVODE(TOUT, T, U, ITASK, IER)
C
         WRITE(6,50) T, IOUT(LNST), IOUT(LQ), ROUT(LH)
 50      FORMAT(/' t = ', E11.3, 3X, 'nst = ', I5,
     &           '  q = ', I2, '  h = ', E14.6)
         WRITE(6,55) U(1,1,1), U(1,5,5), U(1,10,10),
     &               U(2,1,1), U(2,5,5), U(2,10,10)
 55      FORMAT('  c1 (bot.left/middle/top rt.) = ', 3E14.6/
     &          '  c2 (bot.left/middle/top rt.) = ', 3E14.6)
C
         IF (IER .NE. 0) THEN
            WRITE(6,60) IER, IOUT(15)
 60         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     &             '                 Linear Solver returned IER = ', I5)
            CALL FCVFREE
            STOP
         ENDIF
C
         TOUT = TOUT + TWOHR
C
      ENDDO

C     Print final statistics.
      NST = IOUT(LNST)
      NFE = IOUT(LNFE)
      NPSET = IOUT(LNSETUP)
      NPE = IOUT(LNPE)
      NPS = IOUT(LNPS)
      NNI = IOUT(LNNI)
      NLI = IOUT(LNLI)
      AVDIM = DBLE(NLI) / DBLE(NNI)
      NCFN = IOUT(LNCF)
      NCFL = IOUT(LNCFL)
      NETF = IOUT(LNETF)
      WRITE(6,80) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &     NCFL, NETF
  80  FORMAT(//'Final statistics:'//
     &     ' number of steps        = ', I5, 5X,
     &     ' number of f evals.     =', I5/
     &     ' number of prec. setups = ', I5/
     &     ' number of prec. evals. = ', I5, 5X,
     &     ' number of prec. solves = ', I5/
     &     ' number of nonl. iters. = ', I5, 5X,
     &     ' number of lin. iters.  = ', I5/
     &     ' average Krylov subspace dimension (NLI/NNI)  = ', E14.6/
     &     ' number of conv. failures.. nonlinear = ', I3,
     &     ' linear = ', I3/
     &     ' number of error test failures = ', I3)
C
      CALL FCVFREE
C
      STOP
      END

C     ----------------------------------------------------------------

      SUBROUTINE INITKX(MX, MY, U0, IPAR, RPAR)
C     Routine to set problem constants and initial values
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION RPAR(*)
C
      INTEGER*4 MM, JY, JX, P_IPP, P_BD, P_P
      DOUBLE PRECISION U0
      DIMENSION U0(2,MX,MY)
      DOUBLE PRECISION Q1, Q2, Q3, Q4, A3, A4, OM, C3, DY, HDCO
      DOUBLE PRECISION VDCO, HACO, X, Y
      DOUBLE PRECISION CX, CY, DKH, DKV0, DX, HALFDA, PI, VEL
C
      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/,
     1     PI/3.1415926535898D0/
C
C     Problem constants
      MM = MX * MY
      Q1 = 1.63D-16
      Q2 = 4.66D-16
      Q3 = 0.0D0
      Q4 = 0.0D0
      A3 = 22.62D0
      A4 = 7.601D0
      OM = PI / HALFDA
      C3 = 3.7D16
      DX = 20.0D0 / (MX - 1.0D0)
      DY = 20.0D0 / (MY - 1.0D0)
      HDCO = DKH / DX**2
      HACO = VEL / (2.0D0 * DX)
      VDCO = (1.0D0 / DY**2) * DKV0
C
C     Load constants in IPAR and RPAR
      IPAR(1) = MX
      IPAR(2) = MY
      IPAR(3) = MM
C
      RPAR(1)  = Q1
      RPAR(2)  = Q2
      RPAR(3)  = Q3
      RPAR(4)  = Q4
      RPAR(5)  = A3
      RPAR(6)  = A4
      RPAR(7)  = OM
      RPAR(8)  = C3
      RPAR(9)  = DY
      RPAR(10) = HDCO
      RPAR(11) = VDCO
      RPAR(12) = HACO
C
C     Pointers into IPAR and RPAR
      P_IPP = 7
      P_BD  = 13
      P_P   = P_BD + 4*MM
C
      IPAR(4) = P_IPP
      IPAR(5) = P_BD
      IPAR(6) = P_P
C
C     Set initial profiles.
      DO JY = 1, MY
        Y = 30.0D0 + (JY - 1.0D0) * DY
        CY = (0.1D0 * (Y - 40.0D0))**2
        CY = 1.0D0 - CY + 0.5D0 * CY**2
        DO JX = 1, MX
          X = (JX - 1.0D0) * DX
          CX = (0.1D0 * (X - 10.0D0))**2
          CX = 1.0D0 - CX + 0.5D0 * CX**2
          U0(1,JX,JY) = 1.0D6 * CX * CY
          U0(2,JX,JY) = 1.0D12 * CX * CY
         ENDDO
      ENDDO
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVFUN(T, U, UDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function f
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, U(2,*), UDOT(2,*), RPAR(*)
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      INTEGER*4 IER
C
      INTEGER*4 ILEFT, IRIGHT
      INTEGER*4 JX, JY, MX, MY, MM, IBLOK0, IBLOK, IDN, IUP
      DOUBLE PRECISION Q1, Q2, Q3, Q4, A3, A4, OM, C3, DY, HDCO
      DOUBLE PRECISION VDCO, HACO
      DOUBLE PRECISION C1, C2, C1DN, C2DN, C1UP, C2UP, C1LT, C2LT
      DOUBLE PRECISION C1RT, C2RT, CYDN, CYUP, HORD1, HORD2, HORAD1
      DOUBLE PRECISION HORAD2, QQ1, QQ2, QQ3, QQ4, RKIN1, RKIN2, S
      DOUBLE PRECISION VERTD1, VERTD2, YDN, YUP
C
C     Extract constants from IPAR and RPAR
      MX = IPAR(1)
      MY = IPAR(2)
      MM = IPAR(3)
C
      Q1 = RPAR(1)
      Q2 = RPAR(2)
      Q3 = RPAR(3)
      Q4 = RPAR(4)
      A3 = RPAR(5)
      A4 = RPAR(6)
      OM = RPAR(7)
      C3 = RPAR(8)
      DY = RPAR(9)
      HDCO = RPAR(10)
      VDCO = RPAR(11)
      HACO = RPAR(12)
C
C     Set diurnal rate coefficients.
      S = SIN(OM * T)
      IF (S .GT. 0.0D0) THEN
         Q3 = EXP(-A3 / S)
         Q4 = EXP(-A4 / S)
      ELSE
         Q3 = 0.0D0
         Q4 = 0.0D0
      ENDIF
      RPAR(3) = Q3
      RPAR(4) = Q4
C
C     Loop over all grid points.
      DO JY = 1, MY
         YDN = 30.0D0 + (JY - 1.5D0) * DY
         YUP = YDN + DY
         CYDN = VDCO * EXP(0.2D0 * YDN)
         CYUP = VDCO * EXP(0.2D0 * YUP)
         IBLOK0 = (JY - 1) * MX
         IDN = -MX
         IF (JY .EQ. 1) IDN = MX
         IUP = MX
         IF (JY .EQ. MY) IUP = -MX
         DO JX = 1, MX
            IBLOK = IBLOK0 + JX
            C1 = U(1,IBLOK)
            C2 = U(2,IBLOK)
C     Set kinetic rate terms.
            QQ1 = Q1 * C1 * C3
            QQ2 = Q2 * C1 * C2
            QQ3 = Q3 * C3
            QQ4 = Q4 * C2
            RKIN1 = -QQ1 - QQ2 + 2.0D0 * QQ3 + QQ4
            RKIN2 = QQ1 - QQ2 - QQ4
C     Set vertical diffusion terms.
            C1DN = U(1,IBLOK + IDN)
            C2DN = U(2,IBLOK + IDN)
            C1UP = U(1,IBLOK + IUP)
            C2UP = U(2,IBLOK + IUP)
            VERTD1 = CYUP * (C1UP - C1) - CYDN * (C1 - C1DN)
            VERTD2 = CYUP * (C2UP - C2) - CYDN * (C2 - C2DN)
C     Set horizontal diffusion and advection terms.
            ILEFT = -1
            IF (JX .EQ. 1) ILEFT = 1
            IRIGHT = 1
            IF (JX .EQ. MX) IRIGHT = -1
            C1LT = U(1,IBLOK + ILEFT)
            C2LT = U(2,IBLOK + ILEFT)
            C1RT = U(1,IBLOK + IRIGHT)
            C2RT = U(2,IBLOK + IRIGHT)
            HORD1 = HDCO * (C1RT - 2.0D0 * C1 + C1LT)
            HORD2 = HDCO * (C2RT - 2.0D0 * C2 + C2LT)
            HORAD1 = HACO * (C1RT - C1LT)
            HORAD2 = HACO * (C2RT - C2LT)
C     Load all terms into UDOT.
            UDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1
            UDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2
         ENDDO
      ENDDO
C
      IER = 0
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVPSET(T, U, FU, JOK, JCUR, GAMMA, H,
     &                   IPAR, RPAR, IER)
C     Routine to set and preprocess block-diagonal preconditioner.
C     Note: The dimensions in /BDJ/ below assume at most 100 mesh points.
C
      IMPLICIT NONE
C
      INTEGER*4 IER, JOK, JCUR
      DOUBLE PRECISION T, U(2,*), FU(*), GAMMA, H
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION RPAR(*)
C
      INTEGER*4 MX, MY, MM, P_IPP, P_BD, P_P
      DOUBLE PRECISION Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO
C
      IER = 0
C
C     Extract constants from IPAR and RPAR
      MX = IPAR(1)
      MY = IPAR(2)
      MM = IPAR(3)
C
      Q1 = RPAR(1)
      Q2 = RPAR(2)
      Q3 = RPAR(3)
      Q4 = RPAR(4)
      C3 = RPAR(8)
      DY = RPAR(9)
      HDCO = RPAR(10)
      VDCO = RPAR(11)
C
C     Extract pointers into IPAR and RPAR
      P_IPP = IPAR(4)
      P_BD  = IPAR(5)
      P_P   = IPAR(6)
C
C     If needed, recompute BD
C
      IF (JOK .EQ. 1) THEN
C     JOK = 1. Reuse saved BD
        JCUR = 0
      ELSE
C     JOK = 0. Compute diagonal Jacobian blocks.
C     (using q4 value computed on last FCVFUN call).
         CALL PREC_JAC(MX, MY, MM, U, RPAR(P_BD),
     &        Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO)
         JCUR = 1
      ENDIF
C
C     Copy BD to P and scale by -GAMMA
      P_P = -GAMMA*P_BD
C
C     Perform LU decomposition
      CALL PREC_LU(MM, RPAR(P_P), IPAR(P_IPP), IER)
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVPSOL(T, U, FU, R, Z, GAMMA, DELTA, LR,
     &                   IPAR, RPAR, IER)
C     Routine to solve preconditioner linear system.
C
      IMPLICIT NONE
C
      INTEGER*4 IER, LR
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, U(*), FU(*), R(*), Z(2,*)
      DOUBLE PRECISION GAMMA, DELTA, RPAR(*)
C
      INTEGER*4 MM, P_IPP, P_P
C
      IER = 0
C
C     Extract constants from IPAR and RPAR
      MM  = IPAR(3)
C
C     Extract pointers into IPAR and RPAR
      P_IPP = IPAR(4)
      P_P   = IPAR(6)
C
C     Copy RHS into Z
c$$$      CALL DCOPY(2*MM, R, 1, Z, 1)
      Z = R
C
C     Solve the block-diagonal system Px = r using LU factors stored in P
C     and pivot data in IPP, and return the solution in Z.
      CALL PREC_SOL(MM, RPAR(P_P), IPAR(P_IPP), Z)

      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_JAC(MX, MY, MM, U, BD,
     &     Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO)
C     Routine to compute diagonal Jacobian blocks
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY, MM
      DOUBLE PRECISION U(2,*), BD(2,2,MM)
      DOUBLE PRECISION Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO
C
      INTEGER*4 JY, JX, IBLOK, IBLOK0
      DOUBLE PRECISION C1, C2, CYDN, CYUP, DIAG, YDN, YUP
C
      DO JY = 1, MY
         YDN = 30.0D0 + (JY - 1.5D0) * DY
         YUP = YDN + DY
         CYDN = VDCO * EXP(0.2D0 * YDN)
         CYUP = VDCO * EXP(0.2D0 * YUP)
         DIAG = -(CYDN + CYUP + 2.0D0 * HDCO)
         IBLOK0 = (JY - 1) * MX
         DO JX = 1, MX
            IBLOK = IBLOK0 + JX
            C1 = U(1,IBLOK)
            C2 = U(2,IBLOK)
            BD(1,1,IBLOK) = (-Q1 * C3 - Q2 * C2) + DIAG
            BD(1,2,IBLOK) = -Q2 * C1 + Q4
            BD(2,1,IBLOK) =  Q1 * C3 - Q2 * C2
            BD(2,2,IBLOK) = (-Q2 * C1 - Q4) + DIAG
         ENDDO
      ENDDO

      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_LU(MM, P, IPP, IER)
C     Routine to perform LU decomposition on (P+I)
C
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*4 MM
      INTEGER*8 IPP(2,MM)
      DOUBLE PRECISION P(2,2,MM)
      DOUBLE PRECISION p11, p12, p21, p22, det
C
      INTEGER*4 I
C
C     Add identity matrix and do LU decompositions on blocks, in place.
      DO I = 1, MM
         p11 = P(1,1,I) + 1.0D0
         p22 = P(2,2,I) + 1.0D0
         p12 = P(1,2,I)
         p21 = P(1,2,I)
         det = p11*p22 - p12*p21
         IF (det .eq. 0.D0) RETURN

         P(1,1,I) = p22/det
         P(2,2,I) = p11/det
         P(1,2,I) = -p21/det
         P(2,1,I) = -p12/det
      ENDDO
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_SOL(MM, P, IPP, Z)
C     Routine for backsolve
C
      IMPLICIT NONE
C
      INTEGER*4 MM
      INTEGER*8 IPP(2,MM)
      DOUBLE PRECISION P(2,2,MM), Z(2,MM)
      DOUBLE PRECISION z1, z2
C
      INTEGER*4 I
C
      DO I = 1, MM
         z1 = Z(1,I)
         z2 = Z(2,I)
         Z(1,I) = P(1,1,I) * z1 + P(1,2,I) * z2
         Z(2,I) = P(2,1,I) * z1 + P(2,2,I) * z2
      ENDDO

      RETURN
      END
