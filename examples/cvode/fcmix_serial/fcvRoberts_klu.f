!     ----------------------------------------------------------------
!     Programmer(s): Ting Yan @ SMU
!          Based on cvRoberts_klu.c and modified to Fortran 77
!     ----------------------------------------------------------------
!     Copyright (c) 2016, Southern Methodist University.
!     All rights reserved.
!     For details, see the LICENSE file.
!     ----------------------------------------------------------------
!     $Revision: $
!     $Date: $
!     ----------------------------------------------------------------
!     FCVODE Example Problem: Robertson kinetics
!
!     The following is a simple example problem, with the coding
!     needed for its solution by CVODE. The problem is from chemical
!     kinetics, and consists of the following three rate equations:
!
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
!
!     on the interval from t = 0.0 to t = 4.e10, with initial
!     conditions:
!
!     y1 = 1.0, y2 = y3 = 0.
!
!     The problem is stiff. While integrating the system, we also
!     use the root finding feature to find the points at which
!     y1 = 1.e-4 or at which y3 = 0.01. The following coding solves
!     this problem with CVODE, using the Fortran/C interface routine
!     package. This solution uses the BDF method, Newton iteration with
!     the the FCVKLU sparse direct linear solver,and a user-supplied
!     Jacobian routine.
!     It uses a scalar relative tolerance and a vector absolute
!     tolerance. Output is printed in decades from t = .4 to t = 4.e10.
!     Run statistics (optional outputs) are printed at the end.
!     ----------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER IER, I, IERROOT
      INTEGER LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE, LNGE
      INTEGER METH, ITMETH, ITOL, ITASK, JOUT, NOUT
      INTEGER INFO(2)
!     The following declaration specification should match C type
      INTEGER*8 IPAR, NEQ, IOUT(22), IVAL, MXNLI, MXETF
      INTEGER INEQ, NNZ, SPARSETYPE, ORDERING
      DOUBLE PRECISION RTOL, T, T0, TOUT, H0, NLCONV
      DOUBLE PRECISION Y(3), ATOL(3), ROUT(6), RPAR
!
      DATA LNST/3/, LNFE/4/, LNETF/5/, LNCF/6/, LNNI/7/, LNSETUP/8/
      DATA LNGE/12/, LNJE/17/

!     Problem constant
      NEQ = 3                   ! number of equations
      Y(1) = 1.0D0              ! initial y components
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      METH = 2                  ! basic integration method, 2 for BDF
      ITMETH = 2                ! nonlinear iteration method, 2 for Newton iteration
      ITOL = 2                  ! type for absolute tolerance, 2 for array
      RTOL = 1.0D-4             ! scalar relative tolerance
      ATOL(1) = 1.0D-6          ! vector absolute tolerance components
      ATOL(2) = 1.0D-12
      ATOL(3) = 1.0D-4
      T0 = 0.0D0                ! initial time
      TOUT = 0.4D0              ! first output time
      ITASK = 1                 ! Using in call FCVODE, 1 for normal mode
      JOUT = 0
      NOUT = 12                 ! number of output times
!     
      WRITE(6, 10) NEQ
 10   FORMAT('Klu example problem:'//
     1     'Robertson kinetics, NEQ = ', I2//)

!
      
!     Create serial vector of length NEQ for I.C.
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 20) IER
 20      FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
         STOP
      ENDIF
!     

!     Call FCVMALLOC to create the solver memory and specify the 
!     Backward Differentiation Formula and the use of a Newton iteration
      CALL FCVMALLOC(T0, Y, METH, ITMETH, ITOL, RTOL, ATOL,
     1               IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 30) IER
 30      FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
         STOP
      ENDIF
!
!     Set the FCVODE input
      IVAL = 1000         ! max no. of internal steps before t_out
      CALL FCVSETIIN('MAX_NSTEPS', IVAL, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 31) IER
 31      FORMAT(///' SUNDIALS_ERROR: FCVSETIIN returned IER = ', I5)
         STOP
      ENDIF

      MXETF = 20          ! max no. of error test failures
      CALL FCVSETIIN('MAX_ERRFAIL', MXETF, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 31) IER
         STOP
      ENDIF
      
      H0 = 1.0D-4 * RTOL  ! initial step size
      CALL FCVSETRIN('INIT_STEP', H0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 32) IER
 32      FORMAT(///' SUNDIALS_ERROR: FCVSETRIN returned IER = ', I5)
         STOP
      ENDIF
      
      NLCONV = 1.0D-4     ! coefficient in the nonlinear convergence test
      CALL FCVSETRIN('NLCONV_COEF', NLCONV, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 32) IER
         STOP
      ENDIF
!

!     Call FCVROOTINIT to specify the root function g with 2 components
      CALL FCVROOTINIT(2, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 45) IER
 45      FORMAT(///' SUNDIALS_ERROR: FCVROOTINIT returned IER = ' I5)
         CALL FCVFREE
         STOP
      ENDIF
!
      
!     Call FCVKLU to specify the FCVKLU sparse direct linear solver
      INEQ = NEQ              ! convert to 'normal' integer type
      NNZ = INEQ * INEQ       ! maximum number of nonzeros in the Jac matrix
      SPARSETYPE = 0          ! CSC 
      ORDERING = 1            ! matrix ordering desired, 1=COLAMD
      CALL FCVKLU(INEQ, NNZ, SPARSETYPE, ORDERING, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 40) IER
 40      FORMAT(///' SUNDIALS_ERROR: FCVKLU returned IER = ', I5)
         CALL FCVFREE
         STOP
      ENDIF

      CALL FCVSPARSESETJAC(IER)
!     

!     In loop, call FCVODE, print results, and test for error.
      DO WHILE(JOUT .LT. NOUT)
!
         CALL FCVODE(TOUT, T, Y, ITASK, IER)
!     
         WRITE(6, 50) T, Y(1), Y(2), Y(3)
 50      FORMAT('At t = ', E12.4, '   y = ', 3E14.6)
!     
         IF(IER .LT. 0) THEN
            WRITE(6, 60) IER, IOUT(15)
 60         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER =  ', I5, /,
     1             '                 Linear Solver returned IER = ', I5)
            CALL FCVROOTFREE
            CALL FCVFREE
            STOP
         ENDIF
!
         IF(IER .EQ. 2) THEN
            CALL FCVROOTINFO(2, INFO, IERROOT)
            IF (IERROOT .LT. 0) THEN
               WRITE(6, 65) IERROOT
 65            FORMAT(///' SUNDIALS_ERROR: FCVROOTINFO returned IER = ',
     1                I5)
               CALL FCVROOTFREE
               CALL FCVFREE
               STOP
            ENDIF
            WRITE(6, 70) (INFO(I), I = 1, 2)
 70         FORMAT(5X, 'rootfound = ', 2I3)
         ENDIF
!
         IF (IER .EQ. 0) THEN
            TOUT = TOUT * 10.0D0
            JOUT = JOUT + 1
         ENDIF
!
      ENDDO

!     Obtain a derivative of the solution
      CAlL FCVDKY(T, 1, Y, IER)
      IF (IER .NE. 0) THEN
         WRITE(6, 80) IER
 80      FORMAT(///' SUNDIALS_ERROR: FCVDKY returned IER = ' I4)
         CALL FCVROOTFREE
         CALL FCVFREE
         STOP
      ENDIF
      WRITE(6, 85) Y(1), Y(2), Y(3)
 85   FORMAT(/'Final value of ydot = ', 3E14.6)

!     Print out the statistic information
      WRITE(6, 90) IOUT(LNST), IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP),
     1             IOUT(LNNI), IOUT(LNCF), IOUT(LNETF), IOUT(LNGE)
 90   FORMAT(//'Final statistics:'//
     1       ' No. steps = ', I4, '    No. f-s = ', I4
     2       '   No. J-s = ', I4, '    No. LU-s = ', I4/ 
     3       ' No. nonlinear iterations = ', I4/
     4       ' No. nonlinear convergence failures = ', I4/
     5       ' No. error test failures = ', I4/
     6       ' No. root function evals = ', I4)
!
      CALL FCVROOTFREE
      CALL FCVFREE
!
      STOP
      END

!     ----------------------------------------------------------------
      
      SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
!     Fortran routine for right-hand side function      
      IMPLICIT NONE
      
!     The following declaration specification should match C type long int
      INTEGER*8 IPAR(*)
      INTEGER IER
      DOUBLE PRECISION T, Y(*), YDOT(*), RPAR(*)
      DOUBLE PRECISION Y1, Y2, Y3

      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)
      
      YDOT(1) = -0.04D0 * Y1 + 1.0D4 * Y2 * Y3
      YDOT(3) = 3.0D7 * Y2 * Y2
      YDOT(2) = -YDOT(1) - YDOT(3)
!
      IER = 0
!
      RETURN
      END

!     ----------------------------------------------------------------
      
      SUBROUTINE FCVROOTFN(T, Y, G, IPAR, RPAR, IER)
! Fortran routine for root finding
      IMPLICIT NONE

!     The following declaration specification should match C type long int
      INTEGER*8 IPAR(*)
      INTEGER IER
      DOUBLE PRECISION T, Y(*), G(*), RPAR(*)
      DOUBLE PRECISION Y1, Y3
      
      Y1 = Y(1)
      Y3 = Y(3)

      G(1) = Y1 - 1.0D-4
      G(2) = Y3 - 1.0D-2
!
      IER = 0
!
      RETURN
      END

!     ----------------------------------------------------------------
      
      SUBROUTINE FCVSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS,
     1     JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)
!     Fortran routine for user-supplied CSC format KLU Jacobian
      IMPLICIT NONE

!     The following declaration specification should match C type long int
      INTEGER*8  IPAR(*)
      INTEGER N, NNZ, IER
      INTEGER JRVALS(NNZ), JCPTRS(N+1)
      DOUBLE PRECISION T, Y(*), FY(*), H, RPAR(*)
      DOUBLE PRECISION JDATA(NNZ)
      DOUBLE PRECISION WK1(*), WK2(*), WK3(*)
      
      ! Local data     
      DOUBLE PRECISION Y1, Y2, Y3

      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)

      JCPTRS(1) = 0
      JCPTRS(2) = 3
      JCPTRS(3) = 6
      JCPTRS(4) = 9

      JDATA(1) = -0.04D0
      JRVALS(1) = 0
      JDATA(2) = 0.04D0
      JRVALS(2) = 1
      JDATA(3) = 0.0D0
      JRVALS(3) = 2

      JDATA(4) = 1.0D4 * Y3
      JRVALS(4) = 0
      JDATA(5) = -1.0D4 * Y3 - 6.0D7 * Y2
      JRVALS(5) = 1
      JDATA(6) = 6.0D7 * Y2
      JRVALS(6) = 2
      
      JDATA(7) = 1.0D4 * Y2
      JRVALS(7) = 0
      JDATA(8) = -1.0D4 * Y2
      JRVALS(8) = 1
      JDATA(9) = 0.0D0
      JRVALS(9) = 2

      IER = 0
!
      RETURN
      END


!     ------End file------------------
